#include "HybridAStar/planner.hpp"
#include "utils.h"
#include "voxel_map_tool2D.hpp"
#include <geometry_msgs/PoseWithCovarianceStamped.h>

double space_resolution = 0.5;
double map_x_size = 10.0; // m
double map_y_size = 10.0; // m

constexpr double DT = 0.01;

class Task
{
  public:
    Task(ros::NodeHandle &_nh, std::string _map_vis_name = "vis_map", std::string _path_vis_name = "path")
        : traj(_nh, _path_vis_name)
    {
        MapPtr = std::make_shared<voxel_map_tool>(_nh, _map_vis_name);
        // MapPtr = std::make_shared<voxel_map_tool>();
        MapPtr->initGridMap(space_resolution, Eigen::Vector2d{map_x_size, map_y_size}, Eigen::Vector2d::Zero());
        setObs_sub = _nh.subscribe("/initialpose", 10, &Task::setObsCB, this);
        setGoal_sub = _nh.subscribe("/move_base_simple/goal", 10, &Task::setGoalCB, this);
        SimulationLoopTimer = _nh.createTimer(ros::Duration(DT), &Task::SimulationLoop, this);

        Start_P = Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Ones();
        Start_V = Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Zero();
    };

    void SimulationLoop(const ros::TimerEvent &event);
    void setObsCB(geometry_msgs::PoseWithCovarianceStamped::ConstPtr msg_ptr);
    void setGoalCB(geometry_msgs::PoseStamped::ConstPtr msg_ptr);

  private:
    std::shared_ptr<voxel_map_tool> MapPtr;
    ros::Subscriber setObs_sub;
    ros::Subscriber setGoal_sub;
    bool has_goal = false;
    bool has_path = false;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Start_P;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Start_V;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Goal_P;
    HybridAStar::V_A_Traj traj;
    HybridAStar::Planner planner;
    ros::Timer SimulationLoopTimer;
    std::vector<HybridAStar::StateVertex> path;
};

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "main");
    ros::NodeHandle nh;

    Task task(nh);

    ros::spin();

    return 0;
}

int N_planner_frame_passing = 10;

void Task::SimulationLoop(const ros::TimerEvent &event)
{
    MapPtr->vis_map();

    static double t_on_path;
    static int planner_frame_passing = 0;

    if (has_goal == false)
    {
        return;
    }

    /*----------------Planning----------------*/
    if (++planner_frame_passing > N_planner_frame_passing)
    {
        planner_frame_passing = 0;
        Eigen::Matrix<double, HybridAStar::N_coeff * HybridAStar::N_poly, HybridAStar::SpaceDim> C;
        double AE_dT;

        std::clock_t c_start = std::clock();

        bool isSuccess = planner.SearchByHash(MapPtr, Start_P, Start_V, Goal_P);

        std::clock_t c_end = std::clock();
        auto time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
        std::cout << "CPU time used for Hybrid A*: " << time_elapsed_ms << " ms\n";

        if (isSuccess)
        {
            std::cout << "Plannning Succeeded." << std::endl;
            planner.getPath(path, C, AE_dT);
            traj.get(path, C, AE_dT);
            traj.vis_path();
            has_path = true;
            t_on_path = 0.0;
        }
        else
        {
            std::cout << "Planning Failed." << std::endl;
        }
    }

    /*----------------Simulation----------------*/
    if (has_path == true)
    {
        t_on_path += DT;
        Start_P = traj.getPosition(t_on_path);
        Start_V = traj.getVelocity(t_on_path);
    }
}

void Task::setObsCB(geometry_msgs::PoseWithCovarianceStamped::ConstPtr msg_ptr)
{
    double x = msg_ptr->pose.pose.position.x;
    double y = msg_ptr->pose.pose.position.y;

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            MapPtr->setObs(x + i * space_resolution, y + j * space_resolution);
        }
    }
}

void Task::setGoalCB(geometry_msgs::PoseStamped::ConstPtr msg_ptr)
{
    Goal_P(0) = msg_ptr->pose.position.x;
    Goal_P(1) = msg_ptr->pose.position.y;
    has_goal = true;
}