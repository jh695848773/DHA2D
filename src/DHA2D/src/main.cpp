#include "HybridAStar/planner.hpp"
#include "utils.h"
#include "voxel_map_tool2D.hpp"
#include <geometry_msgs/PoseWithCovarianceStamped.h>

double space_resolution = 0.5;
double map_x_size = 10.0; // m
double map_y_size = 10.0; // m

constexpr double DT = 0.001;

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
        CurrentP_pub = _nh.advertise<visualization_msgs::Marker>("robot_position", 10);

        Curr_P = Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Ones();
        Curr_V = Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Zero();
    };

    void SimulationLoop(const ros::TimerEvent &event);
    void setObsCB(geometry_msgs::PoseWithCovarianceStamped::ConstPtr msg_ptr);
    void setGoalCB(geometry_msgs::PoseStamped::ConstPtr msg_ptr);

  private:
    std::shared_ptr<voxel_map_tool> MapPtr;
    ros::Subscriber setObs_sub;
    ros::Subscriber setGoal_sub;
    ros::Publisher CurrentP_pub;
    bool has_goal = false;
    bool has_path = false;
    bool new_goal = false;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Curr_P;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Curr_V;
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

int N_planner_frame_passing = 100;

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
    if (++planner_frame_passing > N_planner_frame_passing &&
        (new_goal || !has_path || planner.checkPathCollision(MapPtr, space_resolution) == false))
    {
        planner_frame_passing = 0;
        Eigen::Matrix<double, HybridAStar::N_coeff * HybridAStar::N_poly, HybridAStar::SpaceDim> C;
        double AE_dT;

        std::clock_t c_start = std::clock();

        bool isSuccess = planner.SearchByHash(MapPtr, Curr_P, Curr_V, Goal_P);

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
            new_goal = false;
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
        Curr_P = traj.getPosition(t_on_path);
        Curr_V = traj.getVelocity(t_on_path);
    }

    /*----------------Visualization----------------*/
    visualization_msgs::Marker point_marker;

    point_marker.header.stamp = ros::Time::now();
    point_marker.type = visualization_msgs::Marker::SPHERE;
    point_marker.header.frame_id = "map";
    point_marker.pose.orientation.w = 1.00;
    point_marker.action = visualization_msgs::Marker::ADD;
    point_marker.color.r = 0.00;
    point_marker.color.g = 1.00;
    point_marker.color.b = 1.00;
    point_marker.color.a = 1.00;
    point_marker.scale.x = 0.4;
    point_marker.scale.y = 0.4;
    point_marker.scale.z = 0.0;

    point_marker.id = 0;
    point_marker.ns = "trajectory";
    point_marker.pose.position.x = Curr_P.x();
    point_marker.pose.position.y = Curr_P.y();
    point_marker.pose.position.z = 0.0;

    CurrentP_pub.publish(point_marker);
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
    new_goal = true;
}