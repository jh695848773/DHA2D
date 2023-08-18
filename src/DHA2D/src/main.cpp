#include "HybridAStar/planner.hpp"
#include "utils.h"
#include "voxel_map_tool2D.hpp"
#include <geometry_msgs/PoseWithCovarianceStamped.h>

double space_resolution = 0.5;
double map_x_size = 10.0; // m
double map_y_size = 10.0; // m

class Task
{
  public:
    Task(ros::NodeHandle &_nh, std::string _map_vis_name = "vis_map", std::string _path_vis_name = "path")
        : path2show(_nh, _path_vis_name)
    {
        MapPtr = std::make_shared<voxel_map_tool>(_nh, _map_vis_name);
        // MapPtr = std::make_shared<voxel_map_tool>();
        MapPtr->initGridMap(space_resolution, Eigen::Vector2d{map_x_size, map_y_size}, Eigen::Vector2d::Zero());
        setObs_sub = _nh.subscribe("/initialpose", 10, &Task::setObsCB, this);
        setGoal_sub = _nh.subscribe("/move_base_simple/goal", 10, &Task::setGoalCB, this);
    };

    void loop();
    void setObsCB(geometry_msgs::PoseWithCovarianceStamped::ConstPtr msg_ptr);
    void setGoalCB(geometry_msgs::PoseStamped::ConstPtr msg_ptr);

  private:
    std::shared_ptr<voxel_map_tool> MapPtr;
    ros::Subscriber setObs_sub;
    ros::Subscriber setGoal_sub;
    bool has_goal = false;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Goal_P;
    HybridAStar::V_A_Traj path2show;
};

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "main");
    ros::NodeHandle nh;

    Task task(nh);

    task.loop();

    return 0;
}

void Task::loop()
{
    HybridAStar::Planner planner;
    double freq = 10.0;
    ros::Rate r(freq);

    while (ros::ok())
    {

        MapPtr->vis_map();

        if (has_goal == true)
        {
            std::vector<HybridAStar::StateVertex> path;
            Eigen::Matrix<double, HybridAStar::N_coeff * HybridAStar::N_poly, HybridAStar::SpaceDim> C;
            double AE_dT;

            std::clock_t c_start = std::clock();

            bool isSuccess =
                planner.SearchByHash(MapPtr, Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Ones(),
                                     Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Zero(), Goal_P);

            std::clock_t c_end = std::clock();
            auto time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
            std::cout << "CPU time used for Hybrid A*: " << time_elapsed_ms << " ms\n";

            if (isSuccess)
            {
                std::cout << "Plannning Succeeded." << std::endl;
                planner.getPath(path, C, AE_dT);
                path2show.get(path, C, AE_dT);
                path2show.vis_path();
            }
            else
            {
                std::cout << "Planning Failed." << std::endl;
            }
        }

        ros::spinOnce();
        r.sleep();
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
            // MapPtr->setObs(x + i * space_resolution, y + j * space_resolution, 1);
        }
    }
}

void Task::setGoalCB(geometry_msgs::PoseStamped::ConstPtr msg_ptr)
{
    Goal_P(0) = msg_ptr->pose.position.x;
    Goal_P(1) = msg_ptr->pose.position.y;
    has_goal = true;
}