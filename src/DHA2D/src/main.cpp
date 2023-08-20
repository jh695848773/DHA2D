#include "HybridAStar/planner.hpp"
#include "trajectory.hpp"
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
        : traj(_nh, _path_vis_name), obstacle_traj(_nh, "obs_traj")
    {
        MapPtr = std::make_shared<voxel_map_tool>(_nh, _map_vis_name);
        // MapPtr = std::make_shared<voxel_map_tool>();
        MapPtr->initGridMap(space_resolution, Eigen::Vector2d{map_x_size, map_y_size}, Eigen::Vector2d::Zero());
        setObs_sub = _nh.subscribe("/initialpose", 10, &Task::setObsCB, this);
        setGoal_sub = _nh.subscribe("/move_base_simple/goal", 10, &Task::setGoalCB, this);
        SimulationLoopTimer = _nh.createTimer(ros::Duration(DT), &Task::SimulationLoop, this);
        CurrentP_pub = _nh.advertise<visualization_msgs::Marker>("robot_position", 10);
        Obs_CurrentP_pub = _nh.advertise<visualization_msgs::Marker>("obs_position", 10);
        GoalP_pub = _nh.advertise<visualization_msgs::Marker>("goal_position", 10);
        Curr_P = Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Ones() * 0.5;
        Curr_V = Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Zero();

        Curr_Obs_P(0) = map_x_size / 2;
        Curr_Obs_P(1) = map_y_size / 2;
    };

    void SimulationLoop(const ros::TimerEvent &event);
    void setObsCB(geometry_msgs::PoseWithCovarianceStamped::ConstPtr msg_ptr);
    void setGoalCB(geometry_msgs::PoseStamped::ConstPtr msg_ptr);
    void setScene1();
    void setScene0();

  private:
    std::shared_ptr<voxel_map_tool> MapPtr;
    ros::Subscriber setObs_sub;
    ros::Subscriber setGoal_sub;
    ros::Publisher CurrentP_pub;
    ros::Publisher GoalP_pub;
    ros::Publisher Obs_CurrentP_pub;
    bool has_goal = false;
    bool has_path = false;
    bool new_goal = false;
    bool collision_flag = false;
    bool has_obs_traj = false;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Curr_P;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Curr_V;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Goal_P;
    HybridAStar::V_A_Traj traj;
    HybridAStar::Planner planner;
    ros::Timer SimulationLoopTimer;
    std::vector<HybridAStar::StateVertex> path;

    minco_trajectory::min_acc<HybridAStar::SpaceDim> obstacle_traj;

    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Curr_Obs_P;
    Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Curr_Obs_V =
        Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Zero();

    double Obs_radius = 1.0;
    double Obs_V_End_max = 1.0;
    double safety_rate = 2.0;

    int SceneFlag = 0;

    void GoalSettingScene1();
};

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "main");
    ros::NodeHandle nh;

    Task task(nh);
    task.setScene1();

    ros::spin();

    return 0;
}

int N_planner_frame_passing = 10;
int N_obs_frame_passing = 150;

void Task::SimulationLoop(const ros::TimerEvent &event)
{
    MapPtr->vis_map();

    static double t_on_path = 0;
    static double obs_t_on_path = 0;
    static int planner_frame_passing = 0;
    static int obs_frame_passing = 0;

    if (SceneFlag == 1)
        GoalSettingScene1();
    else if (SceneFlag == 0)
        if (has_goal == false)
            return;

    /*----------------Planning----------------*/
    if (++planner_frame_passing > N_planner_frame_passing &&
        (planner_frame_passing > 5 * N_planner_frame_passing || new_goal || !has_path ||
         planner.checkPathCollision(MapPtr, std::vector<HybridAStar::SVector>{Curr_Obs_P},
                                    std::vector<HybridAStar::SVector>{Curr_Obs_V},
                                    std::vector<double>{Obs_radius * safety_rate}, space_resolution) == false))
    {
        planner_frame_passing = 0;
        Eigen::Matrix<double, HybridAStar::N_coeff * HybridAStar::N_poly, HybridAStar::SpaceDim> C;
        double AE_dT;

        std::clock_t c_start = std::clock();

        bool isSuccess =
            planner.SearchByHash(MapPtr, Curr_P, Curr_V, Goal_P, HybridAStar::SVector{0.0, 0.0},
                                 std::vector<HybridAStar::SVector>{Curr_Obs_P},
                                 std::vector<HybridAStar::SVector>{Curr_Obs_V}, std::vector<double>{Obs_radius});

        // bool isSuccess =
        //     planner.SearchByHash(MapPtr, Curr_P, Curr_V, Goal_P, HybridAStar::SVector{0.0, 0.0});

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

    /*----------------Update Obs movement----------------*/
    if (++obs_frame_passing > N_obs_frame_passing)
    {
        obs_frame_passing = 0;
        std::vector<double> time_vector = {N_obs_frame_passing * DT * 2.0};
        std::vector<Eigen::Matrix<double, HybridAStar::SpaceDim, 1>> waypoints = {};
        std::vector<Eigen::Matrix<double, HybridAStar::SpaceDim, 1>> start_pv = {Curr_Obs_P, Curr_Obs_V};
        Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Goal_Obs_P;
        Eigen::Matrix<double, HybridAStar::SpaceDim, 1> Goal_Obs_V;

        Goal_Obs_P = Curr_P + Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Random() * 5;
        Goal_Obs_V = Curr_Obs_V + Eigen::Matrix<double, HybridAStar::SpaceDim, 1>::Random() * 10;
        Goal_Obs_P(0) = std::max(std::min(Goal_Obs_P(0), map_x_size - Obs_radius), 0.0 + Obs_radius);
        Goal_Obs_P(1) = std::max(std::min(Goal_Obs_P(1), map_x_size - Obs_radius), 0.0 + Obs_radius);
        Goal_Obs_V(0) = std::max(std::min(Goal_Obs_V(0), Obs_V_End_max), -Obs_V_End_max);
        Goal_Obs_V(1) = std::max(std::min(Goal_Obs_V(1), Obs_V_End_max), -Obs_V_End_max);

        std::vector<Eigen::Matrix<double, HybridAStar::SpaceDim, 1>> end_pv = {Goal_Obs_P, Goal_Obs_V};

        obstacle_traj.calculateTraj(time_vector, waypoints, start_pv, end_pv, 0.0);
        obstacle_traj.vis_traj();
        obs_t_on_path = 0.0;
        has_obs_traj = true;
    }

    /*----------------Simulation----------------*/
    // // Move the Obs
    if (has_obs_traj == true)
    {
        obs_t_on_path += DT;
        Curr_Obs_P = obstacle_traj.getPosition(obs_t_on_path);
        Curr_Obs_V = obstacle_traj.getVelocity(obs_t_on_path);
    }

    // Collision check
    if ((Curr_P - Curr_Obs_P).norm() <= Obs_radius * 1.1)
    {
        collision_flag = true;
    }
    else
    {
        collision_flag = false;
    }

    // Control the path
    if (has_path == true && collision_flag == false)
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
    if (collision_flag == false)
    {
        point_marker.color.r = 0.00;
        point_marker.color.g = 1.00;
        point_marker.color.b = 1.00;
    }
    else
    {
        point_marker.color.r = 1.00;
        point_marker.color.g = 0.00;
        point_marker.color.b = 0.00;
    }

    point_marker.color.a = 1.00;
    point_marker.scale.x = 0.4;
    point_marker.scale.y = 0.4;
    point_marker.scale.z = 0.4;

    point_marker.id = 0;
    point_marker.ns = "trajectory";
    point_marker.pose.position.x = Curr_P.x();
    point_marker.pose.position.y = Curr_P.y();
    point_marker.pose.position.z = 0.0;

    CurrentP_pub.publish(point_marker);

    visualization_msgs::Marker obs_point_marker;

    obs_point_marker.header.stamp = ros::Time::now();
    obs_point_marker.type = visualization_msgs::Marker::SPHERE;
    obs_point_marker.header.frame_id = "map";
    obs_point_marker.pose.orientation.w = 1.00;
    obs_point_marker.action = visualization_msgs::Marker::ADD;

    obs_point_marker.color.r = 1.00;
    obs_point_marker.color.g = 1.00;
    obs_point_marker.color.b = 0.00;

    obs_point_marker.color.a = 1.00;
    obs_point_marker.scale.x = Obs_radius * 2;
    obs_point_marker.scale.y = Obs_radius * 2;
    obs_point_marker.scale.z = 0.0;

    obs_point_marker.id = 0;
    obs_point_marker.ns = "trajectory";
    obs_point_marker.pose.position.x = Curr_Obs_P.x();
    obs_point_marker.pose.position.y = Curr_Obs_P.y();
    obs_point_marker.pose.position.z = 0.0;

    Obs_CurrentP_pub.publish(obs_point_marker);

    {
        visualization_msgs::Marker goal_point_marker;

        goal_point_marker.header.stamp = ros::Time::now();
        goal_point_marker.type = visualization_msgs::Marker::SPHERE;
        goal_point_marker.header.frame_id = "map";
        goal_point_marker.pose.orientation.w = 1.00;
        goal_point_marker.action = visualization_msgs::Marker::ADD;

        goal_point_marker.color.r = 1.00;
        goal_point_marker.color.g = 1.00;
        goal_point_marker.color.b = 0.00;

        goal_point_marker.color.a = 1.00;
        goal_point_marker.scale.x = 0.4;
        goal_point_marker.scale.y = 0.4;
        goal_point_marker.scale.z = 0.4;

        goal_point_marker.id = 0;
        goal_point_marker.ns = "trajectory";
        goal_point_marker.pose.position.x = Goal_P.x();
        goal_point_marker.pose.position.y = Goal_P.y();
        goal_point_marker.pose.position.z = 0.0;

        GoalP_pub.publish(goal_point_marker);
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
    new_goal = true;
}

void Task::setScene0()
{
    SceneFlag = 0;
    return;
}

void Task::setScene1()
{
    auto SetInfla = [&](double x, double y) {
        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                MapPtr->setObs(x + i * space_resolution, y + j * space_resolution);
            }
        }
    };

    /*---------------Set Obs---------------*/
    for (double d = 0; d < map_y_size / 3; d += 1)
        SetInfla(1.5, d);

    for (double d = 0; d < map_x_size / 3; d += 1)
        SetInfla(d, map_y_size - 2.0);

    for (double d = 0; d < map_y_size / 3; d += 1)
        SetInfla(map_x_size - 2.0, map_y_size - d - 0.5);

    for (double d = 0; d < map_x_size / 3; d += 1)
        SetInfla(map_x_size - d - 0.5, 1.5);

    SetInfla(map_x_size / 2, map_y_size / 2);
    SceneFlag = 1;
}

void Task::GoalSettingScene1()
{
    static int goal_flag = 1;

    if (has_goal == false)
    {
        has_goal = true;
        Goal_P = Eigen::Vector2d{0.5, map_y_size - 0.5};
        goal_flag = 2;
    }

    if ((Curr_P - Goal_P).norm() < 0.1)
    {
        if (goal_flag == 1)
        {
            Goal_P = Eigen::Vector2d{0.5, map_y_size - 0.5};
            goal_flag = 2;
        }
        else if (goal_flag == 2)
        {
            Goal_P = Eigen::Vector2d{map_x_size - 0.5, map_y_size - 0.5};
            goal_flag = 3;
        }
        else if (goal_flag == 3)
        {
            Goal_P = Eigen::Vector2d{map_x_size - 0.5, 0.5};
            goal_flag = 4;
        }
        else if (goal_flag == 4)
        {
            Goal_P = Eigen::Vector2d{0.5, 0.5};
            goal_flag = 1;
        }
    }
}