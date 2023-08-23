#include "HybridAStar/planner.hpp"
#include "trajectory.hpp"
#include "utils.h"
#include "voxel_map_tool2D.hpp"
#include <geometry_msgs/PoseWithCovarianceStamped.h>

double space_resolution = 0.5;
double map_x_size = 15.0; // m
double map_y_size = 15.0; // m

constexpr double DT = 0.01;

class Task
{
  public:
    Task(ros::NodeHandle &_nh, std::string _map_vis_name = "vis_map", std::string _path_vis_name = "path")
        : traj(_nh, _path_vis_name),
          obstacle_traj{{_nh, "obs_traj1"}, {_nh, "obs_traj2"}, {_nh, "obs_traj3"}, {_nh, "obs_traj4"},
                        {_nh, "obs_traj4"}, {_nh, "obs_traj4"}, {_nh, "obs_traj5"}},
          Curr_Obs_P{{5, 4},
                     {3, 3},
                     {1.5, map_y_size / 2},
                     {4, map_y_size - 3},
                     {map_y_size / 2, map_y_size - 3},
                     {map_x_size - 3, map_y_size - 3},
                     {map_x_size - 3, 4}},
          Curr_Obs_V{{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
          Obs_radius{1.0, 1.25, 1.25, 1.25, 1.50, 1.50, 2.0}
    //   Curr_Obs_P{{map_x_size / 2, map_y_size / 2}}, Curr_Obs_V{{0, 0}}, Obs_radius{1.5}
    {
        MapPtr = std::make_shared<voxel_map_tool>(_nh, _map_vis_name);
        // MapPtr = std::make_shared<voxel_map_tool>();
        MapPtr->initGridMap(space_resolution, Eigen::Vector2d{map_x_size, map_y_size}, Eigen::Vector2d::Zero());

        MapPtr->encircleBoundary();

        setObs_sub = _nh.subscribe("/initialpose", 10, &Task::setObsCB, this);
        setGoal_sub = _nh.subscribe("/move_base_simple/goal", 10, &Task::setGoalCB, this);
        SimulationLoopTimer = _nh.createTimer(ros::Duration(DT), &Task::SimulationLoop, this);
        CurrentP_pub = _nh.advertise<visualization_msgs::Marker>("robot_position", 10);
        Obs_CurrentP_pub = _nh.advertise<visualization_msgs::MarkerArray>("obs_position", 10);
        GoalP_pub = _nh.advertise<visualization_msgs::Marker>("goal_position", 10);
        Curr_P = HybridAStar::SVector{map_x_size / 2 - 1.5, map_y_size / 2 + 0.25};
        Curr_V = HybridAStar::SVector::Zero();
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
    bool obs_collision_flag = false;
    bool static_collision_flag = false;

    HybridAStar::SVector Curr_P;
    HybridAStar::SVector Curr_V;
    HybridAStar::SVector Goal_P;
    HybridAStar::V_A_Traj traj;
    HybridAStar::Planner planner;

    ros::Timer SimulationLoopTimer;
    std::vector<HybridAStar::StateVertex> path;

    std::vector<minco_trajectory::min_acc<HybridAStar::SpaceDim>> obstacle_traj;

    std::vector<HybridAStar::SVector> Curr_Obs_P;
    std::vector<HybridAStar::SVector> Curr_Obs_V;
    std::vector<double> Obs_radius;

    double Obs_V_End_max = 1.0;
    double safety_rate = 1.25;

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
std::vector<int> N_obs_frame_passing = {150, 150, 200, 200, 250, 250, 250};

void Task::SimulationLoop(const ros::TimerEvent &event)
{
    MapPtr->vis_map();

    static double t_on_path = 0;
    static std::vector<double> obs_t_on_path = {0, 0, 0, 0, 0, 0, 0};
    static int planner_frame_passing = 0;
    static std::vector<int> obs_frame_passing = {0, 0, 0, 0, 0, 0, 0};
    static std::vector<bool> has_obs_traj = {false, false, false, false, false, false, false};

    if (SceneFlag == 1)
        GoalSettingScene1();
    else if (SceneFlag == 0)
        if (has_goal == false)
            return;

    static bool last_collision_flag = false;
    static int collision_time = 0;
    static double time_pass = 0;
    static double energy_used = 0.0;
    static int planning_times = 0;
    static double total_planning_time = 0.0;
    static int loop_times = 0;
    static double total_dist2goal = 0.0;

    loop_times += 1;

    auto Goal_P_ptr = std::make_shared<HybridAStar::SVector>(map_x_size / 2 - 1.5, map_y_size / 2 + 0.25);
    auto Goal_V_ptr = std::make_shared<HybridAStar::SVector>(0, 0);

    if (Goal_P_ptr != nullptr)
        total_dist2goal += (Curr_P - *Goal_P_ptr).norm();

    // std::shared_ptr<HybridAStar::SVector> Goal_P_ptr = nullptr;
    // std::shared_ptr<HybridAStar::SVector> Goal_V_ptr = nullptr;

    auto expandObs = [](std::vector<double> Obs_radius, double safety_rate) {
        std::vector<double> expanded_Obs_radius;
        for (auto a : Obs_radius)

            expanded_Obs_radius.push_back(a * safety_rate);
        return expanded_Obs_radius;
    };

    /*----------------Planning----------------*/
    if (++planner_frame_passing > N_planner_frame_passing &&
        (planner_frame_passing > 5 * N_planner_frame_passing || new_goal || !has_path ||
         planner.checkPathCollision(MapPtr, Curr_Obs_P, Curr_Obs_V, expandObs(Obs_radius, safety_rate),
                                    space_resolution) == false))
    {
        planner_frame_passing = 0;
        Eigen::Matrix<double, HybridAStar::N_coeff * HybridAStar::N_poly, HybridAStar::SpaceDim> C;
        double AE_dT;

        std::clock_t c_start = std::clock();

        HybridAStar::SVector Goal_V = {0, 0};

        bool isSuccess =
            planner.SearchByHash(MapPtr, Curr_P, Curr_V, Goal_P_ptr, Goal_V_ptr, Curr_Obs_P, Curr_Obs_V, Obs_radius);

        planning_times += 1;

        std::clock_t c_end = std::clock();
        auto time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
        std::cout << "CPU time used for Hybrid A*: " << time_elapsed_ms << " ms\n";

        total_planning_time += time_elapsed_ms;

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
    for (int i = 0; i < Curr_Obs_P.size(); ++i)
    {
        if (++obs_frame_passing[i] > N_obs_frame_passing[i])
        {
            obs_frame_passing[i] = 0;
            std::vector<double> time_vector = {N_obs_frame_passing[i] * DT * 2.0};
            std::vector<HybridAStar::SVector> waypoints = {};
            std::vector<HybridAStar::SVector> start_pv = {Curr_Obs_P[i], Curr_Obs_V[i]};
            HybridAStar::SVector Goal_Obs_P;
            HybridAStar::SVector Goal_Obs_V;

            Goal_Obs_P = Curr_Obs_P[i] + HybridAStar::SVector::Random() * 5;
            Goal_Obs_V = Curr_Obs_V[i] + HybridAStar::SVector::Random() * 10;
            Goal_Obs_P(0) = std::max(std::min(Goal_Obs_P(0), map_x_size - Obs_radius[i]), 0.0 + Obs_radius[i]);
            Goal_Obs_P(1) = std::max(std::min(Goal_Obs_P(1), map_x_size - Obs_radius[i]), 0.0 + Obs_radius[i]);
            Goal_Obs_V(0) = std::max(std::min(Goal_Obs_V(0), Obs_V_End_max), -Obs_V_End_max);
            Goal_Obs_V(1) = std::max(std::min(Goal_Obs_V(1), Obs_V_End_max), -Obs_V_End_max);

            std::vector<HybridAStar::SVector> end_pv = {Goal_Obs_P, Goal_Obs_V};

            obstacle_traj[i].calculateTraj(time_vector, waypoints, start_pv, end_pv, 0.0);
            obstacle_traj[i].vis_traj();
            obs_t_on_path[i] = 0;
            has_obs_traj[i] = true;
        }
    }
    /*----------------Simulation----------------*/
    // // Move the Obs
    time_pass += DT;
    std::cout << "time_pass: " << time_pass / 60 << std::endl;

    std::cout << "ave planning time (ms): " << total_planning_time / planning_times << std::endl;

    if (Goal_P_ptr != nullptr)
        std::cout << "ave dist2goal: " << total_dist2goal / loop_times << std::endl;

    if (time_pass / 60 > 5 - DT && time_pass / 60 < 5 + DT)
    {
        std::cout << "Press Enter to Continue";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if (time_pass / 60 > 10 - DT && time_pass / 60 < 10 + DT)
    {
        std::cout << "Press Enter to Continue";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if (time_pass / 60 > 15 - DT && time_pass / 60 < 15 + DT)
    {
        std::cout << "Press Enter to Continue";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if (time_pass / 60 > 20 - DT && time_pass / 60 < 20 + DT)
    {
        std::cout << "Press Enter to Continue";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if (time_pass / 60 > 25 - DT && time_pass / 60 < 25 + DT)
    {
        std::cout << "Press Enter to Continue";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if (time_pass / 60 > 30 - DT && time_pass / 60 < 30 + DT)
    {
        std::cout << "Press Enter to Continue";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    for (int i = 0; i < Curr_Obs_P.size(); ++i)
    {
        if (has_obs_traj[i] == true)
        {
            obs_t_on_path[i] += DT;
            Curr_Obs_P[i] = obstacle_traj[i].getPosition(obs_t_on_path[i]);
            Curr_Obs_V[i] = obstacle_traj[i].getVelocity(obs_t_on_path[i]);
        }
    }

    // Collision check
    obs_collision_flag = false;
    if (last_collision_flag == false)
    {
        for (int i = 0; i < Curr_Obs_P.size(); ++i)
            if ((Curr_P - Curr_Obs_P[i]).norm() <= Obs_radius[i] + 0.1)
                obs_collision_flag = true;
    }
    else
    {
        for (int i = 0; i < Curr_Obs_P.size(); ++i)
            if ((Curr_P - Curr_Obs_P[i]).norm() <= Obs_radius[i] + 0.5)
                obs_collision_flag = true;
    }

    static_collision_flag = false;
    if (!MapPtr->isObsFree(Curr_P))
    {
        static_collision_flag = true;
    }

    if (last_collision_flag == false && (obs_collision_flag == true || static_collision_flag == true))
    {
        collision_time += 1;
    }
    std::cout << "collision_times: " << collision_time << std::endl;

    last_collision_flag = obs_collision_flag || static_collision_flag;

    //  || MapPtr->isObsFree(Curr_P) == false
    // Control the path
    if (has_path == true && obs_collision_flag == false)
    {
        t_on_path += DT;
        Curr_P = traj.getPosition(t_on_path);
        Curr_V = traj.getVelocity(t_on_path);
        auto acc = traj.getAcc(t_on_path);
        energy_used += acc.norm() * DT;
        std::cout << "energy_used: " << energy_used << std::endl;
    }

    /*----------------Visualization----------------*/
    visualization_msgs::Marker point_marker;

    point_marker.header.stamp = ros::Time::now();
    point_marker.type = visualization_msgs::Marker::SPHERE;
    point_marker.header.frame_id = "map";
    point_marker.pose.orientation.w = 1.00;
    point_marker.action = visualization_msgs::Marker::ADD;
    if (obs_collision_flag == false)
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
    point_marker.ns = "self";
    point_marker.pose.position.x = Curr_P.x();
    point_marker.pose.position.y = Curr_P.y();
    point_marker.pose.position.z = 0.0;

    CurrentP_pub.publish(point_marker);

    visualization_msgs::MarkerArray obs_point_marker_array;
    for (int i = 0; i < Curr_Obs_P.size(); ++i)
    {
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
        obs_point_marker.scale.x = Obs_radius[i] * 2;
        obs_point_marker.scale.y = Obs_radius[i] * 2;
        obs_point_marker.scale.z = 0.0;

        obs_point_marker.id = i;
        obs_point_marker.ns = "obs";
        obs_point_marker.pose.position.x = Curr_Obs_P[i].x();
        obs_point_marker.pose.position.y = Curr_Obs_P[i].y();
        obs_point_marker.pose.position.z = 0;
        obs_point_marker_array.markers.push_back(obs_point_marker);
    }
    Obs_CurrentP_pub.publish(obs_point_marker_array);

    if (Goal_P_ptr != nullptr)
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
        goal_point_marker.ns = "goal";
        goal_point_marker.pose.position.x = Goal_P_ptr->x();
        goal_point_marker.pose.position.y = Goal_P_ptr->y();
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
    for (double d = 2.5; d < map_y_size / 4 - 0.5; d += 1)
        SetInfla(2.5, d);

    for (double d = 0; d < map_x_size / 2; d += 1)
        SetInfla(d, map_y_size - 6.0);

    for (double d = 0; d < map_x_size / 2; d += 1)
        SetInfla(d, 6.0);

    SetInfla(5.0, map_y_size - 3);

    for (double d = 2.5; d < map_y_size / 2; d += 1)
        SetInfla(map_x_size - 3.0, map_y_size - d - 0.5);

    SetInfla(map_x_size - 3.5, 3.5);

    SetInfla(map_x_size / 2, 3.0);

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

    if ((Curr_P - Goal_P).norm() < 0.5)
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