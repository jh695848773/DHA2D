#include "HybridAStar/planner.hpp"
#include "utils.h"
#include "voxel_map_tool.hpp"
#include <geometry_msgs/PoseWithCovarianceStamped.h>

double space_resolution = 0.5;
double map_x_size = 10.0; // m
double map_y_size = 10.0; // m

class Task
{
  public:
    Task(ros::NodeHandle &_nh, std::string _map_vis_name = "vis_map") : Map(_nh, _map_vis_name)
    {
        Map.initGridMap(space_resolution, Eigen::Vector2d{map_x_size, map_y_size}, Eigen::Vector2d::Zero());
        setObs_sub = _nh.subscribe("/initialpose", 10, &Task::setObsCB, this);
    };

    void loop();
    void setObsCB(geometry_msgs::PoseWithCovarianceStamped::ConstPtr msg_ptr);

  private:
    voxel_map_tool Map;
    ros::Subscriber setObs_sub;
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
    while (ros::ok())
    {
        double freq = 10.0;
        ros::Rate r(freq);
        Map.vis_map();
        ros::spinOnce();
        r.sleep();
    }
}

void Task::setObsCB(geometry_msgs::PoseWithCovarianceStamped::ConstPtr msg_ptr)
{
    double x = msg_ptr->pose.pose.position.x;
    double y = msg_ptr->pose.pose.position.y;

    Map.setObs(x, y);
}