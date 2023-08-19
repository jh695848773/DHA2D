#ifndef TRAJECTROY_HPP_
#define TRAJECTROY_HPP_

#include "ros/ros.h"
#include <pcl/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>
#include <ros/spinner.h>

#include <geometry_msgs/TransformStamped.h>
#include <pcl/common/transforms.h>

#include <tf/transform_listener.h>

#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

namespace minco_trajectory
{
template <int dim> class min_snap
{
  public:
    std::vector<double> time_table;
    double traj_start_time_sec;

    static constexpr int N_order = 7;
    static constexpr int N_coeff = N_order + 1;

    min_snap(ros::NodeHandle &nh_, std::string vis_name = "/visualizer/trajectory") : nh(nh_)
    {
        trajectoryPub = nh.advertise<visualization_msgs::Marker>(vis_name, 10);
    }

    inline void calculateTraj(const std::vector<double> &time_vector,
                              const std::vector<Eigen::Matrix<double, dim, Eigen::Dynamic>> &waypoints,
                              const std::vector<Eigen::Matrix<double, dim, 1>> &start_pvaj,
                              const std::vector<Eigen::Matrix<double, dim, 1>> &end_pvaj,
                              const double &_traj_start_time_sec)
    {
        const unsigned int N_poly = waypoints.size() + 1;

        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N_poly * N_coeff, N_poly * N_coeff);
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(N_poly * N_coeff, dim);

        M.block<1, N_coeff>(0, 0) = Beta(0, 0);
        M.block<1, N_coeff>(1, 0) = Beta(1, 0);
        M.block<1, N_coeff>(2, 0) = Beta(2, 0);
        M.block<1, N_coeff>(3, 0) = Beta(3, 0);

        D.block<1, dim>(0, 0) = start_pvaj[0];
        D.block<1, dim>(1, 0) = start_pvaj[1];
        D.block<1, dim>(2, 0) = start_pvaj[2];
        D.block<1, dim>(3, 0) = start_pvaj[3];

        for (unsigned int i = 0; i < N_poly - 1; ++i)
        {
            unsigned int j = 0;
            unsigned int k = 0;

            for (; k < waypoints[i].cols(); ++k)
            {
                M.block<1, N_coeff>(4 + N_coeff * i + j, N_coeff * i) = Beta(k, time_vector[i]);
                M.block<1, N_coeff>(4 + N_coeff * i + j, N_coeff * i + N_coeff) = Eigen::MatrixXd::Zero(1, N_coeff);
                D.block<1, dim>(4 + N_coeff * i + j, 0) = waypoints[i].col(k);
                j += 1;
            }

            for (unsigned int l = 0; l + k < 8; ++l)
            {
                M.block<1, N_coeff>(4 + N_coeff * i + j, N_coeff * i) = Beta(l, time_vector[i]);
                M.block<1, N_coeff>(4 + N_coeff * i + j, N_coeff * i + N_coeff) = -Beta(l, 0.0);
                j += 1;
            }
        }

        M.block<1, N_coeff>(N_coeff * N_poly - 4, N_coeff * N_poly - N_coeff) = Beta(0, time_vector.back());
        M.block<1, N_coeff>(N_coeff * N_poly - 3, N_coeff * N_poly - N_coeff) = Beta(1, time_vector.back());
        M.block<1, N_coeff>(N_coeff * N_poly - 2, N_coeff * N_poly - N_coeff) = Beta(2, time_vector.back());
        M.block<1, N_coeff>(N_coeff * N_poly - 1, N_coeff * N_poly - N_coeff) = Beta(3, time_vector.back());

        D.block<1, dim>(N_coeff * N_poly - 4, 0) = end_pvaj[0];
        D.block<1, dim>(N_coeff * N_poly - 3, 0) = end_pvaj[1];
        D.block<1, dim>(N_coeff * N_poly - 2, 0) = end_pvaj[2];
        D.block<1, dim>(N_coeff * N_poly - 1, 0) = end_pvaj[3];

        C = M.inverse() * D;

        time_table.clear();
        time_table.push_back(0);
        for (unsigned long i = 0; i < time_vector.size(); ++i)
        {
            time_table.push_back(time_table.back() + time_vector[i]);
        }
        traj_start_time_sec = _traj_start_time_sec;
    }

    inline Eigen::MatrixXd getPosition(double t)
    {
        if (t > time_table.back())
        {
            t = time_table.back() - 1e-10;
        }

        if (t < time_table.front() + 1e-10)
        {
            t = time_table.front() + 1e-10;
        }

        unsigned long i = 0;
        for (; i < time_table.size() - 1; ++i)
        {
            if (t > time_table.at(i) - 1e-10 && t < time_table.at(i + 1) + 1e-10)
            {
                break;
            }
        }

        Eigen::Matrix<double, dim, 1> current_position;

        for (int j = 0; j < dim; ++j)
            current_position(j) = C.block<N_coeff, 1>(N_coeff * i, j).dot(Beta(0, t - time_table[i]));

        return current_position;
    }

    inline void vis_traj()
    {
        visualization_msgs::Marker trajMarker;

        trajMarker.header.stamp = ros::Time::now();
        trajMarker.type = visualization_msgs::Marker::LINE_LIST;
        trajMarker.header.frame_id = "map";
        trajMarker.pose.orientation.w = 1.00;
        trajMarker.action = visualization_msgs::Marker::ADD;
        trajMarker.id = 0;
        trajMarker.ns = "trajectory";
        trajMarker.color.r = 0.00;
        trajMarker.color.g = 0.50;
        trajMarker.color.b = 1.00;
        trajMarker.color.a = 1.00;
        trajMarker.scale.x = 0.15;
        trajMarker.scale.y = 0.15;
        trajMarker.scale.z = 0.15;

        Eigen::Matrix<double, dim, 1> lastX = getPosition(0);

        double T = 0.05;
        for (double t = T; t < time_table.back(); t += T)
        {
            Eigen::Matrix<double, dim, 1> current_position = getPosition(t);

            geometry_msgs::Point point;

            point.x = lastX(0);
            point.y = lastX(1);
            if constexpr (dim > 2)
                point.z = lastX(2);
            trajMarker.points.push_back(point);
            point.x = current_position(0);
            point.y = current_position(1);
            if constexpr (dim > 2)
                point.z = current_position(2);
            trajMarker.points.push_back(point);
            lastX = current_position;
        }

        trajectoryPub.publish(trajMarker);
    }

  private:
    ros::NodeHandle nh;
    Eigen::MatrixXd C;

    ros::Publisher trajectoryPub;

    inline Eigen::Matrix<double, 8, 1> Beta(const unsigned int d, const double t)
    {
        Eigen::Matrix<double, 8, 1> Beta_t;

        if (d == 0)
        {
            Beta_t(0) = 1;
            Beta_t(1) = t;
            Beta_t(2) = t * t;
            Beta_t(3) = t * t * t;
            Beta_t(4) = t * t * t * t;
            Beta_t(5) = t * t * t * t * t;
            Beta_t(6) = t * t * t * t * t * t;
            Beta_t(7) = t * t * t * t * t * t * t;
        }
        else if (d == 1)
        {
            Beta_t(0) = 0;
            Beta_t(1) = 1;
            Beta_t(2) = 2 * t;
            Beta_t(3) = 3 * t * t;
            Beta_t(4) = 4 * t * t * t;
            Beta_t(5) = 5 * t * t * t * t;
            Beta_t(6) = 6 * t * t * t * t * t;
            Beta_t(7) = 7 * t * t * t * t * t * t;
        }
        else if (d == 2)
        {
            Beta_t(0) = 0;
            Beta_t(1) = 0;
            Beta_t(2) = 2;
            Beta_t(3) = 6 * t;
            Beta_t(4) = 12 * t * t;
            Beta_t(5) = 20 * t * t * t;
            Beta_t(6) = 30 * t * t * t * t;
            Beta_t(7) = 42 * t * t * t * t * t;
        }
        else if (d == 3)
        {

            Beta_t(0) = 0;
            Beta_t(1) = 0;
            Beta_t(2) = 0;
            Beta_t(3) = 6;
            Beta_t(4) = 24 * t;
            Beta_t(5) = 60 * t * t;
            Beta_t(6) = 120 * t * t * t;
            Beta_t(7) = 210 * t * t * t * t;
        }
        else if (d == 4)
        {
            Beta_t(0) = 0;
            Beta_t(1) = 0;
            Beta_t(2) = 0;
            Beta_t(3) = 0;
            Beta_t(4) = 24;
            Beta_t(5) = 120 * t;
            Beta_t(6) = 360 * t * t;
            Beta_t(7) = 840 * t * t * t;
        }
        else if (d == 5)
        {
            Beta_t(0) = 0;
            Beta_t(1) = 0;
            Beta_t(2) = 0;
            Beta_t(3) = 0;
            Beta_t(4) = 0;
            Beta_t(5) = 120;
            Beta_t(6) = 720 * t;
            Beta_t(7) = 2520 * t * t;
        }
        else if (d == 6)
        {
            Beta_t(0) = 0;
            Beta_t(1) = 0;
            Beta_t(2) = 0;
            Beta_t(3) = 0;
            Beta_t(4) = 0;
            Beta_t(5) = 0;
            Beta_t(6) = 720;
            Beta_t(7) = 5040 * t;
        }
        else if (d == 7)
        {
            Beta_t(0) = 0;
            Beta_t(1) = 0;
            Beta_t(2) = 0;
            Beta_t(3) = 0;
            Beta_t(4) = 0;
            Beta_t(5) = 0;
            Beta_t(6) = 0;
            Beta_t(7) = 5040;
        }
        else
        {
            Beta_t(0) = 0;
            Beta_t(1) = 0;
            Beta_t(2) = 0;
            Beta_t(3) = 0;
            Beta_t(4) = 0;
            Beta_t(5) = 0;
            Beta_t(6) = 0;
            Beta_t(7) = 0;
        }
        return Beta_t;
    }
};

template <int dim> class min_acc
{
  public:
    std::vector<double> time_table;
    double traj_start_time_sec;

    static constexpr int N_order = 3;
    static constexpr int N_coeff = N_order + 1;

    min_acc(ros::NodeHandle &nh_, std::string vis_name) : nh(nh_)
    {
        trajectoryPub = nh.advertise<visualization_msgs::Marker>(vis_name, 10);
    }

    inline void calculateTraj(const std::vector<double> &time_vector,
                              const std::vector<Eigen::Matrix<double, dim, 1>> &waypoints,
                              const std::vector<Eigen::Matrix<double, dim, 1>> &start_pv,
                              const std::vector<Eigen::Matrix<double, dim, 1>> &end_pv,
                              const double &_traj_start_time_sec)
    {
        const unsigned int N_poly = waypoints.size() + 1;

        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N_poly * N_coeff, N_poly * N_coeff);
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(N_poly * N_coeff, dim);

        M.block<1, N_coeff>(0, 0) = Beta(0, 0);
        M.block<1, N_coeff>(1, 0) = Beta(1, 0);

        D.block<1, dim>(0, 0) = start_pv[0];
        D.block<1, dim>(1, 0) = start_pv[1];

        for (unsigned int i = 0; i < N_poly - 1; ++i)
        {
            M.block<1, N_coeff>(2 + N_coeff * i + 0, N_coeff * i) = Beta(0, time_vector[i]);
            M.block<1, N_coeff>(2 + N_coeff * i + 1, N_coeff * i) = Beta(0, time_vector[i]);
            M.block<1, N_coeff>(2 + N_coeff * i + 2, N_coeff * i) = Beta(1, time_vector[i]);
            M.block<1, N_coeff>(2 + N_coeff * i + 3, N_coeff * i) = Beta(2, time_vector[i]);

            M.block<1, N_coeff>(2 + N_coeff * i + 0, N_coeff * i + N_coeff) = Eigen::Vector4d::Zero();
            M.block<1, N_coeff>(2 + N_coeff * i + 1, N_coeff * i + N_coeff) = -Beta(0, 0.0);
            M.block<1, N_coeff>(2 + N_coeff * i + 2, N_coeff * i + N_coeff) = -Beta(1, 0.0);
            M.block<1, N_coeff>(2 + N_coeff * i + 3, N_coeff * i + N_coeff) = -Beta(2, 0.0);

            D.block<1, dim>(2 + N_coeff * i + 0, 0) = waypoints[i];
        }

        M.block<1, N_coeff>(N_coeff * N_poly - 2, N_coeff * N_poly - N_coeff) = Beta(0, time_vector.back());
        M.block<1, N_coeff>(N_coeff * N_poly - 1, N_coeff * N_poly - N_coeff) = Beta(1, time_vector.back());

        D.block<1, dim>(N_coeff * N_poly - 2, 0) = end_pv[0];
        D.block<1, dim>(N_coeff * N_poly - 1, 0) = end_pv[1];

        C = M.inverse() * D;

        time_table.clear();
        time_table.push_back(0);
        for (unsigned long i = 0; i < time_vector.size(); ++i)
        {
            time_table.push_back(time_table.back() + time_vector[i]);
        }
        traj_start_time_sec = _traj_start_time_sec;
    }

    inline Eigen::MatrixXd getPosition(double t)
    {
        if (t > time_table.back())
        {
            t = time_table.back() - 1e-10;
        }

        if (t < time_table.front() + 1e-10)
        {
            t = time_table.front() + 1e-10;
        }

        unsigned long i = 0;
        for (; i < time_table.size() - 1; ++i)
        {
            if (t > time_table.at(i) - 1e-10 && t < time_table.at(i + 1) + 1e-10)
            {
                break;
            }
        }

        Eigen::Matrix<double, dim, 1> current_position;

        for (int j = 0; j < dim; ++j)
            current_position(j) = C.block<N_coeff, 1>(N_coeff * i, j).dot(Beta(0, t - time_table[i]));

        return current_position;
    }

    inline Eigen::MatrixXd getVelocity(double t)
    {
        if (t > time_table.back())
        {
            t = time_table.back() - 1e-10;
        }

        unsigned long i = 0;
        for (; i < time_table.size() - 1; ++i)
        {
            if (t >= time_table[i] && t < time_table[i + 1])
            {
                break;
            }
        }

        Eigen::Matrix<double, dim, 1> current_position;

        for (int j = 0; j < dim; ++j)
            current_position(j) = C.block<N_coeff, 1>(N_coeff * i, j).dot(Beta(1, t - time_table[i]));

        return current_position;
    }

    inline void vis_traj()
    {
        visualization_msgs::Marker trajMarker;

        trajMarker.header.stamp = ros::Time::now();
        trajMarker.type = visualization_msgs::Marker::LINE_LIST;
        trajMarker.header.frame_id = "map";
        trajMarker.pose.orientation.w = 1.00;
        trajMarker.action = visualization_msgs::Marker::ADD;
        trajMarker.id = 0;
        trajMarker.ns = "trajectory";
        trajMarker.color.r = 0.00;
        trajMarker.color.g = 0.50;
        trajMarker.color.b = 1.00;
        trajMarker.color.a = 1.00;
        trajMarker.scale.x = 0.15;
        trajMarker.scale.y = 0.15;
        trajMarker.scale.z = 0.15;

        Eigen::Matrix<double, dim, 1> lastX = getPosition(0);

        double T = 0.0005;
        for (double t = T; t < time_table.back(); t += T)
        {
            Eigen::Matrix<double, dim, 1> current_position = getPosition(t);

            geometry_msgs::Point point;

            point.x = lastX(0);
            point.y = lastX(1);
            if constexpr (dim > 2)
                point.z = lastX(2);
            trajMarker.points.push_back(point);
            point.x = current_position(0);
            point.y = current_position(1);
            if constexpr (dim > 2)
                point.z = current_position(2);
            trajMarker.points.push_back(point);
            lastX = current_position;
        }

        trajectoryPub.publish(trajMarker);
    }

  private:
    ros::NodeHandle nh;
    Eigen::MatrixXd C;

    ros::Publisher trajectoryPub;

    Eigen::Vector4d Beta(const unsigned int d, const double t)
    {
        if (d == 0)
        {
            return Eigen::Vector4d{1, t, t * t, t * t * t};
        }
        else if (d == 1)
        {
            return Eigen::Vector4d{0, 1, 2 * t, 3 * t * t};
        }
        else if (d == 2)
        {
            return Eigen::Vector4d{0, 0, 2, 6 * t};
        }
        else if (d == 3)
        {
            return Eigen::Vector4d{0, 0, 0, 6};
        }
        else
        {
            return Eigen::Vector4d{0, 0, 0, 0};
        }
    }
};
} // namespace minco_trajectory

#endif