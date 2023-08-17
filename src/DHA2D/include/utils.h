#ifndef _UTILS_H_
#define _UTILS_H_

#include <tf2_geometry_msgs/tf2_geometry_msgs.h>

#include <Eigen/Core>

#include <Eigen/Geometry>

#include <cmath>

#include <tf/transform_datatypes.h>

inline constexpr double Pi = 3.1415926;

namespace utils
{

inline double mod_angle_2Pi(const double &input_angle)
{
    double result = fmod(input_angle, 2.0 * Pi);

    if (result < 0)
    {
        result += 2.0 * Pi;
    }

    return result;
}

inline double angle_diff(const double &angle1, const double &angle2)
{
    double diff = fmod(mod_angle_2Pi(angle2) - mod_angle_2Pi(angle1), 2.0 * Pi);

    if (diff > Pi)
    {
        diff -= 2.0 * Pi;
    }

    if (diff < -Pi)
    {
        diff += 2.0 * Pi;
    }

    return diff;
}

inline double get_yaw_from_quat_msg(const geometry_msgs::Quaternion &quat_msg)
{
    tf2::Quaternion quat_tf;
    double roll, pitch, yaw;
    tf2::fromMsg(quat_msg, quat_tf);
    tf2::Matrix3x3(quat_tf).getRPY(roll, pitch, yaw);
    return mod_angle_2Pi(yaw);
}

inline Eigen::Matrix4d get_Eigen_T_from_vector3_quat_msg(const geometry_msgs::Vector3 &Vector3_msg,
                                                         const geometry_msgs::Quaternion &quat_msg)
{
    Eigen::Quaterniond q;
    q.w() = quat_msg.w;
    q.x() = quat_msg.x;
    q.y() = quat_msg.y;
    q.z() = quat_msg.z;

    auto R = q.toRotationMatrix();
    Eigen::Vector3d t(Vector3_msg.x, Vector3_msg.y, Vector3_msg.z);

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.block<3, 3>(0, 0) = R;
    T.block<3, 1>(0, 3) = t;

    return T;
}

inline geometry_msgs::Quaternion get_quat_msg_from_yaw(double yaw)
{
    // Normalize yaw to be within [0, 2*pi]
    yaw = mod_angle_2Pi(yaw);

    // Convert yaw angle (in radians) to quaternion
    tf::Quaternion quat = tf::createQuaternionFromYaw(yaw);

    // Convert tf::Quaternion to geometry_msgs::Quaternion and return
    geometry_msgs::Quaternion quat_msg;
    quaternionTFToMsg(quat, quat_msg);
    return quat_msg;
}

} // namespace utils

#endif