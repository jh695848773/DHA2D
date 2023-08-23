#ifndef _VOXEL_MAP_TOOL_H_
#define _VOXEL_MAP_TOOL_H_

#include "math.h"
#include <Eigen/Eigen>
#include <iostream>
#include <nav_msgs/OccupancyGrid.h>
#include <ros/console.h>
#include <ros/ros.h>

class voxel_map_tool
{
  public:
    voxel_map_tool(ros::NodeHandle &nh, const std::string &_vis_name)
    {
        vis_pub = nh.advertise<nav_msgs::OccupancyGrid>(_vis_name, 10);
    };

    voxel_map_tool(){};

    ~voxel_map_tool()
    {
        if (data != nullptr)
        {
            delete[] data;
        }
    };

    bool initGridMap(const double &_resolution, const Eigen::Vector2d &map_size_xz, const Eigen::Vector2d &global_xy_l);
    bool setObs(const double coord_x, const double coord_y);

    // All points out of boundary are considered Obstacles
    bool isObsFree(const double coord_x, const double coord_y) const;
    bool isObsFree(int array_idx) const;
    bool isObsFree(const Eigen::Vector2d &pt) const;
    uint8_t getCost(const double coord_x, const double coord_y) const;
    uint8_t getCost(const Eigen::Vector2d &pt) const;
    void encircleBoundary();

    void vis_map() const;

    double getResolution() const
    {
        return resolution;
    }

    Eigen::Vector2i Coord2GridIdx(const Eigen::Vector2d &pt) const
    {
        Eigen::Vector2i idx;

        idx(0) = std::min(std::max(int((pt(0) - x_offest) * inv_resolution), 0), X_SIZE - 1);
        idx(1) = std::min(std::max(int((pt(1) - y_offest) * inv_resolution), 0), Y_SIZE - 1);

        return idx;
    }

    Eigen::Vector2d GridIdx2Coord(const Eigen::Vector2i &index) const
    {
        Eigen::Vector2d pt;

        pt(0) = ((double)index(0) + 0.5) * resolution + x_offest;
        pt(1) = ((double)index(1) + 0.5) * resolution + y_offest;

        return pt;
    }

    Eigen::Vector3d GridIdx2Coord(const int px, const int py) const
    {
        Eigen::Vector3d m_coord;
        m_coord(0) = (static_cast<double>(px) + 0.5) * resolution + x_offest;
        m_coord(1) = (static_cast<double>(py) + 0.5) * resolution + y_offest;
        return m_coord;
    }

    int GridIdx2Array(int idx_x, int idx_y) const
    {
        return idx_y * X_SIZE + idx_x;
    }

    Eigen::Vector3d Array2Coord(const int idx) const
    {
        // To GridIdx
        int px = idx % X_SIZE; // idx == py * X_SIZE + px
        int py = (idx - px) / X_SIZE;

        Eigen::Vector3d m_coord = GridIdx2Coord(px, py);

        return m_coord;
    }

    Eigen::Vector2d coordRounding(const Eigen::Vector2d &coord) const
    {
        return GridIdx2Coord(Coord2GridIdx(coord));
    }

    int getSize() const
    {
        return XY_SIZE;
    }

    std::string frame_name;

  protected:
    uint8_t *data = nullptr;

    int X_SIZE, Y_SIZE;
    int XY_SIZE;

    double resolution, inv_resolution;
    double x_offest, y_offest;

    static constexpr uint8_t Unoccupied = 0;
    static constexpr uint8_t Occupied = 99;
    ros::Publisher vis_pub;

  private:
};

inline bool voxel_map_tool::initGridMap(const double &_resolution, const Eigen::Vector2d &map_size_xy,
                                        const Eigen::Vector2d &global_xy_l)
{
    if (_resolution <= 0.0 || map_size_xy(0) < 0.0 || map_size_xy(1) < 0.0)
        return false;

    resolution = _resolution;
    inv_resolution = 1.0 / resolution;

    X_SIZE = map_size_xy(0) * inv_resolution;
    Y_SIZE = map_size_xy(1) * inv_resolution;
    XY_SIZE = X_SIZE * Y_SIZE;

    data = new uint8_t[XY_SIZE];
    memset(data, Unoccupied, XY_SIZE * sizeof(uint8_t));

    x_offest = global_xy_l(0);
    y_offest = global_xy_l(1);

    return true;
}

inline bool voxel_map_tool::setObs(const double coord_x, const double coord_y)
{
    int idx_x = static_cast<int>((coord_x - x_offest) * inv_resolution);
    int idx_y = static_cast<int>((coord_y - y_offest) * inv_resolution);

    if (idx_x < 0 || idx_x >= X_SIZE || idx_y < 0 || idx_y >= Y_SIZE)
        return false;

    data[GridIdx2Array(idx_x, idx_y)] = Occupied;

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            if (idx_x + i < 0 || idx_y + j >= X_SIZE || idx_y < 0 || idx_y >= Y_SIZE)
                continue;
            data[GridIdx2Array(idx_x + i, idx_y + j)] =
                30 > data[GridIdx2Array(idx_x + i, idx_y + j)] ? 30 : data[GridIdx2Array(idx_x + i, idx_y + j)];
        }
    }

    return true;
}

inline void voxel_map_tool::encircleBoundary()
{
    int idx_x = 0;
    int idx_y = 0;
    idx_x = 0;
    for (idx_y = 0; idx_y < Y_SIZE; ++idx_y)
    {
        data[GridIdx2Array(idx_x, idx_y)] = Occupied;
        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                if (idx_x + i < 0 || idx_x + i >= X_SIZE || idx_y + j < 0 || idx_y + j >= Y_SIZE)
                    continue;
                data[GridIdx2Array(idx_x + i, idx_y + j)] =
                    30 > data[GridIdx2Array(idx_x + i, idx_y + j)] ? 30 : data[GridIdx2Array(idx_x + i, idx_y + j)];
            }
        }
    }

    idx_y = Y_SIZE - 1;
    for (idx_x = 0; idx_x < X_SIZE; ++idx_x)
    {
        data[GridIdx2Array(idx_x, idx_y)] = Occupied;
        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                if (idx_x + i < 0 || idx_x + i >= X_SIZE || idx_y + j < 0 || idx_y + j >= Y_SIZE)
                    continue;
                data[GridIdx2Array(idx_x + i, idx_y + j)] =
                    30 > data[GridIdx2Array(idx_x + i, idx_y + j)] ? 30 : data[GridIdx2Array(idx_x + i, idx_y + j)];
            }
        }
    }

    idx_x = X_SIZE - 1;
    for (idx_y = Y_SIZE - 1; idx_y >= 0; --idx_y)
    {
        data[GridIdx2Array(idx_x, idx_y)] = Occupied;
        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                if (idx_x + i < 0 || idx_x + i >= X_SIZE || idx_y + j < 0 || idx_y + j >= Y_SIZE)
                    continue;
                data[GridIdx2Array(idx_x + i, idx_y + j)] =
                    30 > data[GridIdx2Array(idx_x + i, idx_y + j)] ? 30 : data[GridIdx2Array(idx_x + i, idx_y + j)];
            }
        }
    }

    idx_y = 0;
    for (idx_x = X_SIZE - 1; idx_x >= 0; --idx_x)
    {
        data[GridIdx2Array(idx_x, idx_y)] = Occupied;
        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                if (idx_x + i < 0 || idx_x + i >= X_SIZE || idx_y + j < 0 || idx_y + j >= Y_SIZE)
                    continue;
                data[GridIdx2Array(idx_x + i, idx_y + j)] =
                    30 > data[GridIdx2Array(idx_x + i, idx_y + j)] ? 30 : data[GridIdx2Array(idx_x + i, idx_y + j)];
            }
        }
    }
}

inline bool voxel_map_tool::isObsFree(const double coord_x, const double coord_y) const
{
    int idx_x = (coord_x - x_offest) * inv_resolution;
    int idx_y = (coord_y - y_offest) * inv_resolution;

    return (idx_x >= 0 && idx_x < X_SIZE && idx_y >= 0 && idx_y < Y_SIZE &&
            (data[GridIdx2Array(idx_x, idx_y)] != Occupied));
}

inline bool voxel_map_tool::isObsFree(int array_idx) const
{
    return (array_idx >= 0 && array_idx < XY_SIZE) && (data[array_idx] != Occupied);
}

inline bool voxel_map_tool::isObsFree(const Eigen::Vector2d &pt) const
{
    int idx_x = (pt(0) - x_offest) * inv_resolution;
    int idx_y = (pt(1) - y_offest) * inv_resolution;

    return (idx_x >= 0 && idx_x < X_SIZE && idx_y >= 0 && idx_y < Y_SIZE &&
            (data[GridIdx2Array(idx_x, idx_y)] != Occupied));
}

inline uint8_t voxel_map_tool::getCost(const double coord_x, const double coord_y) const
{
    int idx_x = (coord_x - x_offest) * inv_resolution;
    int idx_y = (coord_y - y_offest) * inv_resolution;

    if (idx_x >= 0 && idx_x < X_SIZE && idx_y >= 0 && idx_y < Y_SIZE)
        return data[GridIdx2Array(idx_x, idx_y)];
    else
        return Occupied;
}

inline uint8_t voxel_map_tool::getCost(const Eigen::Vector2d &pt) const
{
    int idx_x = (pt(0) - x_offest) * inv_resolution;
    int idx_y = (pt(1) - y_offest) * inv_resolution;

    if (idx_x >= 0 && idx_x < X_SIZE && idx_y >= 0 && idx_y < Y_SIZE)
        return data[GridIdx2Array(idx_x, idx_y)];
    else
        return Occupied;
}

inline void voxel_map_tool::vis_map() const
{
    static int seq = 0;
    nav_msgs::OccupancyGrid msg2pub;
    msg2pub.header.frame_id = "map";
    msg2pub.header.stamp = ros::Time::now();
    msg2pub.header.seq = ++seq;

    msg2pub.info.resolution = resolution;
    msg2pub.info.width = X_SIZE;
    msg2pub.info.height = Y_SIZE;
    msg2pub.info.origin.position.x = x_offest;
    msg2pub.info.origin.position.y = y_offest;
    msg2pub.info.origin.position.z = 0.0;
    msg2pub.info.origin.orientation.w = 1.0;
    msg2pub.info.origin.orientation.x = 0.0;
    msg2pub.info.origin.orientation.y = 0.0;
    msg2pub.info.origin.orientation.z = 0.0;

    for (int j = 0; j < Y_SIZE; ++j)
    {
        for (int i = 0; i < X_SIZE; ++i)
        {
            // if (data[GridIdx2Array(i, j)] == Occupied)
            // {
            //     msg2pub.data.push_back(99);
            // }
            // else
            // {
            //     msg2pub.data.push_back(0);
            // }

            msg2pub.data.push_back(data[GridIdx2Array(i, j)]);
        }
    }

    vis_pub.publish(msg2pub);
}

#endif