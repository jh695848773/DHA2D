#ifndef _VOXEL_MAP_TOOL_H_
#define _VOXEL_MAP_TOOL_H_

#include "math.h"
#include <Eigen/Eigen>
#include <iostream>
#include <ros/console.h>
#include <ros/ros.h>
#include <nav_msgs/OccupancyGrid.h>

class voxel_map_tool
{
  public:
    voxel_map_tool(){};
    ~voxel_map_tool()
    {
        if (data != nullptr)
        {
            delete[] data;
        }
    };

    bool initGridMap(const double &_resolution, const Eigen::Vector3d &map_size_xyz,
                     const Eigen::Vector3d &global_xyz_l);
    bool setObs(const double coord_x, const double coord_y, const double coord_z);
    bool setObs(const double coord_x, const double coord_y, const double coord_z, const int N_infla_grid);
    // All points out of boundary are considered Obstacles
    bool isObsFree(const double coord_x, const double coord_y, const double coord_z) const;
    bool isObsFree(int array_idx) const;
    bool isObsFree(Eigen::Vector3d pt) const;

    double getResolution() const
    {
        return resolution;
    }

    Eigen::Vector3i Coord2GridIdx(const Eigen::Vector3d &pt) const
    {
        Eigen::Vector3i idx;

        idx(0) = std::min(std::max(int((pt(0) - x_offest) * inv_resolution), 0), X_SIZE - 1),
        idx(1) = std::min(std::max(int((pt(1) - y_offest) * inv_resolution), 0), Y_SIZE - 1),
        idx(2) = std::min(std::max(int((pt(2) - z_offest) * inv_resolution), 0), Z_SIZE - 1);

        return idx;
    }

    Eigen::Vector3d GridIdx2Coord(const Eigen::Vector3i &index) const
    {
        Eigen::Vector3d pt;

        pt(0) = ((double)index(0) + 0.5) * resolution + x_offest;
        pt(1) = ((double)index(1) + 0.5) * resolution + y_offest;
        pt(2) = ((double)index(2) + 0.5) * resolution + z_offest;

        return pt;
    }

    Eigen::Vector3d GridIdx2Coord(const int px, const int py, const int pz) const
    {
        Eigen::Vector3d m_coord;
        m_coord(2) = (static_cast<double>(pz) + 0.5) * resolution + z_offest;
        m_coord(1) = (static_cast<double>(py) + 0.5) * resolution + y_offest;
        m_coord(0) = (static_cast<double>(px) + 0.5) * resolution + x_offest;
        return m_coord;
    }

    int GridIdx2Array(int idx_x, int idx_y, int idx_z) const
    {
        return idx_z * XY_SIZE + idx_y * X_SIZE + idx_x;
    }

    Eigen::Vector3d Array2Coord(const int idx) const
    {
        // To GridIdx
        int pz = idx / (XY_SIZE);
        int py = (idx - pz * XY_SIZE) / X_SIZE;
        int px = idx - pz * XY_SIZE - py * X_SIZE;

        Eigen::Vector3d m_coord = GridIdx2Coord(px, py, pz);

        return m_coord;
    }

    Eigen::Vector3d coordRounding(const Eigen::Vector3d &coord) const
    {
        return GridIdx2Coord(Coord2GridIdx(coord));
    }

    int getSize() const
    {
        return XYZ_SIZE;
    }

    std::string frame_name;

  protected:
    uint8_t *data = nullptr;

    int X_SIZE, Y_SIZE, Z_SIZE;
    int XYZ_SIZE, XY_SIZE;

    double resolution, inv_resolution;
    double x_offest, y_offest, z_offest;

    static constexpr uint8_t Unoccupied = 0;
    static constexpr uint8_t Occupied = 1;

  private:
};

inline bool voxel_map_tool::initGridMap(const double &_resolution, const Eigen::Vector3d &map_size_xyz,
                                        const Eigen::Vector3d &global_xyz_l)
{
    if (_resolution <= 0.0 || map_size_xyz(0) < 0.0 || map_size_xyz(1) < 0.0 || map_size_xyz(2) < 0.0)
        return false;

    resolution = _resolution;
    inv_resolution = 1.0 / resolution;

    X_SIZE = map_size_xyz(0) * inv_resolution;
    Y_SIZE = map_size_xyz(1) * inv_resolution;
    Z_SIZE = map_size_xyz(2) * inv_resolution;
    XY_SIZE = X_SIZE * Y_SIZE;
    XYZ_SIZE = Z_SIZE * XY_SIZE;

    data = new uint8_t[XYZ_SIZE];
    memset(data, Unoccupied, XYZ_SIZE * sizeof(uint8_t));

    x_offest = global_xyz_l(0);
    y_offest = global_xyz_l(1);
    z_offest = global_xyz_l(2);

    return true;
}

inline bool voxel_map_tool::setObs(const double coord_x, const double coord_y, const double coord_z)
{

    int idx_x = static_cast<int>((coord_x - x_offest) * inv_resolution);
    int idx_y = static_cast<int>((coord_y - y_offest) * inv_resolution);
    int idx_z = static_cast<int>((coord_z - z_offest) * inv_resolution);

    if (idx_x < 0 || idx_x >= X_SIZE || idx_y < 0 || idx_y >= Y_SIZE || idx_z < 0 || idx_z >= Z_SIZE)
        return false;

    data[GridIdx2Array(idx_x, idx_y, idx_z)] = Occupied;

    return true;
}

// If use this, the points on fringe of the map will be ignore (within N_infla_grid cells from the boundary)
inline bool voxel_map_tool::setObs(const double coord_x, const double coord_y, const double coord_z,
                                   const int N_infla_grid)
{

    int idx_x = static_cast<int>((coord_x - x_offest) * inv_resolution);
    int idx_y = static_cast<int>((coord_y - y_offest) * inv_resolution);
    int idx_z = static_cast<int>((coord_z - z_offest) * inv_resolution);

    if (idx_x < N_infla_grid || idx_x >= (X_SIZE - N_infla_grid) || idx_y < N_infla_grid ||
        idx_y >= (Y_SIZE - N_infla_grid) || idx_z < N_infla_grid || idx_z >= (Z_SIZE - N_infla_grid))
        return false;

    for (int i = -N_infla_grid; i <= N_infla_grid; ++i)
    {
        for (int j = -N_infla_grid; j <= N_infla_grid; ++j)
        {
            for (int k = -N_infla_grid; k <= N_infla_grid; ++k)
            {

                data[GridIdx2Array(idx_x + i, idx_y + j, idx_z + k)] = Occupied;
            }
        }
    }

    return true;
}

inline bool voxel_map_tool::isObsFree(const double coord_x, const double coord_y, const double coord_z) const
{

    int idx_x = (coord_x - x_offest) * inv_resolution;
    int idx_y = (coord_y - y_offest) * inv_resolution;
    int idx_z = (coord_z - z_offest) * inv_resolution;

    return (idx_x >= 0 && idx_x < X_SIZE && idx_y >= 0 && idx_y < Y_SIZE && idx_z >= 0 && idx_z < Z_SIZE &&
            (data[GridIdx2Array(idx_x, idx_y, idx_z)] == Unoccupied));
}

inline bool voxel_map_tool::isObsFree(int array_idx) const
{
    return (array_idx >= 0 && array_idx < XYZ_SIZE) && (data[array_idx] == Unoccupied);
}

inline bool voxel_map_tool::isObsFree(Eigen::Vector3d pt) const
{
    int idx_x = (pt(0) - x_offest) * inv_resolution;
    int idx_y = (pt(1) - y_offest) * inv_resolution;
    int idx_z = (pt(2) - z_offest) * inv_resolution;

    return (idx_x >= 0 && idx_x < X_SIZE && idx_y >= 0 && idx_y < Y_SIZE && idx_z >= 0 && idx_z < Z_SIZE &&
            (data[GridIdx2Array(idx_x, idx_y, idx_z)] == Unoccupied));
}

#endif