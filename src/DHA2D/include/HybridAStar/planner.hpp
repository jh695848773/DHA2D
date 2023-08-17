#ifndef _HA_PLANNER_HPP_
#define _HA_PLANNER_HPP_

#include "HybridAStar/CoreTools.hpp"

namespace HybridAStar
{
class Planner
{
  public:
    bool SearchByHash(std::shared_ptr<const voxel_map_tool> map_ptr, Eigen::Matrix<double, SpaceDim, 1> Start_P,
                      Eigen::Matrix<double, SpaceDim, 1> Start_V, Eigen::Matrix<double, SpaceDim, 1> Goal_P,
                      Eigen::Matrix<double, SpaceDim, 1> Goal_V = {0.0, 0.0}, const double &space_resolution = 0.5,
                      const double &time_resolution = 0.5, const double &PlanningRadius = 20.0,
                      const double time_horizon = 10.0)
    {
        /*-------------------Init the hash grid map-------------------*/
        StateVertexPtrGrid_hash state_vertex_ptr_grid;
        const Eigen::Vector3d map_size = {2.0 * PlanningRadius, 2.0 * PlanningRadius, time_horizon};
        Eigen::Vector3d m_offset;
        m_offset.head(SpaceDim) = Start_P - PlanningRadius * Eigen::Vector2d::Ones();
        m_offset(SpaceDim) = 0.0;
        state_vertex_ptr_grid.initGrid(space_resolution, time_resolution, m_offset);

        /*-------------------Init the priority queue-------------------*/
        std::multimap<double, StateVertex *> OpenSet;

        // Push the first element
        StateVertex *v_ptr = new StateVertex;
        v_ptr->setState(Start_P, Start_V, map_ptr, 0.0, Goal_P, Goal_V);
        state_vertex_ptr_grid.insert(v_ptr->P, 0.0, v_ptr);
        v_ptr->MultimapIt = OpenSet.insert(std::make_pair(v_ptr->f, v_ptr));
        v_ptr->status = InOpenSet;

        while (OpenSet.empty() == false)
        {
            // Pop the first element
            v_ptr = OpenSet.begin()->second;
            v_ptr->status = InClosedSet;
            OpenSet.erase(OpenSet.begin());

            // Check if we can directly connect the goal from the current point
            const double Dist2Goal = (v_ptr->P - Goal_P).norm();
            if (Dist2Goal < MaxAE_Dist)
            {
                AE_dT = v_ptr->time2go;
                Eigen::Matrix<double, N_coeff * N_poly, SpaceDim> TempC;
                if (AnalyticExpansion(v_ptr->P, v_ptr->V, Goal_P, Goal_V, TempC, map_ptr, AE_dT) == true)
                {
                    C = TempC;
                    v_ptr->getPathFromHere(path_points);

                    has_path = true;
                    return true;
                }
            }

            // Expand
            std::vector<Eigen::Vector2d> a_list;
            std::vector<double> dT_list;
            getInput(a_list, dT_list);

            for (auto &a : a_list)
            {
                for (auto &used_dT : dT_list)
                {
                    StateVertex *new_v_ptr = new StateVertex;

                    new_v_ptr->simForward(*v_ptr, a, used_dT, Goal_P, Goal_V);

                    // If it's out of the planning radius or out of the time horizon or hit obstacle
                    if (((new_v_ptr->P - Start_P).norm() > PlanningRadius) || (new_v_ptr->time_stamp > time_horizon) ||
                        (new_v_ptr->checkCollision(map_ptr, space_resolution) == false))
                    {
                        delete new_v_ptr;
                    }
                    else
                    {
                        std::unordered_map<Eigen::Vector3i, StateVertex *, matrix_hash<Eigen::Vector3i>>::iterator
                            hitted_grid_iter;
                        bool isExist =
                            state_vertex_ptr_grid.find(new_v_ptr->P, new_v_ptr->time_stamp, hitted_grid_iter);
                        if (isExist == false)
                        {
                            state_vertex_ptr_grid.insert(new_v_ptr->P, new_v_ptr->time_stamp, new_v_ptr);
                            new_v_ptr->MultimapIt = OpenSet.insert(std::make_pair(new_v_ptr->f, new_v_ptr));
                            new_v_ptr->status = InOpenSet;
                        }
                        else if (hitted_grid_iter->second->status == InOpenSet)
                        {
                            if (new_v_ptr->cost2come < hitted_grid_iter->second->cost2come)
                            {
                                OpenSet.erase(hitted_grid_iter->second->MultimapIt);
                                delete hitted_grid_iter->second;
                                hitted_grid_iter->second = new_v_ptr;

                                new_v_ptr->MultimapIt = OpenSet.insert(std::make_pair(new_v_ptr->f, new_v_ptr));
                                new_v_ptr->status = InOpenSet;
                            }
                            else
                            {
                                delete new_v_ptr;
                            }
                        }
                        else // The hitted grid is in the Closed Set
                        {
                            delete new_v_ptr;
                        }
                    }
                }
            }
        }

        return false;
    }

    bool getPath(std::vector<StateVertex> &_path_points, Eigen::Matrix<double, N_coeff * N_poly, SpaceDim> &_C,
                 double &_AE_dT) // C: the last section of the path)
    {
        if (has_path == false)
        {
            return false;
        }

        _path_points = path_points;
        _C = C;
        _AE_dT = AE_dT;
        return true;
    }

    bool checkPathCollision(std::shared_ptr<const voxel_map_tool> map_ptr, const double &resolution = 0.5)
    {
        if (has_path == false)
        {
            return false;
        }

        for (auto ptr = path_points.begin(); ptr < path_points.cend(); ++ptr)
        {
            if (ptr->checkCollision(map_ptr, resolution) == false)
            {
                return false;
            }
        }

        const double IncreT = AE_dT / AE_CheckSteps;
        for (double t = IncreT; t < AE_dT + 1e-10; t += IncreT)
        {
            Eigen::Vector2d current_position;

            current_position(0) = C.block<N_coeff, 1>(0, 0).dot(Beta(0, t));
            current_position(1) = C.block<N_coeff, 1>(0, 1).dot(Beta(0, t));

            if (map_ptr->isObsFree(current_position) == false)
            {
                return false;
            }
        }

        return true;
    }

  private:
    std::vector<StateVertex> path_points;
    Eigen::Matrix<double, N_coeff * N_poly, SpaceDim> C; // C: the last section of the path
    double AE_dT;

    void getInput(std::vector<Eigen::Vector2d> &a_list, std::vector<double> &dT_list);

    Eigen::Vector4d Beta(const unsigned int d, const double &t);

    bool AnalyticExpansion(const Eigen::Matrix<double, SpaceDim, 1> &Start_P,
                           const Eigen::Matrix<double, SpaceDim, 1> &Start_V,
                           const Eigen::Matrix<double, SpaceDim, 1> &End_P,
                           const Eigen::Matrix<double, SpaceDim, 1> &End_V,
                           Eigen::Matrix<double, N_coeff * N_poly, SpaceDim> &C,
                           std::shared_ptr<const voxel_map_tool> map_ptr, double UseAE_dT);
    bool has_path = false;
};

inline void Planner::getInput(std::vector<Eigen::Vector2d> &a_list, std::vector<double> &dT_list)
{
    Eigen::Vector2d a;
    const Eigen::Vector2d IncreA = 2 * max_a / (a_steps - 1);
    for (a(0) = -max_a(0); a(0) < max_a(0) + 1e-10; a(0) += IncreA(0))
    {
        for (a(1) = -max_a(1); a(1) < max_a(1) + 1e-10; a(1) += IncreA(1))
        {
            a_list.push_back(a);
        }
    }

    const double IncredT = max_dT / dT_steps;
    for (double t = IncredT; t < max_dT + 1e-10; t += IncredT)
        dT_list.push_back(t);
}

inline Eigen::Vector4d Planner::Beta(const unsigned int d, const double &t)
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

inline bool Planner::AnalyticExpansion(const Eigen::Matrix<double, SpaceDim, 1> &Start_P,
                                       const Eigen::Matrix<double, SpaceDim, 1> &Start_V,
                                       const Eigen::Matrix<double, SpaceDim, 1> &End_P,
                                       const Eigen::Matrix<double, SpaceDim, 1> &End_V,
                                       Eigen::Matrix<double, N_coeff * N_poly, SpaceDim> &C,
                                       std::shared_ptr<const voxel_map_tool> map_ptr, double UseAE_dT)
{
    Eigen::Matrix<double, N_coeff * N_poly, N_coeff * N_poly> M;
    Eigen::Matrix<double, N_poly * N_coeff, SpaceDim> D;

    M.block<1, N_coeff>(0, 0) = Beta(0, 0);
    M.block<1, N_coeff>(1, 0) = Beta(1, 0);

    M.block<1, N_coeff>(N_coeff * N_poly - 2, N_coeff * N_poly - N_coeff) = Beta(0, UseAE_dT);
    M.block<1, N_coeff>(N_coeff * N_poly - 1, N_coeff * N_poly - N_coeff) = Beta(1, UseAE_dT);

    D.block<1, SpaceDim>(0, 0) = Start_P;
    D.block<1, SpaceDim>(1, 0) = Start_V;

    D.block<1, SpaceDim>(N_coeff * N_poly - 2, 0) = End_P;
    D.block<1, SpaceDim>(N_coeff * N_poly - 1, 0) = End_V;

    C = M.inverse() * D;

    const double IncreT = UseAE_dT / AE_CheckSteps;
    for (double t = IncreT; t < UseAE_dT + 1e-10; t += IncreT)
    {
        Eigen::Vector2d current_position;

        current_position(0) = C.block<N_coeff, 1>(0, 0).dot(Beta(0, t));
        current_position(1) = C.block<N_coeff, 1>(0, 1).dot(Beta(0, t));

        if (map_ptr->isObsFree(current_position) == false)
        {
            return false;
        }
    }

    return true;
}

} // namespace HybridAStar

#endif