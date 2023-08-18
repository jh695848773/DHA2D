#ifndef _HA_CORE_TOOLS_HPP_
#define _HA_CORE_TOOLS_HPP_

#include <Eigen/Core>

#include "voxel_map_tool.hpp"

#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#include "utils.h"

namespace HybridAStar
{
inline constexpr double max_dT = 1.0;

inline constexpr int AE_CheckSteps = 15;

inline constexpr double MaxAE_Dist = 3.0;

inline constexpr unsigned int TimeDim = 1;

inline constexpr unsigned int SpaceDim = 2;

inline const Eigen::Matrix<double, SpaceDim, 1> V_max = {3.0, 3.0};
inline const Eigen::Matrix<double, SpaceDim, 1> max_a{3.0, 3.0};
inline constexpr double a_steps = 7;
inline constexpr double dT_steps = 2;

static_assert(a_steps >= 2, "a_steps should be >= 2");

inline constexpr double time_weight = 1.0;

inline constexpr char InOpenSet = 1;
inline constexpr char InClosedSet = 2;

inline constexpr unsigned int N_order = 3;
inline constexpr unsigned int N_coeff = N_order + 1;
inline constexpr unsigned int N_poly = 1;

inline constexpr double tie_breaker = 1.0 + 1.0 / 10000;

struct StateVertex
{
    Eigen::Matrix<double, SpaceDim, 1> P;
    Eigen::Matrix<double, SpaceDim, 1> V;

    Eigen::Matrix<double, SpaceDim, 1> at_1;

    struct StateVertex *prev_ptr;

    double dT;

    double time_stamp;

    double cost2go;
    double cost2come;
    double time2go;
    double f;

    bool isCollision = false;

    char status;

    std::multimap<double, StateVertex *>::iterator MultimapIt;

    StateVertex()
    {
    }
    bool simForward(StateVertex &prev, const Eigen::Matrix<double, SpaceDim, 1> &_at_1, const double &_dT,
                    Eigen::Matrix<double, SpaceDim, 1> Goal_P, Eigen::Matrix<double, SpaceDim, 1> Goal_V);

    bool checkCollision(std::shared_ptr<const voxel_map_tool> map_ptr, const double map_resolution);

    bool setState(const Eigen::Matrix<double, SpaceDim, 1> &Curr_P, const Eigen::Matrix<double, SpaceDim, 1> &Curr_V,
                  std::shared_ptr<const voxel_map_tool> map_ptr, const double &time,
                  Eigen::Matrix<double, SpaceDim, 1> Goal_P, Eigen::Matrix<double, SpaceDim, 1> Goal_V);

    void getPathFromHere(std::vector<StateVertex> &path_points);

    double estimateHeuristic(const Eigen::VectorXd &StartState, const Eigen::VectorXd &GoalState, double &optimal_time);

    std::vector<double> cubic(double a, double b, double c, double d);
    std::vector<double> quartic(double a, double b, double c, double d, double e);
};

template <typename T> struct matrix_hash : std::unary_function<T, size_t>
{
    std::size_t operator()(T const &matrix) const
    {
        size_t seed = 0;
        for (size_t i = 0; i < matrix.size(); ++i)
        {
            auto elem = *(matrix.data() + i);
            seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

class StateVertexPtrGrid_hash
{
  public:
    StateVertexPtrGrid_hash()
    {
    }
    ~StateVertexPtrGrid_hash()
    {
        for (auto &ptr2SV : DataPtrTable)
        {
            delete ptr2SV.second;
        }
    }

    bool find(const Eigen::Vector2d &coord, const double &time,
              std::unordered_map<Eigen::Vector3i, StateVertex *, matrix_hash<Eigen::Vector3i>>::iterator &iter);

    void insert(const Eigen::Vector2d &coord, const double &time, StateVertex *ptr); // coord: x, y, time

    bool initGrid(const double &_space_resolution, const double &_time_resolution, const Eigen::Vector3d &global_xyt_l);

  protected:
    std::unordered_map<Eigen::Vector3i, StateVertex *, matrix_hash<Eigen::Vector3i>> DataPtrTable;

    double space_resolution;
    double inv_space_resolution;
    double time_resolution;
    double inv_time_resolution;
    double x_offest, y_offest, t_offest;
};

class V_A_Traj
{
  public:
    std::vector<double> time_table;
    double traj_start_time_sec;

    V_A_Traj(ros::NodeHandle &nh_) : nh(nh_)
    {
    }

    V_A_Traj(ros::NodeHandle &nh_, std::string vis_name) : nh(nh_)
    {
        trajectoryPub = nh.advertise<visualization_msgs::MarkerArray>(vis_name, 10);
    }

    void get(const std::vector<StateVertex> &state_path_points,
             const Eigen::Matrix<double, N_coeff * N_poly, SpaceDim> _C, const double &_C_durr);

    Eigen::Vector2d getPosition(double t);

    bool vis_path();

  private:
    struct State
    {
        Eigen::Vector2d P;
        Eigen::Vector2d V;
        Eigen::Vector2d at_1;
        double dT;
    };

    ros::NodeHandle nh;

    std::vector<State> path_points;
    Eigen::Matrix<double, N_coeff * N_poly, SpaceDim> C;
    double C_dur;

    ros::Publisher trajectoryPub;

    Eigen::Vector4d Beta(const unsigned int d, const double t);
};

inline bool StateVertexPtrGrid_hash::find(
    const Eigen::Vector2d &coord, const double &time,
    std::unordered_map<Eigen::Vector3i, StateVertex *, matrix_hash<Eigen::Vector3i>>::iterator &iter)
{
    Eigen::Vector3i idx;

    idx(0) = (coord(0) - x_offest) * inv_space_resolution;
    idx(1) = (coord(1) - y_offest) * inv_space_resolution;
    idx(2) = (time - t_offest) * inv_time_resolution;

    if ((iter = DataPtrTable.find(idx)) == DataPtrTable.cend())
        return false;
    else
        return true;
}

inline void StateVertexPtrGrid_hash::insert(const Eigen::Vector2d &coord, const double &time, StateVertex *ptr)
{
    Eigen::Vector3i idx;

    idx(0) = (coord(0) - x_offest) * inv_space_resolution;
    idx(1) = (coord(1) - y_offest) * inv_space_resolution;
    idx(2) = (time - t_offest) * inv_time_resolution;

    DataPtrTable.insert(std::make_pair(idx, ptr));
}

inline bool StateVertexPtrGrid_hash::initGrid(const double &_space_resolution, const double &_time_resolution,
                                              const Eigen::Vector3d &global_xyt_l)
{
    space_resolution = _space_resolution;
    time_resolution = _time_resolution;
    inv_space_resolution = 1.0 / space_resolution;
    inv_time_resolution = 1.0 / time_resolution;

    if (TimeDim == 1)
    {
        inv_time_resolution = 0.0;
    }

    x_offest = global_xyt_l(0);
    y_offest = global_xyt_l(1);
    t_offest = global_xyt_l(2);

    return true;
}

inline bool StateVertex::simForward(StateVertex &prev, const Eigen::Matrix<double, SpaceDim, 1> &_at_1,
                                    const double &_dT, Eigen::Matrix<double, SpaceDim, 1> Goal_P,
                                    Eigen::Matrix<double, SpaceDim, 1> Goal_V)
{
    dT = _dT;
    at_1 = _at_1;
    prev_ptr = &prev;
    time_stamp = prev.time_stamp + dT;

    V = prev.V + at_1 * dT;

    // Limit the velocity
    // We only need to consider two cases which is V be out of the limit (bigger or lower).
    for (int i = 0; i < V.size(); ++i)
    {
        if (V(i) > V_max(i))
        {
            if (at_1(i) > 0.0) // V(i) == prev.V(i) + at_1(i) * dT > prev.V(i)
            {
                if (V_max(i) >= prev.V(i)) // V(i) > V_max(i) >= prev.V(i)
                {
                    V(i) = prev.V(i);
                    at_1(i) = 0.0;
                }
                else // V(i) > prev.V(i) > V_max(i)
                {
                    return false;
                }
            }
            else if (at_1(i) == 0.0)
            {
                return false;
            }
        }
        else if (V(i) < -V_max(i))
        {
            if (at_1(i) < 0.0) // V(i) < prev.V(i)
            {
                if (-V_max(i) <= prev.V(i)) // V(i) < -V_max(i) <= prev.V(i)
                {
                    V(i) = prev.V(i);
                    at_1(i) = 0.0;
                }
                else // V(i) < prev.V(i) < -V_max(i)
                {
                    return false;
                }
            }
            else if (at_1(i) == 0.0)
            {
                return false;
            }
        }
    }

    P = prev.P + prev.V * dT + 0.5 * at_1 * dT * dT;

    cost2come = prev.cost2come + (at_1.norm() + time_weight) * dT; // Estimated energy + time wasted

    Eigen::VectorXd StartState(2 * SpaceDim);
    StartState.head(SpaceDim) = P;
    StartState.tail(SpaceDim) = V;
    Eigen::VectorXd GoalState(2 * SpaceDim);
    GoalState.head(SpaceDim) = Goal_P;
    GoalState.tail(SpaceDim) = Goal_V;

    cost2go = estimateHeuristic(StartState, GoalState, time2go);

    f = cost2come + cost2go;

    return true;
}

inline bool StateVertex::checkCollision(std::shared_ptr<const voxel_map_tool> map_ptr, const double map_resolution)
{
    // const double estimated_dist = prev_ptr->V.norm() * dT + 0.5 * at_1.norm() * dT * dT;

    if (prev_ptr == nullptr)
    {
        if (map_ptr->isObsFree(P) == false)
        {
            isCollision = true;
            return false;
        }

        isCollision = false;
        return true;
    }

    const double estimated_dist = std::max(prev_ptr->V.norm(), V.norm()) * dT + 5;

    int check_steps = (int)(estimated_dist / map_resolution);

    if (check_steps < 1)
        check_steps = 1;

    const double IncreT = dT / check_steps;
    for (double t = IncreT; t < dT + 1e-10; t += IncreT)
    {
        Eigen::Vector2d current_position = (prev_ptr->P + prev_ptr->V * t + 0.5 * at_1 * t * t);

        if (map_ptr->isObsFree(current_position) == false)
        {
            isCollision = true;
            return false;
        }
    }

    return true;
}

inline bool StateVertex::setState(const Eigen::Matrix<double, SpaceDim, 1> &Curr_P,
                                  const Eigen::Matrix<double, SpaceDim, 1> &Curr_V,
                                  std::shared_ptr<const voxel_map_tool> map_ptr, const double &_time_stamp,
                                  Eigen::Matrix<double, SpaceDim, 1> Goal_P, Eigen::Matrix<double, SpaceDim, 1> Goal_V)
{
    at_1 = Eigen::Matrix<double, SpaceDim, 1>::Zero();
    prev_ptr = nullptr;
    time_stamp = _time_stamp;

    P = Curr_P;
    V = Curr_V;

    if (map_ptr->isObsFree(P) == false)
    {
        isCollision = true;
        return false;
    }

    cost2come = 0.0; // T + ax^2 + ay^2 + az^2

    Eigen::VectorXd StartState(2 * SpaceDim);
    StartState.head(SpaceDim) = P;
    StartState.tail(SpaceDim) = V;
    Eigen::VectorXd GoalState(2 * SpaceDim);
    GoalState.head(SpaceDim) = Goal_P;
    GoalState.tail(SpaceDim) = Goal_V;

    cost2go = estimateHeuristic(StartState, GoalState, time2go);

    f = cost2come + cost2go;

    dT = 0.0;

    return true;
}

inline void StateVertex::getPathFromHere(std::vector<StateVertex> &path_points)
{
    path_points.clear();
    path_points.push_back(*this);

    auto v_ptr = this->prev_ptr;

    while (v_ptr != nullptr)
    {
        path_points.push_back(*v_ptr);
        v_ptr = v_ptr->prev_ptr;
    }

    std::reverse(path_points.begin(), path_points.end());

    path_points.begin()->prev_ptr = nullptr;
    for (auto ptr = path_points.begin() + 1; ptr < path_points.cend(); ++ptr)
    {
        ptr->prev_ptr = &ptr[-1];
    }
}

inline double StateVertex::estimateHeuristic(const Eigen::VectorXd &StartState, const Eigen::VectorXd &GoalState,
                                             double &optimal_time)
{
    const Eigen::Vector3d dp = GoalState.head(SpaceDim) - StartState.head(SpaceDim);
    const Eigen::Vector3d v0 = StartState.tail(SpaceDim);
    const Eigen::Vector3d v1 = GoalState.tail(SpaceDim);

    double c1 = -36 * dp.dot(dp);
    double c2 = 24 * (v0 + v1).dot(dp);
    double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));
    double c4 = 0;
    double c5 = time_weight;

    std::vector<double> ts = quartic(c5, c4, c3, c2, c1);

    double v_max = V_max.maxCoeff() * 0.5;
    double t_bar = (StartState.head(SpaceDim) - GoalState.head(SpaceDim)).lpNorm<Eigen::Infinity>() / v_max;
    ts.push_back(t_bar);

    double cost = 100000000;
    double t_d = t_bar;

    for (auto t : ts)
    {
        if (t < t_bar)
            continue;
        double c = -c1 / (3 * t * t * t) - c2 / (2 * t * t) - c3 / t + time_weight * t;
        if (c < cost)
        {
            cost = c;
            t_d = t;
        }
    }

    optimal_time = t_d;

    // return 1.0 * (1 + tie_breaker) * cost;
    return 0.0;
}

inline std::vector<double> StateVertex::cubic(double a, double b, double c, double d)
{
    std::vector<double> dts;

    double a2 = b / a;
    double a1 = c / a;
    double a0 = d / a;

    double Q = (3 * a1 - a2 * a2) / 9;
    double R = (9 * a1 * a2 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
    double D = Q * Q * Q + R * R;
    if (D > 0)
    {
        double S = std::cbrt(R + sqrt(D));
        double T = std::cbrt(R - sqrt(D));
        dts.push_back(-a2 / 3 + (S + T));
        return dts;
    }
    else if (D == 0)
    {
        double S = std::cbrt(R);
        dts.push_back(-a2 / 3 + S + S);
        dts.push_back(-a2 / 3 - S);
        return dts;
    }
    else
    {
        double theta = acos(R / sqrt(-Q * Q * Q));
        dts.push_back(2 * sqrt(-Q) * cos(theta / 3) - a2 / 3);
        dts.push_back(2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) - a2 / 3);
        dts.push_back(2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) - a2 / 3);
        return dts;
    }
}

inline std::vector<double> StateVertex::quartic(double a, double b, double c, double d, double e)
{
    std::vector<double> dts;

    double a3 = b / a;
    double a2 = c / a;
    double a1 = d / a;
    double a0 = e / a;

    std::vector<double> ys = cubic(1, -a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);
    double y1 = ys.front();
    double r = a3 * a3 / 4 - a2 + y1;
    if (r < 0)
        return dts;

    double R = sqrt(r);
    double D, E;
    if (R != 0)
    {
        D = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 + 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
        E = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 - 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
    }
    else
    {
        D = sqrt(0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0));
        E = sqrt(0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0));
    }

    if (!std::isnan(D))
    {
        dts.push_back(-a3 / 4 + R / 2 + D / 2);
        dts.push_back(-a3 / 4 + R / 2 - D / 2);
    }
    if (!std::isnan(E))
    {
        dts.push_back(-a3 / 4 - R / 2 + E / 2);
        dts.push_back(-a3 / 4 - R / 2 - E / 2);
    }

    return dts;
}

inline void V_A_Traj::get(const std::vector<StateVertex> &state_path_points,
                          const Eigen::Matrix<double, N_coeff * N_poly, SpaceDim> _C, const double &_C_durr)
{
    path_points.clear();
    path_points.resize(state_path_points.size());
    for (unsigned long i = 0; i < path_points.size(); ++i)
    {
        path_points.at(i).at_1.x() = state_path_points.at(i).at_1.x();
        path_points.at(i).at_1.y() = state_path_points.at(i).at_1.y();

        path_points.at(i).P.x() = state_path_points.at(i).P.x();
        path_points.at(i).P.y() = state_path_points.at(i).P.y();

        path_points.at(i).V.x() = state_path_points.at(i).V.x();
        path_points.at(i).V.y() = state_path_points.at(i).V.y();

        path_points.at(i).dT = state_path_points.at(i).dT;
    }

    C = _C;

    traj_start_time_sec = state_path_points.at(0).time_stamp;

    time_table.clear();
    time_table.push_back(0.0);
    for (unsigned int i = 1; i < path_points.size(); ++i)
    {
        time_table.push_back(time_table.back() + path_points.at(i).dT);
    }

    C_dur = _C_durr;
}

inline Eigen::Vector2d V_A_Traj::getPosition(double t)
{
    Eigen::Vector2d current_position;
    if (t < time_table.front() + 1e-10)
    {
        t = time_table.front() + 1e-10;
    }

    if (t > time_table.back() + C_dur - 1e-10)
    {
        t = time_table.back() + C_dur - 1e-10;
    }

    if (t < time_table.back())
    {
        unsigned long i = 0;
        while (i < time_table.size() - 1)
        {
            if (t > time_table.at(i) - 1e-10 && t < time_table.at(i + 1) + 1e-10)
            {
                break;
            }
            i += 1;
        }
        const double curr_sec_t = t - time_table.at(i);

        return path_points.at(i).P + path_points.at(i).V * curr_sec_t +
               0.5 * path_points.at(i + 1).at_1 * curr_sec_t * curr_sec_t;
    }
    else
    {
        double curr_sec_t = t - time_table.back();
        current_position(0) = C.block<N_coeff, 1>(0, 0).dot(Beta(0, curr_sec_t));
        current_position(1) = C.block<N_coeff, 1>(0, 1).dot(Beta(0, curr_sec_t));

        return current_position;
    }
}

inline bool V_A_Traj::vis_path()
{
    visualization_msgs::MarkerArray trajMarker;
    visualization_msgs::Marker trajMarkerSect;

    trajMarkerSect.header.stamp = ros::Time::now();
    trajMarkerSect.type = visualization_msgs::Marker::LINE_LIST;
    trajMarkerSect.header.frame_id = "map";
    trajMarkerSect.pose.orientation.w = 1.00;
    trajMarkerSect.action = visualization_msgs::Marker::ADD;
    trajMarkerSect.color.r = 0.00;
    trajMarkerSect.color.g = 0.00;
    trajMarkerSect.color.b = 1.00;
    trajMarkerSect.color.a = 1.00;
    trajMarkerSect.scale.x = 0.2;

    trajMarkerSect.id = 0;
    trajMarkerSect.ns = "trajectory";

    Eigen::Vector2d lastX = getPosition(0);

    double T = 0.05;
    double t = T;
    for (; t < time_table.back(); t += T)
    {
        Eigen::Vector2d current_position = getPosition(t);

        geometry_msgs::Point point;

        point.x = lastX(0);
        point.y = lastX(1);
        trajMarkerSect.points.push_back(point);
        point.x = current_position(0);
        point.y = current_position(1);
        trajMarkerSect.points.push_back(point);
        lastX = current_position;
    }

    trajMarker.markers.push_back(trajMarkerSect);

    trajMarkerSect.header.stamp = ros::Time::now();
    trajMarkerSect.type = visualization_msgs::Marker::LINE_LIST;
    trajMarkerSect.header.frame_id = "map";
    trajMarkerSect.pose.orientation.w = 1.00;
    trajMarkerSect.action = visualization_msgs::Marker::ADD;
    trajMarkerSect.color.r = 0.00;
    trajMarkerSect.color.g = 0.00;
    trajMarkerSect.color.b = 1.00;
    trajMarkerSect.color.a = 1.00;
    trajMarkerSect.scale.x = 0.2;

    trajMarkerSect.id = 1;
    trajMarkerSect.ns = "trajectory";

    trajMarkerSect.points.clear();

    for (; t < time_table.back() + C_dur; t += T)
    {
        Eigen::Vector2d current_position = getPosition(t);

        geometry_msgs::Point point;

        point.x = lastX(0);
        point.y = lastX(1);
        trajMarkerSect.points.push_back(point);
        point.x = current_position(0);
        point.y = current_position(1);
        trajMarkerSect.points.push_back(point);
        lastX = current_position;
    }

    trajMarker.markers.push_back(trajMarkerSect);
    trajectoryPub.publish(trajMarker);

    return true;
}

inline Eigen::Vector4d V_A_Traj::Beta(const unsigned int d, const double t)
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
} // namespace HybridAStar

#endif