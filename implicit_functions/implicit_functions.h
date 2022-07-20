#pragma once

#include <string>
#include <array>
#include <Eigen/Core>

bool load_functions(const std::string& filename,
                    const std::vector<std::array<double,3>> &pts,
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcVals);

inline double compute_Euclidean_distance(const std::array<double,3> &p, const std::array<double,3>& q)
{
    return sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]));
}

inline double compute_squared_distance(const std::array<double,3> &p, const std::array<double,3>& q)
{
    return (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]);
}

inline double compute_norm(const std::array<double,3> &v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline void normalize_vector(std::array<double,3> &v)
{
    double norm = compute_norm(v);
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= norm;
}

inline double compute_signed_sphere_distance(const std::array<double, 3>& center, double r, const std::array<double, 3>& p)
{
    return r - compute_Euclidean_distance(center, p);
}

inline double compute_unsigned_sphere_distance(const std::array<double, 3>& center, double r, const std::array<double, 3>& p)
{
    return -abs(r - compute_Euclidean_distance(center, p));
}

inline double compute_sphere_distance(const std::array<double, 3>& center, double r, const std::array<double, 3>& p)
{
    if (r >= 0) {
        return compute_signed_sphere_distance(center, r, p);
    } else {
        return compute_unsigned_sphere_distance(center, -r, p);
    }
}



inline double compute_dot(const std::array<double,3> &a, const std::array<double,3> &b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double compute_cone_distance(const std::array<double,3>& apex,
                                    const std::array<double,3>& axis_unit_vector, double apex_angle,
                                    const std::array<double,3>& p)
{
    return compute_dot(axis_unit_vector, {p[0]-apex[0], p[1]-apex[1], p[2]-apex[2]})
           - cos(apex_angle) * compute_Euclidean_distance(p,apex);
}

inline double compute_cylinder_distance(const std::array<double,3>& axis_point,
                                        const std::array<double,3>& axis_unit_vector, double radius,
                                        const std::array<double,3>& p)
{
    std::array<double,3> vec {p[0]-axis_point[0], p[1]-axis_point[1], p[2]-axis_point[2]};
    double d = compute_dot(axis_unit_vector, vec);
    vec[0] = vec[0] - d * axis_unit_vector[0];
    vec[1] = vec[1] - d * axis_unit_vector[1];
    vec[2] = vec[2] - d * axis_unit_vector[2];
    return radius - compute_norm(vec);
}

inline double compute_plane_distance(const std::array<double,3> &point,
                                     const std::array<double,3> &normal,
                                     const std::array<double,3> &p)
{
    return compute_dot(normal, {p[0]-point[0], p[1]-point[1], p[2]-point[2]});
}

inline double compute_torus_distance(const std::array<double,3>& center,
                                     const std::array<double,3>& axis_unit_vector,
                                     double major_radius, double minor_radius,
                                     const std::array<double,3>& p)
{
    std::array<double,3> vec {p[0] - center[0], p[1]-center[1], p[2]-center[2]};
    double d = compute_dot(vec, axis_unit_vector);
    std::array<double,3> vec_para = {d * axis_unit_vector[0],
                                     d*axis_unit_vector[1], d*axis_unit_vector[2]};
    std::array<double,3> vec_perp = {vec[0]-vec_para[0], vec[1]-vec_para[1], vec[2]-vec_para[2]};
    double vec_perp_norm = compute_norm(vec_perp);
    if (vec_perp_norm == 0) { // point p lies on torus axis
        return minor_radius - sqrt(compute_dot(vec_para, vec_para) + major_radius * major_radius);
    } else {
        return minor_radius - sqrt(major_radius * (major_radius - 2 * vec_perp_norm)
                                   + compute_dot(vec, vec));
    }
}

//
inline double compute_line_distance(const std::array<double,3>& point, const std::array<double,3>& unit_vector,
                                    const std::array<double,3> &p)
{
    std::array<double,3> vec {p[0] - point[0], p[1] - point[1], p[2] - point[2]};
    double d = compute_dot(unit_vector, vec);
    std::array<double,3> vec_para {d * unit_vector[0], d * unit_vector[1], d * unit_vector[2]};
    std::array<double,3> vec_perp {vec[0] - vec_para[0], vec[1] - vec_para[1], vec[2] - vec_para[2]};
    return -compute_norm(vec_perp);
}

inline double compute_circle_distance(const std::array<double,3>& center,
                                      const std::array<double,3> &axis_unit_vector, double radius,
                                      const std::array<double,3> &p)
{
    std::array<double,3> vec {p[0] - center[0], p[1]-center[1], p[2]-center[2]};
    double d = compute_dot(vec, axis_unit_vector);
    std::array<double,3> vec_para = {d * axis_unit_vector[0],
                                     d*axis_unit_vector[1], d*axis_unit_vector[2]};
    std::array<double,3> vec_perp = {vec[0]-vec_para[0], vec[1]-vec_para[1], vec[2]-vec_para[2]};
    double vec_perp_norm = compute_norm(vec_perp);
    if (vec_perp_norm == 0) { // point p lies on torus axis
        return -sqrt(compute_dot(vec_para, vec_para) + radius * radius);
    } else {
        return -sqrt(radius * (radius - 2 * vec_perp_norm) + compute_dot(vec, vec));
    }
}

inline int sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}