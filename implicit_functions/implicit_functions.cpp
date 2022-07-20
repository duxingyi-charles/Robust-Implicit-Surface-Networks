#include "implicit_functions.h"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

bool load_functions(const std::string& filename,
                    const std::vector<std::array<double, 3>>& pts,
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& funcVals)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "function file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    size_t n_pts = pts.size();
    size_t n_func = data.size();
    funcVals.resize(n_pts, n_func);
    for (int j = 0; j < n_func; ++j) {
        std::string type = data[j]["type"].get<std::string>();
        if (type == "plane") {
            std::array<double,3> point;
            for (int i = 0; i < 3; ++i) {
                point[i] = data[j]["point"][i].get<double>();
            }
            std::array<double,3> normal;
            for (int i = 0; i < 3; ++i) {
                normal[i] = data[j]["normal"][i].get<double>();
            }
            //
            for (int i = 0; i < n_pts; i++) {
                funcVals(i,j) = compute_plane_distance(point, normal, pts[i]);
            }
        }
        else if (type == "line") {
            std::array<double,3> point;
            for (int i = 0; i < 3; ++i) {
                point[i] = data[j]["point"][i].get<double>();
            }
            std::array<double,3> unit_vector;
            for (int i = 0; i < 3; ++i) {
                unit_vector[i] = data[j]["unit_vector"][i].get<double>();
            }
            //
            for (int i = 0; i < n_pts; i++) {
                funcVals(i,j) = compute_line_distance(point, unit_vector, pts[i]);
            }
        }
        else if (type == "cylinder") {
            std::array<double,3> axis_point;
            for (int i = 0; i < 3; ++i) {
                axis_point[i] = data[j]["axis_point"][i].get<double>();
            }
            std::array<double,3> axis_unit_vector;
            for (int i = 0; i < 3; ++i) {
                axis_unit_vector[i] = data[j]["axis_vector"][i].get<double>();
            }
            double radius = data[j]["radius"].get<double>();
            //
            for (int i = 0; i < n_pts; i++) {
                funcVals(i,j) = compute_cylinder_distance(axis_point,
                                                          axis_unit_vector, radius,
                                                          pts[i]);
            }
        }
        else if (type == "sphere") {
            std::array<double,3> center;
            for (int i = 0; i < 3; ++i) {
                center[i] = data[j]["center"][i].get<double>();
            }
            double radius = data[j]["radius"].get<double>();
            //
            for (int i = 0; i < n_pts; i++) {
                funcVals(i,j) = compute_sphere_distance(center,
                                                        radius,pts[i]);
            }
        }
        else if (type == "torus") {
            std::array<double,3> center;
            for (int i = 0; i < 3; ++i) {
                center[i] = data[j]["center"][i].get<double>();
            }
            std::array<double,3> axis_unit_vector;
            for (int i = 0; i < 3; ++i) {
                axis_unit_vector[i] = data[j]["axis_vector"][i].get<double>();
            }
            double major_radius = data[j]["major_radius"].get<double>();
            double minor_radius = data[j]["minor_radius"].get<double>();
            //
            for (int i = 0; i < n_pts; i++) {
                funcVals(i,j) = compute_torus_distance(center, axis_unit_vector,
                                                       major_radius, minor_radius, pts[i]);
            }
        }
        else if (type == "circle") {
            std::array<double,3> center;
            for (int i = 0; i < 3; ++i) {
                center[i] = data[j]["center"][i].get<double>();
            }
            std::array<double,3> axis_unit_vector;
            for (int i = 0; i < 3; ++i) {
                axis_unit_vector[i] = data[j]["axis_vector"][i].get<double>();
            }
            double radius = data[j]["radius"].get<double>();
            //
            for (int i = 0; i < n_pts; i++) {
                funcVals(i,j) = compute_circle_distance(center, axis_unit_vector,
                                                        radius, pts[i]);
            }
        }
        else if (type == "cone") {
            std::array<double, 3> apex;
            for (int i = 0; i < 3; ++i) {
                apex[i] = data[j]["apex"][i].get<double>();
            }
            std::array<double, 3> axis_unit_vector;
            for (int i = 0; i < 3; ++i) {
                axis_unit_vector[i] = data[j]["axis_vector"][i].get<double>();
            }
            double apex_angle = data[j]["apex_angle"].get<double>();
            //
            for (int i = 0; i < n_pts; i++) {
                funcVals(i, j) = compute_cone_distance(apex, axis_unit_vector, apex_angle, pts[i]);
            }
        }
        else if (type == "zero") {
            for (int i = 0; i < n_pts; i++) {
                funcVals(i, j) = 0;
            }
        }
        else {
            std::cout << "undefined type: " << type << std::endl;
            return false;
        }
    }
    return true;
}