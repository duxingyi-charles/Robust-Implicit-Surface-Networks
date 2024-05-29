//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_IO_H
#define ROBUST_IMPLICIT_NETWORKS_IO_H

#include <string>
#include <nlohmann/json.hpp>
#include "mesh.h"

struct Config {
    // Input
    std::string tet_mesh_file;
    std::string func_file;
    std::string output_dir;

    // Lookup setting.
    bool use_lookup;
    bool use_secondary_lookup;
    bool use_topo_ray_shooting;

    // Parameter for tet grid generation.
    // (Only used if tet_mesh_file is empty.)
    size_t tet_mesh_resolution;
    std::array<double, 3> tet_mesh_bbox_min;
    std::array<double, 3> tet_mesh_bbox_max;
};

Config parse_config_file(const std::string &filename);

bool load_tet_mesh(const std::string &filename,
                   std::vector<std::array<double, 3>> &pts,
                   std::vector<std::array<size_t, 4>> &tets);

bool generate_tet_mesh(size_t resolution,
                   const std::array<double, 3>& bbox_min,
                   const std::array<double, 3>& bbox_max,
                   std::vector<std::array<double, 3>> &pts,
                   std::vector<std::array<size_t, 4>> &tets);

bool save_result(const std::string& filename,
                 const std::vector<std::array<double, 3>>& mesh_pts,
                 const std::vector<PolygonFace>& mesh_faces,
                 const std::vector<std::vector<size_t>>& patches,
                 const std::vector<size_t>& patch_function_label,
                 const std::vector<Edge>& edges,
                 const std::vector<std::vector<size_t>>& chains,
                 const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                 const std::vector<std::vector<size_t>>& shells,
                 const std::vector<std::vector<size_t>>& cells,
                 const std::vector<std::vector<bool>>& cell_function_label);

bool save_result_MI(const std::string& filename,
                 const std::vector<std::array<double, 3>>& mesh_pts,
                 const std::vector<PolygonFace>& mesh_faces,
                 const std::vector<std::vector<size_t>>& patches,
                 const std::vector<std::pair<size_t, size_t>>& patch_function_label,
                 const std::vector<Edge>& edges,
                 const std::vector<std::vector<size_t>>& chains,
                 const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                 const std::vector<std::vector<size_t>>& shells,
                 const std::vector<std::vector<size_t>>& cells,
                 const std::vector<size_t>& cell_function_label);

bool save_result_msh(const std::string& filename,
                     const std::vector<std::array<double, 3>>& mesh_pts,
                     const std::vector<PolygonFace>& mesh_faces,
                     const std::vector<std::vector<size_t>>& patches,
                     const std::vector<Edge>& edges,
                     const std::vector<std::vector<size_t>>& chains,
                     const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                     const std::vector<std::vector<size_t>>& shells,
                     const std::vector<std::vector<size_t>>& cells);

bool save_timings(const std::string& filename,
                  const std::vector<std::string> &timing_labels,
                  const std::vector<double> &timings);

bool save_statistics(const std::string& filename,
                     const std::vector<std::string> &stats_labels,
                     const std::vector<size_t> &stats);

#endif //ROBUST_IMPLICIT_NETWORKS_IO_H


