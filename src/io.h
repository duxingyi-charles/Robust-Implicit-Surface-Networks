//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_IO_H
#define ROBUST_IMPLICIT_NETWORKS_IO_H

#include <string>
#include <nlohmann/json.hpp>
#include "mesh.h"

// parse config file for Implicit Arrangement
bool parse_config_file(const std::string &filename,
                       std::string& tet_mesh_file,
                       std::string& func_file,
                       std::string& output_dir,
                       bool& use_lookup,
                       bool& use_2func_lookup,
                       bool& use_topo_ray_shooting);

// parse config file for Material Interface
bool parse_config_file_MI(const std::string& filename,
                          std::string& tet_mesh_file,
                          std::string& material_file,
                          std::string& output_dir,
                          bool& use_lookup,
                          bool& use_3func_lookup,
                          bool& use_topo_ray_shooting);

bool load_tet_mesh(const std::string &filename,
                   std::vector<std::array<double, 3>> &pts,
                   std::vector<std::array<size_t, 4>> &tets);

bool save_result(const std::string& filename,
                 const std::vector<std::array<double, 3>>& mesh_pts,
                 const std::vector<PolygonFace>& mesh_faces,
                 const std::vector<std::vector<size_t>>& patches,
                 const std::vector<Edge>& edges,
                 const std::vector<std::vector<size_t>>& chains,
                 const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                 const std::vector<std::vector<std::pair<std::pair<size_t, int>,std::pair<size_t, int>>>>& half_patch_pair_list,
                 const std::vector<std::vector<size_t>>& shells,
                 const std::vector<std::vector<size_t>>& components,
                 const std::vector<std::vector<size_t>>& cells);

bool save_result_msh(const std::string& filename,
                     const std::vector<std::array<double, 3>>& mesh_pts,
                     const std::vector<PolygonFace>& mesh_faces,
                     const std::vector<std::vector<size_t>>& patches,
                     const std::vector<Edge>& edges,
                     const std::vector<std::vector<size_t>>& chains,
                     const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                     const std::vector<std::vector<size_t>>& shells,
                     const std::vector<std::vector<size_t>>& components,
                     const std::vector<std::vector<size_t>>& cells);

bool save_timings(const std::string& filename,
                  const std::vector<std::string> &timing_labels,
                  const std::vector<double> &timings);

bool save_statistics(const std::string& filename,
                     const std::vector<std::string> &stats_labels,
                     const std::vector<size_t> &stats);

#endif //ROBUST_IMPLICIT_NETWORKS_IO_H


