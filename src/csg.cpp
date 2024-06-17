//
//  csg.cpp
//  robust_implicit_networks
//
//  Created by Yiwen Ju on 6/3/24.
//

#include "csg.h"
#include <queue>

bool csg(bool robust_test,
         bool use_lookup,
         bool use_secondary_lookup,
         bool use_topo_ray_shooting,
         bool positive_inside,
         //
         const std::vector<std::array<double, 3>>& pts,
         const std::vector<std::array<size_t, 4>>& tets,
         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& funcVals,
         const std::function<bool(std::vector<bool>)>& lambda,
         //
         std::vector<std::array<double, 3>>& mesh_pts,
         std::vector<PolygonFace>& mesh_faces,
         std::vector<std::vector<size_t>>& patches,
         std::vector<size_t>& patch_function_label,
         std::vector<bool>& patch_sign_label,
         std::vector<Edge>& edges,
         std::vector<std::vector<size_t>>& chains,
         std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
         std::vector<std::vector<size_t>>& shells,
         std::vector<std::vector<size_t>>& cells,
         std::vector<std::vector<bool>>& cell_function_label,
         std::vector<std::string>& timing_labels,
         std::vector<double>& timings,
         std::vector<std::string>& stats_labels,
         std::vector<size_t>& stats){
    if (!positive_inside) funcVals = funcVals * -1;
    if (!implicit_arrangement(
                              robust_test,
                              use_lookup,
                              use_secondary_lookup,
                              use_topo_ray_shooting,
                              //
                              pts, tets, funcVals,
                              //
                              mesh_pts,mesh_faces,patches, patch_function_label,
                              edges,chains,
                              non_manifold_edges_of_vert,
                              shells,cells,cell_function_label,
                              timing_labels,timings,
                              stats_labels,stats)) {
                                  return -1;
                              }
    if (robust_test) return 0;
    std::vector<bool> cells_label(cells.size(), false);
    for (size_t i = 0; i < cells_label.size(); i++){
        cells_label[i] = lambda(cell_function_label[i]);
    }
    std::vector<std::vector<size_t>> pruned_patches;
    std::vector<std::vector<size_t>> pruned_chains;
    pruned_patches.reserve(patches.size());
    pruned_chains.reserve(chains.size());
    //for each patch, find two cells from each side
    std::vector<std::pair<size_t, size_t>> patch_to_cell(patches.size(), {Mesh_None, Mesh_None});
    for (size_t cell_index = 0; cell_index < cells.size(); cell_index ++){
        for (auto shell_index : cells[cell_index]){
            for (auto half_patch : shells[shell_index]){
                size_t patch_id = half_patch / 2;
                if (half_patch % 2){
                    patch_to_cell[patch_id].second = cell_index;
                }else{
                    patch_to_cell[patch_id].first = cell_index;
                }
            }
        }
    }
    std::vector<size_t> vertex_to_patch_count(mesh_pts.size(), 0);
    for (size_t patch_index = 0; patch_index < patches.size(); patch_index++){
        if (cells_label[patch_to_cell[patch_index].first] != cells_label[patch_to_cell[patch_index].second]){
            pruned_patches.emplace_back(patches[patch_index]);
            if (cells_label[patch_to_cell[patch_index].first]){
                patch_sign_label.emplace_back(true);
            }else{
                patch_sign_label.emplace_back(false);
            }
            std::vector<size_t> vertices_on_patch;
            vertices_on_patch.reserve(3 * patches[patch_index].size());
            for (auto face : patches[patch_index]){
                for (auto vertex : mesh_faces[face].vert_indices){
                    vertices_on_patch.emplace_back(vertex);
                }
            }
            for (auto& vertex_index : vertices_on_patch) {
                vertex_to_patch_count[vertex_index] ++;
            }
        }
    }
    //prune chains
    for (auto chain : chains){
        bool active = true;
        for (auto edge_index : chain){
            if (vertex_to_patch_count[edges[edge_index].v1] < 2||vertex_to_patch_count[edges[edge_index].v2] < 2){
                active = false;
                break;
            }
        }
        if (active){
            pruned_chains.emplace_back(chain);
        }
    }
    //prune corners
    for (size_t corner_index = 0; corner_index < non_manifold_edges_of_vert.size(); corner_index ++ ){
        if (vertex_to_patch_count[corner_index] <= 2){
            non_manifold_edges_of_vert[corner_index] = {};
        }
    }
    patches = pruned_patches;
    chains = pruned_chains;
    
    return true;
}
