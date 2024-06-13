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
         //
         const std::vector<std::array<double, 3>>& pts,
         const std::vector<std::array<size_t, 4>>& tets,
         const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& funcVals,
         const std::function<bool(std::vector<bool>)>& lambda,
         //
         std::vector<std::array<double, 3>>& mesh_pts,
         std::vector<PolygonFace>& mesh_faces,
         std::vector<std::vector<size_t>>& patches,
         std::vector<size_t>& patch_function_label,
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
    std::vector<std::vector<size_t>> pruned_shells;
    pruned_shells.reserve(shells.size());
    pruned_patches.reserve(patches.size());
    pruned_chains.reserve(chains.size());
    std::unordered_map<size_t, size_t> patch_to_patch;
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
    std::vector<std::vector<bool>> merged_cell(cells.size(), std::vector<bool>(cells.size(), false));
    for (size_t patch_index = 0; patch_index < patches.size(); patch_index++){
        if (cells_label[patch_to_cell[patch_index].first] != cells_label[patch_to_cell[patch_index].second]){
            pruned_patches.emplace_back(patches[patch_index]);
            patch_to_patch[patch_index] = pruned_patches.size() - 1;
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
        }else{
            merged_cell[patch_to_cell[patch_index].first][patch_to_cell[patch_index].second] = true;
            merged_cell[patch_to_cell[patch_index].second][patch_to_cell[patch_index].first] = true;
        }
    }
    //merge cells together
    std::vector<std::vector<bool>> merged_cell_copy = merged_cell;
    std::vector<std::vector<size_t>> csg_cells;
    for (size_t i = 0; i < merged_cell.size(); i++){
        bool isolate = true;
        for (size_t j = 0; j < merged_cell[i].size(); j++){
            if (merged_cell_copy[i][j]){
                isolate = false;
            }
            if (merged_cell[i][j]){
                std::vector<size_t> visited_cells(cells.size(), false);
                std::vector<size_t> connected_cells = {};
                visited_cells[i] = true;
                visited_cells[j] = true;
                std::queue<size_t> Q;
                Q.push(i);
                Q.push(j);
                merged_cell[i][j] = false;
                merged_cell[j][i] = false;
                while (!Q.empty()){
                    size_t other_cell = Q.front();
                    connected_cells.emplace_back(other_cell);
                    Q.pop();
                    for (size_t cell_index = 0; cell_index < merged_cell[other_cell].size(); cell_index++){
                        if (merged_cell[other_cell][cell_index]){
                            merged_cell[other_cell][cell_index] = false;
                            merged_cell[cell_index][other_cell] = false;
                            if (!visited_cells[cell_index]){
                                visited_cells[cell_index] = true;
                                Q.push(cell_index);
                            }
                        }
                    }
                }
                csg_cells.emplace_back(connected_cells);
            }
        }
        if (isolate)
            csg_cells.emplace_back(std::vector<size_t>({i}));
    }
    
    //collect one-sided patches to form shells from merged csg cells
    std::vector<std::vector<size_t>> patch_to_shell(patches.size(), std::vector<size_t>({}));
    for (size_t csg_cell_index = 0; csg_cell_index < csg_cells.size(); csg_cell_index ++){
        auto csg_cell = csg_cells[csg_cell_index];
        std::vector<size_t> shell_candidate;
        std::vector<bool> patch_active(patches.size(), false);
        for (auto cell_elements : csg_cell){
            for (auto shell_index : cells[cell_elements]){
                for (auto patch_index : shells[shell_index]){
                    patch_active[patch_index / 2] = !patch_active[patch_index / 2];
                }
            }
        }
        for (size_t patch_index = 0; patch_index < patch_active.size(); patch_index ++){
            if (patch_active[patch_index]){
                shell_candidate.emplace_back(patch_index);
                patch_to_shell[patch_index].emplace_back(csg_cell_index);
            }
        }
        pruned_shells.emplace_back(shell_candidate);
    }
    //build shell connectivity
    std::vector<bool> valid_shell_index (pruned_shells.size(), false);
    for (auto patch_map : patch_to_shell){
        size_t min_length = Mesh_None;
        size_t min_index = Mesh_None;
        for (auto shell_index : patch_map){
            if (pruned_shells[shell_index].size() < min_length){
                min_length = pruned_shells[shell_index].size();
                min_index = shell_index;
            }
        }
        if (min_index != Mesh_None){
            valid_shell_index[min_index] = true;
        }
    }
    shells = {};
    for (size_t csg_shell_index = 0; csg_shell_index < valid_shell_index.size(); csg_shell_index ++){
        if (valid_shell_index[csg_shell_index]){
            std::vector<size_t> new_shell;
            new_shell.reserve(pruned_shells[csg_shell_index].size());
            for (auto patch_index : pruned_shells[csg_shell_index]){
                new_shell.emplace_back(patch_to_patch[patch_index]);
            }
            shells.emplace_back(new_shell);
        }
    }
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
    for (size_t corner_index = 0; corner_index < non_manifold_edges_of_vert.size(); corner_index ++ ){
        if (vertex_to_patch_count[corner_index] <= 2){
            non_manifold_edges_of_vert[corner_index] = {};
        }
    }
    patches = pruned_patches;
    chains = pruned_chains;
    
    return true;
}
