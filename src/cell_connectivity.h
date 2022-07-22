//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_CELL_CONNECTIVITY_H
#define ROBUST_IMPLICIT_NETWORKS_CELL_CONNECTIVITY_H

#include <vector>
#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/material_interface.h>

// build simplicial cell adjacency graph (implicit arrangement)
void build_simplicial_cell_adjacency(
        const std::vector<std::array<size_t, 4>> &tets,
        const std::vector<simplicial_arrangement::Arrangement<3>> &cut_results,
        const std::vector<size_t> &cut_result_index,
        const std::vector<long long> &global_vId_of_tet_vert,
        const std::vector<size_t> &global_vId_start_index_of_tet,
        const std::vector<size_t> &iso_fId_of_tet_face,
        const std::vector<size_t> &iso_fId_start_index_of_tet,
        const std::vector<size_t> &patch_of_face,
        const std::vector<size_t> &shell_of_half_patch,
        std::vector<std::pair<size_t, size_t>> &tet_cell_of_simp_cell,
        std::vector<long long> &simp_half_face_info,
        std::vector<size_t> &simp_hFace_start_index
);

// build simplicial cell adjacency graph (material interface)
void build_simplicial_cell_adjacency(
        const std::vector<std::array<size_t, 4>> &tets,
        const std::vector<simplicial_arrangement::MaterialInterface<3>> &cut_results,
        const std::vector<size_t> &cut_result_index,
        const std::vector<long long> &global_vId_of_tet_vert,
        const std::vector<size_t> &global_vId_start_index_of_tet,
        const std::vector<size_t> &MI_fId_of_tet_face,
        const std::vector<size_t> &MI_fId_start_index_of_tet,
        const std::vector<size_t> &patch_of_face,
        const std::vector<size_t> &shell_of_half_patch,
        std::vector<std::pair<size_t, size_t>> &tet_cell_of_simp_cell,
        std::vector<long long> &simp_half_face_info,
        std::vector<size_t> &simp_hFace_start_index
);

void compute_simplicial_cell_connected_components(
        const std::vector<std::pair<size_t, size_t>> &tet_cell_of_simp_cell,
        const std::vector<long long> &simp_half_face_info,
        const std::vector<size_t> &simp_hFace_start_index,
        std::vector<std::vector<size_t>> &arrangement_cells);

#endif //ROBUST_IMPLICIT_NETWORKS_CELL_CONNECTIVITY_H
