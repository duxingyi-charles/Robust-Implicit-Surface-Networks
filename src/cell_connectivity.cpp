//
// Created by Charles Du on 7/20/22.
//
#include "cell_connectivity.h"
#include "mesh.h"
#include "extract_mesh.h"

#include <queue>
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>




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
)
{
//    timing_labels.emplace_back("arrCells(build simpCell graph)");
//    ScopedTimer<> timer("arrangement cells: build simplicial cells graph");
    size_t num_simplicial_cells = 0;
    size_t num_simp_half_faces = 0;
    for (size_t i = 0; i < tets.size(); i++) {
        if (cut_result_index[i] != Mesh_None) {
            const auto& cells = cut_results[cut_result_index[i]].cells;
            num_simplicial_cells += cells.size();
            for (const auto& cell : cells) {
                num_simp_half_faces += cell.faces.size();
            }
        } else { // empty tet
            ++num_simplicial_cells;
            num_simp_half_faces += 4;
        }
    }
    tet_cell_of_simp_cell.reserve(num_simplicial_cells);
    simp_half_face_info.reserve(num_simp_half_faces);
    simp_hFace_start_index.reserve(num_simplicial_cells + 1);
    simp_hFace_start_index.push_back(0);
    // hash table: face key --> (simplicial cell index, face index in simplicial
    // cell)
    absl::flat_hash_map<std::array<long long, 3>, std::pair<size_t, size_t>>
            incident_cell_of_face;
    //                    incident_cell_of_face.reserve();
    std::vector<long long> face_verts;
    face_verts.reserve(4);
    std::array<long long, 3> key;
    std::vector<std::array<long long, 3>> four_face_verts(4);
    for (size_t i = 0; i < tets.size(); i++) {
        if (cut_result_index[i] != Mesh_None) {
            size_t iso_fId_start_index = iso_fId_start_index_of_tet[i];
            size_t global_vId_start_index = global_vId_start_index_of_tet[i];
            const auto& arrangement = cut_results[cut_result_index[i]];
            const auto& faces = arrangement.faces;
            const auto& cells = arrangement.cells;
            for (size_t j = 0; j < cells.size(); j++) {
                // create a new simplicial cell
                const auto& cell = cells[j];
                size_t cur_simp_cell_id = tet_cell_of_simp_cell.size();
                tet_cell_of_simp_cell.emplace_back(i, j);
                for (size_t k = 0; k < cell.faces.size(); k++) {
                    size_t fId = cell.faces[k];
                    size_t iso_fId = iso_fId_of_tet_face[iso_fId_start_index + fId];
                    if (iso_fId !=
                        Mesh_None) { // face k is on iso-surface
                        size_t patch_id = patch_of_face[iso_fId];
                        size_t half_patch_id = (faces[fId].positive_cell == j)
                                               ? (2 * patch_id)
                                               : (2 * patch_id + 1);
                        size_t shell_id = shell_of_half_patch[half_patch_id];
                        simp_half_face_info.push_back(-shell_id - 1);
                    } else { // face k is not on iso-surface
                        face_verts.clear();
                        const auto& faceVertices = faces[fId].vertices;
                        // convert from local vert index to global vert index
                        for (unsigned long vId : faceVertices) {
                            face_verts.push_back(
                                    global_vId_of_tet_vert[global_vId_start_index + vId]);
                        }
                        //
                        compute_iso_face_key(face_verts, key);
                        auto iter_inserted = incident_cell_of_face.try_emplace(
                                key, std::make_pair(cur_simp_cell_id, k));
                        if (!iter_inserted
                                .second) { // same face has been inserted before
                            size_t opposite_simp_cell_id =
                                    iter_inserted.first->second.first;
                            size_t opposite_cell_face_id =
                                    iter_inserted.first->second.second;
                            // make neighbor
                            simp_half_face_info.push_back(opposite_simp_cell_id);
                            simp_half_face_info
                            [simp_hFace_start_index[opposite_simp_cell_id] +
                             opposite_cell_face_id] = cur_simp_cell_id;
                            // delete face in hash table since a face can only be
                            // shared by two cells
                            incident_cell_of_face.erase(iter_inserted.first);
                        } else { // the face is inserted for the first time
                            simp_half_face_info.push_back(
                                    std::numeric_limits<long long>::max());
                        }
                    }
                }
                simp_hFace_start_index.push_back(simp_half_face_info.size());
            }
        } else { // tet i has no isosurface
            // create a new simplicial cell
            size_t cur_simp_cell_id = tet_cell_of_simp_cell.size();
            tet_cell_of_simp_cell.emplace_back(i, 0);
            // global index of four tet vertices
            long long global_vId0 = -tets[i][0] - 1;
            long long global_vId1 = -tets[i][1] - 1;
            long long global_vId2 = -tets[i][2] - 1;
            long long global_vId3 = -tets[i][3] - 1;
            // face 0
            four_face_verts[0] = {global_vId1, global_vId2, global_vId3};
            std::sort(four_face_verts[0].begin(), four_face_verts[0].end());
            // face 1
            four_face_verts[1] = {global_vId2, global_vId3, global_vId0};
            std::sort(four_face_verts[1].begin(), four_face_verts[1].end());
            // face 2
            four_face_verts[2] = {global_vId3, global_vId0, global_vId1};
            std::sort(four_face_verts[2].begin(), four_face_verts[2].end());
            // face 3
            four_face_verts[3] = {global_vId0, global_vId1, global_vId2};
            std::sort(four_face_verts[3].begin(), four_face_verts[3].end());
            //
            for (size_t j = 0; j < 4; j++) {
                auto iter_inserted = incident_cell_of_face.try_emplace(
                        four_face_verts[j], std::make_pair(cur_simp_cell_id, j));
                if (!iter_inserted.second) { // same face has been inserted before
                    size_t opposite_simp_cell_id =
                            iter_inserted.first->second.first;
                    size_t opposite_cell_face_id =
                            iter_inserted.first->second.second;
                    // make neighbor
                    simp_half_face_info.push_back(opposite_simp_cell_id);
                    simp_half_face_info
                    [simp_hFace_start_index[opposite_simp_cell_id] +
                     opposite_cell_face_id] = cur_simp_cell_id;
                    // delete face in hash table since a face can only be shared by
                    // two cells
                    incident_cell_of_face.erase(iter_inserted.first);
                } else { // face inserted for the first time
                    simp_half_face_info.push_back(
                            std::numeric_limits<long long>::max());
                }
            }
            simp_hFace_start_index.push_back(simp_half_face_info.size());
        }
    }
//    timings.push_back(timer.toc());
}

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
)
{
    size_t num_simplicial_cells = 0;
    size_t num_simp_half_faces = 0;
    for (size_t i = 0; i < tets.size(); i++) {
        if (cut_result_index[i] != Mesh_None) {
            const auto& cells = cut_results[cut_result_index[i]].cells;
            num_simplicial_cells += cells.size();
            for (const auto& cell : cells) {
                num_simp_half_faces += cell.faces.size();
            }
        } else { // empty tet
            ++num_simplicial_cells;
            num_simp_half_faces += 4;
        }
    }
    tet_cell_of_simp_cell.reserve(num_simplicial_cells);
    simp_half_face_info.reserve(num_simp_half_faces);
    simp_hFace_start_index.reserve(num_simplicial_cells + 1);
    simp_hFace_start_index.push_back(0);
    // hash table: face key --> (simplicial cell index, face index in simplicial
    // cell)
    absl::flat_hash_map<std::array<long long, 3>, std::pair<size_t, size_t>>
            incident_cell_of_face;
    //                    incident_cell_of_face.reserve();
    std::vector<long long> face_verts;
    face_verts.reserve(4);
    std::array<long long, 3> key;
    std::vector<std::array<long long, 3>> four_face_verts(4);
    for (size_t i = 0; i < tets.size(); i++) {
        if (cut_result_index[i] != Mesh_None) {
            size_t MI_fId_start_index = MI_fId_start_index_of_tet[i];
            size_t global_vId_start_index = global_vId_start_index_of_tet[i];
            const auto& matInterface = cut_results[cut_result_index[i]];
            const auto& faces = matInterface.faces;
            const auto& cells = matInterface.cells;
            for (size_t j = 0; j < cells.size(); j++) {
                // create a new simplicial cell
                const auto& cell = cells[j];
                size_t cur_simp_cell_id = tet_cell_of_simp_cell.size();
                tet_cell_of_simp_cell.emplace_back(i, j);
                for (size_t k = 0; k < cell.faces.size(); k++) {
                    size_t fId = cell.faces[k];
                    size_t MI_fId = MI_fId_of_tet_face[MI_fId_start_index + fId];
                    if (MI_fId !=
                        Mesh_None) { // face k is on iso-surface
                        size_t patch_id = patch_of_face[MI_fId];
                        size_t half_patch_id = (faces[fId].positive_material_label == cell.material_label)
                                               ? (2 * patch_id)
                                               : (2 * patch_id + 1);
                        size_t shell_id = shell_of_half_patch[half_patch_id];
                        simp_half_face_info.push_back(-shell_id - 1);
                    } else { // face k is not on material interface
                        face_verts.clear();
                        const auto& faceVertices = faces[fId].vertices;
                        // convert from local vert index to global vert index
                        for (unsigned long vId : faceVertices) {
                            face_verts.push_back(
                                    global_vId_of_tet_vert[global_vId_start_index + vId]);
                        }
                        //
                        compute_iso_face_key(face_verts, key);
                        auto iter_inserted = incident_cell_of_face.try_emplace(
                                key, std::make_pair(cur_simp_cell_id, k));
                        if (!iter_inserted
                                .second) { // same face has been inserted before
                            size_t opposite_simp_cell_id =
                                    iter_inserted.first->second.first;
                            size_t opposite_cell_face_id =
                                    iter_inserted.first->second.second;
                            // make neighbor
                            simp_half_face_info.push_back(opposite_simp_cell_id);
                            simp_half_face_info
                            [simp_hFace_start_index[opposite_simp_cell_id] +
                             opposite_cell_face_id] = cur_simp_cell_id;
                            // delete face in hash table since a face can only be
                            // shared by two cells
                            incident_cell_of_face.erase(iter_inserted.first);
                        } else { // the face is inserted for the first time
                            simp_half_face_info.push_back(
                                    std::numeric_limits<long long>::max());
                        }
                    }
                }
                simp_hFace_start_index.push_back(simp_half_face_info.size());
            }
        } else { // tet i has no material interface
            // create a new simplicial cell
            size_t cur_simp_cell_id = tet_cell_of_simp_cell.size();
            tet_cell_of_simp_cell.emplace_back(i, 0);
            // global index of four tet vertices
            long long global_vId0 = -tets[i][0] - 1;
            long long global_vId1 = -tets[i][1] - 1;
            long long global_vId2 = -tets[i][2] - 1;
            long long global_vId3 = -tets[i][3] - 1;
            // face 0
            four_face_verts[0] = {global_vId1, global_vId2, global_vId3};
            std::sort(four_face_verts[0].begin(), four_face_verts[0].end());
            // face 1
            four_face_verts[1] = {global_vId2, global_vId3, global_vId0};
            std::sort(four_face_verts[1].begin(), four_face_verts[1].end());
            // face 2
            four_face_verts[2] = {global_vId3, global_vId0, global_vId1};
            std::sort(four_face_verts[2].begin(), four_face_verts[2].end());
            // face 3
            four_face_verts[3] = {global_vId0, global_vId1, global_vId2};
            std::sort(four_face_verts[3].begin(), four_face_verts[3].end());
            //
            for (size_t j = 0; j < 4; j++) {
                auto iter_inserted = incident_cell_of_face.try_emplace(
                        four_face_verts[j], std::make_pair(cur_simp_cell_id, j));
                if (!iter_inserted.second) { // same face has been inserted before
                    size_t opposite_simp_cell_id =
                            iter_inserted.first->second.first;
                    size_t opposite_cell_face_id =
                            iter_inserted.first->second.second;
                    // make neighbor
                    simp_half_face_info.push_back(opposite_simp_cell_id);
                    simp_half_face_info
                    [simp_hFace_start_index[opposite_simp_cell_id] +
                     opposite_cell_face_id] = cur_simp_cell_id;
                    // delete face in hash table since a face can only be shared by
                    // two cells
                    incident_cell_of_face.erase(iter_inserted.first);
                } else { // face inserted for the first time
                    simp_half_face_info.push_back(
                            std::numeric_limits<long long>::max());
                }
            }
            simp_hFace_start_index.push_back(simp_half_face_info.size());
        }
    }
}

void compute_simplicial_cell_connected_components(
        const std::vector<std::pair<size_t, size_t>> &tet_cell_of_simp_cell,
        const std::vector<long long> &simp_half_face_info,
        const std::vector<size_t> &simp_hFace_start_index,
        std::vector<std::vector<size_t>> &arrangement_cells)
{
    // --------------- group simplicial cells into arrangement cells
    // ---------------
    size_t num_simp_cells = tet_cell_of_simp_cell.size();
    std::vector<bool> visited_simp_cell(num_simp_cells, false);
    std::vector<absl::flat_hash_set<size_t>> arrangement_cell_incident_shells;
    for (size_t i = 0; i < num_simp_cells; i++) {
        if (!visited_simp_cell[i]) {
            // new arrangement cell
            arrangement_cell_incident_shells.emplace_back();
            auto& incident_shells = arrangement_cell_incident_shells.back();
            //
            std::queue<size_t> Q;
            Q.push(i);
            visited_simp_cell[i] = true;
            while (!Q.empty()) {
                size_t simp_cell_id = Q.front();
                Q.pop();
                size_t start_index = simp_hFace_start_index[simp_cell_id];
                size_t face_info_size =
                        simp_hFace_start_index[simp_cell_id + 1] - start_index;
                for (size_t j = 0; j < face_info_size; j++) {
                    if (simp_half_face_info[start_index + j] < 0) {
                        // the half-face is on isosurface
                        size_t shell_id = -simp_half_face_info[start_index + j] - 1;
                        incident_shells.insert(shell_id);
                    } else {
                        // the half-face is not on isosurface
                        long long opposite_simp_cell_id =
                                simp_half_face_info[start_index + j];
                        if ((opposite_simp_cell_id !=
                             std::numeric_limits<long long>::max()) &&
                            !visited_simp_cell[opposite_simp_cell_id]) {
                            // the opposite simplicial cell exists and is not
                            // visited
                            Q.push(opposite_simp_cell_id);
                            visited_simp_cell[opposite_simp_cell_id] = true;
                        }
                    }
                }
            }
        }
    }
    // convert vector-of-set to vector-of-vector
    arrangement_cells.reserve(arrangement_cell_incident_shells.size());
    for (auto & arrangement_cell_incident_shell : arrangement_cell_incident_shells) {
        arrangement_cells.emplace_back(arrangement_cell_incident_shell.begin(),
                                       arrangement_cell_incident_shell.end());
    }
}
