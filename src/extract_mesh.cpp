//
// Created by Charles Du on 7/20/22.
//

#include "extract_mesh.h"
#include <absl/container/flat_hash_map.h>

using namespace simplicial_arrangement;

void extract_iso_mesh(size_t num_1_func,
                           size_t num_2_func,
                           size_t num_more_func,
                           const std::vector<Arrangement<3>>& cut_results,
                           const std::vector<size_t>& cut_result_index,
                           const std::vector<size_t>& func_in_tet,
                           const std::vector<size_t>& start_index_of_tet,
                           const std::vector<std::array<size_t, 4>>& tets,
                           std::vector<IsoVert>& iso_verts,
                           std::vector<PolygonFace>& iso_faces)
{
    size_t n_tets = tets.size();
    // estimate number of iso-verts and iso-faces
    size_t max_num_face = num_1_func + 4 * num_2_func + 8 * num_more_func;
    size_t max_num_vert = max_num_face;
    iso_verts.reserve(max_num_vert);
    iso_faces.reserve(max_num_face);
    // hash table for vertices on the boundary of tetrahedron
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    absl::flat_hash_map<std::array<size_t, 3>, size_t> vert_on_tetEdge;
    vert_on_tetEdge.reserve(num_1_func + 3 * num_2_func + num_more_func);
    absl::flat_hash_map<std::array<size_t, 5>, size_t> vert_on_tetFace;
    vert_on_tetFace.reserve(num_2_func + 6 * num_more_func);
    // hash table for faces on the boundary of tetrahedron
    absl::flat_hash_map<std::array<size_t, 3>, size_t> face_on_tetFace;
    //
    std::vector<bool> is_iso_vert;
    is_iso_vert.reserve(8);
    std::vector<bool> is_iso_face;
    is_iso_face.reserve(9);
    std::vector<size_t> iso_vId_of_vert;
    iso_vId_of_vert.reserve(8);
    std::vector<size_t> face_verts;
    face_verts.reserve(4);
    std::array<size_t, 3> key3;
    std::array<size_t, 5> key5;
    std::array<bool, 4> used_pId;
    std::array<size_t, 2> vIds2;
    std::array<size_t, 3> vIds3;
    std::array<size_t, 3> implicit_pIds;
    std::array<size_t, 3> bndry_pIds;
    //
    for (size_t i = 0; i < n_tets; i++) {
        if (cut_result_index[i] != Arrangement<3>::None) {
            const auto& arrangement = cut_results[cut_result_index[i]];
            const auto& vertices = arrangement.vertices;
            const auto& faces = arrangement.faces;
            auto start_index = start_index_of_tet[i];
            auto num_func = start_index_of_tet[i + 1] - start_index;
            // find vertices and faces on isosurface
            is_iso_vert.clear();
            for (int j = 0; j < vertices.size(); ++j) {
                is_iso_vert.push_back(false);
            }
            is_iso_face.clear();
            if (arrangement.unique_planes.empty()) { // all planes are unique
                for (const auto& face : faces) {
                    is_iso_face.push_back(false);
                    if (face.supporting_plane > 3) { // plane 0,1,2,3 are tet boundaries
                        is_iso_face.back() = true;
                        for (const auto& vid : face.vertices) {
                            is_iso_vert[vid] = true;
                        }
                    }
                }
            } else {
                for (const auto& face : faces) {
                    is_iso_face.push_back(false);
                    auto pid = face.supporting_plane;
                    auto uid = arrangement.unique_plane_indices[pid];
                    for (const auto& plane_id : arrangement.unique_planes[uid]) {
                        if (plane_id > 3) { // plane 0,1,2,3 are tet boundaries
                            //                        is_iso_face[j] = true;
                            is_iso_face.back() = true;
                            for (const auto& vid : face.vertices) {
                                is_iso_vert[vid] = true;
                            }
                            break;
                        }
                    }
                }
            }
            // map: local vert index --> iso-vert index
            iso_vId_of_vert.clear();
            // create iso-vertices
            for (size_t j = 0; j < vertices.size(); j++) {
                iso_vId_of_vert.push_back(Arrangement<3>::None);
                if (is_iso_vert[j]) {
                    size_t num_bndry_planes = 0;
                    size_t num_impl_planes = 0;
                    const auto& vertex = vertices[j];
                    // vertex.size() == 3
                    for (size_t k = 0; k < 3; k++) {
                        if (vertex[k] > 3) { // plane 0,1,2,3 are tet boundaries
                            implicit_pIds[num_impl_planes] =
                                    func_in_tet[vertex[k] - 4 + start_index];
                            ++num_impl_planes;
                        } else {
                            bndry_pIds[num_bndry_planes] = vertex[k];
                            ++num_bndry_planes;
                        }
                    }
                    switch (num_bndry_planes) {
                        case 2: // on tet edge
                        {
                            used_pId[0] = false;
                            used_pId[1] = false;
                            used_pId[2] = false;
                            used_pId[3] = false;
                            used_pId[bndry_pIds[0]] = true;
                            used_pId[bndry_pIds[1]] = true;
                            //                        std::array<size_t, 2> vIds;
                            size_t num_vIds = 0;
                            for (size_t k = 0; k < 4; k++) {
                                if (!used_pId[k]) {
                                    vIds2[num_vIds] = tets[i][k];
                                    ++num_vIds;
                                }
                            }
                            size_t vId1 = vIds2[0];
                            size_t vId2 = vIds2[1];
                            if (vId1 > vId2) {
                                size_t tmp = vId1;
                                vId1 = vId2;
                                vId2 = tmp;
                            }
                            key3[0] = vId1;
                            key3[1] = vId2;
                            key3[2] = implicit_pIds[0];
                            auto iter_inserted = vert_on_tetEdge.try_emplace(key3, iso_verts.size());
                            if (iter_inserted.second) {
                                iso_verts.emplace_back();
                                auto& iso_vert = iso_verts.back();
                                iso_vert.tet_index = i;
                                iso_vert.tet_vert_index = j;
                                iso_vert.simplex_size = 2;
                                iso_vert.simplex_vert_indices[0] = vId1;
                                iso_vert.simplex_vert_indices[1] = vId2;
                                iso_vert.func_indices[0] = implicit_pIds[0];
                            }
                            iso_vId_of_vert.back() = iter_inserted.first->second;
                            break;
                        }
                        case 1: // on tet face
                        {
                            size_t pId = bndry_pIds[0];
                            size_t num_vIds = 0;
                            for (size_t k = 0; k < 4; k++) {
                                if (k != pId) {
                                    vIds3[num_vIds] = tets[i][k];
                                    ++num_vIds;
                                }
                            }
                            std::sort(vIds3.begin(), vIds3.end());
                            key5[0] = vIds3[0];
                            key5[1] = vIds3[1];
                            key5[2] = vIds3[2];
                            key5[3] = implicit_pIds[0];
                            key5[4] = implicit_pIds[1];
                            auto iter_inserted = vert_on_tetFace.try_emplace(key5, iso_verts.size());
                            if (iter_inserted.second) {
                                iso_verts.emplace_back();
                                auto& iso_vert = iso_verts.back();
                                iso_vert.tet_index = i;
                                iso_vert.tet_vert_index = j;
                                iso_vert.simplex_size = 3;
                                iso_vert.simplex_vert_indices[0] = vIds3[0];
                                iso_vert.simplex_vert_indices[1] = vIds3[1];
                                iso_vert.simplex_vert_indices[2] = vIds3[2];
                                iso_vert.func_indices[0] = implicit_pIds[0];
                                iso_vert.func_indices[1] = implicit_pIds[1];
                            }
                            iso_vId_of_vert.back() = iter_inserted.first->second;
                            break;
                        }
                        case 0: // in tet cell
                        {
                            iso_vId_of_vert.back() = iso_verts.size();
                            iso_verts.emplace_back();
                            auto& iso_vert = iso_verts.back();
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 4;
                            iso_vert.simplex_vert_indices = tets[i];
                            iso_vert.func_indices = implicit_pIds;
                            break;
                        }
                        case 3: // on tet vertex
                        {
                            used_pId[0] = false;
                            used_pId[1] = false;
                            used_pId[2] = false;
                            used_pId[3] = false;
                            for (const auto& pId : bndry_pIds) {
                                used_pId[pId] = true;
                            }
                            size_t vId;
                            for (size_t k = 0; k < 4; k++) {
                                if (!used_pId[k]) {
                                    vId = k;
                                    break;
                                }
                            }
                            auto key = tets[i][vId];
                            auto iter_inserted = vert_on_tetVert.try_emplace(key, iso_verts.size());
                            if (iter_inserted.second) {
                                iso_verts.emplace_back();
                                auto& iso_vert = iso_verts.back();
                                iso_vert.tet_index = i;
                                iso_vert.tet_vert_index = j;
                                iso_vert.simplex_size = 1;
                                iso_vert.simplex_vert_indices[0] = tets[i][vId];
                            }
                            //                        iso_vId_of_vert[j] = iter_inserted.first->second;
                            iso_vId_of_vert.back() = iter_inserted.first->second;
                            break;
                        }
                        default: break;
                    }
                }
            }
            // create iso-faces
            for (size_t j = 0; j < faces.size(); j++) {
                if (is_iso_face[j]) {
                    face_verts.clear();
                    for (unsigned long vId : faces[j].vertices) {
                        face_verts.push_back(iso_vId_of_vert[vId]);
                    }
                    //
                    // face is on tet boundary if face.negative_cell is NONE
                    bool face_on_tet_boundary = (faces[j].negative_cell == Arrangement<3>::None);
                    //
                    if (face_on_tet_boundary) {
                        compute_iso_face_key(face_verts, key3);
                        auto iter_inserted = face_on_tetFace.try_emplace(key3, iso_faces.size());
                        if (iter_inserted.second) {
                            iso_faces.emplace_back();
                            iso_faces.back().vert_indices = face_verts;
                            iso_faces.back().tet_face_indices.emplace_back(i, j);
                            iso_faces.back().func_index = func_in_tet[faces[j].supporting_plane - 4 + start_index];
                        } else { // iso_face inserted before
                            size_t iso_face_id = (iter_inserted.first)->second;
                            iso_faces[iso_face_id].tet_face_indices.emplace_back(i, j);
                        }
                    } else { // face not on tet boeundary
                        iso_faces.emplace_back();
                        iso_faces.back().vert_indices = face_verts;
                        iso_faces.back().tet_face_indices.emplace_back(i, j);
                        iso_faces.back().func_index = func_in_tet[faces[j].supporting_plane - 4 + start_index];
                    }
                }
            }
        }
    }
    //
}


void extract_iso_mesh(size_t num_1_func,
                      size_t num_2_func,
                      size_t num_more_func,
                      const std::vector<Arrangement<3>>& cut_results,
                      const std::vector<size_t>& cut_result_index,
                      const std::vector<size_t>& func_in_tet,
                      const std::vector<size_t>& start_index_of_tet,
                      const std::vector<std::array<size_t, 4>>& tets,
                      std::vector<IsoVert>& iso_verts,
                      std::vector<PolygonFace>& iso_faces,
                      std::vector<long long>& global_vId_of_tet_vert,
                      std::vector<size_t>& global_vId_start_index_of_tet,
                      std::vector<size_t>& iso_fId_of_tet_face,
                      std::vector<size_t>& iso_fId_start_index_of_tet)
{
    size_t n_tets = tets.size();
    // get total number of verts and faces
    size_t total_num_vert = 0;
    size_t total_num_face = 0;
    for (const auto& arrangement : cut_results) {
        total_num_vert += arrangement.vertices.size();
        total_num_face += arrangement.faces.size();
    }
    global_vId_of_tet_vert.reserve(total_num_vert);
    global_vId_start_index_of_tet.reserve(n_tets + 1);
    global_vId_start_index_of_tet.push_back(0);
    iso_fId_of_tet_face.reserve(total_num_face);
    iso_fId_start_index_of_tet.reserve(n_tets + 1);
    iso_fId_start_index_of_tet.push_back(0);
    // estimate number of iso-verts and iso-faces
    size_t max_num_face = num_1_func + 4 * num_2_func + 8 * num_more_func;
    size_t max_num_vert = max_num_face;
    iso_verts.reserve(max_num_vert);
    iso_faces.reserve(max_num_face);
    // hash table for vertices on the boundary of tetrahedron
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    absl::flat_hash_map<std::array<size_t, 3>, size_t> vert_on_tetEdge;
    vert_on_tetEdge.reserve(num_1_func + 3 * num_2_func + num_more_func);
    absl::flat_hash_map<std::array<size_t, 5>, size_t> vert_on_tetFace;
    vert_on_tetFace.reserve(num_2_func + 6 * num_more_func);
    // hash table for faces on the boundary of tetrahedron
    absl::flat_hash_map<std::array<size_t, 3>, size_t> face_on_tetFace;
    //
    std::vector<bool> is_iso_vert;
    is_iso_vert.reserve(8);
    std::vector<bool> is_iso_face;
    is_iso_face.reserve(9);
    std::vector<size_t> iso_vId_of_vert;
    iso_vId_of_vert.reserve(8);
    std::vector<size_t> face_verts;
    face_verts.reserve(4);
    std::array<size_t, 3> key3;
    std::array<size_t, 5> key5;
    std::array<bool, 4> used_pId;
    std::array<size_t, 2> vIds2;
    std::array<size_t, 3> vIds3;
    std::array<size_t, 3> implicit_pIds;
    std::array<size_t, 3> bndry_pIds;
    //
    for (size_t i = 0; i < n_tets; i++) {
        if (cut_result_index[i] == Arrangement<3>::None) {
            global_vId_start_index_of_tet.push_back(global_vId_of_tet_vert.size());
            iso_fId_start_index_of_tet.push_back(iso_fId_of_tet_face.size());
        } else {
            const auto& arrangement = cut_results[cut_result_index[i]];
            const auto& vertices = arrangement.vertices;
            const auto& faces = arrangement.faces;
            auto start_index = start_index_of_tet[i];
            auto num_func = start_index_of_tet[i + 1] - start_index;
            // find vertices and faces on isosurface
            is_iso_vert.clear();
            for (int j = 0; j < vertices.size(); ++j) {
                is_iso_vert.push_back(false);
            }
            is_iso_face.clear();
            if (arrangement.unique_planes.empty()) { // all planes are unique
                for (const auto& face : faces) {
                    is_iso_face.push_back(false);
                    if (face.supporting_plane > 3) { // plane 0,1,2,3 are tet boundaries
                        is_iso_face.back() = true;
                        for (const auto& vid : face.vertices) {
                            is_iso_vert[vid] = true;
                        }
                    }
                }
            } else {
                for (const auto& face : faces) {
                    is_iso_face.push_back(false);
                    auto uid = arrangement.unique_plane_indices[face.supporting_plane];
                    for (const auto& plane_id : arrangement.unique_planes[uid]) {
                        if (plane_id > 3) { // plane 0,1,2,3 are tet boundaries
                            //                        is_iso_face[j] = true;
                            is_iso_face.back() = true;
                            for (const auto& vid : face.vertices) {
                                is_iso_vert[vid] = true;
                            }
                            break;
                        }
                    }
                }
            }
            // map: local vert index --> iso-vert index
            iso_vId_of_vert.clear();
            // create iso-vertices
            for (size_t j = 0; j < vertices.size(); j++) {
                size_t num_bndry_planes = 0;
                size_t num_impl_planes = 0;
                const auto& vertex = vertices[j];
                // vertex.size() == 3
                for (size_t k = 0; k < 3; k++) {
                    if (vertex[k] > 3) { // plane 0,1,2,3 are tet boundaries
                        implicit_pIds[num_impl_planes] = func_in_tet[vertex[k] - 4 + start_index];
                        ++num_impl_planes;
                    } else {
                        bndry_pIds[num_bndry_planes] = vertex[k];
                        ++num_bndry_planes;
                    }
                }
                //
                if (!is_iso_vert[j]) { // vert[j] must be tet vertex
                    used_pId[0] = false;
                    used_pId[1] = false;
                    used_pId[2] = false;
                    used_pId[3] = false;
                    for (const auto& pId : bndry_pIds) {
                        used_pId[pId] = true;
                    }
                    size_t vId;
                    for (size_t k = 0; k < 4; k++) {
                        if (!used_pId[k]) {
                            vId = k;
                            break;
                        }
                    }
                    global_vId_of_tet_vert.push_back(-tets[i][vId] - 1);
                    iso_vId_of_vert.push_back(Arrangement<3>::None);
                } else { // iso-vertex
                    switch (num_bndry_planes) {
                        case 2: // on tet edge
                        {
                            used_pId[0] = false;
                            used_pId[1] = false;
                            used_pId[2] = false;
                            used_pId[3] = false;
                            used_pId[bndry_pIds[0]] = true;
                            used_pId[bndry_pIds[1]] = true;
                            size_t num_vIds = 0;
                            for (size_t k = 0; k < 4; k++) {
                                if (!used_pId[k]) {
                                    vIds2[num_vIds] = tets[i][k];
                                    ++num_vIds;
                                }
                            }
                            size_t vId1 = vIds2[0];
                            size_t vId2 = vIds2[1];
                            if (vId1 > vId2) {
                                size_t tmp = vId1;
                                vId1 = vId2;
                                vId2 = tmp;
                            }
                            key3[0] = vId1;
                            key3[1] = vId2;
                            key3[2] = implicit_pIds[0];
                            auto iter_inserted = vert_on_tetEdge.try_emplace(key3, iso_verts.size());
                            if (iter_inserted.second) {
                                iso_verts.emplace_back();
                                auto& iso_vert = iso_verts.back();
                                iso_vert.tet_index = i;
                                iso_vert.tet_vert_index = j;
                                iso_vert.simplex_size = 2;
                                iso_vert.simplex_vert_indices[0] = vId1;
                                iso_vert.simplex_vert_indices[1] = vId2;
                                iso_vert.func_indices[0] = implicit_pIds[0];
                            }
                            global_vId_of_tet_vert.push_back(iter_inserted.first->second);
                            iso_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        case 1: // on tet face
                        {
                            size_t pId = bndry_pIds[0];
                            size_t num_vIds = 0;
                            for (size_t k = 0; k < 4; k++) {
                                if (k != pId) {
                                    vIds3[num_vIds] = tets[i][k];
                                    ++num_vIds;
                                }
                            }
                            std::sort(vIds3.begin(), vIds3.end());
                            key5[0] = vIds3[0];
                            key5[1] = vIds3[1];
                            key5[2] = vIds3[2];
                            key5[3] = implicit_pIds[0];
                            key5[4] = implicit_pIds[1];
                            auto iter_inserted = vert_on_tetFace.try_emplace(key5, iso_verts.size());
                            if (iter_inserted.second) {
                                iso_verts.emplace_back();
                                auto& iso_vert = iso_verts.back();
                                iso_vert.tet_index = i;
                                iso_vert.tet_vert_index = j;
                                iso_vert.simplex_size = 3;
                                iso_vert.simplex_vert_indices[0] = vIds3[0];
                                iso_vert.simplex_vert_indices[1] = vIds3[1];
                                iso_vert.simplex_vert_indices[2] = vIds3[2];
                                iso_vert.func_indices[0] = implicit_pIds[0];
                                iso_vert.func_indices[1] = implicit_pIds[1];
                            }
                            global_vId_of_tet_vert.push_back(iter_inserted.first->second);
                            iso_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        case 0: // in tet cell
                        {
                            global_vId_of_tet_vert.push_back(iso_verts.size());
                            iso_vId_of_vert.push_back(iso_verts.size());
                            iso_verts.emplace_back();
                            auto& iso_vert = iso_verts.back();
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 4;
                            iso_vert.simplex_vert_indices = tets[i];
                            iso_vert.func_indices = implicit_pIds;
                            break;
                        }
                        case 3: // on tet vertex
                        {
                            used_pId[0] = false;
                            used_pId[1] = false;
                            used_pId[2] = false;
                            used_pId[3] = false;
                            for (const auto& pId : bndry_pIds) {
                                used_pId[pId] = true;
                            }
                            size_t vId;
                            for (size_t k = 0; k < 4; k++) {
                                if (!used_pId[k]) {
                                    vId = k;
                                    break;
                                }
                            }
                            auto key = tets[i][vId];
                            auto iter_inserted = vert_on_tetVert.try_emplace(key, iso_verts.size());
                            if (iter_inserted.second) {
                                iso_verts.emplace_back();
                                auto& iso_vert = iso_verts.back();
                                iso_vert.tet_index = i;
                                iso_vert.tet_vert_index = j;
                                iso_vert.simplex_size = 1;
                                iso_vert.simplex_vert_indices[0] = tets[i][vId];
                            }
                            global_vId_of_tet_vert.push_back(-key - 1);
                            iso_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        default: break;
                    }
                }
            }
            global_vId_start_index_of_tet.push_back(global_vId_of_tet_vert.size());
            // create iso-faces
            for (size_t j = 0; j < faces.size(); j++) {
                iso_fId_of_tet_face.push_back(Arrangement<3>::None);
                if (is_iso_face[j]) {
                    face_verts.clear();
                    for (unsigned long vId : faces[j].vertices) {
                        face_verts.push_back(iso_vId_of_vert[vId]);
                    }
                    //
                    // face is on tet boundary if face.negative_cell is NONE
                    bool face_on_tet_boundary = (faces[j].negative_cell == Arrangement<3>::None);
                    //
                    if (face_on_tet_boundary) {
                        compute_iso_face_key(face_verts, key3);
                        auto iter_inserted = face_on_tetFace.try_emplace(key3, iso_faces.size());
                        if (iter_inserted.second) {
                            iso_fId_of_tet_face.back() = iso_faces.size();
                            iso_faces.emplace_back();
                            iso_faces.back().vert_indices = face_verts;
                            iso_faces.back().tet_face_indices.emplace_back(i, j);
                        } else { // iso_face inserted before
                            size_t iso_face_id = (iter_inserted.first)->second;
                            iso_fId_of_tet_face.back() = iso_face_id;
                            iso_faces[iso_face_id].tet_face_indices.emplace_back(i, j);
                        }
                    } else { // face not on tet boundary
                        iso_fId_of_tet_face.back() = iso_faces.size();
                        iso_faces.emplace_back();
                        iso_faces.back().vert_indices = face_verts;
                        iso_faces.back().tet_face_indices.emplace_back(i, j);
                    }
                }
            }
            iso_fId_start_index_of_tet.push_back(iso_fId_of_tet_face.size());
        }
    }
    //
}


void extract_MI_mesh(size_t num_2_func,
                          size_t num_3_func,
                          size_t num_more_func,
                          const std::vector<MaterialInterface<3>>& cut_results,
                          const std::vector<size_t>& cut_result_index,
                          const std::vector<size_t>& material_in_tet,
                          const std::vector<size_t>& start_index_of_tet,
                          const std::vector<std::array<size_t, 4>>& tets,
                          std::vector<MI_Vert>& MI_verts,
                          std::vector<PolygonFace>& MI_faces)
{
    size_t n_tets = tets.size();
    // estimate number of MI-verts and MI-faces
    // C(2,2) = 1, C(3,2) = 3, C(4,2) = 6
    size_t max_num_face = num_2_func + 3 * num_3_func + 6 * num_more_func;
    // for triangle mesh, 3|F| is about 2|F|, |V|-|E|+|F| is about 0, so |V| is about 0.5|F|.
    // Since most non-empty tets has 2 functions (i.e. 1 material interface),
    // most polygon's in the mesh are triangle or quad.
    // so number of triangles is roughly between |P| and 2|P| (|P| = number of polygons)
    // so, 0.5 * 2|P| = |P| is an estimate of |V|'s upper-bound
    size_t max_num_vert = max_num_face;
    MI_verts.reserve(max_num_vert);
    MI_faces.reserve(max_num_face);
    // hash table for vertices on the boundary of tetrahedron
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    // key is (edge v1, edge v2, material 1, material 2)
    absl::flat_hash_map<std::array<size_t, 4>, size_t> vert_on_tetEdge;
    vert_on_tetEdge.reserve(num_2_func + 2 * num_3_func + 4 * num_more_func);
    // key is (tri v1, tri v2, tri v3, material 1, material 2, material 3)
    absl::flat_hash_map<std::array<size_t, 6>, size_t> vert_on_tetFace;
    // 2func: 0 face vert, 3func: 2 face vert, 4func: 4 face vert
    // num face vert estimate = (2 * n_3func + 4 * n_moreFunc) /2
    vert_on_tetFace.reserve(num_3_func + 2 * num_more_func);
    // map: face on tet boundary -> material label
    absl::flat_hash_map<std::array<long long, 3>, size_t> material_of_bndry_face;
    // map: face on tet boundary -> material labels
    absl::flat_hash_map<std::array<long long, 3>, std::vector<size_t>> materials_of_bndry_face;
    std::vector<size_t> materials;
    materials.reserve(2);
    //
    std::vector<bool> is_MI_vert;
    is_MI_vert.reserve(8);
    std::vector<bool> is_MI_face;
    is_MI_face.reserve(9);
    std::vector<long long> MI_vId_of_vert;
    MI_vId_of_vert.reserve(8);
    std::vector<size_t> face_verts;
    face_verts.reserve(4);
    std::vector<long long> bndry_face_verts;
    bndry_face_verts.reserve(4);
    std::array<long long, 3> key3;
    std::array<size_t, 4> key4;
    std::array<size_t, 6> key6;
    std::array<bool, 4> used_mId;
    std::array<size_t, 2> vIds2;
    std::array<size_t, 3> vIds3;
    std::array<size_t, 4> func_mIds;
    std::array<size_t, 3> bndry_mIds;
    //
    for (size_t i = 0; i < n_tets; i++) {
        if (cut_result_index[i] != MaterialInterface<3>::None) {
            const auto& mInterface = cut_results[cut_result_index[i]];
            const auto& vertices = mInterface.vertices;
            const auto& faces = mInterface.faces;
            auto start_index = start_index_of_tet[i];
            auto num_func = start_index_of_tet[i + 1] - start_index;
            // find vertices and faces on material interface
            is_MI_vert.clear();
            for (int j = 0; j < vertices.size(); ++j) {
                is_MI_vert.push_back(false);
            }
            is_MI_face.clear();
            for (const auto& face : faces) {
                // material label 0,1,2,3 represents tet boundary
                if (face.positive_material_label > 3) {
                    is_MI_face.push_back(true);
                    for (auto vId : face.vertices) {
                        is_MI_vert[vId] = true;
                    }
                } else {
                    is_MI_face.push_back(false);
                }
            }
            // map: local vert index --> MI-vert index
            MI_vId_of_vert.clear();
            // create MI-vertices
            for (size_t j = 0; j < vertices.size(); j++) {
                size_t num_bndry_materials = 0;
                size_t num_func_materials = 0;
                const auto& vertex = vertices[j];
                // vertex.size() == 4
                for (size_t k = 0; k < 4; k++) {
                    if (vertex[k] > 3) { // material 0,1,2,3 are tet boundaries
                        func_mIds[num_func_materials] =
                                material_in_tet[vertex[k] - 4 + start_index];
                        ++num_func_materials;
                    } else {
                        bndry_mIds[num_bndry_materials] = vertex[k];
                        ++num_bndry_materials;
                    }
                }
                if (!is_MI_vert[j]) {
                    // vert not on any interior face, therefore must be a tet vertex
                    used_mId[0] = false;
                    used_mId[1] = false;
                    used_mId[2] = false;
                    used_mId[3] = false;
                    for (const auto& mId : bndry_mIds) {
                        used_mId[mId] = true;
                    }
                    size_t vId;
                    for (size_t k = 0; k < 4; k++) {
                        if (!used_mId[k]) {
                            vId = k;
                            break;
                        }
                    }
                    // tet vert i is mapped to (-i-1)
                    MI_vId_of_vert.push_back(-tets[i][vId] - 1);
                } else { // vert on an interior face
                    switch (num_bndry_materials) {
                        case 2: // on tet edge
                        {
                            used_mId[0] = false;
                            used_mId[1] = false;
                            used_mId[2] = false;
                            used_mId[3] = false;
                            used_mId[bndry_mIds[0]] = true;
                            used_mId[bndry_mIds[1]] = true;
                            size_t num_vIds = 0;
                            for (size_t k = 0; k < 4; k++) {
                                if (!used_mId[k]) {
                                    vIds2[num_vIds] = tets[i][k];
                                    ++num_vIds;
                                }
                            }
                            // sort {vId1, vId2}
                            size_t vId1 = vIds2[0];
                            size_t vId2 = vIds2[1];
                            if (vId1 > vId2) {
                                size_t tmp = vId1;
                                vId1 = vId2;
                                vId2 = tmp;
                            }
                            // sort {mId1, mId2}
                            size_t mId1 = func_mIds[0];
                            size_t mId2 = func_mIds[1];
                            if (mId1 > mId2) {
                                size_t tmp = mId1;
                                mId1 = mId2;
                                mId2 = tmp;
                            }
                            key4[0] = vId1;
                            key4[1] = vId2;
                            key4[2] = mId1;
                            key4[3] = mId2;
                            auto iter_inserted = vert_on_tetEdge.try_emplace(key4, MI_verts.size());
                            if (iter_inserted.second) {
                                MI_verts.emplace_back();
                                auto& MI_vert = MI_verts.back();
                                MI_vert.tet_index = i;
                                MI_vert.tet_vert_index = j;
                                MI_vert.simplex_size = 2;
                                MI_vert.simplex_vert_indices[0] = vId1;
                                MI_vert.simplex_vert_indices[1] = vId2;
                                MI_vert.material_indices[0] = mId1;
                                MI_vert.material_indices[1] = mId2;
                            }
                            MI_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        case 1: // on tet face
                        {
                            size_t mId = bndry_mIds[0];
                            size_t num_vIds = 0;
                            for (size_t k = 0; k < 4; k++) {
                                if (k != mId) {
                                    vIds3[num_vIds] = tets[i][k];
                                    ++num_vIds;
                                }
                            }
                            std::sort(vIds3.begin(), vIds3.end());
                            std::sort(func_mIds.begin(), func_mIds.begin() + 3);
                            key6[0] = vIds3[0];
                            key6[1] = vIds3[1];
                            key6[2] = vIds3[2];
                            key6[3] = func_mIds[0];
                            key6[4] = func_mIds[1];
                            key6[5] = func_mIds[2];
                            auto iter_inserted = vert_on_tetFace.try_emplace(key6, MI_verts.size());
                            if (iter_inserted.second) {
                                MI_verts.emplace_back();
                                auto& MI_vert = MI_verts.back();
                                MI_vert.tet_index = i;
                                MI_vert.tet_vert_index = j;
                                MI_vert.simplex_size = 3;
                                MI_vert.simplex_vert_indices[0] = vIds3[0];
                                MI_vert.simplex_vert_indices[1] = vIds3[1];
                                MI_vert.simplex_vert_indices[2] = vIds3[2];
                                MI_vert.material_indices[0] = func_mIds[0];
                                MI_vert.material_indices[1] = func_mIds[1];
                                MI_vert.material_indices[2] = func_mIds[2];
                            }
                            MI_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        case 0: // in tet cell
                        {
                            MI_vId_of_vert.push_back(MI_verts.size());
                            MI_verts.emplace_back();
                            auto& MI_vert = MI_verts.back();
                            MI_vert.tet_index = i;
                            MI_vert.tet_vert_index = j;
                            MI_vert.simplex_size = 4;
                            MI_vert.simplex_vert_indices = tets[i];
                            MI_vert.material_indices = func_mIds;
                            break;
                        }
                        case 3: // on tet vertex
                        {
                            used_mId[0] = false;
                            used_mId[1] = false;
                            used_mId[2] = false;
                            used_mId[3] = false;
                            for (const auto& mId : bndry_mIds) {
                                used_mId[mId] = true;
                            }
                            size_t vId;
                            for (size_t k = 0; k < 4; k++) {
                                if (!used_mId[k]) {
                                    vId = k;
                                    break;
                                }
                            }
                            auto key = tets[i][vId];
                            auto iter_inserted = vert_on_tetVert.try_emplace(key, MI_verts.size());
                            if (iter_inserted.second) {
                                MI_verts.emplace_back();
                                auto& MI_vert = MI_verts.back();
                                MI_vert.tet_index = i;
                                MI_vert.tet_vert_index = j;
                                MI_vert.simplex_size = 1;
                                MI_vert.simplex_vert_indices[0] = tets[i][vId];
                            }
                            MI_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        default: break;
                    }
                }
            }
            // create MI-faces
            for (size_t j = 0; j < faces.size(); j++) {
                const auto& face = faces[j];
                if (is_MI_face[j]) { // face j is in tet interior
                    face_verts.clear();
                    for (unsigned long vId : face.vertices) {
                        face_verts.push_back(MI_vId_of_vert[vId]);
                    }
                    MI_faces.emplace_back();
                    MI_faces.back().vert_indices = face_verts;
                    MI_faces.back().tet_face_indices.emplace_back(i, j);
                } else { // face j is on tet boundary
                    bndry_face_verts.clear();
                    for (unsigned long vId : face.vertices) {
                        if (MI_vId_of_vert[vId] >= 0 &&
                            MI_verts[MI_vId_of_vert[vId]].simplex_size == 1) {
                            // intersection on tet vertex
                            bndry_face_verts.push_back(
                                    -MI_verts[MI_vId_of_vert[vId]].simplex_vert_indices[0] - 1);
                        } else {
                            bndry_face_verts.push_back(MI_vId_of_vert[vId]);
                        }
                    }
                    compute_iso_face_key(bndry_face_verts, key3);
                    // get material label(s) of the cell incident to the face in the current tet
                    size_t m = MaterialInterface<3>::None;
                    if (mInterface.unique_materials.empty()) {
                        // no duplicate materials
                        m = material_in_tet[face.negative_material_label - 4 + start_index];
                    }
                    else {
                        size_t uId = mInterface.unique_material_indices[face.negative_material_label];
                        if (mInterface.unique_materials[uId].size() == 1) {
                            // face's incident cell has a unique material label
                            m = material_in_tet[face.negative_material_label - 4 + start_index];
                        } else {
                            // face's incident cell has multiple material labels
                            materials.clear();
                            for (auto mId : mInterface.unique_materials[uId]) {
                                materials.push_back(material_in_tet[mId - 4 + start_index]);
                            }
                        }
                    }
                    //
                    bool is_material_interface = false;
                    if (m != MaterialInterface<3>::None) {
                        // face's incident cell has unique material label
                        auto iter_inserted = material_of_bndry_face.try_emplace(key3, m);
                        if (!iter_inserted.second) { // inserted before
                            if (iter_inserted.first->second == m) {
                                // material labels on both sides of the face are the same
                                // the face is not on material interface
                                material_of_bndry_face.erase(iter_inserted.first);
                            } else {
                                // different material labels on different sides
                                // the face is on material interface
                                is_material_interface = true;
                            }
                        }
                        else {
                            auto iter = materials_of_bndry_face.find(key3);
                            if (iter != materials_of_bndry_face.end()) {
                                // the cell on the other side has multiple material labels
                                bool has_m = false;
                                for (auto mId : iter->second) {
                                    if (mId == m) {
                                        has_m = true;
                                        break;
                                    }
                                }
                                if (has_m) {
                                    // cells on two sides have a common material label
                                    // the face is not on material interface
                                    material_of_bndry_face.erase(iter_inserted.first);
                                    materials_of_bndry_face.erase(iter);
                                } else {
                                    // cells on two sides have no material label in common
                                    // the face is on material interface
                                    is_material_interface = true;
                                }
                            }
                        }
                    }
                    else {
                        // face's incident cell has multiple material labels
                        auto iter_inserted = materials_of_bndry_face.try_emplace(key3, materials);
                        if (!iter_inserted.second) { // inserted before
                            bool has_common_material_label = false;
                            for (auto mi : iter_inserted.first->second) {
                                for (auto mj : materials) {
                                    if (mi == mj) {
                                        has_common_material_label = true;
                                        break;
                                    }
                                }
                            }
                            if (has_common_material_label) {
                                materials_of_bndry_face.erase(iter_inserted.first);
                            } else {
                                // different material labels on different sides
                                // the face is on material interface
                                is_material_interface = true;
                            }
                        }
                        else {
                            auto iter = material_of_bndry_face.find(key3);
                            if (iter != material_of_bndry_face.end()) {
                                // the cell on the other side has a unique material label
                                m = iter->second;
                                bool has_m = false;
                                for (auto mId : materials) {
                                    if (mId == m) {
                                        has_m = true;
                                        break;
                                    }
                                }
                                if (has_m) {
                                    // cells on two sides have a common material label
                                    // the face is not on material interface
                                    material_of_bndry_face.erase(iter);
                                    materials_of_bndry_face.erase(iter_inserted.first);
                                } else {
                                    // cells on two sides have no material label in common
                                    // the face is on material interface
                                    is_material_interface = true;
                                }
                            }
                        }
                    }
                    //
                    if (is_material_interface) {
                        face_verts.clear();
                        for (size_t k = 0; k < face.vertices.size(); ++k) {
                            if (bndry_face_verts[k] < 0) {
                                // try to create new vert on tet vertex
                                size_t key = -bndry_face_verts[k] - 1;
                                auto iterInserted =
                                        vert_on_tetVert.try_emplace(key, MI_verts.size());
                                if (iterInserted.second) {
                                    MI_verts.emplace_back();
                                    auto& MI_vert = MI_verts.back();
                                    MI_vert.tet_index = i;
                                    //                                        MI_vert.tet_vert_index = j;
                                    MI_vert.tet_vert_index = face.vertices[k];
                                    MI_vert.simplex_size = 1;
                                    MI_vert.simplex_vert_indices[0] = key;
                                }
                                //                                    MI_vId_of_vert.push_back(iterInserted.first->second);
                                face_verts.push_back(iterInserted.first->second);
                            } else {
                                face_verts.push_back(MI_vId_of_vert[face.vertices[k]]);
                            }
                        }
                        // insert the face to MI_faces
                        MI_faces.emplace_back();
                        MI_faces.back().vert_indices = face_verts;
                        MI_faces.back().tet_face_indices.emplace_back(i, j);
                    }
                }
            }
        }
    }
}

void extract_MI_mesh(size_t num_2_func,
                     size_t num_3_func,
                     size_t num_more_func,
                     const std::vector<MaterialInterface<3>>& cut_results,
                     const std::vector<size_t>& cut_result_index,
                     const std::vector<size_t>& material_in_tet,
                     const std::vector<size_t>& start_index_of_tet,
                     const std::vector<std::array<size_t, 4>>& tets,
                     std::vector<MI_Vert>& MI_verts,
                     std::vector<PolygonFace>& MI_faces,
                     std::vector<long long int>& global_vId_of_tet_vert,
                     std::vector<size_t>& global_vId_start_index_of_tet,
                     std::vector<size_t>& MI_fId_of_tet_face,
                     std::vector<size_t>& MI_fId_start_index_of_tet)
{
    size_t n_tets = tets.size();
    // get total number of verts and faces
    size_t total_num_vert = 0;
    size_t total_num_face = 0;
    for (const auto& result : cut_results) {
        total_num_vert += result.vertices.size();
        total_num_face += result.faces.size();
    }
    global_vId_of_tet_vert.reserve(total_num_vert);
    global_vId_start_index_of_tet.reserve(n_tets + 1);
    global_vId_start_index_of_tet.push_back(0);
    MI_fId_of_tet_face.reserve(total_num_face);
    MI_fId_start_index_of_tet.reserve(n_tets + 1);
    MI_fId_start_index_of_tet.push_back(0);
    // estimate number of MI-verts and MI-faces
    // C(2,2) = 1, C(3,2) = 3, C(4,2) = 6
    size_t max_num_face = num_2_func + 3 * num_3_func + 6 * num_more_func;
    // for triangle mesh, 3|F| is about 2|F|, |V|-|E|+|F| is about 0, so |V| is about 0.5|F|.
    // Since most non-empty tets has 2 functions (i.e. 1 material interface),
    // most polygon's in the mesh are triangle or quad.
    // so number of triangles is roughly between |P| and 2|P| (|P| = number of polygons)
    // so, 0.5 * 2|P| = |P| is an estimate of |V|'s upper-bound
    size_t max_num_vert = max_num_face;
    MI_verts.reserve(max_num_vert);
    MI_faces.reserve(max_num_face);
    // hash table for vertices on the boundary of tetrahedron
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    // key is (edge v1, edge v2, material 1, material 2)
    absl::flat_hash_map<std::array<size_t, 4>, size_t> vert_on_tetEdge;
    vert_on_tetEdge.reserve(num_2_func + 2 * num_3_func + 4 * num_more_func);
    // key is (tri v1, tri v2, tri v3, material 1, material 2, material 3)
    absl::flat_hash_map<std::array<size_t, 6>, size_t> vert_on_tetFace;
    // 2func: 0 face vert, 3func: 2 face vert, 4func: 4 face vert
    // num face vert estimate = (2 * n_3func + 4 * n_moreFunc) /2
    vert_on_tetFace.reserve(num_3_func + 2 * num_more_func);
    // map: face on tet boundary -> material label
    absl::flat_hash_map<std::array<long long, 3>, size_t> material_of_bndry_face;
    // map: face on tet boundary -> material labels
    absl::flat_hash_map<std::array<long long, 3>, std::vector<size_t>> materials_of_bndry_face;
    std::vector<size_t> materials;
    materials.reserve(2);
//    absl::flat_hash_map<std::array<long long, 3>, std::array<size_t,3>> mat_tet_face_of_bndry_face;
    // map: face on tet boundary -> (tetId, faceId in tet)
    absl::flat_hash_map<std::array<long long, 3>, std::pair<size_t, size_t>> tet_face_of_bndry_face;
    //
    std::vector<bool> is_MI_vert;
    is_MI_vert.reserve(8);
    std::vector<bool> is_MI_face;
    is_MI_face.reserve(9);
    std::vector<long long> MI_vId_of_vert;
    MI_vId_of_vert.reserve(8);
    std::vector<size_t> face_verts;
    face_verts.reserve(4);
    std::vector<long long> bndry_face_verts;
    bndry_face_verts.reserve(4);
    std::array<long long, 3> key3;
    std::array<size_t, 4> key4;
    std::array<size_t, 6> key6;
    std::array<bool, 4> used_mId;
    std::array<size_t, 2> vIds2;
    std::array<size_t, 3> vIds3;
    std::array<size_t, 4> func_mIds;
    std::array<size_t, 3> bndry_mIds;
    //
    for (size_t i = 0; i < n_tets; i++) {
        if (cut_result_index[i] == MaterialInterface<3>::None) {
            global_vId_start_index_of_tet.push_back(global_vId_of_tet_vert.size());
            MI_fId_start_index_of_tet.push_back(MI_fId_of_tet_face.size());
        } else {
            const auto& mInterface = cut_results[cut_result_index[i]];
            const auto& vertices = mInterface.vertices;
            const auto& faces = mInterface.faces;
            auto start_index = start_index_of_tet[i];
            auto num_func = start_index_of_tet[i + 1] - start_index;
            // find vertices and faces on material interface
            is_MI_vert.clear();
            for (int j = 0; j < vertices.size(); ++j) {
                is_MI_vert.push_back(false);
            }
            is_MI_face.clear();
            for (const auto& face : faces) {
                // material label 0,1,2,3 represents tet boundary
                if (face.positive_material_label > 3) {
                    is_MI_face.push_back(true);
                    for (auto vId : face.vertices) {
                        is_MI_vert[vId] = true;
                    }
                } else {
                    is_MI_face.push_back(false);
                }
            }
            // map: local vert index --> MI-vert index
            MI_vId_of_vert.clear();
            // create MI-vertices
            for (size_t j = 0; j < vertices.size(); j++) {
                size_t num_bndry_materials = 0;
                size_t num_func_materials = 0;
                const auto& vertex = vertices[j];
                // vertex.size() == 4
                for (size_t k = 0; k < 4; k++) {
                    if (vertex[k] > 3) { // material 0,1,2,3 are tet boundaries
                        func_mIds[num_func_materials] =
                                material_in_tet[vertex[k] - 4 + start_index];
                        ++num_func_materials;
                    } else {
                        bndry_mIds[num_bndry_materials] = vertex[k];
                        ++num_bndry_materials;
                    }
                }
                if (!is_MI_vert[j]) {
                    // vert not on any interior face, therefore must be a tet vertex
                    used_mId[0] = false;
                    used_mId[1] = false;
                    used_mId[2] = false;
                    used_mId[3] = false;
                    for (const auto& mId : bndry_mIds) {
                        used_mId[mId] = true;
                    }
                    size_t vId;
                    for (size_t k = 0; k < 4; k++) {
                        if (!used_mId[k]) {
                            vId = k;
                            break;
                        }
                    }
                    // tet vert i is mapped to (-i-1)
                    MI_vId_of_vert.push_back(-tets[i][vId] - 1);
                    global_vId_of_tet_vert.push_back(-tets[i][vId] - 1);
                } else { // vert on an interior face
                    switch (num_bndry_materials) {
                        case 2: // on tet edge
                        {
                            used_mId[0] = false;
                            used_mId[1] = false;
                            used_mId[2] = false;
                            used_mId[3] = false;
                            used_mId[bndry_mIds[0]] = true;
                            used_mId[bndry_mIds[1]] = true;
                            size_t num_vIds = 0;
                            for (size_t k = 0; k < 4; k++) {
                                if (!used_mId[k]) {
                                    vIds2[num_vIds] = tets[i][k];
                                    ++num_vIds;
                                }
                            }
                            // sort {vId1, vId2}
                            size_t vId1 = vIds2[0];
                            size_t vId2 = vIds2[1];
                            if (vId1 > vId2) {
                                size_t tmp = vId1;
                                vId1 = vId2;
                                vId2 = tmp;
                            }
                            // sort {mId1, mId2}
                            size_t mId1 = func_mIds[0];
                            size_t mId2 = func_mIds[1];
                            if (mId1 > mId2) {
                                size_t tmp = mId1;
                                mId1 = mId2;
                                mId2 = tmp;
                            }
                            key4[0] = vId1;
                            key4[1] = vId2;
                            key4[2] = mId1;
                            key4[3] = mId2;
                            auto iter_inserted = vert_on_tetEdge.try_emplace(key4, MI_verts.size());
                            if (iter_inserted.second) {
                                MI_verts.emplace_back();
                                auto& MI_vert = MI_verts.back();
                                MI_vert.tet_index = i;
                                MI_vert.tet_vert_index = j;
                                MI_vert.simplex_size = 2;
                                MI_vert.simplex_vert_indices[0] = vId1;
                                MI_vert.simplex_vert_indices[1] = vId2;
                                MI_vert.material_indices[0] = mId1;
                                MI_vert.material_indices[1] = mId2;
                            }
                            global_vId_of_tet_vert.push_back(iter_inserted.first->second);
                            MI_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        case 1: // on tet face
                        {
                            size_t mId = bndry_mIds[0];
                            size_t num_vIds = 0;
                            for (size_t k = 0; k < 4; k++) {
                                if (k != mId) {
                                    vIds3[num_vIds] = tets[i][k];
                                    ++num_vIds;
                                }
                            }
                            std::sort(vIds3.begin(), vIds3.end());
                            std::sort(func_mIds.begin(), func_mIds.begin() + 3);
                            key6[0] = vIds3[0];
                            key6[1] = vIds3[1];
                            key6[2] = vIds3[2];
                            key6[3] = func_mIds[0];
                            key6[4] = func_mIds[1];
                            key6[5] = func_mIds[2];
                            auto iter_inserted = vert_on_tetFace.try_emplace(key6, MI_verts.size());
                            if (iter_inserted.second) {
                                MI_verts.emplace_back();
                                auto& MI_vert = MI_verts.back();
                                MI_vert.tet_index = i;
                                MI_vert.tet_vert_index = j;
                                MI_vert.simplex_size = 3;
                                MI_vert.simplex_vert_indices[0] = vIds3[0];
                                MI_vert.simplex_vert_indices[1] = vIds3[1];
                                MI_vert.simplex_vert_indices[2] = vIds3[2];
                                MI_vert.material_indices[0] = func_mIds[0];
                                MI_vert.material_indices[1] = func_mIds[1];
                                MI_vert.material_indices[2] = func_mIds[2];
                            }
                            global_vId_of_tet_vert.push_back(iter_inserted.first->second);
                            MI_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        case 0: // in tet cell
                        {
                            global_vId_of_tet_vert.push_back(MI_verts.size());
                            MI_vId_of_vert.push_back(MI_verts.size());
                            MI_verts.emplace_back();
                            auto& MI_vert = MI_verts.back();
                            MI_vert.tet_index = i;
                            MI_vert.tet_vert_index = j;
                            MI_vert.simplex_size = 4;
                            MI_vert.simplex_vert_indices = tets[i];
                            MI_vert.material_indices = func_mIds;
                            break;
                        }
                        case 3: // on tet vertex
                        {
                            used_mId[0] = false;
                            used_mId[1] = false;
                            used_mId[2] = false;
                            used_mId[3] = false;
                            for (const auto& mId : bndry_mIds) {
                                used_mId[mId] = true;
                            }
                            size_t vId;
                            for (size_t k = 0; k < 4; k++) {
                                if (!used_mId[k]) {
                                    vId = k;
                                    break;
                                }
                            }
                            auto key = tets[i][vId];
                            auto iter_inserted = vert_on_tetVert.try_emplace(key, MI_verts.size());
                            if (iter_inserted.second) {
                                MI_verts.emplace_back();
                                auto& MI_vert = MI_verts.back();
                                MI_vert.tet_index = i;
                                MI_vert.tet_vert_index = j;
                                MI_vert.simplex_size = 1;
                                MI_vert.simplex_vert_indices[0] = tets[i][vId];
                            }
                            global_vId_of_tet_vert.push_back(-key - 1);
                            MI_vId_of_vert.push_back(iter_inserted.first->second);
                            break;
                        }
                        default: break;
                    }
                }
            }
            global_vId_start_index_of_tet.push_back(global_vId_of_tet_vert.size());
            // create MI-faces
            for (size_t j = 0; j < faces.size(); j++) {
                MI_fId_of_tet_face.push_back(MaterialInterface<3>::None);
                const auto& face = faces[j];
                if (is_MI_face[j]) { // face j is in tet interior
                    face_verts.clear();
                    for (unsigned long vId : face.vertices) {
                        face_verts.push_back(MI_vId_of_vert[vId]);
                    }
                    MI_fId_of_tet_face.back() = MI_faces.size();
                    MI_faces.emplace_back();
                    MI_faces.back().vert_indices = face_verts;
                    MI_faces.back().tet_face_indices.emplace_back(i, j);
                } else { // face j is on tet boundary
                    bndry_face_verts.clear();
                    for (unsigned long vId : face.vertices) {
                        if (MI_vId_of_vert[vId] >= 0 &&
                            MI_verts[MI_vId_of_vert[vId]].simplex_size == 1) {
                            // intersection on tet vertex
                            bndry_face_verts.push_back(
                                    -MI_verts[MI_vId_of_vert[vId]].simplex_vert_indices[0] - 1);
                        } else {
                            bndry_face_verts.push_back(MI_vId_of_vert[vId]);
                        }
                    }
                    compute_iso_face_key(bndry_face_verts, key3);
                    // get material label(s) of the cell incident to the face in the current tet
                    size_t m = MaterialInterface<3>::None;
                    if (mInterface.unique_materials.empty()) {
                        // no duplicate materials
                        m = material_in_tet[face.negative_material_label - 4 + start_index];
                    }
                    else {
                        size_t uId = mInterface.unique_material_indices[face.negative_material_label];
                        if (mInterface.unique_materials[uId].size() == 1) {
                            // face's incident cell has a unique material label
                            m = material_in_tet[face.negative_material_label - 4 + start_index];
                        } else {
                            // face's incident cell has multiple material labels
                            materials.clear();
                            for (auto mId : mInterface.unique_materials[uId]) {
                                materials.push_back(material_in_tet[mId - 4 + start_index]);
                            }
                        }
                    }
                    //
                    bool is_material_interface = false;
                    if (m != MaterialInterface<3>::None) {
                        // face's incident cell has unique material label
                        auto iter_inserted = material_of_bndry_face.try_emplace(key3, m);
                        if (!iter_inserted.second) { // inserted before
                            if (iter_inserted.first->second == m) {
                                // material labels on both sides of the face are the same
                                // the face is not on material interface
                                material_of_bndry_face.erase(iter_inserted.first);
                            } else {
                                // different material labels on different sides
                                // the face is on material interface
                                is_material_interface = true;
                            }
                        }
                        else {
                            auto iter = materials_of_bndry_face.find(key3);
                            if (iter != materials_of_bndry_face.end()) {
                                // the cell on the other side has multiple material labels
                                bool has_m = false;
                                for (auto mId : iter->second) {
                                    if (mId == m) {
                                        has_m = true;
                                        break;
                                    }
                                }
                                if (has_m) {
                                    // cells on two sides have a common material label
                                    // the face is not on material interface
                                    material_of_bndry_face.erase(iter_inserted.first);
                                    materials_of_bndry_face.erase(iter);
                                } else {
                                    // cells on two sides have no material label in common
                                    // the face is on material interface
                                    is_material_interface = true;
                                }
                            }
                        }
                    }
                    else {
                        // face's incident cell has multiple material labels
                        auto iter_inserted = materials_of_bndry_face.try_emplace(key3, materials);
                        if (!iter_inserted.second) { // inserted before
                            bool has_common_material_label = false;
                            for (auto mi : iter_inserted.first->second) {
                                for (auto mj : materials) {
                                    if (mi == mj) {
                                        has_common_material_label = true;
                                        break;
                                    }
                                }
                            }
                            if (has_common_material_label) {
                                materials_of_bndry_face.erase(iter_inserted.first);
                            } else {
                                // different material labels on different sides
                                // the face is on material interface
                                is_material_interface = true;
                            }
                        }
                        else {
                            auto iter = material_of_bndry_face.find(key3);
                            if (iter != material_of_bndry_face.end()) {
                                // the cell on the other side has a unique material label
                                m = iter->second;
                                bool has_m = false;
                                for (auto mId : materials) {
                                    if (mId == m) {
                                        has_m = true;
                                        break;
                                    }
                                }
                                if (has_m) {
                                    // cells on two sides have a common material label
                                    // the face is not on material interface
                                    material_of_bndry_face.erase(iter);
                                    materials_of_bndry_face.erase(iter_inserted.first);
                                } else {
                                    // cells on two sides have no material label in common
                                    // the face is on material interface
                                    is_material_interface = true;
                                }
                            }
                        }
                    }
                    //
                    auto iter_inserted = tet_face_of_bndry_face.try_emplace(key3, std::make_pair(i,j));
                    if (is_material_interface) {
                        // assert(!iter_inserted.second)
                        auto opposite_tetId = iter_inserted.first->second.first;
                        auto opposite_faceId = iter_inserted.first->second.second;
                        face_verts.clear();
                        for (size_t k = 0; k < face.vertices.size(); ++k) {
                            if (bndry_face_verts[k] < 0) {
                                // try to create new vert on tet vertex
                                size_t key = -bndry_face_verts[k] - 1;
                                auto iterInserted =
                                        vert_on_tetVert.try_emplace(key, MI_verts.size());
                                if (iterInserted.second) {
                                    MI_verts.emplace_back();
                                    auto& MI_vert = MI_verts.back();
                                    MI_vert.tet_index = i;
                                    //                                        MI_vert.tet_vert_index = j;
                                    MI_vert.tet_vert_index = face.vertices[k];
                                    MI_vert.simplex_size = 1;
                                    MI_vert.simplex_vert_indices[0] = key;
                                }
                                //MI_vId_of_vert.push_back(iterInserted.first->second);
                                face_verts.push_back(iterInserted.first->second);
                            } else {
                                face_verts.push_back(MI_vId_of_vert[face.vertices[k]]);
                            }
                        }
                        // insert the face to MI_faces
                        MI_fId_of_tet_face.back() = MI_faces.size();
                        MI_fId_of_tet_face[MI_fId_start_index_of_tet[opposite_tetId] + opposite_faceId] = MI_faces.size();
                        MI_faces.emplace_back();
                        MI_faces.back().vert_indices = face_verts;
                        MI_faces.back().tet_face_indices.emplace_back(i, j);
                    }
                }
            }
            MI_fId_start_index_of_tet.push_back(MI_fId_of_tet_face.size());
        }
    }
}


void compute_iso_vert_xyz(const std::vector<IsoVert>& iso_verts,
                          const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcVals,
                          const std::vector<std::array<double, 3>>& pts,
                          std::vector<std::array<double, 3>>& iso_pts)
{
    iso_pts.resize(iso_verts.size());
    std::array<double, 2> b2;
    std::array<double, 3> f1s3;
    std::array<double, 3> f2s3;
    std::array<double, 3> b3;
    std::array<double, 4> f1s4;
    std::array<double, 4> f2s4;
    std::array<double, 4> f3s4;
    std::array<double, 4> b4;
    for (size_t i = 0; i < iso_verts.size(); i++) {
        const auto& iso_vert = iso_verts[i];
        switch (iso_vert.simplex_size) {
            case 2: // on tet edge
            {
                auto vId1 = iso_vert.simplex_vert_indices[0];
                auto vId2 = iso_vert.simplex_vert_indices[1];
                auto fId = iso_vert.func_indices[0];
                auto f1 = funcVals(vId1, fId);
                auto f2 = funcVals(vId2, fId);
                //
                compute_barycentric_coords(f1, f2, b2);
                iso_pts[i][0] = b2[0] * pts[vId1][0] + b2[1] * pts[vId2][0];
                iso_pts[i][1] = b2[0] * pts[vId1][1] + b2[1] * pts[vId2][1];
                iso_pts[i][2] = b2[0] * pts[vId1][2] + b2[1] * pts[vId2][2];
                break;
            }
            case 3: // on tet face
            {
                auto vId1 = iso_vert.simplex_vert_indices[0];
                auto vId2 = iso_vert.simplex_vert_indices[1];
                auto vId3 = iso_vert.simplex_vert_indices[2];
                auto fId1 = iso_vert.func_indices[0];
                auto fId2 = iso_vert.func_indices[1];
                //
                f1s3[0] = funcVals(vId1, fId1);
                f1s3[1] = funcVals(vId2, fId1);
                f1s3[2] = funcVals(vId3, fId1);
                //
                f2s3[0] = funcVals(vId1, fId2);
                f2s3[1] = funcVals(vId2, fId2);
                f2s3[2] = funcVals(vId3, fId2);
                //
                compute_barycentric_coords(f1s3, f2s3, b3);
                iso_pts[i][0] = b3[0] * pts[vId1][0] + b3[1] * pts[vId2][0] + b3[2] * pts[vId3][0];
                iso_pts[i][1] = b3[0] * pts[vId1][1] + b3[1] * pts[vId2][1] + b3[2] * pts[vId3][1];
                iso_pts[i][2] = b3[0] * pts[vId1][2] + b3[1] * pts[vId2][2] + b3[2] * pts[vId3][2];
                break;
            }
            case 4: // in tet cell
            {
                auto vId1 = iso_vert.simplex_vert_indices[0];
                auto vId2 = iso_vert.simplex_vert_indices[1];
                auto vId3 = iso_vert.simplex_vert_indices[2];
                auto vId4 = iso_vert.simplex_vert_indices[3];
                auto fId1 = iso_vert.func_indices[0];
                auto fId2 = iso_vert.func_indices[1];
                auto fId3 = iso_vert.func_indices[2];
                f1s4[0] = funcVals(vId1, fId1);
                f1s4[1] = funcVals(vId2, fId1);
                f1s4[2] = funcVals(vId3, fId1);
                f1s4[3] = funcVals(vId4, fId1);
                //
                f2s4[0] = funcVals(vId1, fId2);
                f2s4[1] = funcVals(vId2, fId2);
                f2s4[2] = funcVals(vId3, fId2);
                f2s4[3] = funcVals(vId4, fId2);
                //
                f3s4[0] = funcVals(vId1, fId3);
                f3s4[1] = funcVals(vId2, fId3);
                f3s4[2] = funcVals(vId3, fId3);
                f3s4[3] = funcVals(vId4, fId3);
                //
                compute_barycentric_coords(f1s4, f2s4, f3s4, b4);
                iso_pts[i][0] = b4[0] * pts[vId1][0] + b4[1] * pts[vId2][0] + b4[2] * pts[vId3][0] +
                                b4[3] * pts[vId4][0];
                iso_pts[i][1] = b4[0] * pts[vId1][1] + b4[1] * pts[vId2][1] + b4[2] * pts[vId3][1] +
                                b4[3] * pts[vId4][1];
                iso_pts[i][2] = b4[0] * pts[vId1][2] + b4[1] * pts[vId2][2] + b4[2] * pts[vId3][2] +
                                b4[3] * pts[vId4][2];
                break;
            }
            case 1: // on tet vertex
                iso_pts[i] = pts[iso_vert.simplex_vert_indices[0]];
                break;
            default: break;
        }
    }
}


void compute_MI_vert_xyz(const std::vector<MI_Vert>& MI_verts,
                         const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcVals,
                         const std::vector<std::array<double, 3>>& pts,
                         std::vector<std::array<double, 3>>& MI_pts)
{
    MI_pts.resize(MI_verts.size());
    std::array<double, 2> b2;
    std::array<double, 3> f1s3;
    std::array<double, 3> f2s3;
    std::array<double, 3> b3;
    std::array<double, 4> f1s4;
    std::array<double, 4> f2s4;
    std::array<double, 4> f3s4;
    std::array<double, 4> b4;
    for (size_t i = 0; i < MI_verts.size(); i++) {
        const auto& MI_vert = MI_verts[i];
        switch (MI_vert.simplex_size) {
            case 2: // on tet edge
            {
                auto vId1 = MI_vert.simplex_vert_indices[0];
                auto vId2 = MI_vert.simplex_vert_indices[1];
                auto fId1 = MI_vert.material_indices[0];
                auto fId2 = MI_vert.material_indices[1];
                auto f1 = funcVals(vId1,fId1) - funcVals(vId1,fId2);
                auto f2 = funcVals(vId2,fId1) - funcVals(vId2,fId2);
                //
                compute_barycentric_coords(f1, f2, b2);
                MI_pts[i][0] = b2[0] * pts[vId1][0] + b2[1] * pts[vId2][0];
                MI_pts[i][1] = b2[0] * pts[vId1][1] + b2[1] * pts[vId2][1];
                MI_pts[i][2] = b2[0] * pts[vId1][2] + b2[1] * pts[vId2][2];
                break;
            }
            case 3: // on tet face
            {
                auto vId1 = MI_vert.simplex_vert_indices[0];
                auto vId2 = MI_vert.simplex_vert_indices[1];
                auto vId3 = MI_vert.simplex_vert_indices[2];
                auto fId1 = MI_vert.material_indices[0];
                auto fId2 = MI_vert.material_indices[1];
                auto fId3 = MI_vert.material_indices[2];
                // f1 - f2
                f1s3[0] = funcVals(vId1,fId1) - funcVals(vId1,fId2);
                f1s3[1] = funcVals(vId2,fId1) - funcVals(vId2,fId2);
                f1s3[2] = funcVals(vId3,fId1) - funcVals(vId3,fId2);
                // f2 - f3
                f2s3[0] = funcVals(vId1,fId2) - funcVals(vId1,fId3);
                f2s3[1] = funcVals(vId2,fId2) - funcVals(vId2,fId3);
                f2s3[2] = funcVals(vId3,fId2) - funcVals(vId3,fId3);
                //
                compute_barycentric_coords(f1s3, f2s3, b3);
                MI_pts[i][0] = b3[0] * pts[vId1][0] + b3[1] * pts[vId2][0] + b3[2] * pts[vId3][0];
                MI_pts[i][1] = b3[0] * pts[vId1][1] + b3[1] * pts[vId2][1] + b3[2] * pts[vId3][1];
                MI_pts[i][2] = b3[0] * pts[vId1][2] + b3[1] * pts[vId2][2] + b3[2] * pts[vId3][2];
                break;
            }
            case 4: // in tet cell
            {
                auto vId1 = MI_vert.simplex_vert_indices[0];
                auto vId2 = MI_vert.simplex_vert_indices[1];
                auto vId3 = MI_vert.simplex_vert_indices[2];
                auto vId4 = MI_vert.simplex_vert_indices[3];
                auto fId1 = MI_vert.material_indices[0];
                auto fId2 = MI_vert.material_indices[1];
                auto fId3 = MI_vert.material_indices[2];
                auto fId4 = MI_vert.material_indices[3];
                // f1 - f2
                f1s4[0] = funcVals(vId1,fId1) - funcVals(vId1,fId2);
                f1s4[1] = funcVals(vId2,fId1) - funcVals(vId2,fId2);
                f1s4[2] = funcVals(vId3,fId1) - funcVals(vId3,fId2);
                f1s4[3] = funcVals(vId4,fId1) - funcVals(vId4,fId2);
                // f2 - f3
                f2s4[0] = funcVals(vId1,fId2) - funcVals(vId1,fId3);
                f2s4[1] = funcVals(vId2,fId2) - funcVals(vId2,fId3);
                f2s4[2] = funcVals(vId3,fId2) - funcVals(vId3,fId3);
                f2s4[3] = funcVals(vId4,fId2) - funcVals(vId4,fId3);
                // f3 - f4
                f3s4[0] = funcVals(vId1,fId3) - funcVals(vId1,fId4);
                f3s4[1] = funcVals(vId2,fId3) - funcVals(vId2,fId4);
                f3s4[2] = funcVals(vId3,fId3) - funcVals(vId3,fId4);
                f3s4[3] = funcVals(vId4,fId3) - funcVals(vId4,fId4);
                //
                compute_barycentric_coords(f1s4, f2s4, f3s4, b4);
                MI_pts[i][0] = b4[0] * pts[vId1][0] + b4[1] * pts[vId2][0] + b4[2] * pts[vId3][0] +
                               b4[3] * pts[vId4][0];
                MI_pts[i][1] = b4[0] * pts[vId1][1] + b4[1] * pts[vId2][1] + b4[2] * pts[vId3][1] +
                               b4[3] * pts[vId4][1];
                MI_pts[i][2] = b4[0] * pts[vId1][2] + b4[1] * pts[vId2][2] + b4[2] * pts[vId3][2] +
                               b4[3] * pts[vId4][2];
                break;
            }
            case 1: // on tet vertex
                MI_pts[i] = pts[MI_vert.simplex_vert_indices[0]];
                break;
            default: break;
        }
    }
}
