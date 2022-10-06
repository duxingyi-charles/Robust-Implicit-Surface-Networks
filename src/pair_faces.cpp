//
// Created by Charles Du on 7/20/22.
//

#include "pair_faces.h"
#include "extract_mesh.h"

using namespace simplicial_arrangement;

class Degenerate_Iso_Edge_Exception : public std::exception {
    virtual const char *what() const throw() {
        return "More than one tets containing the iso-faces incident to the iso-edge.";
    }
} degenerate_iso_edge_exception;

void compute_face_order(const Edge &iso_edge, const std::vector<std::array<size_t, 4>> &tets,
                        const std::vector<IsoVert> &iso_verts,
                        const std::vector<PolygonFace> &iso_faces,
                        const std::vector<Arrangement<3>> &cut_results,
                        const std::vector<size_t> &cut_result_index,
                        const std::vector<size_t> &func_in_tet,
                        const std::vector<size_t> &start_index_of_tet,
                        const absl::flat_hash_map<size_t, std::vector<size_t>> &incident_tets,
                        std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>> &ordered_face_pairs) {
// collect all iso-faces incident to the iso-edge
    const auto &face_edge_indices = iso_edge.face_edge_indices;
// collect all tets containing the incident iso-faces
    std::set < size_t > containing_tets;
    for (const auto &face_edge: face_edge_indices) {
        const auto &tetrahedrons = iso_faces[face_edge.first].tet_face_indices;
        for (const auto &t: tetrahedrons) {
            containing_tets.insert(t.first);
        }
    }
// pair faces
    if (containing_tets.size() == 1) {
//        std::cout << ">>>>>>>> iso-edge in tet" << std::endl;
        auto tet_id = *containing_tets.begin();
        pair_faces_in_one_tet(
                cut_results[cut_result_index[tet_id]], iso_faces, iso_edge,
                ordered_face_pairs);
    } else {
        size_t v1 = iso_edge.v1;
        size_t v2 = iso_edge.v2;
        const auto &iso_v1 = iso_verts[v1];
        const auto &iso_v2 = iso_verts[v2];
        bool on_tet_edge = false;
        size_t vId1, vId2, vId3;
        if (iso_v1.simplex_size == 1 && iso_v2.simplex_size == 1) {
            on_tet_edge = true;
            vId1 = iso_v1.simplex_vert_indices[0];
            vId2 = iso_v2.simplex_vert_indices[0];
        } else if (iso_v1.simplex_size == 2 && iso_v2.simplex_size == 1) {
            vId1 = iso_v1.simplex_vert_indices[0];
            vId2 = iso_v1.simplex_vert_indices[1];
            vId3 = iso_v2.simplex_vert_indices[0];
            on_tet_edge = (vId3 == vId1 || vId3 == vId2);
        } else if (iso_v1.simplex_size == 1 && iso_v2.simplex_size == 2) {
            vId1 = iso_v2.simplex_vert_indices[0];
            vId2 = iso_v2.simplex_vert_indices[1];
            vId3 = iso_v1.simplex_vert_indices[0];
            on_tet_edge = (vId3 == vId1 || vId3 == vId2);
        } else if (iso_v1.simplex_size == 2 && iso_v2.simplex_size == 2) {
            vId1 = iso_v1.simplex_vert_indices[0];
            vId2 = iso_v1.simplex_vert_indices[1];
            // assume: indices in Iso_Vert::simplex_vert_indices are sorted
            on_tet_edge = (vId1 == iso_v2.simplex_vert_indices[0] &&
                           vId2 == iso_v2.simplex_vert_indices[1]);
        }
        if (on_tet_edge) {
//            std::cout << ">>>>>>>> iso-edge on tet edge" << std::endl;
            // iso_edge lies on tet edge (vId1, vId2)
            // tet vertices vId1, vId2 must be degenerate vertices
            // find all tets incident to edge (vId1, vId2)
            std::unordered_set < size_t > incident_tets1(incident_tets.at(vId1).begin(), incident_tets.at(vId1).end());
            std::vector<size_t> common_tetIds;
            for (size_t tId: incident_tets.at(vId2)) {
                if (incident_tets1.find(tId) != incident_tets1.end()) {
                    common_tetIds.push_back(tId);
                }
            }
            // pair half-faces
            pair_faces_in_tets(iso_edge, {vId1, vId2}, common_tetIds,
                               tets, iso_faces, cut_results, cut_result_index,
                               func_in_tet, start_index_of_tet,
                               ordered_face_pairs);
        } else {
//            std::cout << ">>>>>>>> iso-edge on tet face" << std::endl;
            // iso_edge lies on a tet boundary face
            // assert: containing_tets.size() == 2
            size_t tet_id1 = *containing_tets.begin();
            size_t tet_id2 = *containing_tets.rbegin();
            std::unordered_set < size_t > tet1_vIds(tets[tet_id1].begin(), tets[tet_id1].end());
            std::vector<size_t> common_vIds;
            common_vIds.reserve(3);
            for (size_t vId: tets[tet_id2]) {
                if (tet1_vIds.find(vId) != tet1_vIds.end()) {
                    common_vIds.push_back(vId);
                }
            }
            // pair half-faces
            pair_faces_in_tets(iso_edge, common_vIds, {tet_id1, tet_id2},
                               tets, iso_faces, cut_results, cut_result_index,
                               func_in_tet, start_index_of_tet,
                               ordered_face_pairs);
        }
    }
}

void compute_face_order(const Edge &MI_edge, const std::vector<PolygonFace> &MI_faces,
                        const std::vector<simplicial_arrangement::MaterialInterface<3>> &cut_results,
                        const std::vector<size_t> &cut_result_index,
                        const absl::flat_hash_map<size_t, std::vector<size_t>> &incident_tets,
                        std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>> &ordered_face_pairs) {
// collect all MI-faces incident to the MI-edge
    const auto &face_edge_indices = MI_edge.face_edge_indices;
// collect all tets containing the incident iso-faces
    std::set < size_t > containing_tets;
    for (const auto &face_edge: face_edge_indices) {
        const auto &tets = MI_faces[face_edge.first].tet_face_indices;
        for (const auto &t: tets) {
            containing_tets.insert(t.first);
        }
    }
// pair faces
    if (containing_tets.size() == 1) {
        auto tet_id = *containing_tets.begin();
        pair_faces_in_one_tet(
                cut_results[cut_result_index[tet_id]], MI_faces, MI_edge,
                ordered_face_pairs);
    } else {
        std::cout << containing_tets.size() << " containing tets." << std::endl;
        throw degenerate_iso_edge_exception;
    }
}


void pair_faces_in_one_tet(const Arrangement<3> &tet_cut_result,
                           const std::vector<PolygonFace> &iso_faces,
                           const Edge &iso_edge,
                           std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>> &ordered_face_pairs) {
    // find tet faces that are incident to the iso_edge
    std::vector<bool> is_incident_faces(tet_cut_result.faces.size(), false);
    // map: tet face id --> iso face id
    std::vector<size_t> iso_face_Id_of_face(tet_cut_result.faces.size(), Arrangement<3>::None);
    for (const auto &fId_eId_pair: iso_edge.face_edge_indices) {
        size_t iso_face_id = fId_eId_pair.first;
        size_t face_id = iso_faces[iso_face_id].tet_face_indices[0].second;
        is_incident_faces[face_id] = true;
        iso_face_Id_of_face[face_id] = iso_face_id;
    }

    // travel around the edge
    std::vector<bool> visited_cell(tet_cut_result.cells.size(), false);
    // start from the first incident face, and travel to its positive cell
    size_t iso_face_id1 = iso_edge.face_edge_indices[0].first;
    size_t face_id1 = iso_faces[iso_face_id1].tet_face_indices[0].second;
    int face_sign_1 = 1;
    size_t cell_id = tet_cut_result.faces[face_id1].positive_cell;
    //
    size_t iso_face_id2 = Arrangement<3>::None;
    size_t face_id2 = Arrangement<3>::None;
    int face_sign_2 = 0; // unknown
    while (cell_id != Arrangement<3>::None && !visited_cell[cell_id]) {
        visited_cell[cell_id] = true;
        // find next face
        for (const auto &fId: tet_cut_result.cells[cell_id].faces) {
            if (is_incident_faces[fId] && fId != face_id1) {
                face_id2 = fId;
            }
        }
        if (face_id2 == Arrangement<3>::None) {
            // next face is not found
            break;
        } else {
            // get sign of face2 and find next cell
            if (tet_cut_result.faces[face_id2].positive_cell == cell_id) {
                cell_id = tet_cut_result.faces[face_id2].negative_cell;
                face_sign_2 = 1;
            } else {
                cell_id = tet_cut_result.faces[face_id2].positive_cell;
                face_sign_2 = -1;
            }
            // add (face1, face2) to the list of face pairs
            iso_face_id2 = iso_face_Id_of_face[face_id2];
            ordered_face_pairs.emplace_back(std::make_pair(iso_face_id1, face_sign_1),
                                            std::make_pair(iso_face_id2, face_sign_2));
            // update face1 and clear face2
            face_id1 = face_id2;
            iso_face_id1 = iso_face_id2;
            face_sign_1 = -face_sign_2;
            face_id2 = Arrangement<3>::None;
            iso_face_id2 = Arrangement<3>::None;
            face_sign_2 = 0;
        }
    }
    // travel in a different direction
    iso_face_id1 = iso_edge.face_edge_indices[0].first;
    face_id1 = iso_faces[iso_face_id1].tet_face_indices[0].second;
    face_sign_1 = -1;
    cell_id = tet_cut_result.faces[face_id1].negative_cell;
    iso_face_id2 = Arrangement<3>::None;
    face_id2 = Arrangement<3>::None;
    face_sign_2 = 0;
    while (cell_id != Arrangement<3>::None && !visited_cell[cell_id]) {
        visited_cell[cell_id] = true;
        // find next face
        for (const auto &fId: tet_cut_result.cells[cell_id].faces) {
            if (is_incident_faces[fId] && fId != face_id1) {
                face_id2 = fId;
            }
        }
        if (face_id2 == Arrangement<3>::None) {
            // next face is not found
            break;
        } else {
            // get sign of face2 and find next cell
            if (tet_cut_result.faces[face_id2].positive_cell == cell_id) {
                cell_id = tet_cut_result.faces[face_id2].negative_cell;
                face_sign_2 = 1;
            } else {
                cell_id = tet_cut_result.faces[face_id2].positive_cell;
                face_sign_2 = -1;
            }
            // add (face1, face2) to the list of face pairs
            iso_face_id2 = iso_face_Id_of_face[face_id2];
            ordered_face_pairs.emplace_back(std::make_pair(iso_face_id1, face_sign_1),
                                            std::make_pair(iso_face_id2, face_sign_2));
            // update face1 and clear face2
            face_id1 = face_id2;
            iso_face_id1 = iso_face_id2;
            face_sign_1 = -face_sign_2;
            face_id2 = Arrangement<3>::None;
            iso_face_id2 = Arrangement<3>::None;
            face_sign_2 = 0;
        }
    }
}

void pair_faces_in_one_tet(const simplicial_arrangement::MaterialInterface<3> &tet_cut_result,
                           const std::vector<PolygonFace> &MI_faces,
                           const Edge &edge,
                           std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>> &ordered_face_pairs) {
    // find tet faces that are incident to the edge
    std::vector<bool> is_incident_faces(tet_cut_result.faces.size(), false);
    // map: tet face id --> MI-face id
    std::vector<size_t> MI_face_Id_of_face(tet_cut_result.faces.size(), MaterialInterface<3>::None);
    for (const auto &fId_eId_pair: edge.face_edge_indices) {
        size_t MI_face_id = fId_eId_pair.first;
        size_t face_id = MI_faces[MI_face_id].tet_face_indices[0].second;
        is_incident_faces[face_id] = true;
        MI_face_Id_of_face[face_id] = MI_face_id;
    }

    // map: material label -> cell id
    size_t max_material_label = 0;
    for (const auto &cell: tet_cut_result.cells) {
        if (cell.material_label > max_material_label) {
            max_material_label = cell.material_label;
        }
    }
    std::vector<size_t> cell_of_material(max_material_label + 1, MaterialInterface<3>::None);
    for (size_t i = 0; i < tet_cut_result.cells.size(); ++i) {
        cell_of_material[tet_cut_result.cells[i].material_label] = i;
    }

    // travel around the edge
    std::vector<bool> visited_cell(tet_cut_result.cells.size(), false);
    // start from the first incident face, and travel to its positive cell
    size_t MI_face_id1 = edge.face_edge_indices[0].first;
    size_t face_id1 = MI_faces[MI_face_id1].tet_face_indices[0].second;
    int face_sign_1 = 1;
    size_t cell_id = cell_of_material[tet_cut_result.faces[face_id1].positive_material_label];
    //
    size_t MI_face_id2 = MaterialInterface<3>::None;
    size_t face_id2 = MaterialInterface<3>::None;
    int face_sign_2 = 0; //unknown
    while (cell_id != MaterialInterface<3>::None && !visited_cell[cell_id]) {
        visited_cell[cell_id] = true;
        // find next face
        for (const auto &fId: tet_cut_result.cells[cell_id].faces) {
            if (is_incident_faces[fId] && fId != face_id1) {
                face_id2 = fId;
            }
        }
        if (face_id2 == MaterialInterface<3>::None) {
            // next face is not found
            break;
        } else {
            // get sign of face2 and find next cell
            if (cell_of_material[tet_cut_result.faces[face_id2].positive_material_label] == cell_id) {
                cell_id = cell_of_material[tet_cut_result.faces[face_id2].negative_material_label];
                face_sign_2 = 1;
            } else {
                cell_id = cell_of_material[tet_cut_result.faces[face_id2].positive_material_label];
                face_sign_2 = -1;
            }
            // add (face1, face2) to the list of face pairs
            MI_face_id2 = MI_face_Id_of_face[face_id2];
            ordered_face_pairs.emplace_back(std::make_pair(MI_face_id1, face_sign_1),
                                            std::make_pair(MI_face_id2, face_sign_2));
            // update face1 and clear face2
            face_id1 = face_id2;
            MI_face_id1 = MI_face_id2;
            face_sign_1 = -face_sign_2;
            face_id2 = MaterialInterface<3>::None;
            MI_face_id2 = MaterialInterface<3>::None;
            face_sign_2 = 0;
        }
    }
    // visit in a different direction
    MI_face_id1 = edge.face_edge_indices[0].first;
    face_id1 = MI_faces[MI_face_id1].tet_face_indices[0].second;
    face_sign_1 = -1;
    cell_id = cell_of_material[tet_cut_result.faces[face_id1].negative_material_label];
    MI_face_id2 = MaterialInterface<3>::None;
    face_id2 = MaterialInterface<3>::None;
    face_sign_2 = 0; //unknown
    while (cell_id != MaterialInterface<3>::None && !visited_cell[cell_id]) {
        visited_cell[cell_id] = true;
        // find next face
        for (const auto &fId: tet_cut_result.cells[cell_id].faces) {
            if (is_incident_faces[fId] && fId != face_id1) {
                face_id2 = fId;
            }
        }
        if (face_id2 == MaterialInterface<3>::None) {
            // next face is not found
            break;
        } else {
            // get sign of face2 and find next cell
            if (cell_of_material[tet_cut_result.faces[face_id2].positive_material_label] == cell_id) {
                cell_id = cell_of_material[tet_cut_result.faces[face_id2].negative_material_label];
                face_sign_2 = 1;
            } else {
                cell_id = cell_of_material[tet_cut_result.faces[face_id2].positive_material_label];
                face_sign_2 = -1;
            }
            // add (face1, face2) to the list of face pairs
            MI_face_id2 = MI_face_Id_of_face[face_id2];
            ordered_face_pairs.emplace_back(std::make_pair(MI_face_id1, face_sign_1),
                                            std::make_pair(MI_face_id2, face_sign_2));
            // update face1 and clear face2
            face_id1 = face_id2;
            MI_face_id1 = MI_face_id2;
            face_sign_1 = -face_sign_2;
            face_id2 = MaterialInterface<3>::None;
            MI_face_id2 = MaterialInterface<3>::None;
            face_sign_2 = 0;
        }
    }
}

void pair_faces_in_tets(const Edge &iso_edge, const std::vector<size_t> &containing_simplex,
                        const std::vector<size_t> &containing_tetIds, const std::vector<std::array<size_t, 4>> &tets,
                        const std::vector<PolygonFace> &iso_faces,
                        const std::vector<simplicial_arrangement::Arrangement<3>> &cut_results,
                        const std::vector<size_t> &cut_result_index,
                        const std::vector<size_t> &func_in_tet,
                        const std::vector<size_t> &start_index_of_tet,
                        std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>> &ordered_face_pairs) {
    //// pre-processing
    // collect all iso-faces incident to the iso-edge
    const auto &face_edge_indices = iso_edge.face_edge_indices;
    // map: (tet_id, tet_face_id) -> iso_face_id
    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> iso_face_Id_of_face;
    for (const auto &face_edge: face_edge_indices) {
        size_t iso_face_id = face_edge.first;
        const auto &tetrahedrons = iso_faces[iso_face_id].tet_face_indices;
        for (const auto &t: tetrahedrons) {
            iso_face_Id_of_face[t] = iso_face_id;
        }
    }
    // find identical tet boundary planes incident to iso-edge
    // map: (tet_Id, tet_plane_Id) -> (oppo_tet_Id, oppo_tet_plane_Id)
    absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> identical_tet_planes;
    if (containing_simplex.size() == 3) {
        // iso_edge contained in a tet boundary triangle
        size_t vId1 = containing_simplex[0];
        size_t vId2 = containing_simplex[1];
        size_t vId3 = containing_simplex[2];
        size_t tId1 = containing_tetIds[0];
        size_t tId2 = containing_tetIds[1];
        const auto &tet1 = tets[tId1];
        const auto &tet2 = tets[tId2];
        size_t pId1, pId2;
        size_t vId;
        for (size_t i = 0; i < 4; ++i) {
            vId = tet1[i];
            if (vId != vId1 && vId != vId2 && vId != vId3) {
                pId1 = i;
                break;
            }
        }
        for (size_t i = 0; i < 4; ++i) {
            vId = tet2[i];
            if (vId != vId1 && vId != vId2 && vId != vId3) {
                pId2 = i;
                break;
            }
        }
        // (tId1, pId1) and (tId2, pId2)
        identical_tet_planes[{tId1, pId1}] = {tId2, pId2};
        identical_tet_planes[{tId2, pId2}] = {tId1, pId1};
    } else {
        // iso_edge contained in a tet edge
        size_t vId1 = containing_simplex[0];
        size_t vId2 = containing_simplex[1];
        absl::flat_hash_map<std::array<size_t, 3>, std::pair<size_t, size_t>> tet_plane_of_tri;
        size_t vId;
        std::array<size_t, 3> tri;
        for (size_t tId: containing_tetIds) {
            for (size_t i = 0; i < 4; ++i) {
                vId = tets[tId][i];
                if (vId != vId1 && vId != vId2) {
                    tri = {tets[tId][(i + 1) % 4], tets[tId][(i + 2) % 4], tets[tId][(i + 3) % 4]};
                    std::sort(tri.begin(), tri.end());
                    auto iter_inserted = tet_plane_of_tri.try_emplace(tri, std::make_pair(tId, i));
                    if (!iter_inserted.second) {
                        identical_tet_planes[{tId, i}] = iter_inserted.first->second;
                        identical_tet_planes[iter_inserted.first->second] = {tId, i};
                    }
                }
            }
        }
    }
    // find identical faces (tet_id, tet_face_id) on tet boundary planes incident to iso-edge
    // map: (tet_Id, tet_face_Id) -> (oppo_tet_Id, oppo_tet_face_Id)
    absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> opposite_face;
    // map: (tet_Id, tet_vert_Id) -> boundary vert Id
    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> bndry_vert_id;
    size_t num_bndry_vert = 0;
    // hash table for vertices on the boundary of tetrahedron
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    absl::flat_hash_map<std::array<size_t, 3>, size_t> vert_on_tetEdge;
    absl::flat_hash_map<std::array<size_t, 5>, size_t> vert_on_tetFace;
    // hash table for faces on the boundary of tetrahedron
    // map: (i,j,k) -> (tet_Id, tet_face_Id)
    absl::flat_hash_map<std::array<size_t, 3>, std::pair<size_t, size_t>> face_on_tetFace;
    // auxiliary data
    std::vector<bool> is_bndry_vert;
    std::vector<bool> is_bndry_face;
    std::vector<size_t> bndry_vId_of_vert;
    std::vector<size_t> face_verts;
    std::array<size_t, 3> key3;
    std::array<size_t, 5> key5;
    std::array<bool, 4> used_pId;
    std::array<size_t, 2> vIds2;
    std::array<size_t, 3> vIds3;
    std::array<size_t, 3> implicit_pIds;
    std::array<size_t, 3> bndry_pIds;
    for (size_t i: containing_tetIds) {
        if (cut_result_index[i] == Arrangement<3>::None) {
            // empty tet i
            for (size_t fi = 0; fi < 4; ++i) {
                if (identical_tet_planes.find(std::make_pair(i, fi)) != identical_tet_planes.end()) {
                    // face fi lies on a tet boundary incident to iso-edge
                    std::array<size_t, 3> bndry_face_verts;
                    for (size_t j = 0; j < 3; ++j) {
                        auto key = tets[i][(fi + j + 1) % 4];
                        auto iter_inserted = vert_on_tetVert.try_emplace(key, num_bndry_vert);
                        if (iter_inserted.second) {
                            num_bndry_vert++;
                        }
                        bndry_face_verts[j] = iter_inserted.first->second;
                    }
                    std::sort(bndry_face_verts.begin(), bndry_face_verts.end());
                    auto iter_inserted = face_on_tetFace.try_emplace(bndry_face_verts, std::make_pair(i, fi));
                    if (!iter_inserted.second) {
                        // face inserted before
                        opposite_face[iter_inserted.first->second] = {i, fi};
                        opposite_face[{i, fi}] = iter_inserted.first->second;
                    }
                }
            }
        } else {
            // non-empty tet i
            const auto &arrangement = cut_results[cut_result_index[i]];
            const auto &vertices = arrangement.vertices;
            const auto &faces = arrangement.faces;
            auto start_index = start_index_of_tet[i];
            auto num_func = start_index_of_tet[i + 1] - start_index;
            // find vertices and faces on tet boundary incident to iso-edge
            is_bndry_vert.resize(vertices.size());
            std::fill(is_bndry_vert.begin(), is_bndry_vert.end(), false);
            is_bndry_face.clear();
            for (const auto &face: faces) {
                is_bndry_face.push_back(false);
                if (face.supporting_plane < 4 &&
                    identical_tet_planes.find({i, face.supporting_plane}) != identical_tet_planes.end()) {
                    is_bndry_face.back() = true;
                    for (const auto &vid: face.vertices) {
                        is_bndry_vert[vid] = true;
                    }
                }
            }
            // map: local vert index --> boundary vert index
            bndry_vId_of_vert.clear();
            // create boundary vertices
            for (size_t j = 0; j < vertices.size(); j++) {
                bndry_vId_of_vert.push_back(Arrangement<3>::None);
                if (is_bndry_vert[j]) {
                    size_t num_bndry_planes = 0;
                    size_t num_impl_planes = 0;
                    const auto &vertex = vertices[j];
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
                            auto iter_inserted = vert_on_tetEdge.try_emplace(key3, num_bndry_vert);
                            if (iter_inserted.second) {
                                num_bndry_vert++;
                            }
                            bndry_vId_of_vert.back() = iter_inserted.first->second;
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
                            auto iter_inserted = vert_on_tetFace.try_emplace(key5, num_bndry_vert);
                            if (iter_inserted.second) {
                                num_bndry_vert++;
                            }
                            bndry_vId_of_vert.back() = iter_inserted.first->second;
                            break;
                        }
                        case 3: // on tet vertex
                        {
                            used_pId[0] = false;
                            used_pId[1] = false;
                            used_pId[2] = false;
                            used_pId[3] = false;
                            for (const auto &pId: bndry_pIds) {
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
                            auto iter_inserted = vert_on_tetVert.try_emplace(key, num_bndry_vert);
                            if (iter_inserted.second) {
                                num_bndry_vert++;
                            }
                            bndry_vId_of_vert.back() = iter_inserted.first->second;
                            break;
                        }
                        default:
                            break;
                    }
                }
            }
            // pair boundary faces
            for (size_t j = 0; j < faces.size(); j++) {
                if (is_bndry_face[j]) {
                    face_verts.clear();
                    for (auto vId: faces[j].vertices) {
                        face_verts.push_back(bndry_vId_of_vert[vId]);
                    }
                    compute_iso_face_key(face_verts, key3);
                    auto iter_inserted = face_on_tetFace.try_emplace(key3, std::make_pair(i, j));
                    if (!iter_inserted.second) {
                        // face inserted before, pair the two faces
                        opposite_face[iter_inserted.first->second] = {i, j};
                        opposite_face[{i, j}] = iter_inserted.first->second;
                    }
                }
            }
        }
    }

    // find the half-iso-face of the local face in tet with given orientation
    // the orientation of an iso-face is defined by the smallest-index implicit function passing the iso-face
    auto get_half_iso_face = [&iso_face_Id_of_face, &cut_results, &cut_result_index]
            (std::pair<size_t, size_t> tet_face, int orient, size_t &iso_face_id, int &iso_orient) {
        iso_face_id = iso_face_Id_of_face[tet_face];
        const auto &cell_complex = cut_results[cut_result_index[tet_face.first]];
        const auto &faces = cell_complex.faces;
        auto supp_pId = faces[tet_face.second].supporting_plane;
        if (supp_pId > 3) { // plane 0,1,2,3 are tet boundary planes
            // supporting plane is not a tet boundary plane
            // assume: supporting plane is the iso-surface of smallest-index implicit function
            iso_orient = orient;
        } else {
            // supporting plane is a tet boundary plane, must be duplicate planes
            auto uid = cell_complex.unique_plane_indices[supp_pId];
            // find the smallest-index non-boundary plane
            size_t min_pId = std::numeric_limits<size_t>::max();
            for (auto pId: cell_complex.unique_planes[uid]) {
                if (pId > 3 && pId < min_pId) {
                    min_pId = pId;
                }
            }
            // orient is the orientation of supporting plane
            // flip orientation if smallest-index non-boundary plane has different orientation
            if (cell_complex.unique_plane_orientations[min_pId] != cell_complex.unique_plane_orientations[supp_pId]) {
                iso_orient = -orient;
            } else {
                iso_orient = orient;
            }
        }
    };

    // pair (tet_id, tet_face_id)
    auto find_next =
            [&cut_results, &cut_result_index, &opposite_face, &iso_face_Id_of_face]
                    (std::pair<size_t, size_t> face, int orient,
                     std::pair<size_t, size_t> &face_next, int &orient_next, auto &&find_next) {
                size_t tet_id = face.first;
                size_t tet_face_id = face.second;
                if (cut_result_index[tet_id] == Arrangement<3>::None) {
                    // empty tet
                    if (orient == 1) { // Positive side: the tet
                        for (size_t fi = 0; fi < 4; ++fi) {
                            if (fi != tet_face_id && opposite_face.find({tet_id, fi}) != opposite_face.end()) {
                                face_next = {tet_id, fi};
                                orient_next = 1;
                                return;
                            }
                        }
                    } else { // Negative side: None
                        find_next(opposite_face[face], 1, face_next, orient_next, find_next);
                    }
                } else {
                    // non-empty tet
                    const auto &cell_complex = cut_results[cut_result_index[tet_id]];
                    size_t cell_id = (orient == 1 ? cell_complex.faces[tet_face_id].positive_cell :
                                      cell_complex.faces[tet_face_id].negative_cell);
                    if (cell_id != Arrangement<3>::None) {
                        for (auto fi: cell_complex.cells[cell_id].faces) {
                            if (fi != tet_face_id && (
                                    iso_face_Id_of_face.find({tet_id, fi}) != iso_face_Id_of_face.end() ||
                                    opposite_face.find({tet_id, fi}) != opposite_face.end())) {
                                face_next = {tet_id, fi};
                                break;
                            }
                        }
                        orient_next = cell_complex.faces[face_next.second].positive_cell == cell_id ? 1 : -1;
                    } else {
                        // cell is None, so the face lies on tet boundary
                        find_next(opposite_face[face], 1, face_next, orient_next, find_next);
                    }
                }
            };


    //// face pairing algorithm
    size_t iso_face_0 = face_edge_indices[0].first;
    // pair (tet_id, tet_face_id)
    std::pair<size_t, size_t> face_curr = iso_faces[iso_face_0].tet_face_indices[0];
    int orient_curr = 1; // orientation: 1 for Positive, -1 for Negative
    std::pair<size_t, size_t> face_next;
    int orient_next, iso_orient_curr, iso_orient_next;
    size_t iso_face_curr = iso_face_0;
    size_t iso_face_next = Arrangement<3>::None;
    while (iso_face_next != iso_face_0) {
        find_next(face_curr, orient_curr, face_next, orient_next, find_next);
        while (iso_face_Id_of_face.find(face_next) == iso_face_Id_of_face.end()) {
            find_next(face_next, -orient_next, face_next, orient_next, find_next);
        }
        get_half_iso_face(face_curr, orient_curr, iso_face_curr, iso_orient_curr);
        get_half_iso_face(face_next, orient_next, iso_face_next, iso_orient_next);
        ordered_face_pairs.emplace_back(std::make_pair(iso_face_curr, iso_orient_curr),
                                        std::make_pair(iso_face_next, iso_orient_next));
        face_curr = face_next;
        orient_curr = -orient_next;
    }
}
