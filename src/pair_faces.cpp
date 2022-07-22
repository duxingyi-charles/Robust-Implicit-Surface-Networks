//
// Created by Charles Du on 7/20/22.
//

#include "pair_faces.h"

using namespace simplicial_arrangement;

class Degenerate_Iso_Edge_Exception : public std::exception {
    virtual const char *what() const throw() {
        return "More than one tets containing the iso-faces incident to the iso-edge.";
    }
} degenerate_iso_edge_exception;

void compute_face_order(const Edge &iso_edge, const std::vector<PolygonFace> &iso_faces,
                        const std::vector<Arrangement<3>> &cut_results,
                        const std::vector<size_t> &cut_result_index,
                        const absl::flat_hash_map<size_t, std::vector<size_t>> &incident_tets,
                        std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>> &ordered_face_pairs) {
// collect all iso-faces incident to the iso-edge
    const auto &face_edge_indices = iso_edge.face_edge_indices;
// collect all tets containing the incident iso-faces
    std::set<size_t> containing_tets;
    for (const auto &face_edge: face_edge_indices) {
        const auto &tets = iso_faces[face_edge.first].tet_face_indices;
        for (const auto &t: tets) {
            containing_tets.insert(t.first);
        }
    }
// pair faces
    if (containing_tets.size() == 1) {
        auto tet_id = *containing_tets.begin();
        compute_face_order_in_one_tet(
                cut_results[cut_result_index[tet_id]], iso_faces, iso_edge,
//            ordered_faces);
                ordered_face_pairs);
    } else {
        std::cout << containing_tets.size() << " containing tets." << std::endl;
        throw degenerate_iso_edge_exception;
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
    std::set<size_t> containing_tets;
    for (const auto &face_edge: face_edge_indices) {
        const auto &tets = MI_faces[face_edge.first].tet_face_indices;
        for (const auto &t: tets) {
            containing_tets.insert(t.first);
        }
    }
// pair faces
    if (containing_tets.size() == 1) {
        auto tet_id = *containing_tets.begin();
        compute_face_order_in_one_tet(
                cut_results[cut_result_index[tet_id]], MI_faces, MI_edge,
                ordered_face_pairs);
    } else {
        std::cout << containing_tets.size() << " containing tets." << std::endl;
        throw degenerate_iso_edge_exception;
    }
}


void compute_face_order_in_one_tet(const Arrangement<3> &tet_cut_result,
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
            ordered_face_pairs.push_back(std::make_pair(std::make_pair(iso_face_id1, face_sign_1),
                                                        std::make_pair(iso_face_id2, face_sign_2)));
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
            ordered_face_pairs.push_back(std::make_pair(std::make_pair(iso_face_id1, face_sign_1),
                                                        std::make_pair(iso_face_id2, face_sign_2)));
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

void compute_face_order_in_one_tet(const simplicial_arrangement::MaterialInterface<3>& tet_cut_result,
                                   const std::vector<PolygonFace>& MI_faces,
                                   const Edge& edge,
                                   std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>>& ordered_face_pairs) {
    // find tet faces that are incident to the edge
    std::vector<bool> is_incident_faces(tet_cut_result.faces.size(), false);
    // map: tet face id --> MI-face id
    std::vector<size_t> MI_face_Id_of_face(tet_cut_result.faces.size(), MaterialInterface<3>::None);
    for (const auto& fId_eId_pair : edge.face_edge_indices) {
        size_t MI_face_id = fId_eId_pair.first;
        size_t face_id = MI_faces[MI_face_id].tet_face_indices[0].second;
        is_incident_faces[face_id] = true;
        MI_face_Id_of_face[face_id] = MI_face_id;
    }

    // map: material label -> cell id
    size_t max_material_label = 0;
    for (const auto& cell : tet_cut_result.cells) {
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
    while(cell_id != MaterialInterface<3>::None && !visited_cell[cell_id]) {
        visited_cell[cell_id] = true;
        // find next face
        for (const auto& fId : tet_cut_result.cells[cell_id].faces) {
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
            ordered_face_pairs.push_back(std::make_pair(std::make_pair(MI_face_id1, face_sign_1),
                                                        std::make_pair(MI_face_id2, face_sign_2)));
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
    while(cell_id != MaterialInterface<3>::None && !visited_cell[cell_id]) {
        visited_cell[cell_id] = true;
        // find next face
        for (const auto& fId : tet_cut_result.cells[cell_id].faces) {
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
            ordered_face_pairs.push_back(std::make_pair(std::make_pair(MI_face_id1, face_sign_1),
                                                        std::make_pair(MI_face_id2, face_sign_2)));
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