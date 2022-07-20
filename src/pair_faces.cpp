//
// Created by Charles Du on 7/20/22.
//

#include "pair_faces.h"

using namespace simplicial_arrangement;

class Degenerate_Iso_Edge_Exception : public std::exception
{
    virtual const char* what() const throw()
    {
        return "More than one tets containing the iso-faces incident to the iso-edge.";
    }
} degenerate_iso_edge_exception;

void compute_face_order(const Edge& iso_edge, const std::vector<PolygonFace>& iso_faces,
                        const std::vector<Arrangement<3>>& cut_results,
const std::vector<size_t>& cut_result_index,
const absl::flat_hash_map<size_t, std::vector<size_t>>& incident_tets,
//    std::vector<std::pair<size_t, int>>& ordered_faces)
std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>>& ordered_face_pairs)
{
// collect all iso-faces incident to the iso-edge
const auto& face_edge_indices = iso_edge.face_edge_indices;
// collect all tets containing the incident iso-faces
std::set<size_t> containing_tets;
for (const auto& face_edge : face_edge_indices) {
const auto& tets = iso_faces[face_edge.first].tet_face_indices;
for (const auto& t : tets) {
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


void compute_face_order_in_one_tet(const Arrangement<3>& tet_cut_result,
                                   const std::vector<PolygonFace>& iso_faces,
                                   const Edge& iso_edge,
                                   std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>>& ordered_face_pairs)
{
    // in a single tetrahedron, num face pairs == num iso-faces incident to the iso-edge
    size_t num_incident_faces = iso_edge.face_edge_indices.size();
    //
    std::vector<bool> is_incident_faces(tet_cut_result.faces.size(), false);
    // map: tet face id --> iso face id
    std::vector<size_t> iso_face_Id_of_face(tet_cut_result.faces.size(), Arrangement<3>::None);
    for (const auto& fId_eId_pair : iso_edge.face_edge_indices) {
        size_t iso_face_id = fId_eId_pair.first;
        size_t face_id = iso_faces[iso_face_id].tet_face_indices[0].second;
        is_incident_faces[face_id] = true;
        iso_face_Id_of_face[face_id] = iso_face_id;
    }

    // travel around the edge
    std::vector<bool> visited_cell(tet_cut_result.cells.size(), false);
    //
    size_t iso_face_id1 = iso_edge.face_edge_indices[0].first;
    size_t face_id1 = iso_faces[iso_face_id1].tet_face_indices[0].second;
    int face_sign_1 = 1;
    size_t cell_id = tet_cut_result.faces[face_id1].positive_cell;
    size_t iso_face_id2 = Arrangement<3>::None;
    size_t face_id2 = Arrangement<3>::None;
    int face_sign_2 = 0;
    while(cell_id != Arrangement<3>::None && !visited_cell[cell_id]) {
        //ordered_faces[i] = std::make_pair(iso_face_id1, face_sign);
        visited_cell[cell_id] = true;
        // find next face
        for (const auto& fId : tet_cut_result.cells[cell_id].faces) {
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
            ordered_face_pairs.push_back(std::make_pair(std::make_pair(iso_face_id1,face_sign_1),
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
    // visit in a different direction
    iso_face_id1 = iso_edge.face_edge_indices[0].first;
    face_id1 = iso_faces[iso_face_id1].tet_face_indices[0].second;
    face_sign_1 = -1;
    cell_id = tet_cut_result.faces[face_id1].negative_cell;
    iso_face_id2 = Arrangement<3>::None;
    face_id2 = Arrangement<3>::None;
    face_sign_2 = 0;
    while(cell_id != Arrangement<3>::None && !visited_cell[cell_id]) {
        //ordered_faces[i] = std::make_pair(iso_face_id1, face_sign);
        visited_cell[cell_id] = true;
        // find next face
        for (const auto& fId : tet_cut_result.cells[cell_id].faces) {
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
            ordered_face_pairs.push_back(std::make_pair(std::make_pair(iso_face_id1,face_sign_1),
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