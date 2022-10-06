//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_PAIR_FACES_H
#define ROBUST_IMPLICIT_NETWORKS_PAIR_FACES_H

#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/material_interface.h>
#include <absl/container/flat_hash_map.h>

#include "mesh.h"

// compute neighboring pair of half-faces around an iso-edge
// output:
// pair<size_t, int> : pair (iso-face index, iso-face orientation)
void compute_face_order(const Edge &iso_edge, const std::vector<std::array<size_t, 4>> &tets,
                        const std::vector<IsoVert> &iso_verts,
                        const std::vector<PolygonFace> &iso_faces,
                        const std::vector<simplicial_arrangement::Arrangement<3>> &cut_results,
                        const std::vector<size_t> &cut_result_index,
                        const std::vector<size_t>& func_in_tet,
                        const std::vector<size_t>& start_index_of_tet,
                        const absl::flat_hash_map<size_t, std::vector<size_t>> &incident_tets,
                        std::vector<std::pair<std::pair<size_t, int>,
                        std::pair<size_t, int>>> &ordered_face_pairs);

// compute neighboring pair of MI-faces around an MI-edge
// output:
// pair<size_t, int> : pair (MI-face index, MI-face orientation)
void compute_face_order(const Edge& MI_edge, const std::vector<PolygonFace>& MI_faces,
                        const std::vector<simplicial_arrangement::MaterialInterface<3>>& cut_results,
                        const std::vector<size_t>& cut_result_index,
                        const absl::flat_hash_map<size_t, std::vector<size_t>>& incident_tets,
                        std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>>& ordered_face_pairs);


// compute neighboring pair of half-faces around an iso-edge in a tetrahedron
// pair<size_t, int> : pair (iso-face index, iso-face orientation)
void pair_faces_in_one_tet(const simplicial_arrangement::Arrangement<3>& tet_cut_result,
                           const std::vector<PolygonFace>& iso_faces,
                           const Edge& iso_edge,
                           std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>>& ordered_face_pairs);

// compute neighboring pair of half-faces around an iso-edge in multiple tetrahedrons
// pair<size_t, int> : pair (iso-face index, iso-face orientation)
void pair_faces_in_tets(const Edge& iso_edge, const std::vector<size_t> &containing_simplex,
                        const std::vector<size_t> &containing_tetIds,
                        const std::vector<std::array<size_t, 4>> &tets,
                        const std::vector<PolygonFace> &iso_faces,
                        const std::vector<simplicial_arrangement::Arrangement<3>> &cut_results,
                        const std::vector<size_t> &cut_result_index,
                        const std::vector<size_t>& func_in_tet,
                        const std::vector<size_t>& start_index_of_tet,
                        std::vector<std::pair<std::pair<size_t, int>,
                                std::pair<size_t, int>>> &ordered_face_pairs);

// compute neighboring pair of half-faces around an MI-edge in a tetrahedron
// pair<size_t, int> : pair (MI-face index, MI-face orientation)
void pair_faces_in_one_tet(const simplicial_arrangement::MaterialInterface<3>& tet_cut_result,
                           const std::vector<PolygonFace>& MI_faces,
                           const Edge& edge,
                           std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>>& ordered_face_pairs);

#endif //ROBUST_IMPLICIT_NETWORKS_PAIR_FACES_H
