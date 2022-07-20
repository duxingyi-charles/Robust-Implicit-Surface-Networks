//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_EXTRACT_MESH_H
#define ROBUST_IMPLICIT_NETWORKS_EXTRACT_MESH_H

#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/material_interface.h>
#include <Eigen/Core>
#include "mesh.h"

// extract iso-mesh (topology only)
void extract_iso_mesh(
        size_t num_1_func, size_t num_2_func, size_t num_more_func,
        const std::vector<simplicial_arrangement::Arrangement<3>>& cut_results,
        const std::vector<size_t>& cut_result_index,
        const std::vector<size_t>& func_in_tet,
        const std::vector<size_t>& start_index_of_tet,
        const std::vector<std::array<size_t, 4>>& tets,
        std::vector<IsoVert>& iso_verts,
        std::vector<PolygonFace>& iso_faces);

// extract iso-mesh (topology only) and create map: local index --> global index
void extract_iso_mesh(
        size_t num_1_func, size_t num_2_func, size_t num_more_func,
        const std::vector<simplicial_arrangement::Arrangement<3>>& cut_results,
        const std::vector<size_t>& cut_result_index,
        const std::vector<size_t>& func_in_tet,
        const std::vector<size_t>& start_index_of_tet,
        const std::vector<std::array<size_t, 4>>& tets,
        std::vector<IsoVert>& iso_verts,
        std::vector<PolygonFace>& iso_faces,
        std::vector<long long>& global_vId_of_tet_vert,
        std::vector<size_t>& global_vId_start_index_of_tet,
        std::vector<size_t>& iso_fId_of_tet_face,
        std::vector<size_t>& iso_fId_start_index_of_tet);

// extract material interface mesh (topology only)
void extract_MI_mesh(
        size_t num_2_func, size_t num_3_func, size_t num_more_func,
        const std::vector<simplicial_arrangement::MaterialInterface<3>>& cut_results,
        const std::vector<size_t>& cut_result_index,
        const std::vector<size_t>& material_in_tet,
        const std::vector<size_t>& start_index_of_tet,
        const std::vector<std::array<size_t, 4>>& tets,
        std::vector<MI_Vert>& MI_verts,
        std::vector<PolygonFace>& MI_faces);

// extract material interface mesh (topology only) and create map: local index --> global index
void extract_MI_mesh(
        size_t num_2_func, size_t num_3_func, size_t num_more_func,
        const std::vector<simplicial_arrangement::MaterialInterface<3>>& cut_results,
        const std::vector<size_t>& cut_result_index,
        const std::vector<size_t>& material_in_tet,
        const std::vector<size_t>& start_index_of_tet,
        const std::vector<std::array<size_t, 4>>& tets,
        std::vector<MI_Vert>& MI_verts,
        std::vector<PolygonFace>& MI_faces,
        std::vector<long long>& global_vId_of_tet_vert,
        std::vector<size_t>& global_vId_start_index_of_tet,
        std::vector<size_t>& MI_fId_of_tet_face,
        std::vector<size_t>& MI_fId_start_index_of_tet);

// given the list of vertex indices of a face, return the unique key of the face: (the smallest vert Id,
// second-smallest vert Id, the largest vert Id) assume: face_verts is a list of non-duplicate natural
// numbers, with at least three elements.
template <typename IndexType>
void compute_iso_face_key(const std::vector<IndexType>& face_verts, std::array<IndexType, 3>& key)
{
    IndexType min_vert = face_verts[0];
    size_t min_pos = 0;
    IndexType max_vert = face_verts[0];
    for (size_t i = 1; i < face_verts.size(); i++) {
        if (face_verts[i] < min_vert) {
            min_vert = face_verts[i];
            min_pos = i;
        } else if (face_verts[i] > max_vert) {
            max_vert = face_verts[i];
        }
    }
    IndexType second_min_vert = max_vert + 1;
    for (size_t i = 0; i < face_verts.size(); i++) {
        if (i != min_pos && face_verts[i] < second_min_vert) {
            second_min_vert = face_verts[i];
        }
    }
    //
    key[0] = min_vert;
    key[1] = second_min_vert;
    key[2] = max_vert;
}

// compute xyz coordinates of iso-vertices
void compute_iso_vert_xyz(
        const std::vector<IsoVert> &iso_verts,
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcVals,
        const std::vector<std::array<double, 3>> &pts,
        std::vector<std::array<double, 3>>& iso_pts);


// compute xyz coordinates of material interface vertices
void compute_MI_vert_xyz(
        const std::vector<MI_Vert> &MI_verts,
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcVals,
        const std::vector<std::array<double,3>> &pts,
        std::vector<std::array<double,3>>& MI_pts);

// compute barycentric coordinate of Point (intersection of three planes)
// Point in tet cell
template <typename Scalar>
inline void compute_barycentric_coords(const std::array<Scalar, 4>& plane1,
                                       const std::array<Scalar, 4>& plane2,
                                       const std::array<Scalar, 4>& plane3,
                                       std::array<Scalar, 4>& bary_coords)
{
    Scalar n1 = plane1[3] * (plane2[2] * plane3[1] - plane2[1] * plane3[2]) +
                plane1[2] * (plane2[1] * plane3[3] - plane2[3] * plane3[1]) +
                plane1[1] * (plane2[3] * plane3[2] - plane2[2] * plane3[3]);
    Scalar n2 = plane1[3] * (plane2[0] * plane3[2] - plane2[2] * plane3[0]) +
                plane1[2] * (plane2[3] * plane3[0] - plane2[0] * plane3[3]) +
                plane1[0] * (plane2[2] * plane3[3] - plane2[3] * plane3[2]);
    Scalar n3 = plane1[3] * (plane2[1] * plane3[0] - plane2[0] * plane3[1]) +
                plane1[1] * (plane2[0] * plane3[3] - plane2[3] * plane3[0]) +
                plane1[0] * (plane2[3] * plane3[1] - plane2[1] * plane3[3]);
    Scalar n4 = plane1[2] * (plane2[0] * plane3[1] - plane2[1] * plane3[0]) +
                plane1[1] * (plane2[2] * plane3[0] - plane2[0] * plane3[2]) +
                plane1[0] * (plane2[1] * plane3[2] - plane2[2] * plane3[1]);
    Scalar d = n1 + n2 + n3 + n4;
    //
    bary_coords[0] = n1 / d;
    bary_coords[1] = n2 / d;
    bary_coords[2] = n3 / d;
    bary_coords[3] = n4 / d;
}

// Point on tet face
template <typename Scalar>
inline void compute_barycentric_coords(const std::array<Scalar, 3>& plane1,
                                       const std::array<Scalar, 3>& plane2,
                                       std::array<Scalar, 3>& bary_coords)
{
    Scalar n1 = plane1[2] * plane2[1] - plane1[1] * plane2[2];
    Scalar n2 = plane1[0] * plane2[2] - plane1[2] * plane2[0];
    Scalar n3 = plane1[1] * plane2[0] - plane1[0] * plane2[1];
    Scalar d = n1 + n2 + n3;
    //
    bary_coords[0] = n1 / d;
    bary_coords[1] = n2 / d;
    bary_coords[2] = n3 / d;
}

// Point on tet edge
template <typename Scalar>
inline void compute_barycentric_coords(Scalar f1, Scalar f2, std::array<Scalar, 2>& bary_coords)
{
    bary_coords[0] = f2 / (f2 - f1);
    bary_coords[1] = 1 - bary_coords[0];
}


#endif //ROBUST_IMPLICIT_NETWORKS_EXTRACT_MESH_H
