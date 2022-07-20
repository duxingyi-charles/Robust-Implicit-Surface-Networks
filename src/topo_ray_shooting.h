//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_TOPO_RAY_SHOOTING_H
#define ROBUST_IMPLICIT_NETWORKS_TOPO_RAY_SHOOTING_H

#include <vector>
#include <simplicial_arrangement/simplicial_arrangement.h>
#include <absl/container/flat_hash_map.h>
#include "mesh.h"


void topo_ray_shooting(
        const std::vector<std::array<double, 3>> &pts,
        const std::vector<std::array<size_t, 4>> &tets,
        const std::vector<simplicial_arrangement::Arrangement<3>> &cut_results,
        const std::vector<size_t> &cut_result_index,
        const std::vector<IsoVert> &iso_verts,
        const std::vector<PolygonFace> &iso_faces,
        const std::vector<std::vector<size_t>> &patches,
        const std::vector<size_t> &patch_of_face,
        const std::vector<std::vector<size_t>> &shells,
        const std::vector<size_t> &shell_of_half_patch,
        const std::vector<std::vector<size_t>> &components,
        const std::vector<size_t> &component_of_patch,
        std::vector<std::vector<size_t>> &arrangement_cells
);

// Given tet mesh,
// build the map: v-->v_next, where v_next has lower order than v
void build_next_vert(const std::vector<std::array<double,3>> &pts,
                     const std::vector<std::array<size_t,4>> &tets,
                     std::vector<size_t> &next_vert);

// find extremal edge for each component
void find_extremal_edges(
        const std::vector<std::array<double, 3>> &pts,
        const std::vector<IsoVert> &iso_verts,
        const std::vector<PolygonFace> &iso_faces,
        const std::vector<std::vector<size_t>> &patches,
        const std::vector<std::vector<size_t>> &components,
        const std::vector<size_t> &component_of_patch,
        const std::vector<size_t> &next_vert,
        // extremal edge of component i is stored at position [2*i], [2*i+1]
        std::vector<size_t> &extremal_edge_of_component,
        // store an iso-vert index on edge (v, v_next), None means there is no such iso-vert
        std::vector<size_t> &iso_vert_on_v_v_next,
        // map: (tet_id, tet_face_id) --> iso_face_id
        absl::flat_hash_map<std::pair<size_t, size_t>, size_t> &iso_face_id_of_tet_face,
        // map: (tet_id, tet_vert_id) --> (iso_vert_id, component_id)
        absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
            &iso_vId_compId_of_tet_vert
        );

// compute the order of iso-vertices on a tet edge v->u, v,u in {0,1,2,3}
// return a list of sorted vertex indices {v_id, i1, i2, ..., u_id}
void compute_edge_intersection_order(const simplicial_arrangement::Arrangement<3>& tet_cut_result, size_t v, size_t u,
                                     std::vector<size_t> &vert_indices);

// find the two faces passing v1 and v2, v1->v2 is part of a tet edge
void compute_passing_face_pair(const simplicial_arrangement::Arrangement<3>& tet_cut_result,
                               size_t v1, size_t v2,
                               std::pair<size_t,int> &face_orient1,
                               std::pair<size_t,int> &face_orient2);

// find the face passing v, v->u is part of a tet edge, and u is a tet vertex
void compute_passing_face(const simplicial_arrangement::Arrangement<3>& tet_cut_result,
                          size_t v, size_t u,
                          std::pair<size_t,int> &face_orient);

// point (x,y,z): dictionary order
inline bool point_xyz_less(const std::array<double, 3>& p, const std::array<double, 3>& q)
{
    if (p[0] == q[0]) {
        if (p[1] == q[1]) {
            return p[2] < q[2];
        } else {
            return p[1] < q[1];
        }
    } else {
        return p[0] < q[0];
    }
}

#endif //ROBUST_IMPLICIT_NETWORKS_TOPO_RAY_SHOOTING_H

