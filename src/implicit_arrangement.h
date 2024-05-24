//
// Created by Charles Du on 8/18/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_IMPLICIT_ARRANGEMENT_H
#define ROBUST_IMPLICIT_NETWORKS_IMPLICIT_ARRANGEMENT_H

#include "robust_implicit_networks.h"
///
/// The body of computing the implicit arrangement; used in `src/implicit_arrangement.cpp`
///
///  @param[in] robust_test         Toggle robust test: run twice and see if both results are consistent
///  @param[in] use_lookup           Toggle to use a look-up table: Use look-up table to accelerate
///  @param[in] use_secondary_lookup            Toggle to use look-up tables for tetrahedral with two active functions
///  @param[in] use_topo_ray_shooting           Toggle to use topological ray shooting to compute the spatial decomposition induced by the arrangement
///  @param[in] pts         Vertices of the background grid
///  @param[in] tets            Tetetrahedra corners' vertex indices
///  @param[in] funcVals            2D matrix of function values at all vertices
///
///
///  @param[out] iso_pts            Vertices at the surface network mesh
///  @param[out] iso_faces          Polygonal faces at the surface network mesh
///  @param[out] patches            A connected component of faces bounded by non-manifold edges
///  @param[out] patch_function_label           index of the zero-valued function that creates the patch
///  @param[out] iso_edges          Edges at the surface network mesh
///  @param[out] chains         Chains of non-manifold edges
///  @param[out] non_manifold_edges_of_vert         Indices of non-manifold vertices
///  @param[out] shells         A connected component of the boundary partitioned by the surface network; encoded by patch indices: positive patch i -> 2i, negative patch i->2i+1
///  @param[out] arrangement_cells          A 3D region partitioned by the surface network; encoded by a vector of shell indices
///  @param[out] cell_function_label            A 2D boolean array for the signs of each pair of a function and an arrangement cell
///  @param[out] timing_labels          Labels for timing
///  @param[out] timings            Timing results
///  @param[out] stats_labels           Labels for geometry metrics
///  @param[out] stats          Geometry metrics results
///
///
///  @see           `PolygonFace` and `Edge` in `mesh.h`

bool implicit_arrangement(
        bool robust_test,
        bool use_lookup,
        bool use_secondary_lookup,
        bool use_topo_ray_shooting,
        //
        const std::vector<std::array<double, 3>>& pts,
        const std::vector<std::array<size_t, 4>>& tets,
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& funcVals,
        //
        std::vector<std::array<double, 3>>& iso_pts,
        std::vector<PolygonFace>& iso_faces,
        std::vector<std::vector<size_t>>& patches,
        std::vector<size_t>& patch_function_label,
        std::vector<Edge>& iso_edges,
        std::vector<std::vector<size_t>>& chains,
        std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
        std::vector<std::vector<size_t>>& shells,
        std::vector<std::vector<size_t>>& arrangement_cells,
        std::vector<std::vector<size_t>>& cell_function_label,
        std::vector<std::string>& timing_labels,
        std::vector<double>& timings,
        std::vector<std::string>& stats_labels,
        std::vector<size_t>& stats);

#endif //ROBUST_IMPLICIT_NETWORKS_IMPLICIT_ARRANGEMENT_H
