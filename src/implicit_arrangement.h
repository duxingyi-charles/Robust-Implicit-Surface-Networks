//
// Created by Charles Du on 8/18/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_IMPLICIT_ARRANGEMENT_H
#define ROBUST_IMPLICIT_NETWORKS_IMPLICIT_ARRANGEMENT_H

#include "robust_implicit_networks.h"

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
        std::vector<Edge>& iso_edges,
        std::vector<std::vector<size_t>>& chains,
        std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
        std::vector<std::vector<size_t>>& shells,
        std::vector<std::vector<size_t>>& arrangement_cells,
        std::vector<std::string>& timing_labels,
        std::vector<double>& timings,
        std::vector<std::string>& stats_labels,
        std::vector<size_t>& stats);

#endif //ROBUST_IMPLICIT_NETWORKS_IMPLICIT_ARRANGEMENT_H
