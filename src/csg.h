//
//  csg.h
//  robust_implicit_networks
//
//  Created by Yiwen Ju on 6/3/24.
//

#pragma once
#include <string>
#include <Eigen/Core>

#include "implicit_arrangement.h"
#include "robust_implicit_networks.h"


using namespace simplicial_arrangement;

///
/// The pruning of `patches`, `chains`, `non_manifold_edges_of_vert`
///
///
/// @param[in] robust_test         Toggle robust test: run twice and see if both results are consistent
/// @param[in] use_lookup           Toggle to use a look-up table: Use look-up table to accelerate
/// @param[in] use_secondary_lookup            Toggle to use look-up tables for tetrahedral with two active functions
/// @param[in] use_topo_ray_shooting           Toggle to use topological ray shooting to compute the spatial decomposition induced by the arrangement
/// @param[in] positive_inside          Defines the positive side of the shape. `true` means positive inside, and `false` means negative inside.
/// @param[in] pts         Vertices of the background grid
/// @param[in] tets            Tetetrahedra corners' vertex indices
/// @param[in] funcVals            2D matrix of function values at all vertices
/// @param[in] lambda           The indicator function takes a 2D vector, which each row represents a list of functions signs (one per implicit function) within a cell, and outputs a 1D vector cell_label, which indicates each cell is within or outside of the shape.
///
///
/// @param[out] mesh_pts            Vertices at the surface network mesh
/// @param[out] mesh_faces          Polygonal faces at the surface network mesh
/// @param[out] patch_function_label           a pair of two indices of function that are domiating the neighboring cells of this patch
/// @param[out] edges          Edges at the surface network mesh
/// @param[out] patches            A connected component of faces bounded by non-manifold edges
/// @param[out] chains         Chains of non-manifold edges
/// @param[out] non_manifold_edges_of_vert         Indices of non-manifold vertices
/// @param[out] shells         An array of shells. Each shell is a connected component consist of patches. There is no mapping involved.
///
/// @param[out] patch_function_label           index of the zero-valued function that creates the patch. It serves as the intermidiate variable.
/// @param[out] cells          A 3D region partitioned by the implicit arrangement (without CSG); encoded by a vector of shell indices. It serves as the intermidiate variable.
/// @param[out] cell_function_label            A 2D boolean array for the signs of each pair of a function and an arrangement cell. It serves as the intermidiate variable.
bool csg(
         bool robust_test,
         bool use_lookup,
         bool use_secondary_lookup,
         bool use_topo_ray_shooting,
         bool positive_inside,
         //
         const std::vector<std::array<double, 3>>& pts,
         const std::vector<std::array<size_t, 4>>& tets,
         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& funcVals,
         const std::function<bool(std::vector<bool>)>& lambda,
         //
         std::vector<std::array<double, 3>>& mesh_pts,
         std::vector<PolygonFace>& mesh_faces,
         std::vector<std::vector<size_t>>& patches,
         std::vector<size_t>& patch_function_label,
         std::vector<Edge>& edges,
         std::vector<std::vector<size_t>>& chains,
         std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
         std::vector<std::vector<size_t>>& shells,
         std::vector<std::vector<size_t>>& cells,
         std::vector<std::vector<bool>>& cell_function_label,
         std::vector<std::string>& timing_labels,
         std::vector<double>& timings,
         std::vector<std::string>& stats_labels,
         std::vector<size_t>& stats);

