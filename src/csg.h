//
//  csg.h
//  robust_implicit_networks
//
//  Created by Yiwen Ju on 6/3/24.
//

#pragma once
#include <string>


#include "implicit_arrangement.h"

///
/// The pruning of `patches`, `chains`, `non_manifold_edges_of_vert`
///
/// @param[in] 
/// @param[in] mesh_pts            Vertices at the surface network mesh
/// @param[in] mesh_faces          Polygonal faces at the surface network mesh
/// @param[in] patch_function_label           a pair of two indices of function that are domiating the neighboring cells of this patch
/// @param[in] edges          Edges at the surface network mesh
/// @param[in] cells          A 3D region partitioned by the surface network; encoded by a vector of shell indices
/// @param[in] lambda           The indicator function takes a 2D vector, which each row represents a list of functions signs (one per implicit function) within a cell, and outputs a 1D vector cell_label, which indicates each cell is within or outside of the shape.
///
///
/// @param[out] patches            A connected component of faces bounded by non-manifold edges
/// @param[out] chains         Chains of non-manifold edges
/// @param[out] non_manifold_edges_of_vert         Indices of non-manifold vertices
/// @param[out] shells         An array of shells. Each shell is a connected component consist of patches. There is no mapping involved.
///
void csg(const std::vector<std::array<double, 3>>& mesh_pts,
                const std::vector<PolygonFace>& mesh_faces,
                std::vector<std::vector<size_t>>& patches,
                const std::vector<size_t>& patch_function_label,
                const std::vector<Edge>& edges,
                std::vector<std::vector<size_t>>& chains,
                std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                std::vector<std::vector<size_t>>& shells,
                const std::vector<std::vector<size_t>>& cells,
                const std::vector<std::vector<bool>>& cell_function_label,
                const std::function<std::vector<bool>(std::vector<std::vector<bool>>)>& lambda);

