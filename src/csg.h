//
//  csg.h
//  robust_implicit_networks
//
//  Created by Yiwen Ju on 6/3/24.
//

#ifndef csg_h
#define csg_h
#include <fstream>
#include <string>
#include <nlohmann/json.hpp>


#include "implicit_arrangement.h"

///
///The structure of a csg unit
///A nested data structure where `operation` takes in a boolean operation, and the elements contain either a function index or later csg unit, where netaives represent the indices, and positives reprsent the index in the CSG tree.
///e.g. {Intersection, -1, 2} represents an intersection operation between the first function and the second csg unit.
///
struct csg_unit{
    int operation;
    std::array<int, 2> elements;
};

///
///Defines the enumeration of CSG operations. In the data structure, Intersection is 0, Union is 1, and Negation is 2.
enum csg_operations{
    Intersection,
    Union,
    Negation
};

///load the csg file
///@param[in] filename          The name of the input CSG tree
///@param[out] tree         The loaded tree structure
///
bool load_csgTree(const std::string filename, std::vector<csg_unit>& tree);

/// an iterative algortihm that traverses through the csg tree
///
///@param[in] csgTree           The CSG structure: a list of csg units
///@param[in] curNode           The current index in the csg structure
///@param[in] funcInt           The intervals of all the functions
///
///@param[out] std::pair
std::pair<std::array<double, 2>, std::vector<int>> iterTree(const std::vector<csg_unit>csgTree,const int curNode,const std::vector<std::array<double , 2>> funcInt);

///
/// The pruning of `patches`, `chains`, `non_manifold_edges_of_vert`
///
/// @param[in] 
/// @param[in] mesh_pts            Vertices at the surface network mesh
/// @param[in] mesh_faces          Polygonal faces at the surface network mesh
/// @param[in] patch_function_label           a pair of two indices of function that are domiating the neighboring cells of this patch
/// @param[in] edges          Edges at the surface network mesh
/// @param[in] shells         An array of shells. Each shell is a connected component consist of patches.
/// @param[in] cells          A 3D region partitioned by the surface network; encoded by a vector of shell indices
/// @param[in] cells_label            a 1D vector of boolens for each cell that represents inside or outside of the CSG shape.
///
///
/// @param[out] patches            A connected component of faces bounded by non-manifold edges
/// @param[out] chains         Chains of non-manifold edges
/// @param[out] non_manifold_edges_of_vert         Indices of non-manifold vertices
///
void prune_data(const std::vector<std::array<double, 3>>& mesh_pts,
                const std::vector<PolygonFace>& mesh_faces,
                std::vector<std::vector<size_t>>& patches,
                const std::vector<size_t>& patch_function_label,
                const std::vector<Edge>& edges,
                std::vector<std::vector<size_t>>& chains,
                std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                const std::vector<std::vector<size_t>>& shells,
                const std::vector<std::vector<size_t>>& cells,
                const std::vector<bool>& cells_label);

#endif /* csg_h */
