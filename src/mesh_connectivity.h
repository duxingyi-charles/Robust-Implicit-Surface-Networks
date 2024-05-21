//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_MESH_CONNECTIVITY_H
#define ROBUST_IMPLICIT_NETWORKS_MESH_CONNECTIVITY_H

#include "mesh.h"

/// Compute iso-edges and edge-face connectivity
///
///@param[in] mesh_faces            Mesh faces taken from the kernel
///@param[out] edges_of_face            A list of face edges containing global edge index
///
///@param[out] mesh_edges           Edges of the mesh; containing two vertices, global face index, and local edge index
void compute_mesh_edges(const std::vector<PolygonFace>& mesh_faces,
                        std::vector<std::vector<size_t>> & edges_of_face,
                        std::vector<Edge>& mesh_edges);

/// Group `iso_faces` into patches for IA:
///
///
/// @param[in] edges_of_faces          the edges' indices of each face with the same orientation from the kernel.
/// @param[in] medg_edges            the list of edges on the mesh.
/// @param[in] mesh_faces            the list of faces on the mesh. Each face contains an implicit function label.
///
///
/// @param[out] patches           output - the list of patches which contain a list of faces' indices.
/// @param[out] patch_function_label             output - the list of patch to implicit funciton label.
void compute_patches(const std::vector<std::vector<size_t>> & edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<PolygonFace>& mesh_faces,
                     std::vector<std::vector<size_t>>& patches,
                     std::vector<size_t>& patch_function_label);
///A validation of patch to function label through taking the majority vote of function labels on all the patch vertices
///
///@param[in] iso_verts             the list of vertices on the mesh. Each vertex contains 1-3 implicit function labesl depending on its connectivity. This purely serves as a verification for the face labels above.
///
///@param[out] Bool            True - passed the patch label checks / False - failed the check
bool check_patch_label(const std::vector<std::vector<size_t>> & edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<PolygonFace>& mesh_faces,
                     const std::vector<IsoVert> iso_verts,
                     std::vector<std::vector<size_t>>& patches,
                     std::vector<size_t>& patch_function_label);

/// Group iso-faces into patches for MI: currently, this support for patch labels is currently missing.
void compute_patches(const std::vector<std::vector<size_t>> & edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<MI_Vert> MI_verts,
                     std::vector<std::vector<size_t>>& patches);

/// Group non-manifold iso-edges into chains
///
/// @param[in] mesh_edges           As above
/// @param[out] non_manifold_edges_of_vert          Indices of non-manifold vertices
/// @param[out] chains          Chains of non-manifold edges; encoded by a vector of edge indices.
void compute_chains(const std::vector<Edge>& mesh_edges,
                    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                    std::vector<std::vector<size_t>>& chains);

/// compute shells and connected components of isosurfaces
/// each shell is a list of half-patches
/// each component is a list of patches
/// we also build maps: half-patch --> shell,  patch --> component
void compute_shells_and_components(size_t num_patch,
                                   const std::vector<std::vector<std::pair<std::pair<size_t, int>,std::pair<size_t, int>>>>& half_patch_pair_list,
                                   std::vector<std::vector<size_t>>& shells,
                                   std::vector<size_t>& shell_of_half_patch,
                                   std::vector<std::vector<size_t>>& components,
                                   std::vector<size_t>& component_of_patch
);

/// group shells into arrangement cells
void compute_arrangement_cells(size_t num_shell,
                               const std::vector<std::pair<size_t,size_t>> &shell_links,
                               std::vector<std::vector<size_t>>& arrangement_cells);


#endif //ROBUST_IMPLICIT_NETWORKS_MESH_CONNECTIVITY_H
