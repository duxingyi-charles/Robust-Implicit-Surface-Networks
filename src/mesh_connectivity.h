//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_MESH_CONNECTIVITY_H
#define ROBUST_IMPLICIT_NETWORKS_MESH_CONNECTIVITY_H

#include "mesh.h"

// compute iso-edges and edge-face connectivity
void compute_mesh_edges(const std::vector<PolygonFace>& mesh_faces,
                        std::vector<std::vector<size_t>> & edges_of_face,
                        std::vector<Edge>& mesh_edges);

// group iso-faces into patches for IA:
// edges_of_faces: the edges' indices of each face with the same orientation from the kernel.
// medg_edges: the list of edges on the mesh.
// mesh_faces: the list of faces on the mesh. Each face contains an implicit function label.
// patches: output - the list of patches which contain a list of faces' indices.
// patches: output - the list of patch to implicit funciton label.
void compute_patches(const std::vector<std::vector<size_t>> & edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<PolygonFace>& mesh_faces,
                     std::vector<std::vector<size_t>>& patches,
                     std::vector<size_t>& patch_function_label);
// iso_verts: the list of vertices on the mesh. Each vertex contains 1-3 implicit function labesl depending on its connectivity. This purely serves as a verification for the face labels above.
// output: True - passed the patch label checks / False - failed the check
bool check_patch_label(const std::vector<std::vector<size_t>> & edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<PolygonFace>& mesh_faces,
                     const std::vector<IsoVert> iso_verts,
                     std::vector<std::vector<size_t>>& patches,
                     std::vector<size_t>& patch_function_label);

// group iso-faces into patches for MI: currently, this support for patch labels is missing.
void compute_patches(const std::vector<std::vector<size_t>> & edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<MI_Vert> MI_verts,
                     std::vector<std::vector<size_t>>& patches);

// group non-manifold iso-edges into chains
void compute_chains(const std::vector<Edge>& mesh_edges,
                    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                    std::vector<std::vector<size_t>>& chains);

// compute shells and connected components of isosurfaces
// each shell is a list of half-patches
// each component is a list of patches
// we also build maps: half-patch --> shell,  patch --> component
void compute_shells_and_components(size_t num_patch,
                                   const std::vector<std::vector<std::pair<std::pair<size_t, int>,std::pair<size_t, int>>>>& half_patch_pair_list,
                                   std::vector<std::vector<size_t>>& shells,
                                   std::vector<size_t>& shell_of_half_patch,
                                   std::vector<std::vector<size_t>>& components,
                                   std::vector<size_t>& component_of_patch
);

// group shells into arrangement cells
void compute_arrangement_cells(size_t num_shell,
                               const std::vector<std::pair<size_t,size_t>> &shell_links,
                               std::vector<std::vector<size_t>>& arrangement_cells);


#endif //ROBUST_IMPLICIT_NETWORKS_MESH_CONNECTIVITY_H
