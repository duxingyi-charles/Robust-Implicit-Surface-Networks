//
// Created by Charles Du on 7/20/22.
//

#ifndef ROBUST_IMPLICIT_NETWORKS_MESH_H
#define ROBUST_IMPLICIT_NETWORKS_MESH_H

#include <vector>
#include <array>

static constexpr size_t Mesh_None = std::numeric_limits<size_t>::max();


// a polygon in a polygonal mesh
struct PolygonFace
{
    // a list of polygon's vertex indices (index into some global list of vertices)
    std::vector<size_t> vert_indices;
    // the local index of this polygon in all the tets that contains it
    // each pair is (tet_Id, tet_face_Id)
    std::vector<std::pair<size_t, size_t>> tet_face_indices;
};

// vertex of isosurface
struct IsoVert
{
    // the tet containing the IsoVert
    size_t tet_index;
    // the index of IsoVert in tet.vertices
    size_t tet_vert_index;
    // minimal simplex that contains the IsoVert
    size_t simplex_size; // 1: point, 2: edge, 3: triangle, 4: tetrahedron
    // index into a list of tet vertices
    std::array<size_t, 4> simplex_vert_indices;
    // list of implicit functions whose isosurfaces pass IsoVert (indexed into a global list of
    // implicit functions)
    std::array<size_t, 3> func_indices;
};


// vertex of material interface
struct MI_Vert
{
    // the tet containing the MI_Vert
    size_t tet_index;
    // the index of MI_Vert in tet.vertices
    size_t tet_vert_index;
    // minimal simplex that contains the MI_Vert
    size_t simplex_size; // 1: point, 2: edge, 3: triangle, 4: tetrahedron
    // index into a list of tet vertices
    std::array<size_t, 4> simplex_vert_indices;
    // list of materials whose values are equal at MI_Vert (indexed into a global list of
    // material functions)
    std::array<size_t, 4> material_indices;
};


struct Edge
{
    size_t v1;
    size_t v2;
    // each pair is (face_Id, edge_face_Id)
    // face_Id: face index in the global list of faces
    // edge_face_Id: edge index in the face
    std::vector<std::pair<size_t, size_t>> face_edge_indices;
};

#endif //ROBUST_IMPLICIT_NETWORKS_MESH_H
