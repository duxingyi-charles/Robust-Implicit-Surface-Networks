//
// Created by Charles Du on 7/20/22.
//

#include "topo_ray_shooting.h"
#include "mesh_connectivity.h"


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
)
{
    // map: tet vert index --> index of next vert (with smaller (x,y,z))
    std::vector<size_t> next_vert;
    build_next_vert(pts, tets, next_vert);

    // find extremal edge for each component
    // extremal edge of component i is stored at position [2*i], [2*i+1]
    std::vector<size_t> extremal_edge_of_component;
    // store an iso-vert index on edge (v, v_next), Mesh_None means there is no such iso-vert
    std::vector<size_t> iso_vert_on_v_v_next;
    // map: (tet_id, tet_face_id) --> iso_face_id
    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> iso_face_id_of_tet_face;
    // map: (tet_id, tet_vert_id) --> (iso_vert_id, component_id)
    absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
            iso_vId_compId_of_tet_vert;
    find_extremal_edges(pts, iso_verts, iso_faces, patches, components, component_of_patch, next_vert,
    extremal_edge_of_component, iso_vert_on_v_v_next, iso_face_id_of_tet_face, iso_vId_compId_of_tet_vert);


    // topological ray shooting
    std::vector<std::pair<size_t, size_t>> shell_links;
    {
//        timing_labels.emplace_back("arrCells(ray shooting)");
//        ScopedTimer<> timer("arrangement cells: topo ray shooting");
        shell_links.reserve(components.size());
        std::vector<size_t> sorted_vert_indices_on_edge;
        sorted_vert_indices_on_edge.reserve(3);
        for (size_t i = 0; i < components.size(); ++i) {
            // extremal edge: v1 -> v2
            auto extreme_v1 = extremal_edge_of_component[2 * i];
            auto extreme_v2 = extremal_edge_of_component[2 * i + 1];
            auto iso_vId = iso_vert_on_v_v_next[extreme_v1];
            auto tetId = iso_verts[iso_vId].tet_index;
            const auto& tet_cut_result = cut_results[cut_result_index[tetId]];
            // get local index of v1 and v2 in the tet
            size_t local_v1, local_v2;
            for (size_t j = 0; j < 4; ++j) {
                if (tets[tetId][j] == extreme_v1) {
                    local_v1 = j;
                } else if (tets[tetId][j] == extreme_v2) {
                    local_v2 = j;
                }
            }
            // get an ordered list of vertices on edge v1 -> v2
            sorted_vert_indices_on_edge.clear();
            compute_edge_intersection_order(
                    tet_cut_result, local_v1, local_v2, sorted_vert_indices_on_edge);
            // find the vertex v_start on v1->v2
            // 1. on current component
            // 2. nearest to v2
            size_t j_start;
            for (size_t j = 0; j + 1 < sorted_vert_indices_on_edge.size(); ++j) {
                const auto& iso_vId_compId = iso_vId_compId_of_tet_vert[std::make_pair(
                        tetId, sorted_vert_indices_on_edge[j])];
                if (iso_vId_compId.second == i) {
                    j_start = j;
                }
            }
            if (j_start + 2 < sorted_vert_indices_on_edge.size()) {
                // there is a vert from another component between v_start -> v2
                std::pair<size_t, int> face_orient1, face_orient2;
                compute_passing_face_pair(tet_cut_result,
                                          sorted_vert_indices_on_edge[j_start],
                                          sorted_vert_indices_on_edge[j_start + 1],
                                          face_orient1,
                                          face_orient2);
                size_t iso_fId1 =
                        iso_face_id_of_tet_face[std::make_pair(tetId, face_orient1.first)];
                size_t iso_fId2 =
                        iso_face_id_of_tet_face[std::make_pair(tetId, face_orient2.first)];
                size_t shell1 =
                        (face_orient1.second == 1)
                        ? shell_of_half_patch[2 * patch_of_face[iso_fId1]]
                        : shell_of_half_patch[2 * patch_of_face[iso_fId1] + 1];
                size_t shell2 =
                        (face_orient2.second == 1)
                        ? shell_of_half_patch[2 * patch_of_face[iso_fId2]]
                        : shell_of_half_patch[2 * patch_of_face[iso_fId2] + 1];
                // link shell1 with shell2
                shell_links.emplace_back(shell1, shell2);
            } else {
                // there is no vert between v_start -> v2
                std::pair<size_t, int> face_orient;
                compute_passing_face(tet_cut_result,
                                     sorted_vert_indices_on_edge[j_start],
                                     sorted_vert_indices_on_edge.back(),
                                     face_orient);
                size_t iso_fId =
                        iso_face_id_of_tet_face[std::make_pair(tetId, face_orient.first)];
                size_t shell_start =
                        (face_orient.second == 1)
                        ? shell_of_half_patch[2 * patch_of_face[iso_fId]]
                        : shell_of_half_patch[2 * patch_of_face[iso_fId] + 1];
                // follow the ray till another iso-vertex or the sink
                auto v_curr = extreme_v2;
                while (next_vert[v_curr] != Mesh_None &&
                       iso_vert_on_v_v_next[v_curr] == Mesh_None) {
                    v_curr = next_vert[v_curr];
                }
                if (iso_vert_on_v_v_next[v_curr] != Mesh_None) {
                    // reached iso-vert at end of the ray
                    auto iso_vId_end = iso_vert_on_v_v_next[v_curr];
                    auto end_tetId = iso_verts[iso_vId_end].tet_index;
                    const auto& end_tet_cut_result =
                            cut_results[cut_result_index[end_tetId]];
                    auto v_next = next_vert[v_curr];
                    // find local vertex indices in the end tetrahedron
                    for (size_t j = 0; j < 4; ++j) {
                        if (tets[end_tetId][j] == v_curr) {
                            local_v1 = j;
                        } else if (tets[end_tetId][j] == v_next) {
                            local_v2 = j;
                        }
                    }
                    // get an ordered list of vertices on edge v_curr -> v_next
                    sorted_vert_indices_on_edge.clear();
                    compute_edge_intersection_order(end_tet_cut_result,
                                                    local_v1,
                                                    local_v2,
                                                    sorted_vert_indices_on_edge);
                    // find the end shell
                    compute_passing_face(end_tet_cut_result,
                                         sorted_vert_indices_on_edge[1],
                                         sorted_vert_indices_on_edge.front(),
                                         face_orient);
                    iso_fId = iso_face_id_of_tet_face[std::make_pair(
                            end_tetId, face_orient.first)];
                    size_t shell_end =
                            (face_orient.second == 1)
                            ? shell_of_half_patch[2 * patch_of_face[iso_fId]]
                            : shell_of_half_patch[2 * patch_of_face[iso_fId] + 1];
                    // link start shell with end shell
                    shell_links.emplace_back(shell_start, shell_end);
                } else {
                    // next_vert[v_curr] is Mesh_None, v_curr is the sink vertex
                    // link shell_start with the sink
                    shell_links.emplace_back(shell_start, Mesh_None);
                }
            }
        }
//        timings.push_back(timer.toc());
    }

    // group shells into arrangement cells
    compute_arrangement_cells(shells.size(), shell_links, arrangement_cells);
}

void topo_ray_shooting(
        const std::vector<std::array<double, 3>> &pts,
        const std::vector<std::array<size_t, 4>> &tets,
        const std::vector<simplicial_arrangement::MaterialInterface<3>> &cut_results,
        const std::vector<size_t> &cut_result_index,
        const std::vector<MI_Vert> &MI_verts,
        const std::vector<PolygonFace> &MI_faces,
        const std::vector<std::vector<size_t>> &patches,
        const std::vector<size_t> &patch_of_face,
        const std::vector<std::vector<size_t>> &shells,
        const std::vector<size_t> &shell_of_half_patch,
        const std::vector<std::vector<size_t>> &components,
        const std::vector<size_t> &component_of_patch,
        std::vector<std::vector<size_t>> &material_cells
) {
// map: tet vert index --> index of next vert (with smaller (x,y,z))
    std::vector<size_t> next_vert;
    build_next_vert(pts, tets, next_vert);

    // find extremal edge for each component
    // extremal edge of component i is stored at position [2*i], [2*i+1]
    std::vector<size_t> extremal_edge_of_component;
    // store an MI-vert index on edge (v, v_next), Mesh_None means there is no such MI-vert
    std::vector<size_t> MI_vert_on_v_v_next;
    // map: (tet_id, tet_face_id) --> MI_face_id
    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> MI_face_id_of_tet_face;
    // map: (tet_id, tet_vert_id) --> (MI_vert_id, component_id)
    absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
            MI_vId_compId_of_tet_vert;
    find_extremal_edges(pts, MI_verts, MI_faces, patches, components, component_of_patch, next_vert,
                        extremal_edge_of_component, MI_vert_on_v_v_next, MI_face_id_of_tet_face, MI_vId_compId_of_tet_vert);


    // topological ray shooting
    std::vector<std::pair<size_t, size_t>> shell_links;
    {
        shell_links.reserve(components.size());
        std::vector<size_t> sorted_vert_indices_on_edge;
        sorted_vert_indices_on_edge.reserve(3);
        for (size_t i = 0; i < components.size(); ++i) {
            // extremal edge: v1 -> v2
            auto extreme_v1 = extremal_edge_of_component[2 * i];
            auto extreme_v2 = extremal_edge_of_component[2 * i + 1];
            auto MI_vId = MI_vert_on_v_v_next[extreme_v1];
            auto tetId = MI_verts[MI_vId].tet_index;
            const auto& tet_cut_result = cut_results[cut_result_index[tetId]];
            // get local index of v1 and v2 in the tet
            size_t local_v1, local_v2;
            for (size_t j = 0; j < 4; ++j) {
                if (tets[tetId][j] == extreme_v1) {
                    local_v1 = j;
                } else if (tets[tetId][j] == extreme_v2) {
                    local_v2 = j;
                }
            }
            //                    std::cout << "local_v1 = " << local_v1 << std::endl;
            //                    std::cout << "local_v2 = " << local_v2 << std::endl;
            // get an ordered list of vertices on edge v1 -> v2
            sorted_vert_indices_on_edge.clear();
            compute_edge_intersection_order(
                    tet_cut_result, local_v1, local_v2, sorted_vert_indices_on_edge);
            
            // find the vertex v_start on v1->v2
            // 1. on current component
            // 2. nearest to v2
            size_t j_start;
            for (size_t j = 0; j + 1 < sorted_vert_indices_on_edge.size(); ++j) {
                const auto& MI_vId_compId = MI_vId_compId_of_tet_vert[std::make_pair(
                        tetId, sorted_vert_indices_on_edge[j])];
                if (MI_vId_compId.second == i) {
                    j_start = j;
                }
            }
            if (j_start + 2 < sorted_vert_indices_on_edge.size()) {
                // there is a vert from another component between v_start -> v2
                std::pair<size_t, int> face_orient1, face_orient2;
                compute_passing_face_pair(tet_cut_result,
                                             sorted_vert_indices_on_edge[j_start],
                                             sorted_vert_indices_on_edge[j_start + 1],
                                             face_orient1,
                                             face_orient2);
                size_t MI_fId1 =
                        MI_face_id_of_tet_face[std::make_pair(tetId, face_orient1.first)];
                size_t MI_fId2 =
                        MI_face_id_of_tet_face[std::make_pair(tetId, face_orient2.first)];
                size_t shell1 =
                        (face_orient1.second == 1)
                        ? shell_of_half_patch[2 * patch_of_face[MI_fId1]]
                        : shell_of_half_patch[2 * patch_of_face[MI_fId1] + 1];
                size_t shell2 =
                        (face_orient2.second == 1)
                        ? shell_of_half_patch[2 * patch_of_face[MI_fId2]]
                        : shell_of_half_patch[2 * patch_of_face[MI_fId2] + 1];
                // link shell1 with shell2
                shell_links.emplace_back(shell1, shell2);
            } else {
                // there is no vert between v_start -> v2
                std::pair<size_t, int> face_orient;
                compute_passing_face(tet_cut_result,
                                        sorted_vert_indices_on_edge[j_start],
                                        sorted_vert_indices_on_edge.back(),
                                        face_orient);
                size_t MI_fId =
                        MI_face_id_of_tet_face[std::make_pair(tetId, face_orient.first)];
                size_t shell_start =
                        (face_orient.second == 1)
                        ? shell_of_half_patch[2 * patch_of_face[MI_fId]]
                        : shell_of_half_patch[2 * patch_of_face[MI_fId] + 1];
                // follow the ray till another iso-vertex or the sink
                auto v_curr = extreme_v2;
                while (next_vert[v_curr] != Mesh_None &&
                       MI_vert_on_v_v_next[v_curr] == Mesh_None) {
                    v_curr = next_vert[v_curr];
                }
                if (MI_vert_on_v_v_next[v_curr] != Mesh_None) {
                    // reached iso-vert at end of the ray
                    auto MI_vId_end = MI_vert_on_v_v_next[v_curr];
                    auto end_tetId = MI_verts[MI_vId_end].tet_index;
                    const auto& end_tet_cut_result =
                            cut_results[cut_result_index[end_tetId]];
                    auto v_next = next_vert[v_curr];
                    // find local vertex indices in the end tetrahedron
                    for (size_t j = 0; j < 4; ++j) {
                        if (tets[end_tetId][j] == v_curr) {
                            local_v1 = j;
                        } else if (tets[end_tetId][j] == v_next) {
                            local_v2 = j;
                        }
                    }
                    // get an ordered list of vertices on edge v_curr -> v_next
                    sorted_vert_indices_on_edge.clear();
                    compute_edge_intersection_order(end_tet_cut_result,
                                                       local_v1,
                                                       local_v2,
                                                       sorted_vert_indices_on_edge);
                    // find the end shell
                    compute_passing_face(end_tet_cut_result,
                                            sorted_vert_indices_on_edge[1],
                                            sorted_vert_indices_on_edge.front(),
                                            face_orient);
                    MI_fId = MI_face_id_of_tet_face[std::make_pair(
                            end_tetId, face_orient.first)];
                    size_t shell_end =
                            (face_orient.second == 1)
                            ? shell_of_half_patch[2 * patch_of_face[MI_fId]]
                            : shell_of_half_patch[2 * patch_of_face[MI_fId] + 1];
                    // link start shell with end shell
                    shell_links.emplace_back(shell_start, shell_end);
                } else {
                    // next_vert[v_curr] is None, v_curr is the sink vertex
                    // link shell_start with the sink
                    shell_links.emplace_back(shell_start, Mesh_None);
                }
            }
        }
    }

    // group shells into arrangement cells
    compute_arrangement_cells(shells.size(), shell_links, material_cells);
}

void build_next_vert(const std::vector<std::array<double, 3>>& pts,
                     const std::vector<std::array<size_t, 4>>& tets,
                     std::vector<size_t>& next_vert)
{
    next_vert.resize(pts.size(), std::numeric_limits<size_t>::max());
    for (const auto& tet : tets) {
        // find the smallest vertex of tet
        size_t min_id = 0;
        for (size_t i = 1; i < 4; ++i) {
            if (point_xyz_less(pts[tet[i]], pts[tet[min_id]])) {
                min_id = i;
            }
        }
        size_t min_vId = tet[min_id];
        //
        for (size_t i = 0; i < 4; ++i) {
            if (i != min_id) {
                next_vert[tet[i]] = min_vId;
            }
        }
    }
}

void find_extremal_edges(
        const std::vector<std::array<double, 3>> &pts,
        const std::vector<IsoVert> &iso_verts,
        const std::vector<PolygonFace> &iso_faces,
        const std::vector<std::vector<size_t>> &patches,
        const std::vector<std::vector<size_t>> &components,
        const std::vector<size_t> &component_of_patch,
        const std::vector<size_t> &next_vert,
        std::vector<size_t> &extremal_edge_of_component,
        std::vector<size_t> &iso_vert_on_v_v_next,
        absl::flat_hash_map<std::pair<size_t, size_t>, size_t> &iso_face_id_of_tet_face,
        absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
        &iso_vId_compId_of_tet_vert)
{
    extremal_edge_of_component.resize(2 * components.size(), Mesh_None);
    iso_vert_on_v_v_next.resize(pts.size(), Mesh_None);
    iso_face_id_of_tet_face.reserve(iso_faces.size());
    iso_vId_compId_of_tet_vert.reserve(iso_faces.size() / 2);
    //
    std::vector<bool> is_iso_vert_visited(iso_verts.size(), false);
    for (size_t i = 0; i < patches.size(); ++i) {
        size_t component_id = component_of_patch[i];
        auto& u1 = extremal_edge_of_component[2 * component_id];
        auto& u2 = extremal_edge_of_component[2 * component_id + 1];
        for (auto fId : patches[i]) {
            for (const auto& tet_face : iso_faces[fId].tet_face_indices) {
                iso_face_id_of_tet_face.try_emplace(tet_face, fId);
            }
            for (auto vId : iso_faces[fId].vert_indices) {
                if (!is_iso_vert_visited[vId]) {
                    is_iso_vert_visited[vId] = true;
                    const auto& vert = iso_verts[vId];
                    if (vert.simplex_size == 2) { // edge iso-vertex
                        auto v1 = vert.simplex_vert_indices[0];
                        auto v2 = vert.simplex_vert_indices[1];
                        if (next_vert[v1] == v2) { // on tree edge v1 -> v2
                            // update extremal edge
                            if (u1 == Mesh_None) {
                                u1 = v1;
                                u2 = v2;
                            } else {
                                if (v2 == u2) {
                                    if (point_xyz_less(pts[v1], pts[u1])) {
                                        u1 = v1;
                                    }
                                } else if (point_xyz_less(pts[v2], pts[u2])) {
                                    u1 = v1;
                                    u2 = v2;
                                }
                            }
                            // record an iso-vert on edge v1 -> v2
                            iso_vert_on_v_v_next[v1] = vId;
                            // fill map
                            iso_vId_compId_of_tet_vert.try_emplace(
                                    std::make_pair(vert.tet_index, vert.tet_vert_index),
                                    std::make_pair(vId, component_id));
                        } else if (next_vert[v2] == v1) { // on tree edge v2 -> v1
                            // update extremal edge
                            if (u1 == Mesh_None) {
                                u1 = v2;
                                u2 = v1;
                            } else {
                                if (v1 == u2) {
                                    if (point_xyz_less(pts[v2], pts[u1])) {
                                        u1 = v2;
                                    }
                                } else if (point_xyz_less(pts[v1], pts[u2])) {
                                    u1 = v2;
                                    u2 = v1;
                                }
                            }
                            // record an iso-vert on v2 -> v1
                            iso_vert_on_v_v_next[v2] = vId;
                            // fill map
                            iso_vId_compId_of_tet_vert.try_emplace(
                                    std::make_pair(vert.tet_index, vert.tet_vert_index),
                                    std::make_pair(vId, component_id));
                        }
                    }
                }
            }
        }
    }
}

void find_extremal_edges(
        const std::vector<std::array<double, 3>> &pts,
        const std::vector<MI_Vert> &MI_verts,
        const std::vector<PolygonFace> &MI_faces,
        const std::vector<std::vector<size_t>> &patches,
        const std::vector<std::vector<size_t>> &components,
        const std::vector<size_t> &component_of_patch,
        const std::vector<size_t> &next_vert,
        // extremal edge of component i is stored at position [2*i], [2*i+1]
        std::vector<size_t> &extremal_edge_of_component,
        // store an MI-vert index on edge (v, v_next), None means there is no such MI-vert
        std::vector<size_t> &MI_vert_on_v_v_next,
        // map: (tet_id, tet_face_id) --> MI_face_id
        absl::flat_hash_map<std::pair<size_t, size_t>, size_t> &MI_face_id_of_tet_face,
        // map: (tet_id, tet_vert_id) --> (MI_vert_id, component_id)
        absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
        &MI_vId_compId_of_tet_vert) {
    extremal_edge_of_component.resize(
            2 * components.size(), Mesh_None);
    MI_vert_on_v_v_next.resize(pts.size(), Mesh_None);
    MI_face_id_of_tet_face.reserve(MI_faces.size());
    MI_vId_compId_of_tet_vert.reserve(MI_faces.size() / 2);
    //
    std::vector<bool> is_MI_vert_visited(MI_verts.size(), false);
    for (size_t i = 0; i < patches.size(); ++i) {
        size_t component_id = component_of_patch[i];
        auto& u1 = extremal_edge_of_component[2 * component_id];
        auto& u2 = extremal_edge_of_component[2 * component_id + 1];
        for (auto fId : patches[i]) {
            for (const auto& tet_face : MI_faces[fId].tet_face_indices) {
                MI_face_id_of_tet_face.try_emplace(tet_face, fId);
            }
            for (auto vId : MI_faces[fId].vert_indices) {
                if (!is_MI_vert_visited[vId]) {
                    is_MI_vert_visited[vId] = true;
                    const auto& vert = MI_verts[vId];
                    if (vert.simplex_size == 2) { // edge MI-vertex
                        auto v1 = vert.simplex_vert_indices[0];
                        auto v2 = vert.simplex_vert_indices[1];
                        if (next_vert[v1] == v2) { // on tree edge v1 -> v2
                            // update extremal edge
                            if (u1 == Mesh_None) {
                                u1 = v1;
                                u2 = v2;
                            } else {
                                if (v2 == u2) {
                                    if (point_xyz_less(pts[v1], pts[u1])) {
                                        u1 = v1;
                                    }
                                } else if (point_xyz_less(pts[v2], pts[u2])) {
                                    u1 = v1;
                                    u2 = v2;
                                }
                            }
                            // record an iso-vert on edge v1 -> v2
                            MI_vert_on_v_v_next[v1] = vId;
                            // fill map
                            MI_vId_compId_of_tet_vert.try_emplace(
                                    std::make_pair(vert.tet_index, vert.tet_vert_index),
                                    std::make_pair(vId, component_id));
                        } else if (next_vert[v2] == v1) { // on tree edge v2 -> v1
                            // update extremal edge
                            if (u1 == Mesh_None) {
                                u1 = v2;
                                u2 = v1;
                            } else {
                                if (v1 == u2) {
                                    if (point_xyz_less(pts[v2], pts[u1])) {
                                        u1 = v2;
                                    }
                                } else if (point_xyz_less(pts[v1], pts[u2])) {
                                    u1 = v2;
                                    u2 = v1;
                                }
                            }
                            // record an iso-vert on v2 -> v1
                            MI_vert_on_v_v_next[v2] = vId;
                            // fill map
                            MI_vId_compId_of_tet_vert.try_emplace(
                                    std::make_pair(vert.tet_index, vert.tet_vert_index),
                                    std::make_pair(vId, component_id));
                        }
                    }
                }
            }
        }
    }
}

void compute_edge_intersection_order(
        const simplicial_arrangement::Arrangement<3>& tet_cut_result, size_t v, size_t u, std::vector<size_t>& vert_indices)
{
    const auto& vertices = tet_cut_result.vertices;
    const auto& faces = tet_cut_result.faces;
    //
    std::array<bool, 4> edge_flag{true, true, true, true};
    edge_flag[v] = false;
    edge_flag[u] = false;

    // find vertices on edge v->u, and index of v and u in tet_cut_result.vertices
    size_t v_id, u_id;
    std::vector<bool> is_on_edge_vu(vertices.size(), false);
    size_t num_vert_in_edge_vu = 0;
    size_t flag_count;
    size_t other_plane;
    for (size_t i = 0; i < vertices.size(); ++i) {
        flag_count = 0;
        other_plane = Mesh_None;
        const auto& vert = vertices[i];
        for (int j = 0; j < 3; ++j) {
            if (vert[j] < 4) { // 0,1,2,3 are tet boundary planes
                if (edge_flag[vert[j]]) {
                    ++flag_count;
                } else {
                    other_plane = vert[j];
                }
            }
        }
        if (flag_count == 2) {
            is_on_edge_vu[i] = true;
            if (other_plane == u) { // current vertex is v
                v_id = i;
            } else if (other_plane == v) { // current vertex is u
                u_id = i;
            } else { // current vertex is in interior of edge v->u
                ++num_vert_in_edge_vu;
            }
        }
    }
    if (num_vert_in_edge_vu == 0) { // no intersection points in edge v->u
        vert_indices.push_back(v_id);
        vert_indices.push_back(u_id);
        return;
    }

    // find all faces on triangle v->u->w, w is a third tet vertex
    size_t pId; // plane index of v->u->w, pick pId != v and pId != u
    for (size_t i = 0; i < 4; ++i) {
        if (edge_flag[i]) {
            pId = i;
            break;
        }
    }
    //    std::vector<bool> is_on_tri_vuw(faces.size(), false);
    std::vector<size_t> faces_on_tri_vuw;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        if (face.negative_cell == Mesh_None) { // face is on tet boundary
            if (face.supporting_plane == pId) {
                faces_on_tri_vuw.push_back(i);
            }
        }
    }

    // build edge-face connectivity on triangle v->u->w
    // map: edge (v1, v2) --> incident faces (f1, f2), f2 is None if (v1,v2) is on triangle boundary
    absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> faces_of_edge;
    for (auto fId : faces_on_tri_vuw) {
        const auto& face = faces[fId];
        size_t num_vert = face.vertices.size();
        for (size_t i = 0; i < num_vert; ++i) {
            size_t i_next = (i + 1) % num_vert;
            size_t vi = face.vertices[i];
            size_t vi_next = face.vertices[i_next];
            // add fId to edge (vi, vi_next)
            auto iter_inserted = faces_of_edge.try_emplace(
                    std::make_pair(vi, vi_next), std::make_pair(fId, Mesh_None));
            if (!iter_inserted.second) { // inserted before
                iter_inserted.first->second.second = fId;
            }
            // add fId to edge (vi_next, vi)
            iter_inserted = faces_of_edge.try_emplace(
                    std::make_pair(vi_next, vi), std::make_pair(fId, Mesh_None));
            if (!iter_inserted.second) { // inserted before
                iter_inserted.first->second.second = fId;
            }
        }
    }

    // find the face on triangle v->u->w:
    // 1. has vertex v
    // 2. has an edge on v->u
    size_t f_start;
    for (auto fId : faces_on_tri_vuw) {
        const auto& face = faces[fId];
        bool find_v = false;
        size_t count = 0;
        for (auto vi : face.vertices) {
            if (vi == v_id) {
                find_v = true;
            }
            if (is_on_edge_vu[vi]) {
                ++count;
            }
        }
        if (find_v && count == 2) {
            f_start = fId;
            break;
        }
    }

    // trace edge v->u
    vert_indices.reserve(num_vert_in_edge_vu + 2);
    vert_indices.push_back(v_id);
    std::vector<bool> visited_face(faces.size(), false);
    size_t v_curr = v_id;
    size_t f_curr = f_start;
    std::pair<size_t, size_t> edge_prev;
    std::pair<size_t, size_t> edge_next;
    std::pair<size_t, size_t> edge_on_vu;
    while (v_curr != u_id) {
        // clear edge_on_vu
        edge_on_vu.first = Mesh_None;
        edge_on_vu.second = Mesh_None;
        //
        const auto& face = faces[f_curr];
        size_t num_vert = face.vertices.size();
        // visit all edges of face, find edge_prev, edge_next and edge_on_vu
        for (size_t i = 0; i < num_vert; ++i) {
            size_t i_next = (i + 1) % num_vert;
            size_t vi = face.vertices[i];
            size_t vi_next = face.vertices[i_next];
            if (is_on_edge_vu[vi] && !is_on_edge_vu[vi_next]) {
                auto& two_faces = faces_of_edge[std::make_pair(vi, vi_next)];
                auto other_face = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
                if (vi == v_id ||
                    (other_face != Mesh_None && visited_face[other_face])) {
                    edge_prev.first = vi;
                    edge_prev.second = vi_next;
                } else {
                    edge_next.first = vi;
                    edge_next.second = vi_next;
                }
            } else if (is_on_edge_vu[vi_next] && !is_on_edge_vu[vi]) {
                auto& two_faces = faces_of_edge[std::make_pair(vi, vi_next)];
                auto other_face = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
                if (vi_next == v_id ||
                    (other_face != Mesh_None && visited_face[other_face])) {
                    edge_prev.first = vi;
                    edge_prev.second = vi_next;
                } else {
                    edge_next.first = vi;
                    edge_next.second = vi_next;
                }
            } else if (is_on_edge_vu[vi] && is_on_edge_vu[vi_next]) {
                edge_on_vu.first = vi;
                edge_on_vu.second = vi_next;
            }
        }
        //
        if (edge_on_vu.first == Mesh_None) {
            // no edge of the face is on v->u
            // keep v_curr, update f_curr
            visited_face[f_curr] = true;
            auto& two_faces = faces_of_edge[edge_next];
            f_curr = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
        } else {
            // there is an edge of the face on v->u
            // update v_curr to be the next vert on edge v->u
            v_curr = (edge_on_vu.first == v_curr) ? edge_on_vu.second : edge_on_vu.first;
            vert_indices.push_back(v_curr);
            // update f_curr
            visited_face[f_curr] = true;
            auto& two_faces = faces_of_edge[edge_next];
            f_curr = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
        }
    }
}

void compute_edge_intersection_order(const simplicial_arrangement::MaterialInterface<3>& tet_cut_result,
                                     size_t v,
                                     size_t u,
                                     std::vector<size_t>& vert_indices)
{
    const auto& vertices = tet_cut_result.vertices;
    const auto& faces = tet_cut_result.faces;
    //
    std::array<bool, 4> edge_flag{true, true, true, true};
    edge_flag[v] = false;
    edge_flag[u] = false;

    // find vertices on edge v->u, and index of v and u in tet_cut_result.vertices
    size_t v_id, u_id;
    std::vector<bool> is_on_edge_vu(vertices.size(), false);
    size_t num_vert_in_edge_vu = 0;
    size_t flag_count;
    size_t other_plane;
    for (size_t i = 0; i < vertices.size(); ++i) {
        flag_count = 0;
        other_plane = Mesh_None;
        const auto& vert = vertices[i];
        // vert.size() == 4
        for (int j = 0; j < 4; ++j) {
            if (vert[j] < 4) { // 0,1,2,3 are tet boundary
                if (edge_flag[vert[j]]) {
                    ++flag_count;
                } else {
                    other_plane = vert[j];
                }
            }
        }
        if (flag_count == 2) {
            is_on_edge_vu[i] = true;
            if (other_plane == u) { // current vertex is v
                v_id = i;
            } else if (other_plane == v) { // current vertex is u
                u_id = i;
            } else { // current vertex is in interior of edge v->u
                ++num_vert_in_edge_vu;
            }
        }
    }
    if (num_vert_in_edge_vu == 0) { // no intersection points in edge v->u
        vert_indices.push_back(v_id);
        vert_indices.push_back(u_id);
        return;
    }

    // find all faces on triangle v->u->w, w is a third tet vertex
    size_t pId; // plane index of v->u->w, pick pId != v and pId != u
    for (size_t i = 0; i < 4; ++i) {
        if (edge_flag[i]) {
            pId = i;
            break;
        }
    }
    std::vector<size_t> faces_on_tri_vuw;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        if (face.positive_material_label == pId) { // face is on tet boundary pId
            faces_on_tri_vuw.push_back(i);
        }
    }

    // build edge-face connectivity on triangle v->u->w
    // map: edge (v1, v2) --> incident faces (f1, f2), f2 is None if (v1,v2) is on triangle boundary
    absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> faces_of_edge;
    for (auto fId : faces_on_tri_vuw) {
        const auto& face = faces[fId];
        size_t num_vert = face.vertices.size();
        for (size_t i = 0; i < num_vert; ++i) {
            size_t i_next = (i + 1) % num_vert;
            size_t vi = face.vertices[i];
            size_t vi_next = face.vertices[i_next];
            // add fId to edge (vi, vi_next)
            auto iter_inserted = faces_of_edge.try_emplace(
                    std::make_pair(vi, vi_next), std::make_pair(fId, Mesh_None));
            if (!iter_inserted.second) { // inserted before
                iter_inserted.first->second.second = fId;
            }
            // add fId to edge (vi_next, vi)
            iter_inserted = faces_of_edge.try_emplace(
                    std::make_pair(vi_next, vi), std::make_pair(fId, Mesh_None));
            if (!iter_inserted.second) { // inserted before
                iter_inserted.first->second.second = fId;
            }
        }
    }

    // find the face on triangle v->u->w:
    // 1. has vertex v
    // 2. has an edge on v->u
    size_t f_start;
    for (auto fId : faces_on_tri_vuw) {
        const auto& face = faces[fId];
        bool find_v = false;
        size_t count = 0;
        for (auto vi : face.vertices) {
            if (vi == v_id) {
                find_v = true;
            }
            if (is_on_edge_vu[vi]) {
                ++count;
            }
        }
        if (find_v && count == 2) {
            f_start = fId;
            break;
        }
    }

    // trace edge v->u
    vert_indices.reserve(num_vert_in_edge_vu + 2);
    vert_indices.push_back(v_id);
    std::vector<bool> visited_face(faces.size(), false);
    size_t v_curr = v_id;
    size_t f_curr = f_start;
    std::pair<size_t, size_t> edge_prev;
    std::pair<size_t, size_t> edge_next;
    std::pair<size_t, size_t> edge_on_vu;
    while (v_curr != u_id) {
        // clear edge_on_vu
        edge_on_vu.first = Mesh_None;
        edge_on_vu.second = Mesh_None;
        //
        const auto& face = faces[f_curr];
        size_t num_vert = face.vertices.size();
        // visit all edges of face, find edge_prev, edge_next and edge_on_vu
        for (size_t i = 0; i < num_vert; ++i) {
            size_t i_next = (i + 1) % num_vert;
            size_t vi = face.vertices[i];
            size_t vi_next = face.vertices[i_next];
            if (is_on_edge_vu[vi] && !is_on_edge_vu[vi_next]) {
                auto& two_faces = faces_of_edge[std::make_pair(vi, vi_next)];
                auto other_face = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
                if (vi == v_id ||
                    (other_face != Mesh_None && visited_face[other_face])) {
                    edge_prev.first = vi;
                    edge_prev.second = vi_next;
                } else {
                    edge_next.first = vi;
                    edge_next.second = vi_next;
                }
            } else if (is_on_edge_vu[vi_next] && !is_on_edge_vu[vi]) {
                auto& two_faces = faces_of_edge[std::make_pair(vi, vi_next)];
                auto other_face = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
                if (vi_next == v_id ||
                    (other_face != Mesh_None && visited_face[other_face])) {
                    edge_prev.first = vi;
                    edge_prev.second = vi_next;
                } else {
                    edge_next.first = vi;
                    edge_next.second = vi_next;
                }
            } else if (is_on_edge_vu[vi] && is_on_edge_vu[vi_next]) {
                edge_on_vu.first = vi;
                edge_on_vu.second = vi_next;
            }
        }
        //
        if (edge_on_vu.first == Mesh_None) {
            // no edge of the face is on v->u
            // keep v_curr, update f_curr
            visited_face[f_curr] = true;
            auto& two_faces = faces_of_edge[edge_next];
            f_curr = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
        } else {
            // there is an edge of the face on v->u
            // update v_curr to be the next vert on edge v->u
            v_curr = (edge_on_vu.first == v_curr) ? edge_on_vu.second : edge_on_vu.first;
            vert_indices.push_back(v_curr);
            // update f_curr
            visited_face[f_curr] = true;
            auto& two_faces = faces_of_edge[edge_next];
            f_curr = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
        }
    }
}

void compute_passing_face_pair(const simplicial_arrangement::Arrangement<3>& tet_cut_result,
                               size_t v1,
                               size_t v2,
                               std::pair<size_t, int>& face_orient1,
                               std::pair<size_t, int>& face_orient2)
{
    // find a face incident to edge v1 -> v2
    const auto& faces = tet_cut_result.faces;
    size_t incident_face_id;
    bool found_incident_face = false;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        size_t num_vert = face.vertices.size();
        for (size_t j = 0; j < num_vert; ++j) {
            size_t vId1 = face.vertices[j];
            size_t vId2 = face.vertices[(j + 1) % num_vert];
            if ((vId1 == v1 && vId2 == v2) || (vId1 == v2 && vId2 == v1)) {
                incident_face_id = i;
                found_incident_face = true;
                break;
            }
        }
        if (found_incident_face) {
            break;
        }
    }
    // assert: found_incident_face == true
    size_t cell_id = faces[incident_face_id].positive_cell;
    const auto& cell = tet_cut_result.cells[cell_id];

    // find the two faces
    // 1. bounding the cell
    // 2. passing vertex v1 or v2
    for (auto fId : cell.faces) {
        const auto& face = faces[fId];
        bool found_v1 = false;
        bool found_v2 = false;
        for (auto vId : face.vertices) {
            if (vId == v1) {
                found_v1 = true;
            } else if (vId == v2) {
                found_v2 = true;
            }
        }
        if (found_v1 && !found_v2) {
            // the current face passes v1 but not v2
            face_orient1.first = fId;
            face_orient1.second = (face.positive_cell == cell_id) ? 1 : -1;
        } else if (!found_v1 && found_v2) {
            // the current face passes v2 but not v1
            face_orient2.first = fId;
            face_orient2.second = (face.positive_cell == cell_id) ? 1 : -1;
        }
    }
}

void compute_passing_face_pair(const simplicial_arrangement::MaterialInterface<3>& tet_cut_result,
                                  size_t v1,
                                  size_t v2,
                                  std::pair<size_t, int>& face_orient1,
                                  std::pair<size_t, int>& face_orient2)
{
    // find a face incident to edge v1 -> v2
    const auto& faces = tet_cut_result.faces;
    size_t incident_face_id;
    bool found_incident_face = false;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        size_t num_vert = face.vertices.size();
        for (size_t j = 0; j < num_vert; ++j) {
            size_t vId1 = face.vertices[j];
            size_t vId2 = face.vertices[(j + 1) % num_vert];
            if ((vId1 == v1 && vId2 == v2) || (vId1 == v2 && vId2 == v1)) {
                incident_face_id = i;
                found_incident_face = true;
                break;
            }
        }
        if (found_incident_face) {
            break;
        }
    }
    // assert: found_incident_face == true

    // map: material label -> cell id
    size_t max_material_label = 0;
    for (const auto& cell : tet_cut_result.cells) {
        if (cell.material_label > max_material_label) {
            max_material_label = cell.material_label;
        }
    }
    std::vector<size_t> cell_of_material(max_material_label + 1, Mesh_None);
    for (size_t i = 0; i < tet_cut_result.cells.size(); ++i) {
        cell_of_material[tet_cut_result.cells[i].material_label] = i;
    }
    // faces[incident_face_id].positive_material_label will be a tet boundary
    size_t cell_id = cell_of_material[faces[incident_face_id].negative_material_label];
    const auto& cell = tet_cut_result.cells[cell_id];

    // find the two faces
    // 1. bounding the cell
    // 2. passing vertex v1 or v2
    for (auto fId : cell.faces) {
        const auto& face = faces[fId];
        bool found_v1 = false;
        bool found_v2 = false;
        for (auto vId : face.vertices) {
            if (vId == v1) {
                found_v1 = true;
            } else if (vId == v2) {
                found_v2 = true;
            }
        }
        if (found_v1 && !found_v2) {
            // the current face passes v1 but not v2
            face_orient1.first = fId;
            face_orient1.second =
                    (cell_of_material[face.positive_material_label] == cell_id) ? 1 : -1;
        } else if (!found_v1 && found_v2) {
            // the current face passes v2 but not v1
            face_orient2.first = fId;
            face_orient2.second =
                    (cell_of_material[face.positive_material_label] == cell_id) ? 1 : -1;
        }
    }
}

void compute_passing_face(
        const simplicial_arrangement::Arrangement<3>& tet_cut_result, size_t v, size_t u, std::pair<size_t, int>& face_orient)
{
    // find a face incident to edge v -> u
    const auto& faces = tet_cut_result.faces;
    size_t incident_face_id;
    bool found_incident_face = false;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        size_t num_vert = face.vertices.size();
        for (size_t j = 0; j < num_vert; ++j) {
            size_t vId1 = face.vertices[j];
            size_t vId2 = face.vertices[(j + 1) % num_vert];
            if ((vId1 == v && vId2 == u) || (vId1 == u && vId2 == v)) {
                incident_face_id = i;
                found_incident_face = true;
                break;
            }
        }
        if (found_incident_face) {
            break;
        }
    }
    // assert: found_incident_face == true
    size_t cell_id = faces[incident_face_id].positive_cell;
    const auto& cell = tet_cut_result.cells[cell_id];

    // find the face
    // 1. bounding the cell
    // 2. passing vertex v
    for (auto fId : cell.faces) {
        const auto& face = faces[fId];
        bool found_v = false;
        bool found_u = false;
        for (auto vId : face.vertices) {
            if (vId == v) {
                found_v = true;
            } else if (vId == u) {
                found_u = true;
            }
        }
        if (found_v && !found_u) {
            // the current face passes v but not u
            face_orient.first = fId;
            face_orient.second = (face.positive_cell == cell_id) ? 1 : -1;
            return;
        }
    }
}

void compute_passing_face(const simplicial_arrangement::MaterialInterface<3>& tet_cut_result,
                             size_t v,
                             size_t u,
                             std::pair<size_t, int>& face_orient)
{
    // find a face incident to edge v -> u
    const auto& faces = tet_cut_result.faces;
    size_t incident_face_id;
    bool found_incident_face = false;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        size_t num_vert = face.vertices.size();
        for (size_t j = 0; j < num_vert; ++j) {
            size_t vId1 = face.vertices[j];
            size_t vId2 = face.vertices[(j + 1) % num_vert];
            if ((vId1 == v && vId2 == u) || (vId1 == u && vId2 == v)) {
                incident_face_id = i;
                found_incident_face = true;
                break;
            }
        }
        if (found_incident_face) {
            break;
        }
    }
    // assert: found_incident_face == true

    // map: material label -> cell id
    size_t max_material_label = 0;
    for (const auto& cell : tet_cut_result.cells) {
        if (cell.material_label > max_material_label) {
            max_material_label = cell.material_label;
        }
    }
    std::vector<size_t> cell_of_material(max_material_label + 1, Mesh_None);
    for (size_t i = 0; i < tet_cut_result.cells.size(); ++i) {
        cell_of_material[tet_cut_result.cells[i].material_label] = i;
    }
    // faces[incident_face_id].positive_material_label will be a tet boundary
    size_t cell_id = cell_of_material[faces[incident_face_id].negative_material_label];
    const auto& cell = tet_cut_result.cells[cell_id];

    // find the face
    // 1. bounding the cell
    // 2. passing vertex v
    for (auto fId : cell.faces) {
        const auto& face = faces[fId];
        bool found_v = false;
        bool found_u = false;
        for (auto vId : face.vertices) {
            if (vId == v) {
                found_v = true;
            } else if (vId == u) {
                found_u = true;
            }
        }
        if (found_v && !found_u) {
            // the current face passes v but not u
            face_orient.first = fId;
            face_orient.second =
                    (cell_of_material[face.positive_material_label] == cell_id) ? 1 : -1;
            return;
        }
    }
}