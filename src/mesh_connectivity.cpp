//
// Created by Charles Du on 7/20/22.
//

#include "mesh_connectivity.h"
#include <absl/container/flat_hash_map.h>
#include <queue>


void compute_mesh_edges(const std::vector<PolygonFace>& mesh_faces,
                        std::vector<std::vector<size_t>>& edges_of_face,
                        std::vector<Edge>& mesh_edges)
{
    size_t max_num_edge = 0;
    for (const auto& iso_face : mesh_faces) {
        max_num_edge += iso_face.vert_indices.size();
    }
    mesh_edges.reserve(max_num_edge / 2);
    edges_of_face.reserve(mesh_faces.size());
    size_t num_iso_edge;
    // map: (v1, v2) -> iso-edge index
    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> edge_id;
    for (size_t i = 0; i < mesh_faces.size(); i++) {
        auto& face = mesh_faces[i];
        size_t num_edge = face.vert_indices.size();
        //        face.edge_indices.resize(num_edge);
        edges_of_face.emplace_back(num_edge);
        auto& face_edges = edges_of_face.back();
        for (size_t j = 0; j < num_edge; j++) {
            size_t v1 = face.vert_indices[j];
            size_t v2 = (j + 1 == num_edge) ? face.vert_indices[0] : face.vert_indices[j + 1];
            // swap if v1 > v2
            size_t tmp = v1;
            if (v1 > v2) {
                v1 = v2;
                v2 = tmp;
            }
            //
            num_iso_edge = mesh_edges.size();
            auto iter_inserted = edge_id.try_emplace(std::make_pair(v1, v2), num_iso_edge);
            if (iter_inserted.second) { // new iso-edge
                mesh_edges.emplace_back();
                mesh_edges.back().v1 = v1;
                mesh_edges.back().v2 = v2;
                mesh_edges.back().face_edge_indices.emplace_back(i, j);
                //                face.edge_indices[j] = num_iso_edge;
                face_edges[j] = num_iso_edge;
            } else { // existing iso-edge
                size_t eId = iter_inserted.first->second;
                mesh_edges[eId].face_edge_indices.emplace_back(i, j);
                //                face.edge_indices[j] = eId;
                face_edges[j] = eId;
            }
        }
    }
}

void compute_patches(const std::vector<std::vector<size_t>>& edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<PolygonFace>& mesh_faces,
                     std::vector<std::vector<size_t>>& patches,
                     std::vector<size_t>& patch_function_label)
{
    std::vector<bool> visisted_face(edges_of_face.size(), false);
    for (size_t i = 0; i < edges_of_face.size(); i++) {
        if (!visisted_face[i]) {
            // new patch
            patches.emplace_back();
            auto& patch = patches.back();
            std::queue<size_t> Q;
            Q.push(i);
            patch.emplace_back(i);
            visisted_face[i] = true;
            size_t patch_label = mesh_faces[Q.front()].func_index;
            patch_function_label.emplace_back(patch_label);
            while (!Q.empty()) {
                auto fId = Q.front();
                Q.pop();
                for (size_t eId : edges_of_face[fId]) {
                    if (mesh_edges[eId].face_edge_indices.size() == 2) { // manifold edge
                        size_t other_fId = (mesh_edges[eId].face_edge_indices[0].first == fId)
                                           ? mesh_edges[eId].face_edge_indices[1].first
                                           : mesh_edges[eId].face_edge_indices[0].first;
                        if (!visisted_face[other_fId]) {
                            Q.push(other_fId);
                            patch.emplace_back(other_fId);
                            visisted_face[other_fId] = true;
                        }
                    }
                }
            }
        }
    }
}

bool check_patch_label(const std::vector<std::vector<size_t>>& edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<PolygonFace>& mesh_faces,
                     const std::vector<IsoVert> iso_verts,
                     std::vector<std::vector<size_t>>& patches,
                       std::vector<size_t>& patch_function_label)
{
    bool check_result = true;
    std::vector<bool> visisted_face(edges_of_face.size(), false);
    for (size_t i = 0; i < patches.size(); i++) {
        absl::flat_hash_map<size_t, size_t> func_to_num; //adding a hash map that finds the number of vertices that are associated to each function
        absl::flat_hash_map<size_t, bool> edge_traversed;
        size_t patch_label = mesh_faces[patches[i][0]].func_index;
        for (size_t fId : patches[i]){
            for (size_t eId : edges_of_face[fId]) {
                if (edge_traversed[eId]){
                    continue;
                }else{
                    edge_traversed[eId] = true;
                }
                for (size_t funcIter = 0; funcIter < iso_verts[mesh_edges[eId].v1].func_indices.size(); funcIter++){
                    size_t func_index_v1 = iso_verts[mesh_edges[eId].v1].func_indices[funcIter];
                    size_t func_index_v2 = iso_verts[mesh_edges[eId].v2].func_indices[funcIter];
                    if (func_index_v1 != Mesh_None){
                        if (func_to_num.contains(func_index_v1)){
                            func_to_num[func_index_v1] ++;
                        }else{
                            func_to_num[func_index_v1] = 1;
                        }
                    }
                    if (func_index_v2 != Mesh_None){
                        if (func_to_num.contains(func_index_v2)){
                            func_to_num[func_index_v2] ++;
                        }else{
                            func_to_num[func_index_v2] = 1;
                        }
                    }
                }
            }
            auto func_label = func_to_num.begin();
            for (auto it = func_to_num.begin(); it != func_to_num.end(); ++it) {
                if (it->second > func_label->second)
                    func_label = it;
            }
            for (auto it = func_to_num.begin(); it != func_to_num.end(); ++it) {
                if (it->second == func_label->second && it->first != func_label->first)
                {
                    throw std::runtime_error("ERROR: Same Common Labels");
                    check_result = false;
                }
            }
            if (func_label->first != patch_label){
                throw std::runtime_error("ERROR: Inconsistent Label between Vertices and Faces");
                check_result = false;
            }
        }
    }
    return check_result;
}

void compute_patches(const std::vector<std::vector<size_t>>& edges_of_face,
                     const std::vector<Edge>& mesh_edges,
                     const std::vector<MI_Vert> MI_verts,
                     std::vector<std::vector<size_t>>& patches)
{
    std::vector<bool> visisted_face(edges_of_face.size(), false);
    for (size_t i = 0; i < edges_of_face.size(); i++) {
        if (!visisted_face[i]) {
            // new patch
            patches.emplace_back();
            auto& patch = patches.back();
            std::queue<size_t> Q;
            Q.push(i);
            patch.emplace_back(i);
            visisted_face[i] = true;
            while (!Q.empty()) {
                auto fId = Q.front();
                Q.pop();
                for (size_t eId : edges_of_face[fId]) {
                    if (mesh_edges[eId].face_edge_indices.size() == 2) { // manifold edge
                        size_t other_fId = (mesh_edges[eId].face_edge_indices[0].first == fId)
                                           ? mesh_edges[eId].face_edge_indices[1].first
                                           : mesh_edges[eId].face_edge_indices[0].first;
                        if (!visisted_face[other_fId]) {
                            Q.push(other_fId);
                            patch.emplace_back(other_fId);
                            visisted_face[other_fId] = true;
                        }
                    }
                }
            }
        }
    }
}


// group non-manifold iso-edges into chains
void compute_chains(const std::vector<Edge>& mesh_edges,
                    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                    std::vector<std::vector<size_t>>& chains)
{
    std::vector<bool> visited_edge(mesh_edges.size(), false);
    for (size_t i = 0; i < mesh_edges.size(); i++) {
        if (!visited_edge[i] && mesh_edges[i].face_edge_indices.size() > 2) {
            // unvisited non-manifold iso-edge (not a boundary edge)
            // new chain
            chains.emplace_back();
            auto& chain = chains.back();
            std::queue<size_t> Q;
            Q.push(i);
            chain.emplace_back(i);
            visited_edge[i] = true;
            while (!Q.empty()) {
                auto eId = Q.front();
                Q.pop();
                // v1
                size_t v = mesh_edges[eId].v1;
                if (non_manifold_edges_of_vert[v].size() == 2) {
                    size_t other_eId = (non_manifold_edges_of_vert[v][0] == eId)
                                       ? non_manifold_edges_of_vert[v][1]
                                       : non_manifold_edges_of_vert[v][0];
                    if (!visited_edge[other_eId]) {
                        Q.push(other_eId);
                        chain.emplace_back(other_eId);
                        visited_edge[other_eId] = true;
                    }
                }
                // v2
                v = mesh_edges[eId].v2;
                if (non_manifold_edges_of_vert[v].size() == 2) {
                    size_t other_eId = (non_manifold_edges_of_vert[v][0] == eId)
                                       ? non_manifold_edges_of_vert[v][1]
                                       : non_manifold_edges_of_vert[v][0];
                    if (!visited_edge[other_eId]) {
                        Q.push(other_eId);
                        chain.emplace_back(other_eId);
                        visited_edge[other_eId] = true;
                    }
                }
            }
        }
    }
}

void compute_shells_and_components(size_t num_patch,
                                   const std::vector<std::vector<std::pair<std::pair<size_t, int>, std::pair<size_t, int>>>>&
                                   half_patch_pair_list,
                                   std::vector<std::vector<size_t>>& shells,
                                   std::vector<size_t>& shell_of_half_patch,
                                   std::vector<std::vector<size_t>>& components,
                                   std::vector<size_t>& component_of_patch)
{
    // (patch i, 1) <--> 2i,  (patch i, -1) <--> 2i+1
    // compute half-patch adjacency list
    std::vector<std::vector<size_t>> half_patch_adj_list(2 * num_patch);
    for (const auto& half_patch_pairs : half_patch_pair_list) {
        for (size_t i = 0; i < half_patch_pairs.size(); i++) {
            const auto& hp1 = half_patch_pairs[i].first;
            const auto& hp2 = half_patch_pairs[i].second;
            // half-patch index of hp1
            size_t hp_Id1 = (hp1.second == 1) ? 2 * hp1.first : (2 * hp1.first + 1);
            // half-patch index of hp2
            size_t hp_Id2 = (hp2.second == 1) ? 2 * hp2.first : (2 * hp2.first + 1);
            half_patch_adj_list[hp_Id1].push_back(hp_Id2);
            half_patch_adj_list[hp_Id2].push_back(hp_Id1);
        }
    }
    // find connected component of half-patch adjacency graph
    // each component is a shell
    std::vector<bool> visited_half_patch(2 * num_patch, false);
    shells.clear();
    shell_of_half_patch.resize(2 * num_patch);
    for (size_t i = 0; i < 2 * num_patch; i++) {
        if (!visited_half_patch[i]) {
            // create new component
            size_t shell_Id = shells.size();
            shells.emplace_back();
            auto& shell = shells.back();
            std::queue<size_t> Q;
            Q.push(i);
            shell.push_back(i);
            shell_of_half_patch[i] = shell_Id;
            visited_half_patch[i] = true;
            while (!Q.empty()) {
                auto half_patch = Q.front();
                Q.pop();
                for (auto hp : half_patch_adj_list[half_patch]) {
                    if (!visited_half_patch[hp]) {
                        shell.push_back(hp);
                        shell_of_half_patch[hp] = shell_Id;
                        Q.push(hp);
                        visited_half_patch[hp] = true;
                    }
                }
            }
        }
    }
    // find connected component of patch-adjacency graph
    // each component is an iso-surface component
    std::vector<bool> visited_patch(num_patch, false);
    components.clear();
    component_of_patch.resize(num_patch);
    for (size_t i = 0; i < num_patch; i++) {
        if (!visited_patch[i]) {
            // create new component
            size_t component_Id = components.size();
            components.emplace_back();
            auto& component = components.back();
            std::queue<size_t> Q;
            Q.push(i);
            component.push_back(i);
            component_of_patch[i] = component_Id;
            visited_patch[i] = true;
            while (!Q.empty()) {
                auto patch = Q.front();
                Q.pop();
                // 1-side of patch
                for (auto hp : half_patch_adj_list[2 * patch]) {
                    if (!visited_patch[hp / 2]) {
                        auto p = hp / 2;
                        component.push_back(p);
                        component_of_patch[p] = component_Id;
                        Q.push(p);
                        visited_patch[p] = true;
                    }
                }
                // -1-side of patch
                for (auto hp : half_patch_adj_list[2 * patch + 1]) {
                    if (!visited_patch[hp / 2]) {
                        auto p = hp / 2;
                        component.push_back(p);
                        component_of_patch[p] = component_Id;
                        Q.push(p);
                        visited_patch[p] = true;
                    }
                }
            }
        }
    }
    // get shells as list of patch indices
    //    for (auto& shell : shells) {
    //        for (auto& pId : shell) {
    //            pId /= 2; // get patch index of half-patch
    //        }
    //    }
}


void compute_arrangement_cells(size_t num_shell,
                               const std::vector<std::pair<size_t, size_t>>& shell_links,
                               std::vector<std::vector<size_t>>& arrangement_cells)
{
    // build shell adjacency list
    size_t sink_shell = num_shell;
    absl::flat_hash_map<size_t, std::vector<size_t>> adjacent_shells;
    for (const auto& link : shell_links) {
        if (link.first == Mesh_None) {
            adjacent_shells[sink_shell].push_back(link.second);
            adjacent_shells[link.second].push_back(sink_shell);
        } else if (link.second == Mesh_None) {
            adjacent_shells[sink_shell].push_back(link.first);
            adjacent_shells[link.first].push_back(sink_shell);
        } else {
            adjacent_shells[link.first].push_back(link.second);
            adjacent_shells[link.second].push_back(link.first);
        }
    }

    // find connected components of shell adjacency graph
    // each component is an arrangement cells
    std::vector<bool> visited_shell(num_shell + 1, false);
    //    arrangement_cells.clear();
    for (size_t i = 0; i < num_shell + 1; ++i) {
        if (!visited_shell[i]) {
            // create new component
            arrangement_cells.emplace_back();
            auto& arr_cell = arrangement_cells.back();
            std::queue<size_t> Q;
            Q.push(i);
            arr_cell.push_back(i);
            visited_shell[i] = true;
            while (!Q.empty()) {
                auto shell_id = Q.front();
                Q.pop();
                for (auto s : adjacent_shells[shell_id]) {
                    if (!visited_shell[s]) {
                        arr_cell.push_back(s);
                        Q.push(s);
                        visited_shell[s] = true;
                    }
                }
            }
        }
    }

    // remove sink shell from arrangement cells
    std::vector<size_t> sink_free_shell_list;
    for (auto& arr_cell : arrangement_cells) {
        sink_free_shell_list.clear();
        for (auto s : arr_cell) {
            if (s < num_shell) {
                sink_free_shell_list.push_back(s);
            }
        }
        arr_cell = sink_free_shell_list;
    }
}
