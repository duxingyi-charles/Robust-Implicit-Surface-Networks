//
//  csg.cpp
//  robust_implicit_networks
//
//  Created by Yiwen Ju on 6/3/24.
//


#include "csg.h"

bool load_csgTree(const std::string filename, std::vector<csg_unit>& tree){
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin)
    {
        std::cout << "function file not exist!" << std::endl;
        return false;
    }
    json tree_data;
    fin >> tree_data;
    fin.close();
    //
    size_t n_units = tree_data.size();
    tree.resize(n_units);
    for (size_t j = 0 ; j < n_units; j++){
        std::string type = tree_data[j]["type"].get<std::string>();
        std::array<int, 2> elements;
        for (int i = 0; i < 2; i ++){
            elements[i] = tree_data[j]["elements"][i].get<int>();
        }
        if (type == "Intersection"){
            tree[j] = {Intersection, elements};
        }else if (type == "Union"){
            tree[j] = {Union, elements};
        }else if (type == "Negation"){
            tree[j] = {Negation, elements};
        }
    }
    return true;
}

std::pair<std::array<double, 2>, std::vector<int>> iterTree(const std::vector<csg_unit>csgTree,const int curNode,const std::vector<std::array<double , 2>> funcInt){
    csg_unit curUnit = csgTree[curNode - 1];
    std::array<double, 2> interval, childInt1, childInt2;
    std::vector<int> af(funcInt.size(), 1), childAF1(funcInt.size(), 1), childAF2(funcInt.size(), 1);
    if (curUnit.elements[0] > 0){
        std::pair<std::array<double, 2>, std::vector<int>> child1 = iterTree(csgTree, curUnit.elements[0], funcInt);
        childInt1 = child1.first;
        childAF1 = child1.second;
    }else{
        childInt1 = funcInt[-curUnit.elements[0] - 1];
        childAF1[-curUnit.elements[0] - 1] = 0;
    }
    if (childInt1[0] * childInt1[1]>0){
        for (size_t i = 0; i < childAF1.size(); i++){
            childAF1[i] = 1;
        }
    }
    if (curUnit.operation != Negation){
        if (curUnit.elements[1] > 0){
            std::pair<std::array<double, 2>, std::vector<int>> child2 = iterTree(csgTree, curUnit.elements[1], funcInt);
            childInt2 = child2.first;
            childAF2 = child2.second;
        }else{
            childInt2 = funcInt[-curUnit.elements[1] - 1];
            childAF2[-curUnit.elements[1] - 1] = 0;
        }
    }
    if (childInt2[0] * childInt2[1]>0){
        for (size_t i = 0; i < childAF2.size(); i++){
            childAF2[i] = 1;
        }
    }
    switch (curUnit.operation){
        case Intersection:
            interval = {std::max(childInt1[0], childInt2[0]), std::max(childInt1[1], childInt2[1])};
            if(interval[0]*interval[1] <= 0){
                for (int i = 0; i < funcInt.size(); i++){
                    af[i] = childAF1[i] * childAF2[i];
                }
            }
            break;
        case Union:
            interval = {std::min(childInt1[0], childInt2[0]), std::min(childInt1[1], childInt2[1])};
            if(interval[0]*interval[1] <= 0){
                for (int i = 0; i < funcInt.size(); i++){
                    af[i] = childAF1[i] * childAF2[i];
                }
            }
            break;
        case Negation:
            interval = {std::min(-childInt1[0], -childInt1[1]), std::max(-childInt1[0], -childInt1[1])};
            if(interval[0]*interval[1] <= 0)
                af = childAF1;
            break;
        default:
            std::cout << "not a valid CSG operation" << std::endl;
    }
    return std::pair(interval, af);
}

void prune_data(const std::vector<std::array<double, 3>>& mesh_pts,
                const std::vector<PolygonFace>& mesh_faces,
                std::vector<std::vector<size_t>>& patches,
                const std::vector<size_t>& patch_function_label,
                const std::vector<Edge>& edges,
                std::vector<std::vector<size_t>>& chains,
                std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                const std::vector<std::vector<size_t>>& shells,
                const std::vector<std::vector<size_t>>& cells,
                const std::vector<bool>& cells_label){
    std::vector<std::vector<size_t>> pruned_patches;
    std::vector<std::vector<size_t>> pruned_chains;
    std::vector<std::vector<size_t>> pruned_corners;
    std::vector<std::vector<size_t>> pruned_shells;
    std::vector<std::vector<size_t>> pruned_cells;
    pruned_patches.reserve(patches.size());
    pruned_chains.reserve(chains.size());
    pruned_corners.reserve(non_manifold_edges_of_vert.size());
    
    //for each patch, find two cells from each side
    std::vector<std::pair<size_t, size_t>> patch_to_cell(patches.size(), {Mesh_None, Mesh_None});
    for (size_t cell_index = 0; cell_index < cells.size(); cell_index ++){
        for (auto shell_index : cells[cell_index]){
            for (auto half_patch : shells[shell_index]){
                size_t patch_id = half_patch / 2;
                if (half_patch % 2){
                    patch_to_cell[patch_id].second = cell_index;
                }else{
                    patch_to_cell[patch_id].first = cell_index;
                }
            }
        }
    }
    std::vector<size_t> vertex_to_patch_count(mesh_pts.size(), 0);
    for (size_t patch_index = 0; patch_index < patches.size(); patch_index++){
        if (cells_label[patch_to_cell[patch_index].first] != cells_label[patch_to_cell[patch_index].second]){
            pruned_patches.emplace_back(patches[patch_index]);
            std::vector<size_t> vertices_on_patch;
            vertices_on_patch.reserve(3 * patches[patch_index].size());
            for (auto face : patches[patch_index]){
                vertices_on_patch.emplace_back(mesh_faces[face].vert_indices[0]);
                vertices_on_patch.emplace_back(mesh_faces[face].vert_indices[1]);
                vertices_on_patch.emplace_back(mesh_faces[face].vert_indices[2]);
            }
            std::sort(vertices_on_patch.begin(), vertices_on_patch.end());
            auto last = std::unique(vertices_on_patch.begin(), vertices_on_patch.end());
            vertices_on_patch.erase(last, vertices_on_patch.end());
            for (auto& vertex_index : vertices_on_patch) {
                vertex_to_patch_count[vertex_index] ++;
            }
        }
    }
    for (auto chain : chains){
        bool active = true;
        for (auto edge_index : chain){
            if (vertex_to_patch_count[edges[edge_index].v1] < 2||vertex_to_patch_count[edges[edge_index].v2] < 2){
                active = false;
                break;
            }
        }
        if (active){
            pruned_chains.emplace_back(chain);
        }
    }
    for (size_t corner_index = 0; corner_index < non_manifold_edges_of_vert.size(); corner_index ++ ){
        if (vertex_to_patch_count[corner_index] > 2){
            pruned_corners.emplace_back(non_manifold_edges_of_vert[corner_index]);
        }
    }
    patches = pruned_patches;
    chains = pruned_chains;
    non_manifold_edges_of_vert = pruned_corners;
}
