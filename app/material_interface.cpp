//
// Created by Charles Du on 1/15/22.
//

#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/material_interface.h>

#include <chrono>
#include <iostream>
#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <absl/container/flat_hash_map.h>

#include "robust_implicit_networks.h"
#include "implicit_functions.h"
#include "ScopedTimer.h"

typedef std::chrono::duration<double> Time_duration;

using namespace simplicial_arrangement;


int main(int argc, const char* argv[])
{
    struct
    {
        std::string config_file;
        bool timing_only = false;
        bool robust_test = false;
    } args;
    CLI::App app{"Material Interface Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_option("-T,--timing-only", args.timing_only, "Record timing without saving results");
    app.add_option("-R,--robust-test",args.robust_test, "Perform robustness test");
    CLI11_PARSE(app, argc, argv);

    // parse configure file
    std::string tet_mesh_file;
    std::string material_file;
    std::string output_dir;
    bool use_lookup = true;
    bool use_3func_lookup = true;
    bool use_topo_ray_shooting = true;
    parse_config_file_MI(args.config_file, tet_mesh_file, material_file, output_dir,
                         use_lookup, use_3func_lookup, use_topo_ray_shooting);
    if (use_lookup) {
        // load lookup table
        std::cout << "load table ..." << std::endl;
        bool loaded = load_lookup_table(simplicial_arrangement::MATERIAL_INTERFACE);
        if (loaded) {
            std::cout << "loading finished." << std::endl;
        } else {
            std::cout << "loading failed." << std::endl;
            return -1;
        }
    } else {
        use_3func_lookup = false;
    }

    // record timings
    std::vector<std::string> timing_labels;
    std::vector<double> timings;

    // record stats
    std::vector<std::string> stats_labels;
    std::vector<size_t> stats;

    // load tet mesh
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    size_t n_tets = tets.size();
    size_t n_pts = pts.size();
    std::cout << "tet mesh: " << n_pts << " verts, " << n_tets << " tets." << std::endl;
    stats_labels.emplace_back("num_pts");
    stats.push_back(n_pts);
    stats_labels.emplace_back("num_tets");
    stats.push_back(n_tets);

    // load implicit functions and compute function values at vertices
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
    if (load_functions(material_file, pts, funcVals)) {
        std::cout << "function loading finished." << std::endl;
    } else {
        std::cout << "function loading failed." << std::endl;
        return -2;
    }
    size_t n_func = funcVals.cols();
    std::cout << "n_func = " << n_func << std::endl;

    // highest material at vertices
    std::vector<size_t> highest_material_at_vert;
    // degenerate vertex: more than one highest material, i.e. material interface passes the vertex
    std::vector<bool> is_degenerate_vertex;
    bool found_degenerate_vertex = false;
    absl::flat_hash_map<size_t, std::vector<size_t>> highest_materials_at_vert;
    {
        timing_labels.emplace_back("highest func");
        ScopedTimer<> timer("highest func");
        is_degenerate_vertex.resize(n_pts, false);
        highest_material_at_vert.reserve(n_pts);
        for (Eigen::Index i = 0; i < n_pts; i++) {
            double max = funcVals(i, 0);
            size_t max_id = 0;
            size_t max_count = 1;
            for (Eigen::Index j = 1; j < n_func; j++) {
                if (funcVals(i,j) > max) {
                    max = funcVals(i,j);
                    max_id = j;
                    max_count = 1;
                } else if (funcVals(i,j) == max) {
                    ++max_count;
                }
            }
            highest_material_at_vert.push_back(max_id);
            //
            if (max_count > 1) {
                is_degenerate_vertex[i] = true;
                found_degenerate_vertex = true;
                auto& materials = highest_materials_at_vert[i];
                materials.reserve(max_count);
                for (Eigen::Index j = 0; j < n_func; j++) {
                    if (funcVals(i,j) == max) {
                        materials.push_back(j);
                    }
                }
            }
        }
        timings.push_back(timer.toc());
    }

    // filter active materials in each tet
    // a tet is non-empty if there are material interface in it
    size_t num_intersecting_tet = 0;
    std::vector<size_t> material_in_tet; // active material indices in CRS vector format
    std::vector<size_t> start_index_of_tet;
    {
        timing_labels.emplace_back("filter");
        ScopedTimer<> timer("filter");
        material_in_tet.reserve(n_tets);
        start_index_of_tet.reserve(n_tets + 1);
        start_index_of_tet.push_back(0);
        std::set<size_t> materials;
        std::array<double, 4> min_h;
        for (size_t i = 0; i < n_tets; ++i) {
            const auto& tet = tets[i];
            // find high materials
            materials.clear();
            for (size_t j = 0; j < 4; ++j) {
                if (is_degenerate_vertex[tet[j]]) {
                    const auto& ms = highest_materials_at_vert[tet[j]];
                    materials.insert(ms.begin(), ms.end());
                } else {
                    materials.insert(highest_material_at_vert[tet[j]]);
                }
            }
            // if only one high material, there is no material interface
            if (materials.size() < 2) {  // no material interface
                start_index_of_tet.push_back(material_in_tet.size());
                continue;
            }
            // find min of high materials
            min_h.fill(std::numeric_limits<double>::max());
            for(auto it = materials.begin(); it != materials.end(); it++) {
                for (size_t j = 0; j < 4; ++j) {
                    if (funcVals(tet[j], *it) < min_h[j]) {
                        min_h[j] = funcVals(tet[j], *it);
                    }
                }
            }
            // find materials greater than at least two mins of high materials
            size_t greater_count;
            for (size_t j = 0; j < n_func; ++j) {
                greater_count = 0;
                for (size_t k = 0; k < 4; ++k) {
                    if (funcVals(tet[k],j) > min_h[k]) {
                        ++greater_count;
                    }
                }
                if (greater_count > 1) {
                    materials.insert(j);
                }
            }
            //
            ++num_intersecting_tet;
            material_in_tet.insert(material_in_tet.end(), materials.begin(), materials.end());
            start_index_of_tet.push_back(material_in_tet.size());
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num_intersecting_tet = " << num_intersecting_tet << std::endl;
    stats_labels.emplace_back("num_intersecting_tet");
    stats.push_back(num_intersecting_tet);

    // compute material interface in each tet
    std::vector<MaterialInterface<3>> cut_results;
    std::vector<size_t> cut_result_index;
    //
    Time_duration time_2_func = Time_duration::zero();
    Time_duration time_3_func = Time_duration::zero();
    Time_duration time_more_func = Time_duration::zero();
    size_t num_2_func = 0;
    size_t num_3_func = 0;
    size_t num_more_func = 0;
    // failure types for robustness test
    bool is_type2 = false; // crash in the normal order
    bool is_type3 = false; // succeed in the normal order, crash in the reverse order
    bool is_type1 = false; // no crash for either order, but results are inconsistent
    //
    {
        ScopedTimer<> timer("material interface in tets");
        if (args.robust_test) {
            cut_results.reserve(num_intersecting_tet);
            cut_result_index.reserve(n_tets);
            std::vector<MaterialInterface<3>> cut_results2; // cut_results when reverse function order
            size_t start_index;
            size_t num_func;
            std::vector<Material<double, 3>> materials;
            std::vector<Material<double, 3>> materials_reverse; // materials in reverse order
            materials.reserve(3);
            materials_reverse.reserve(3);
            for (size_t i = 0; i < tets.size(); i++) {
                start_index = start_index_of_tet[i];
                num_func = start_index_of_tet[i + 1] - start_index;
                if (num_func == 0) {
                    cut_result_index.push_back(MaterialInterface<3>::None);
                    continue;
                }
                std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                //
                size_t v1 = tets[i][0];
                size_t v2 = tets[i][1];
                size_t v3 = tets[i][2];
                size_t v4 = tets[i][3];
                materials.clear();
                for (size_t j = 0; j < num_func; j++) {
                    size_t f_id = material_in_tet[start_index + j];
                    materials.emplace_back();
                    auto& material = materials.back();
                    material[0] = funcVals(v1, f_id);
                    material[1] = funcVals(v2, f_id);
                    material[2] = funcVals(v3, f_id);
                    material[3] = funcVals(v4, f_id);
                }
                // reverse material order
                materials_reverse.clear();
                for (size_t j = 0; j < num_func; ++j) {
                    materials_reverse.emplace_back(materials[num_func -1 -j]);
                }
                //
                bool crashed = false;
                if (!use_3func_lookup && num_func == 3) {
                    cut_result_index.push_back(cut_results.size());
                    disable_lookup_table();
                    try {
                        cut_results.emplace_back(std::move(compute_material_interface(materials)));
                    } catch (std::runtime_error& e) {
                        crashed = true;
                        is_type2 = true;
                        break;
                    }
                    if (!crashed) {
                        try {
                            cut_results2.emplace_back(std::move(compute_material_interface(materials_reverse)));
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            is_type3 = true;
                        }
                        if (!crashed) {
                            const auto& last_cut_result = cut_results.back();
                            const auto& last_cut_result2 = cut_results2.back();
                            if (last_cut_result.vertices.size() != last_cut_result2.vertices.size() ||
                                last_cut_result.faces.size() != last_cut_result2.faces.size() ||
                                last_cut_result.cells.size() != last_cut_result2.cells.size()) {
                                // inconsistent results
                                is_type1 = true;
                            }
                        }
                    }
                    enable_lookup_table();
                } else {  // 3-function lookup table enabled
                    cut_result_index.push_back(cut_results.size());
                    try {
                        cut_results.emplace_back(std::move(compute_material_interface(materials)));
                    } catch (std::runtime_error& e) {
                        crashed = true;
                        is_type2 = true;
                        break;
                    }
                    if (!crashed) {
                        try {
                            cut_results2.emplace_back(std::move(compute_material_interface(materials_reverse)));
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            is_type3 = true;
                        }
                        if (!crashed) {
                            const auto& last_cut_result = cut_results.back();
                            const auto& last_cut_result2 = cut_results2.back();
                            if (last_cut_result.vertices.size() != last_cut_result2.vertices.size() ||
                                last_cut_result.faces.size() != last_cut_result2.faces.size() ||
                                last_cut_result.cells.size() != last_cut_result2.cells.size()) {
                                // inconsistent results
                                is_type1 = true;
                            }
                        }
                    }
                }
                //
                std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                switch (num_func) {
                    case 2:
                        time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_2_func;
                        break;
                    case 3:
                        time_3_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_3_func;
                        break;
                    default:
                        time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_more_func;
                        break;
                }
            }
        } else {  // not performing robustness test
            cut_results.reserve(num_intersecting_tet);
            cut_result_index.reserve(n_tets);
            size_t start_index;
            size_t num_func;
            std::vector<Material<double, 3>> materials;
            materials.reserve(3);
            try {
                for (size_t i = 0; i < tets.size(); i++) {
                    start_index = start_index_of_tet[i];
                    num_func = start_index_of_tet[i + 1] - start_index;
                    if (num_func == 0) {
                        cut_result_index.push_back(MaterialInterface<3>::None);
                        continue;
                    }
                    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                    //
                    size_t v1 = tets[i][0];
                    size_t v2 = tets[i][1];
                    size_t v3 = tets[i][2];
                    size_t v4 = tets[i][3];
                    materials.clear();
                    for (size_t j = 0; j < num_func; j++) {
                        size_t f_id = material_in_tet[start_index + j];
                        materials.emplace_back();
                        auto& material = materials.back();
                        material[0] = funcVals(v1, f_id);
                        material[1] = funcVals(v2, f_id);
                        material[2] = funcVals(v3, f_id);
                        material[3] = funcVals(v4, f_id);
                    }
                    //
                    if (!use_3func_lookup && num_func == 3) {
                        cut_result_index.push_back(cut_results.size());
                        disable_lookup_table();
                        cut_results.emplace_back(std::move(compute_material_interface(materials)));
                        enable_lookup_table();
                    } else {
                        cut_result_index.push_back(cut_results.size());
                        cut_results.emplace_back(std::move(compute_material_interface(materials)));
                    }

                    //
                    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                    switch (num_func) {
                        case 2:
                            time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_2_func;
                            break;
                        case 3:
                            time_3_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_3_func;
                            break;
                        default:
                            time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_more_func;
                            break;
                    }
                }
            } catch (std::runtime_error& e) {
                std::cout << e.what() << std::endl;
                return -1;
            }
        }
        //
        timing_labels.emplace_back("MI(other)");
        timings.push_back(
                timer.toc() - time_2_func.count() - time_3_func.count() - time_more_func.count());
    }
    timing_labels.emplace_back("MI(2 func)");
    timings.push_back(time_2_func.count());
    timing_labels.emplace_back("MI(3 func)");
    timings.push_back(time_3_func.count());
    timing_labels.emplace_back("MI(>=4 func)");
    timings.push_back(time_more_func.count());
    std::cout << " -- [MI(2 func)]: " << time_2_func.count() << " s" << std::endl;
    std::cout << " -- [MI(3 func)]: " << time_3_func.count() << " s" << std::endl;
    std::cout << " -- [MI(>=4 func)]: " << time_more_func.count() << " s" << std::endl;
    //
    stats_labels.emplace_back("num_2_func");
    stats.push_back(num_2_func);
    stats_labels.emplace_back("num_3_func");
    stats.push_back(num_3_func);
    stats_labels.emplace_back("num_more_func");
    stats.push_back(num_more_func);

    if (args.robust_test) {
        if (is_type2) {
            std::cout << "type 2 failure (crash in the normal order)." << std::endl;
        } else if (is_type3) {
            std::cout << "type 3 failure (crash in the reverse order)." << std::endl;
        } else if (is_type1) {
            std::cout << "type 1 failure (inconsistency)." << std::endl;
        } else {
            std::cout << "success." << std::endl;
        }
        return 0;
    }


    // extract material interface mesh
    std::vector<MI_Vert> MI_verts;
    std::vector<PolygonFace> MI_faces;
    // the following data are only needed when we use
    // the baseline nesting algorithm
    // (group simplicial cells into material cells)
    std::vector<long long> global_vId_of_tet_vert;
    std::vector<size_t> global_vId_start_index_of_tet;
    std::vector<size_t> MI_fId_of_tet_face;
    std::vector<size_t> MI_fId_start_index_of_tet;
    {
        timing_labels.emplace_back("extract mesh");
        ScopedTimer<> timer("extract mesh");
        if (use_topo_ray_shooting) {
            extract_MI_mesh(num_2_func,
                                 num_3_func,
                                 num_more_func,
                                 cut_results,
                                 cut_result_index,
                                 material_in_tet,
                                 start_index_of_tet,
                                 tets,
                                 MI_verts,
                                 MI_faces);
        } else { // nesting algorithm: group simplicial cells into material cells
            extract_MI_mesh(num_2_func,
                            num_3_func,
                            num_more_func,
                            cut_results,
                            cut_result_index,
                            material_in_tet,
                            start_index_of_tet,
                            tets,
                            MI_verts,
                            MI_faces,
                            global_vId_of_tet_vert,
                            global_vId_start_index_of_tet,
                            MI_fId_of_tet_face,
                            MI_fId_start_index_of_tet);
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num MI-vertices = " << MI_verts.size() << std::endl;
    std::cout << "num MI-faces = " << MI_faces.size() << std::endl;
    stats_labels.emplace_back("num_MI_verts");
    stats.push_back(MI_verts.size());
    stats_labels.emplace_back("num_MI_faces");
    stats.push_back(MI_faces.size());

    // compute xyz of MI-vertices
    std::vector<std::array<double, 3>> MI_pts;
    {
        timing_labels.emplace_back("compute xyz");
        ScopedTimer<> timer("compute xyz");
        compute_MI_vert_xyz(MI_verts, funcVals, pts, MI_pts);
        timings.push_back(timer.toc());
    }

    //  compute iso-edges and edge-face connectivity
    std::vector<Edge> MI_edges;
    std::vector<std::vector<size_t>> edges_of_MI_face;
    {
        timing_labels.emplace_back("edge-face connectivity");
        ScopedTimer<> timer("edge-face connectivity");
        compute_mesh_edges(MI_faces, edges_of_MI_face, MI_edges);
        timings.push_back(timer.toc());
    }
    std::cout << "num MI-edges = " << MI_edges.size() << std::endl;
    stats_labels.emplace_back("num_MI_edges");
    stats.push_back(MI_edges.size());

    // group iso-faces into patches
    std::vector<std::vector<size_t>> patches;
    {
        timing_labels.emplace_back("patches");
        ScopedTimer<> timer("patches");
        compute_patches(edges_of_MI_face, MI_edges, patches);
        timings.push_back(timer.toc());
    }
    std::cout << "num patches = " << patches.size() << std::endl;
    stats_labels.emplace_back("num_patches");
    stats.push_back(patches.size());

    // compute map: MI-face Id --> patch Id
    std::vector<size_t> patch_of_face;
    {
        timing_labels.emplace_back("face-patch map");
        ScopedTimer<> timer("face-patch map");
        patch_of_face.resize(MI_faces.size());
        for (size_t i = 0; i < patches.size(); i++) {
            for (const auto& fId : patches[i]) {
                patch_of_face[fId] = i;
            }
        }
        timings.push_back(timer.toc());
    }

    // group non-manifold iso-edges into chains
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
    std::vector<std::vector<size_t>> chains;
    {
        timing_labels.emplace_back("chains");
        ScopedTimer<> timer("chains");
        non_manifold_edges_of_vert.resize(MI_pts.size());
        // get incident non-manifold edges for MI-vertices
        for (size_t i = 0; i < MI_edges.size(); i++) {
            if (MI_edges[i].face_edge_indices.size() >
                2) { // non-manifold edge (not a boundary edge)
                // there is only one patch incident to a boundary edge,
                // so there is no need to figure out the "order" of patches around a boundary
                // edge
                non_manifold_edges_of_vert[MI_edges[i].v1].push_back(i);
                non_manifold_edges_of_vert[MI_edges[i].v2].push_back(i);
            }
        }
        // group non-manifold iso-edges into chains
        compute_chains(MI_edges, non_manifold_edges_of_vert, chains);
        timings.push_back(timer.toc());
    }
    std::cout << "num chains = " << chains.size() << std::endl;
    stats_labels.emplace_back("num_chains");
    stats.push_back(chains.size());

    // compute incident tets for degenerate vertices
    absl::flat_hash_map<size_t, std::vector<size_t>> incident_tets;
    {
        timing_labels.emplace_back("vert-tet connectivity");
        ScopedTimer<> timer("vert-tet connectivity(degenerate vert only)");
        if (found_degenerate_vertex) {
            for (size_t i = 0; i < n_tets; ++i) {
                const auto& tet = tets[i];
                for (size_t j = 0; j < 4; ++j) {
                    if (is_degenerate_vertex[tet[j]]) {
                        incident_tets[tet[j]].push_back(i);
                    }
                }
            }
        }
        timings.push_back(timer.toc());
    }
    std::cout << "incident_tets.size() = " << incident_tets.size() << std::endl;

    // compute order of patches around chains
    // pair<size_t, int> : pair (iso-face index, iso-face orientation)
//    std::vector<std::vector<std::pair<size_t, int>>> half_faces_list;
//    std::vector<std::vector<std::pair<size_t, int>>> half_patch_list;
//
    // pair<size_t, int> : pair (iso-face index, iso-face orientation)
    std::vector<std::vector<std::pair<std::pair<size_t,int>,std::pair<size_t,int>>>> half_face_pair_list;
    std::vector<std::vector<std::pair<std::pair<size_t,int>,std::pair<size_t,int>>>> half_patch_pair_list;
    {
        timing_labels.emplace_back("order patches around chains");
        ScopedTimer<> timer("order patches around chains");
        half_face_pair_list.resize(chains.size());
        half_patch_pair_list.resize(chains.size());
        // pick representative iso-edge from each chain
        std::vector<size_t> chain_representatives(chains.size());
        for (size_t i = 0; i < chains.size(); i++) {
            chain_representatives[i] = chains[i][0];
        }
        // order iso-faces incident to each representative iso-edge
        for (size_t i = 0; i < chain_representatives.size(); i++) {
            const auto& MI_edge = MI_edges[chain_representatives[i]];
            // with degeneracy handling
            try {
                compute_face_order(MI_edge,
                                   MI_faces,
                                   cut_results,
                                   cut_result_index,
                                   incident_tets,
                                   half_face_pair_list[i]);
            } catch (std::exception& e) {
                std::cout << "order patches failed: " << e.what() << std::endl;
                return -1;
            }
        }
        // replace MI-face indices by patch indices
        for (size_t i = 0; i < half_face_pair_list.size(); i++) {
            half_patch_pair_list[i].resize(half_face_pair_list[i].size());
            for (size_t j = 0; j < half_face_pair_list[i].size(); j++) {
                half_patch_pair_list[i][j] = std::make_pair(
                        std::make_pair(patch_of_face[half_face_pair_list[i][j].first.first],
                                       half_face_pair_list[i][j].first.second),
                        std::make_pair(patch_of_face[half_face_pair_list[i][j].second.first],
                                       half_face_pair_list[i][j].second.second)
                );
            }
        }
        timings.push_back(timer.toc());
    }

    // group patches into shells and components
    // each shell is represented as a list of half-patch indices
    // each component is represented as a list of patch indices
    std::vector<std::vector<size_t>> shells;
    std::vector<size_t> shell_of_half_patch;
    std::vector<std::vector<size_t>> components;
    std::vector<size_t> component_of_patch;
    {
        timing_labels.emplace_back("shells and components");
        ScopedTimer<> timer("shells and components");
        compute_shells_and_components(patches.size(),
                                      half_patch_pair_list,
                                      shells,
                                      shell_of_half_patch,
                                      components,
                                      component_of_patch);
        timings.push_back(timer.toc());
    }
    std::cout << "num shells = " << shells.size() << std::endl;
    std::cout << "num components = " << components.size() << std::endl;
    stats_labels.emplace_back("num_shells");
    stats.push_back(shells.size());
    stats_labels.emplace_back("num_components");
    stats.push_back(components.size());

    // resolve nesting order, compute arrangement cells
    // a material cell is represented by a list of bounding shells
    std::vector<std::vector<size_t>> material_cells;
    {
        ScopedTimer<> timer("material cells");
        if (components.size() < 2) { // no nesting problem, each shell is an material cell
            material_cells.reserve(shells.size());
            for (size_t i = 0; i < shells.size(); ++i) {
                material_cells.emplace_back(1);
                material_cells.back()[0] = i;
            }
        } else { // resolve nesting order
            if (use_topo_ray_shooting) {
                timing_labels.emplace_back("matCells(ray shooting)");
                ScopedTimer<> timer("material cells: topo ray shooting");
                topo_ray_shooting(pts, tets,
                                  cut_results, cut_result_index,
                                  MI_verts, MI_faces,
                                  patches, patch_of_face,
                                  shells, shell_of_half_patch,
                                  components, component_of_patch,
                                  material_cells);
                timings.push_back(timer.toc());
            }
            else { // group simplicial cells into arrangement cells
                // ------------------- build adjacency graph of simplicial cells ---------------
                // pair (tetId, tet_cell_Id)
                std::vector<std::pair<size_t, size_t>> tet_cell_of_simp_cell;
                // simplicial half-face info
                // if the face is on material interface, let s be its shell index, record (-s-1)
                // if the face is not on material interface, record the simp_cell index on the opposite side
                std::vector<long long> simp_half_face_info;
                std::vector<size_t> simp_hFace_start_index;
                {
                    timing_labels.emplace_back("matCells(build simpCell graph)");
                    ScopedTimer<> timer("material cells: build simplicial cells graph");
                    build_simplicial_cell_adjacency(tets,
                                                    cut_results, cut_result_index,
                                                    global_vId_of_tet_vert, global_vId_start_index_of_tet,
                                                    MI_fId_of_tet_face, MI_fId_start_index_of_tet,
                                                    patch_of_face, shell_of_half_patch,
                                                    tet_cell_of_simp_cell,
                                                    simp_half_face_info,
                                                    simp_hFace_start_index
                    );
                    timings.push_back(timer.toc());
                }
                {
                    // ------------------- group simplicial cells into arrangement cells
                    // ---------------
                    timing_labels.emplace_back("matCells(group simpCells into matCells)");
                    ScopedTimer<> timer("material cells: group simplicial cells");
                    compute_simplicial_cell_connected_components(
                            tet_cell_of_simp_cell, simp_half_face_info, simp_hFace_start_index,
                            material_cells);
                    timings.push_back(timer.toc());
                }
            }
        }
        timings.push_back(timer.toc());
        timing_labels.emplace_back("material cells");
    }
    std::cout << "num_cells = " << material_cells.size() << std::endl;
    stats_labels.emplace_back("num_cells");
    stats.push_back(material_cells.size());

    if (components.size() > 1) {
        timing_labels.emplace_back("matCells(other)");
        size_t num_timings = timings.size();
        if (use_topo_ray_shooting) {
            timings.push_back(timings[num_timings - 1] - timings[num_timings - 2]);
        } else {
            // baseline: group simplicial cells into material cells
            timings.push_back(
                    timings[num_timings - 1] - timings[num_timings - 2] - timings[num_timings - 3]);
        }
    }

    // test: export MI_mesh, patches, chains
    if (!args.timing_only && material_cells.size() > 0) {
        save_result(output_dir + "/MI_mesh.json",
                                         MI_pts,
                                         MI_faces,
                                         patches,
                                         MI_edges,
                                         chains,
                                         non_manifold_edges_of_vert,
                                         half_patch_pair_list,
                                         shells,
                                         components,
                                         material_cells);
        //
        save_result_msh(output_dir + "/MI_mesh",
                                             MI_pts,
                                             MI_faces,
                                             patches,
                                             MI_edges,
                                             chains,
                                             non_manifold_edges_of_vert,
                                             shells,
                                             components,
                                             material_cells);
    }
    // save timing records
    save_timings(output_dir + "/timings.json", timing_labels, timings);
    // save statistics
    save_statistics(output_dir + "/stats.json", stats_labels, stats);

    return 0;
}
