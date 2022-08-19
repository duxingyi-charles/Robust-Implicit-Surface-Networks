//
// Created by Charles Du on 8/18/22.
//
#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include "implicit_arrangement.h"

#include "ScopedTimer.h"
typedef std::chrono::duration<double> Time_duration;

bool implicit_arrangement(
        bool robust_test,
        bool use_lookup,
        bool use_secondary_lookup,
        bool use_topo_ray_shooting,
        //
        const std::vector<std::array<double, 3>>& pts,
        const std::vector<std::array<size_t, 4>>& tets,
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& funcVals,
        //
        std::vector<std::array<double, 3>>& iso_pts,
        std::vector<PolygonFace>& iso_faces,
        std::vector<std::vector<size_t>>& patches,
        std::vector<Edge>& iso_edges,
        std::vector<std::vector<size_t>>& chains,
        std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
        std::vector<std::vector<size_t>>& shells,
        std::vector<std::vector<size_t>>& arrangement_cells,
        std::vector<std::string>& timing_labels,
        std::vector<double>& timings,
        std::vector<std::string>& stats_labels,
        std::vector<size_t>& stats)
{
    using namespace simplicial_arrangement;

    size_t n_tets = tets.size();
    size_t n_pts = pts.size();
    std::cout << "tet mesh: " << n_pts << " verts, " << n_tets << " tets." << std::endl;
    stats_labels.emplace_back("num_pts");
    stats.push_back(n_pts);
    stats_labels.emplace_back("num_tets");
    stats.push_back(n_tets);

    size_t n_func = funcVals.cols();

    // compute function signs at vertices
    auto sign = [](double x) -> int {
        return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
    };
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcSigns;
    std::vector<bool> is_degenerate_vertex;
    bool found_degenerate_vertex = false;
    size_t num_degenerate_vertex = 0;
    {
        timing_labels.emplace_back("func signs");
        ScopedTimer<> timer("func signs");
        is_degenerate_vertex.resize(n_pts, false);
        funcSigns.resize(n_pts, n_func);
        for (Eigen::Index i = 0; i < n_pts; i++) {
            for (Eigen::Index j = 0; j < n_func; j++) {
                funcSigns(i, j) = sign(funcVals(i, j));
                if (funcSigns(i, j) == 0) {
                    is_degenerate_vertex[i] = true;
                    found_degenerate_vertex = true;
                    num_degenerate_vertex++;
                }
            }
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num_degenerate_vertex = " << num_degenerate_vertex << std::endl;
    stats_labels.emplace_back("num_degenerate_vertex");
    stats.push_back(num_degenerate_vertex);

    // filter active functions in each tetrahedron
    size_t num_intersecting_tet = 0;
    std::vector<size_t> func_in_tet;  // active function indices in CRS vector format
    std::vector<size_t> start_index_of_tet;
    {
        timing_labels.emplace_back("filter");
        ScopedTimer<> timer("filter");
        func_in_tet.reserve(n_tets);
        start_index_of_tet.reserve(n_tets + 1);
        start_index_of_tet.push_back(0);
        int pos_count;
        int neg_count;
        for (Eigen::Index i = 0; i < n_tets; i++) {
            for (Eigen::Index j = 0; j < n_func; j++) {
                pos_count = 0;
                neg_count = 0;
                for (size_t vId : tets[i]) {
                    if (funcSigns( vId,j) == 1) {
                        pos_count += 1;
                    } else if (funcSigns( vId,j) == -1) {
                        neg_count += 1;
                    }
                }
                // tets[i].size() == 4
                if (pos_count < 4 && neg_count < 4) {
                    func_in_tet.push_back(j);
                }
            }
            if (func_in_tet.size() > start_index_of_tet.back()) {
                ++num_intersecting_tet;
            }
            start_index_of_tet.push_back(func_in_tet.size());
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num_intersecting_tet = " << num_intersecting_tet << std::endl;
    stats_labels.emplace_back("num_intersecting_tet");
    stats.push_back(num_intersecting_tet);

    // compute arrangement in each tet
    std::vector<Arrangement<3>> cut_results;
    std::vector<size_t> cut_result_index;
    //
    Time_duration time_1_func = Time_duration::zero();
    Time_duration time_2_func = Time_duration::zero();
    Time_duration time_more_func = Time_duration::zero();
    size_t num_1_func = 0;
    size_t num_2_func = 0;
    size_t num_more_func = 0;
    // failure types for robustness test
    bool is_type2 = false; // crash in the normal order
    bool is_type3 = false; // succeed in the normal order, crash in the reverse order
    bool is_type1 = false; // no crash for either order, but results are inconsistent
    {
        ScopedTimer<> timer("simp_arr");
        if (robust_test) {
            Arrangement<3> cut_result, cut_result2;
            size_t start_index;
            size_t num_func;
            std::vector<Plane<double, 3>> planes;
            std::vector<Plane<double, 3>> planes_reverse; // planes in reverse order
            planes.reserve(3);
            planes_reverse.reserve(3);
            for (size_t i = 0; i < tets.size(); i++) {
                start_index = start_index_of_tet[i];
                num_func = start_index_of_tet[i + 1] - start_index;
                if (num_func == 0) {
                    continue;
                }
                std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                //
                size_t v1 = tets[i][0];
                size_t v2 = tets[i][1];
                size_t v3 = tets[i][2];
                size_t v4 = tets[i][3];
                planes.clear();
                for (size_t j = 0; j < num_func; j++) {
                    size_t f_id = func_in_tet[start_index + j];
                    planes.emplace_back();
                    auto& plane = planes.back();
                    plane[0] = funcVals(v1, f_id);
                    plane[1] = funcVals(v2, f_id);
                    plane[2] = funcVals(v3, f_id);
                    plane[3] = funcVals(v4, f_id);
                }
                // reverse plane order
                planes_reverse.clear();
                for (size_t j = 0; j < num_func; ++j) {
                    planes_reverse.emplace_back(planes[num_func -1 -j]);
                }
                //
                bool crashed = false;
                if (use_lookup && !use_secondary_lookup && num_func == 2) {
                    cut_result_index.push_back(cut_results.size());
                    disable_lookup_table();
                    try {
                        cut_result = compute_arrangement(planes);
                    } catch (std::runtime_error& e) {
                        crashed = true;
                        is_type2 = true;
                        break;
                    }
                    if (!crashed) {
                        try {
                            cut_result2 = compute_arrangement(planes_reverse);
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            is_type3 = true;
                        }
                        if (!crashed) {
                            if (cut_result.vertices.size() != cut_result2.vertices.size() ||
                                cut_result.faces.size() != cut_result2.faces.size() ||
                                cut_result.cells.size() != cut_result2.cells.size()) {
                                // inconsistent results
                                is_type1 = true;
                            }
                        }
                    }
                    enable_lookup_table();
                } else {
                    try {
                        cut_result = compute_arrangement(planes);
                    } catch (std::runtime_error& e) {
                        crashed = true;
                        is_type2 = true;
                        break;
                    }
                    if (!crashed) {
                        try {
                            cut_result2 = compute_arrangement(planes_reverse);
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            is_type3 = true;
                        }
                        if (!crashed) {
                            if (cut_result.vertices.size() != cut_result2.vertices.size() ||
                                cut_result.faces.size() != cut_result2.faces.size() ||
                                cut_result.cells.size() != cut_result2.cells.size()) {
                                // inconsistent results
                                is_type1 = true;
                            }
                        }
                    }
                }
                // update timing
                std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                switch (num_func) {
                    case 1:
                        time_1_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_1_func;
                        break;
                    case 2:
                        time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_2_func;
                        break;
                    default:
                        time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_more_func;
                        break;
                }
            }
        }
        else { // not performing robustness test
            cut_results.reserve(num_intersecting_tet);
            cut_result_index.reserve(n_tets);
            size_t start_index;
            size_t num_func;
            std::vector<Plane<double, 3>> planes;
            planes.reserve(3);
            try {
                for (size_t i = 0; i < tets.size(); i++) {
                    start_index = start_index_of_tet[i];
                    num_func = start_index_of_tet[i + 1] - start_index;
                    if (num_func == 0) {
                        cut_result_index.push_back(Arrangement<3>::None);
                        continue;
                    }
                    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                    //
                    size_t v1 = tets[i][0];
                    size_t v2 = tets[i][1];
                    size_t v3 = tets[i][2];
                    size_t v4 = tets[i][3];
                    planes.clear();
                    for (size_t j = 0; j < num_func; j++) {
                        size_t f_id = func_in_tet[start_index + j];
                        planes.emplace_back();
                        auto& plane = planes.back();
                        plane[0] = funcVals(v1, f_id);
                        plane[1] = funcVals(v2, f_id);
                        plane[2] = funcVals(v3, f_id);
                        plane[3] = funcVals(v4, f_id);
                    }
                    //
                    if (use_lookup && !use_secondary_lookup && num_func == 2) {
                        cut_result_index.push_back(cut_results.size());
                        disable_lookup_table();
                        cut_results.emplace_back(std::move(compute_arrangement(planes)));
                        enable_lookup_table();
                    } else {
                        cut_result_index.push_back(cut_results.size());
                        cut_results.emplace_back(std::move(compute_arrangement(planes)));
                    }
                    // update timing
                    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                    switch (num_func) {
                        case 1:
                            time_1_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_1_func;
                            break;
                        case 2:
                            time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_2_func;
                            break;
                        default:
                            time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_more_func;
                            break;
                    }
                }
            } catch (std::runtime_error& e) {
                std::cout << e.what() << std::endl;
                return false;
            }
        }
        //
        timing_labels.emplace_back("simp_arr(other)");
        timings.push_back(
                timer.toc() - time_1_func.count() - time_2_func.count() - time_more_func.count());
    }
    timing_labels.emplace_back("simp_arr(1 func)");
    timings.push_back(time_1_func.count());
    timing_labels.emplace_back("simp_arr(2 func)");
    timings.push_back(time_2_func.count());
    timing_labels.emplace_back("simp_arr(>=3 func)");
    timings.push_back(time_more_func.count());
    std::cout << " -- [simp_arr(1 func)]: " << time_1_func.count() << " s" << std::endl;
    std::cout << " -- [simp_arr(2 func)]: " << time_2_func.count() << " s" << std::endl;
    std::cout << " -- [simp_arr(>=3 func)]: " << time_more_func.count() << " s" << std::endl;
    //
    stats_labels.emplace_back("num_1_func");
    stats.push_back(num_1_func);
    stats_labels.emplace_back("num_2_func");
    stats.push_back(num_2_func);
    stats_labels.emplace_back("num_more_func");
    stats.push_back(num_more_func);

    if (robust_test) {
        if (is_type2) {
            std::cout << "type 2 failure (crash in the normal order)." << std::endl;
            return false;
        } else if (is_type3) {
            std::cout << "type 3 failure (crash in the reverse order)." << std::endl;
            return false;
        } else if (is_type1) {
            std::cout << "type 1 failure (inconsistency)." << std::endl;
            return false;
        } else {
            std::cout << "success." << std::endl;
            return true;
        }
    }

    // extract arrangement mesh: combining results from all tets to produce a mesh
    std::vector<IsoVert> iso_verts;
//    std::vector<PolygonFace> iso_faces;
    // the following data are only needed when we use
    // the baseline nesting algorithm
    // (group simplicial cells into arrangement cells)
    std::vector<long long> global_vId_of_tet_vert;
    std::vector<size_t> global_vId_start_index_of_tet;
    std::vector<size_t> iso_fId_of_tet_face;
    std::vector<size_t> iso_fId_start_index_of_tet;
    //
    {
        timing_labels.emplace_back("extract mesh");
        ScopedTimer<> timer("extract mesh");
        if (use_topo_ray_shooting) {
            extract_iso_mesh(num_1_func,
                             num_2_func,
                             num_more_func,
                             cut_results,
                             cut_result_index,
                             func_in_tet,
                             start_index_of_tet,
                             tets,
                             iso_verts,
                             iso_faces);
        } else { // nesting algorithm: group simplicial cells into arrangement cells
            extract_iso_mesh(num_1_func,
                             num_2_func,
                             num_more_func,
                             cut_results,
                             cut_result_index,
                             func_in_tet,
                             start_index_of_tet,
                             tets,
                             iso_verts,
                             iso_faces,
                             global_vId_of_tet_vert,
                             global_vId_start_index_of_tet,
                             iso_fId_of_tet_face,
                             iso_fId_start_index_of_tet);
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num iso-vertices = " << iso_verts.size() << std::endl;
    std::cout << "num iso-faces = " << iso_faces.size() << std::endl;
    stats_labels.emplace_back("num_iso_verts");
    stats.push_back(iso_verts.size());
    stats_labels.emplace_back("num_iso_faces");
    stats.push_back(iso_faces.size());

    // compute xyz coordinates of iso-vertices
//    std::vector<std::array<double, 3>> iso_pts;
    {
        timing_labels.emplace_back("compute xyz");
        ScopedTimer<> timer("compute xyz");
        compute_iso_vert_xyz(iso_verts, funcVals, pts, iso_pts);
        timings.push_back(timer.toc());
    }

    //  compute iso-edges and edge-face connectivity
//    std::vector<Edge> iso_edges;
    std::vector<std::vector<size_t>> edges_of_iso_face;
    {
        timing_labels.emplace_back("isoEdge-face connectivity");
        ScopedTimer<> timer("isoEdge-face connectivity");
        compute_mesh_edges(iso_faces, edges_of_iso_face, iso_edges);
        timings.push_back(timer.toc());
    }
    std::cout << "num iso-edges = " << iso_edges.size() << std::endl;
    stats_labels.emplace_back("num_iso_edges");
    stats.push_back(iso_edges.size());

    // group iso-faces into patches
//    std::vector<std::vector<size_t>> patches;
    {
        timing_labels.emplace_back("patches");
        ScopedTimer<> timer("patches");
        compute_patches(edges_of_iso_face, iso_edges, patches);
        timings.push_back(timer.toc());
    }
    std::cout << "num patches = " << patches.size() << std::endl;
    stats_labels.emplace_back("num_patches");
    stats.push_back(patches.size());

    // compute map: iso-face Id --> patch Id
    std::vector<size_t> patch_of_face;
    {
        timing_labels.emplace_back("face-patch map");
        ScopedTimer<> timer("face-patch map");
        patch_of_face.resize(iso_faces.size());
        for (size_t i = 0; i < patches.size(); i++) {
            for (const auto& fId : patches[i]) {
                patch_of_face[fId] = i;
            }
        }
        timings.push_back(timer.toc());
    }

    // group non-manifold iso-edges into chains
//    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
//    std::vector<std::vector<size_t>> chains;
    {
        timing_labels.emplace_back("chains");
        ScopedTimer<> timer("chains");
        non_manifold_edges_of_vert.resize(iso_pts.size());
        // get incident non-manifold edges for iso-vertices
        for (size_t i = 0; i < iso_edges.size(); i++) {
            if (iso_edges[i].face_edge_indices.size() >
                2) { // non-manifold edge (not a boundary edge)
                // there is only one patch incident to a boundary edge,
                // so there is no need to figure out the "order" of patches around a boundary
                // edge
                non_manifold_edges_of_vert[iso_edges[i].v1].push_back(i);
                non_manifold_edges_of_vert[iso_edges[i].v2].push_back(i);
            }
        }
        // group non-manifold iso-edges into chains
        compute_chains(iso_edges, non_manifold_edges_of_vert, chains);
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
            const auto& iso_edge = iso_edges[chain_representatives[i]];
            // with degeneracy handling
            try {
                compute_face_order(iso_edge,
                                   iso_faces,
                                   cut_results,
                                   cut_result_index,
                                   incident_tets,
                                   half_face_pair_list[i]);
            } catch (std::exception& e) {
                std::cout << "order patches failed: " << e.what() << std::endl;
                return false;
            }
        }
        // replace iso-face indices by patch indices
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
//    std::vector<std::vector<size_t>> shells;
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
    // an arrangement cell is represented by a list of bounding shells
//    std::vector<std::vector<size_t>> arrangement_cells;
    {
        ScopedTimer<> timer("arrangement cells");
        if (components.size() < 2) { // no nesting problem, each shell is an arrangement cell
            arrangement_cells.reserve(shells.size());
            for (size_t i = 0; i < shells.size(); ++i) {
                arrangement_cells.emplace_back(1);
                arrangement_cells.back()[0] = i;
            }
        } else { // resolve nesting order
            if (use_topo_ray_shooting) {
                timing_labels.emplace_back("arrCells(ray shooting)");
                ScopedTimer<> timer("arrangement cells: topo ray shooting");
                topo_ray_shooting(pts, tets,
                                  cut_results, cut_result_index,
                                  iso_verts, iso_faces,
                                  patches, patch_of_face,
                                  shells, shell_of_half_patch,
                                  components, component_of_patch,
                                  arrangement_cells);
                timings.push_back(timer.toc());
            }
            else { // group simplicial cells into arrangement cells
                // ------------------- build adjacency graph of simplicial cells ---------------
                // pair (tetId, tet_cell_Id)
                std::vector<std::pair<size_t, size_t>> tet_cell_of_simp_cell;
                // simplicial half-face info
                // if the face is on isosurface, let s be its shell index, record (-s-1)
                // if the face is not on isosurface, record the simp_cell index on the opposite side
                std::vector<long long> simp_half_face_info;
                std::vector<size_t> simp_hFace_start_index;
                {
                    timing_labels.emplace_back("arrCells(build simpCell graph)");
                    ScopedTimer<> timer("arrangement cells: build simplicial cells graph");
                    build_simplicial_cell_adjacency(tets,
                                                    cut_results, cut_result_index,
                                                    global_vId_of_tet_vert, global_vId_start_index_of_tet,
                                                    iso_fId_of_tet_face, iso_fId_start_index_of_tet,
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
                    timing_labels.emplace_back("arrCells(group simpCells into arrCells)");
                    ScopedTimer<> timer("arrangement cells: group simplicial cells");
                    compute_simplicial_cell_connected_components(
                            tet_cell_of_simp_cell, simp_half_face_info, simp_hFace_start_index,
                            arrangement_cells);
                    timings.push_back(timer.toc());
                }
            }
        }
        timings.push_back(timer.toc());
        timing_labels.emplace_back("arrangement cells");
    }
    std::cout << "num_cells = " << arrangement_cells.size() << std::endl;
    stats_labels.emplace_back("num_cells");
    stats.push_back(arrangement_cells.size());

    if (components.size() > 1) {
        timing_labels.back() = "arrCells(other)";
        size_t num_timings = timings.size();
        if (use_topo_ray_shooting) {
            timings.back() = timings[num_timings - 1] - timings[num_timings - 2];
        } else {
            // baseline: group simplicial cells into arrangement cells
            timings.back() = timings[num_timings - 1] - timings[num_timings - 2] - timings[num_timings - 3];
        }
    }

    return true;
}
