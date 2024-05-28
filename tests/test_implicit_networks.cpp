//
// Created by Charles Du on 8/19/22.
//
#include <simplicial_arrangement/lookup_table.h>
#include "implicit_arrangement.h"
#include "material_interface.h"
#include "implicit_functions.h"

#include <catch2/catch.hpp>

// tests:
// examples with known structure
// consistency:
// - look-up vs non-look-up
// - topo-ray-shooting vs cell-grouping

TEST_CASE("implicit arrangement on known examples", "[IA][examples]") {
    bool robust_test = false;
    bool use_lookup = true;
    bool loaded = simplicial_arrangement::load_lookup_table();
    REQUIRE(loaded);
    bool use_secondary_lookup = true;
    bool use_topo_ray_shooting = true;

    // generate tet grid
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    int tet_mesh_resolution = 101;
    std::array<double,3> tet_mesh_bbox_min {-1,-1,-1};
    std::array<double,3> tet_mesh_bbox_max {1,1,1};
    generate_tet_mesh(tet_mesh_resolution, tet_mesh_bbox_min,
                      tet_mesh_bbox_max, pts, tets);

    // function values
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
    // implicit arrangement result
    std::vector<std::array<double, 3>> iso_pts;
    std::vector<PolygonFace> iso_faces;
    std::vector<std::vector<size_t>> patches;
    std::vector<size_t> patch_function_label;
    std::vector<Edge> iso_edges;
    std::vector<std::vector<size_t>> chains;
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
    std::vector<std::vector<size_t>> shells;
    std::vector<std::vector<size_t>> arrangement_cells;
    std::vector<std::vector<bool>> cell_function_label;
    // record timings
    std::vector<std::string> timing_labels;
    std::vector<double> timings;
    // record stats
    std::vector<std::string> stats_labels;
    std::vector<size_t> stats;

    SECTION("one plane") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 1;
        funcVals.resize(n_pts, n_func);
        
        if (!load_functions(std::string(TEST_FILE) + "/1-plane.json", pts, funcVals)) {
            throw std::runtime_error("ERROR: Failed to load functions.");
        }


        // compute implicit arrangement
        bool success = implicit_arrangement(
                robust_test,
                use_lookup,
                use_secondary_lookup,
                use_topo_ray_shooting,
                //
                pts, tets, funcVals,
                //
                iso_pts,iso_faces,patches,
                patch_function_label,
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,cell_function_label,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 1);
        REQUIRE(patch_function_label.size() == 1);
        REQUIRE(chains.size() == 0);
        REQUIRE(arrangement_cells.size() == 2);
        
        // function label check
        std::vector<size_t> patch_gt = {0};
        REQUIRE(patch_function_label == patch_gt);
        std::vector<std::vector<bool>> cell_gt = {{1}, {0}};
        REQUIRE(cell_function_label == cell_gt);
    }

    SECTION("one sphere") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 1;
        funcVals.resize(n_pts, n_func);
        
        if (!load_functions(std::string(TEST_FILE) + "/1-sphere.json", pts, funcVals)) {
            throw std::runtime_error("ERROR: Failed to load functions.");
        }

        // compute implicit arrangement
        bool success = implicit_arrangement(
                robust_test,
                use_lookup,
                use_secondary_lookup,
                use_topo_ray_shooting,
                //
                pts, tets, funcVals,
                //
                iso_pts,iso_faces,patches,
                patch_function_label,
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,cell_function_label,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 1);
        REQUIRE(patch_function_label.size() == 1);
        REQUIRE(chains.size() == 0);
        REQUIRE(arrangement_cells.size() == 2);
        
        // function label check
        std::vector<size_t> patch_gt = {0};
        REQUIRE(patch_function_label == patch_gt);
        std::vector<std::vector<bool>> cell_gt = {{1}, {0}};
        REQUIRE(cell_function_label == cell_gt);
    }

    SECTION("plane and sphere") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 2;
        funcVals.resize(n_pts, n_func);
        size_t func_id;
        
        if (!load_functions(std::string(TEST_FILE) + "/2-planesphere.json", pts, funcVals)) {
            throw std::runtime_error("ERROR: Failed to load functions.");
        }

        // compute implicit arrangement
        bool success = implicit_arrangement(
                robust_test,
                use_lookup,
                use_secondary_lookup,
                use_topo_ray_shooting,
                //
                pts, tets, funcVals,
                //
                iso_pts,iso_faces,patches,
                patch_function_label,
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,cell_function_label,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 4);
        REQUIRE(patch_function_label.size() == 4);
        REQUIRE(chains.size() == 1);
        REQUIRE(arrangement_cells.size() == 4);
        
        // function label check
        std::vector<size_t> patch_gt = {1,0,0,1};
        REQUIRE(patch_function_label == patch_gt);
        std::vector<std::vector<bool>> cell_gt = {{0,1}, {0,0},{1,0},{1,1}};
        REQUIRE(cell_function_label == cell_gt);
    }

    SECTION("two spheres") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 2;
        funcVals.resize(n_pts, n_func);
        size_t func_id;
        
        if (!load_functions(std::string(TEST_FILE) + "/2-sphere.json", pts, funcVals)) {
            throw std::runtime_error("ERROR: Failed to load functions.");
        }

        // compute implicit arrangement
        bool success = implicit_arrangement(
                robust_test,
                use_lookup,
                use_secondary_lookup,
                use_topo_ray_shooting,
                //
                pts, tets, funcVals,
                //
                iso_pts,iso_faces,patches,
                patch_function_label,
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,cell_function_label,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 4);
        REQUIRE(patch_function_label.size() == 4);
        REQUIRE(chains.size() == 1);
        REQUIRE(arrangement_cells.size() == 4);
        
        // function label check
        std::vector<size_t> patch_gt = {0,1,1,0};
        REQUIRE(patch_function_label == patch_gt);
        std::vector<std::vector<bool>> cell_gt = {{1,0}, {0,0},{1,1},{0,1}};
        REQUIRE(cell_function_label == cell_gt);
    }
    
    SECTION("three spheres config 1") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 3;
        funcVals.resize(n_pts, n_func);
        size_t func_id;
        
        if (!load_functions(std::string(TEST_FILE) + "/3-sphere-1.json", pts, funcVals)) {
            throw std::runtime_error("ERROR: Failed to load functions.");
        }

        // compute implicit arrangement
        bool success = implicit_arrangement(
                robust_test,
                use_lookup,
                use_secondary_lookup,
                use_topo_ray_shooting,
                //
                pts, tets, funcVals,
                //
                iso_pts,iso_faces,patches,
                patch_function_label,
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,cell_function_label,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 3);
        REQUIRE(patch_function_label.size() == 3);
        REQUIRE(chains.size() == 0);
        REQUIRE(arrangement_cells.size() == 4);
        
        // function label check
        std::vector<size_t> patch_gt = {0,2,1};
        REQUIRE(patch_function_label == patch_gt);
        std::vector<std::vector<bool>> cell_gt = {{1,0,0}, {0,0,0},{0,0,1},{0,1,0}};
        REQUIRE(cell_function_label == cell_gt);
    }
    
    SECTION("three spheres config 2") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 3;
        funcVals.resize(n_pts, n_func);
        size_t func_id;
        
        if (!load_functions(std::string(TEST_FILE) + "/3-sphere-2.json", pts, funcVals)) {
            throw std::runtime_error("ERROR: Failed to load functions.");
        }

        // compute implicit arrangement
        bool success = implicit_arrangement(
                robust_test,
                use_lookup,
                use_secondary_lookup,
                use_topo_ray_shooting,
                //
                pts, tets, funcVals,
                //
                iso_pts,iso_faces,patches,
                patch_function_label,
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,cell_function_label,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 7);
        REQUIRE(patch_function_label.size() == 7);
        REQUIRE(chains.size() == 2);
        REQUIRE(arrangement_cells.size() == 6);
        
        // function label check
        std::vector<size_t> patch_gt = {0, 2, 0, 2, 1, 1, 2};
        REQUIRE(patch_function_label == patch_gt);
        std::vector<std::vector<bool>> cell_gt = {{1, 0, 0}, {0, 0, 0}, {1, 0, 1}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}};
        REQUIRE(cell_function_label == cell_gt);
    }
    
    SECTION("three spheres config 3") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 3;
        funcVals.resize(n_pts, n_func);
        size_t func_id;
        
        if (!load_functions(std::string(TEST_FILE) + "/3-sphere-3.json", pts, funcVals)) {
            throw std::runtime_error("ERROR: Failed to load functions.");
        }

        // compute implicit arrangement
        bool success = implicit_arrangement(
                robust_test,
                use_lookup,
                use_secondary_lookup,
                use_topo_ray_shooting,
                //
                pts, tets, funcVals,
                //
                iso_pts,iso_faces,patches,
                patch_function_label,
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,cell_function_label,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 12);
        REQUIRE(patch_function_label.size() == 12);
        REQUIRE(chains.size() == 6);
        REQUIRE(arrangement_cells.size() == 8);
        
        // function label check
        std::vector<size_t> patch_gt = {0, 2, 0, 2, 1, 1, 2, 1, 0, 1, 2, 0};
        REQUIRE(patch_function_label == patch_gt);
        std::vector<std::vector<bool>> cell_gt = {{1, 0, 0},{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}, {0, 1, 1}, {0, 1, 0}};
        REQUIRE(cell_function_label == cell_gt);
    }
    
    SECTION("three spheres config 4") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 3;
        funcVals.resize(n_pts, n_func);
        size_t func_id;
        
        if (!load_functions(std::string(TEST_FILE) + "/3-sphere-4.json", pts, funcVals)) {
            throw std::runtime_error("ERROR: Failed to load functions.");
        }

        // compute implicit arrangement
        bool success = implicit_arrangement(
                robust_test,
                use_lookup,
                use_secondary_lookup,
                use_topo_ray_shooting,
                //
                pts, tets, funcVals,
                //
                iso_pts,iso_faces,patches,
                patch_function_label,
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,cell_function_label,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 3);
        REQUIRE(patch_function_label.size() == 3);
        REQUIRE(chains.size() == 0);
        REQUIRE(arrangement_cells.size() == 4);
        
        // function label check
        std::vector<size_t> patch_gt = {2, 1, 0};
        REQUIRE(patch_function_label == patch_gt);
        std::vector<std::vector<bool>> cell_gt = {{0, 0, 1}, {0, 0, 0}, {0, 1, 1}, {1, 1, 1}};
        REQUIRE(cell_function_label == cell_gt);
    }


}
