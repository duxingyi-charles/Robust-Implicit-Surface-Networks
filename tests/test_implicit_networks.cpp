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
    std::vector<Edge> iso_edges;
    std::vector<std::vector<size_t>> chains;
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
    std::vector<std::vector<size_t>> shells;
    std::vector<std::vector<size_t>> arrangement_cells;
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
        // plane
        size_t func_id = 0;
        std::array<double,3> point {0,0,0};
        std::array<double,3> normal {1,0,0};
        for (int i = 0; i < n_pts; i++) {
            funcVals(i, func_id) = compute_plane_distance(point, normal, pts[i]);
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
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 1);
        REQUIRE(chains.size() == 0);
        REQUIRE(arrangement_cells.size() == 2);
    }

    SECTION("one sphere") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 1;
        funcVals.resize(n_pts, n_func);
        // sphere
        size_t func_id = 0;
        std::array<double,3> center {0,0,0};
        double radius = 0.5;
        for (int i = 0; i < n_pts; i++) {
            funcVals(i, func_id) = compute_sphere_distance(center, radius, pts[i]);
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
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 1);
        REQUIRE(chains.size() == 0);
        REQUIRE(arrangement_cells.size() == 2);
    }

    SECTION("plane and sphere") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 2;
        funcVals.resize(n_pts, n_func);
        size_t func_id;
        // plane
        func_id = 0;
        std::array<double,3> point {0,0,0};
        std::array<double,3> normal {1,0,0};
        for (int i = 0; i < n_pts; i++) {
            funcVals(i, func_id) = compute_plane_distance(point, normal, pts[i]);
        }
        // sphere
        func_id = 1;
        std::array<double,3> center {0,0,0};
        double radius = 0.5;
        for (int i = 0; i < n_pts; i++) {
            funcVals(i, func_id) = compute_sphere_distance(center, radius, pts[i]);
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
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 4);
        REQUIRE(chains.size() == 1);
        REQUIRE(arrangement_cells.size() == 4);
    }

    SECTION("two spheres") {
        // compute function values on tet grid vertices
        size_t n_pts = pts.size();
        size_t n_func = 2;
        funcVals.resize(n_pts, n_func);
        size_t func_id;
        // sphere
        {
            func_id = 0;
            std::array<double, 3> center{-0.3, 0, 0};
            double radius = 0.5;
            for (int i = 0; i < n_pts; i++) {
                funcVals(i, func_id) = compute_sphere_distance(center, radius, pts[i]);
            }
        }
        // sphere
        {
            func_id = 1;
            std::array<double, 3> center{0.3, 0, 0};
            double radius = 0.5;
            for (int i = 0; i < n_pts; i++) {
                funcVals(i, func_id) = compute_sphere_distance(center, radius, pts[i]);
            }
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
                iso_edges,chains,
                non_manifold_edges_of_vert,
                shells,arrangement_cells,
                timing_labels,timings,
                stats_labels,stats);
        REQUIRE(success);

        // check
        REQUIRE(patches.size() == 4);
        REQUIRE(chains.size() == 1);
        REQUIRE(arrangement_cells.size() == 4);
    }


}
