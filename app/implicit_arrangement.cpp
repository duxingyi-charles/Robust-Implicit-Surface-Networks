#include <simplicial_arrangement/lookup_table.h>

#include <iostream>
#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "implicit_arrangement.h"
#include "implicit_functions.h"

using namespace simplicial_arrangement;

int main(int argc, const char* argv[])
{
    struct
    {
        std::string config_file;
        bool timing_only = false;
        bool robust_test = false;
    } args;
    CLI::App app{"Implicit Arrangement Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_flag("-T,--timing-only", args.timing_only, "Record timing without saving results");
    app.add_flag("-R,--robust-test",args.robust_test, "Perform robustness test");
    CLI11_PARSE(app, argc, argv);

    // parse configure file
    Config config = parse_config_file(args.config_file);

    if (config.use_lookup) {
        // load lookup table
        std::cout << "load table ..." << std::endl;
        bool loaded = load_lookup_table();
        if (loaded) {
            std::cout << "loading finished." << std::endl;
        } else {
            std::cout << "loading failed." << std::endl;
            return -1;
        }
    } else {
        disable_lookup_table();
        config.use_secondary_lookup = false;
    }

    // load tet mesh
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    if (config.tet_mesh_file != "") {
        std::cout << "load mesh file " << config.tet_mesh_file << std::endl;
        load_tet_mesh(config.tet_mesh_file, pts, tets);
    } else {
        std::cout << "generating mesh with resolution "
            << config.tet_mesh_resolution << std::endl;
        generate_tet_mesh(config.tet_mesh_resolution, config.tet_mesh_bbox_min,
                config.tet_mesh_bbox_max, pts, tets);
    }

    // load implicit functions and compute function values at vertices
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
    if (load_functions(config.func_file, pts, funcVals)) {
        std::cout << "function loading finished." << std::endl;
    } else {
        std::cout << "function loading failed." << std::endl;
        return -2;
    }

    // compute implicit arrangement
    std::vector<std::array<double, 3>> iso_pts;
    std::vector<PolygonFace> iso_faces;
    std::vector<std::vector<size_t>> patches;
    std::vector<size_t> patch_function_label;
    std::vector<Edge> iso_edges;
    std::vector<std::vector<size_t>> chains;
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
    std::vector<std::vector<size_t>> shells;
    std::vector<std::vector<size_t>> arrangement_cells;
    std::vector<std::vector<size_t>> cell_function_label;
    // record timings
    std::vector<std::string> timing_labels;
    std::vector<double> timings;
    // record stats
    std::vector<std::string> stats_labels;
    std::vector<size_t> stats;

    if (!implicit_arrangement(
            args.robust_test,
            config.use_lookup,
            config.use_secondary_lookup,
            config.use_topo_ray_shooting,
    //
    pts, tets, funcVals,
    //
    iso_pts,iso_faces,patches, patch_function_label,
    iso_edges,chains,
    non_manifold_edges_of_vert,
    shells,arrangement_cells,cell_function_label,
    timing_labels,timings,
    stats_labels,stats)) {
        return -1;
    }
    if (args.robust_test) return 0;

    // save result
    if (!args.timing_only) {
        save_result(config.output_dir + "/mesh.json",
                    iso_pts,
                    iso_faces,
                    patches,
                    patch_function_label,
                    iso_edges,
                    chains,
                    non_manifold_edges_of_vert,
                    shells,
                    arrangement_cells,
                    cell_function_label);
        //
        save_result_msh(config.output_dir + "/mesh",
                        iso_pts,
                        iso_faces,
                        patches,
                        iso_edges,
                        chains,
                        non_manifold_edges_of_vert,
                        shells,
                        arrangement_cells);
    }
    // save timing records
    save_timings(config.output_dir + "/timings.json", timing_labels, timings);
    // save statistics
    save_statistics(config.output_dir + "/stats.json", stats_labels, stats);

    return 0;
}
