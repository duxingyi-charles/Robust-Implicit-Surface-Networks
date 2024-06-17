#include <simplicial_arrangement/lookup_table.h>

#include <iostream>
#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "implicit_arrangement.h"
#include "implicit_functions.h"
#include "csg.h"

using namespace simplicial_arrangement;

///
///The structure of a csg unit
///A nested data structure where `operation` takes in a boolean operation, and the elements contain either a function index or later csg unit, where netaives represent the indices, and positives reprsent the index in the CSG tree.
///e.g. {Intersection, -1, 2} represents an intersection operation between the first function and the second csg unit.
///
struct csg_unit{
    int operation;
    std::array<int, 2> elements;
};

///
///Defines the enumeration of CSG operations. In the data structure, Intersection is 0, Union is 1, and Negation is 2.
enum csg_operations{
    Intersection,
    Union,
    Negation
};

///load the csg file
///@param[in] filename          The name of the input CSG tree
///@param[out] tree         The loaded tree structure
///
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

/// an iterative algortihm that traverses through the csg tree
///
///@param[in] csgTree           The CSG structure: a list of csg units
///@param[in] curNode           The current index in the csg structure
///@param[in] funcInt           The intervals of all the functions
///
///@param[out] std::pair
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
        case Union:
            interval = {std::max(childInt1[0], childInt2[0]), std::max(childInt1[1], childInt2[1])};
            if(interval[0]*interval[1] <= 0){
                for (int i = 0; i < funcInt.size(); i++){
                    af[i] = childAF1[i] * childAF2[i];
                }
            }
            break;
        case Intersection:
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

int main(int argc, const char* argv[])
{
    struct
    {
        std::string config_file;
        std::string csg_file;
        bool timing_only = false;
        bool robust_test = false;
        bool positive_inside = true;
    } args;
    CLI::App app{"Implicit Arrangement Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_flag("-T,--timing-only", args.timing_only, "Record timing without saving results");
    app.add_flag("-R,--robust-test",args.robust_test, "Perform robustness test");
    app.add_flag("-P,--positive-inside",args.positive_inside, "Using positives as the inside");
    app.add_option("--tree", args.csg_file, "CSG Tree file");
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
    std::vector<bool> patch_sign_label;
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
    
    auto lambda = [&](std::vector<bool> cells_label){
        bool cell_label;
        if (args.csg_file == ""){
            for (auto sign : cells_label){
                if (sign == 1){
                    cell_label = true;
                    break;
                }
            }
            //throw std::runtime_error("ERROR: no csg file provided");
            return cell_label;
        }else{
            std::vector<csg_unit> csgTree;
            bool loaded = load_csgTree(args.csg_file, csgTree);
            if (!loaded){
                throw std::runtime_error("ERROR: reading csg file failed");
            }
                std::vector<std::array<double, 2>> funcInt;
                //funcInt.reserve(label.size());
                for (auto sign : cells_label){
                    funcInt.emplace_back((sign > 0) ? std::array<double, 2>({ 1, 2 }) : std::array<double, 2>({ -1, -2 }));
                }
                //Currently, csg tree traversal is using intervals; will change to signs later.
                std::pair<std::array<double, 2>, std::vector<int>> csgResult = iterTree(csgTree, 1, funcInt);
                cell_label = (csgResult.first[0] > 0) ? true : false;
            return cell_label;
        }
    };
    
    if(!csg(
            args.robust_test,
            config.use_lookup,
            config.use_secondary_lookup,
            config.use_topo_ray_shooting,
            args.positive_inside,
            //
            pts, tets, funcVals, lambda,
            //
            iso_pts,
            iso_faces,
            patches,
            patch_function_label,
            patch_sign_label,
            iso_edges,
            chains,
            non_manifold_edges_of_vert,
            shells,
            arrangement_cells,
            cell_function_label,
            timing_labels,timings,
            stats_labels,stats)) {
                return -1;
            }
       if (args.robust_test) return 0;
       
       // save result
       if (!args.timing_only) {
        save_result_CSG(config.output_dir + "/mesh.json",
                    iso_pts,
                    iso_faces,
                    patches,
                    patch_sign_label,
                    iso_edges,
                    chains,
                    non_manifold_edges_of_vert);
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
