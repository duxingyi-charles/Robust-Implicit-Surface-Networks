//
// Created by Charles Du on 7/20/22.
//

#include "io.h"
#include "msh_io.h"

#include <fstream>
#include <iostream>

bool parse_config_file(const std::string& filename,
                       std::string& tet_mesh_file,
                       std::string& func_file,
                       std::string& output_dir,
                       bool& use_lookup,
                       bool& use_2func_lookup,
                       bool& use_topo_ray_shooting)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "configure file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    tet_mesh_file = data["tetMeshFile"];
    func_file = data["funcFile"];
    output_dir = data["outputDir"];
    use_lookup = data["useLookup"];
    use_2func_lookup = data["use2funcLookup"];
    use_topo_ray_shooting = data["useTopoRayShooting"];
    return true;
}

bool parse_config_file_MI(const std::string& filename,
                          std::string& tet_mesh_file,
                          std::string& material_file,
                          std::string& output_dir,
                          bool& use_lookup,
                          bool& use_3func_lookup,
                          bool& use_topo_ray_shooting)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "configure file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    tet_mesh_file = data["tetMeshFile"];
    material_file = data["materialFile"];
    output_dir = data["outputDir"];
    use_lookup = data["useLookup"];
    use_3func_lookup = data["use3funcLookup"];
    use_topo_ray_shooting = data["useTopoRayShooting"];
    return true;
}

bool load_tet_mesh(const std::string& filename,
                   std::vector<std::array<double, 3>>& pts,
                   std::vector<std::array<size_t, 4>>& tets)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "tet mesh file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    pts.resize(data[0].size());
    for (size_t j = 0; j < pts.size(); j++) {
        for (size_t k = 0; k < 3; k++) {
            pts[j][k] = data[0][j][k].get<double>();
        }
    }
    //
    tets.resize(data[1].size());
    for (size_t j = 0; j < tets.size(); j++) {
        for (size_t k = 0; k < 4; k++) {
            tets[j][k] = data[1][j][k].get<size_t>();
        }
    }
    return true;
}

bool save_result(const std::string& filename,
                 const std::vector<std::array<double, 3>>& mesh_pts,
                 const std::vector<PolygonFace>& mesh_faces,
                 const std::vector<std::vector<size_t>>& patches,
                 const std::vector<Edge>& edges,
                 const std::vector<std::vector<size_t>>& chains,
                 const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                 const std::vector<std::vector<std::pair<std::pair<size_t, int>,std::pair<size_t, int>>>>& half_patch_pair_list,
                 const std::vector<std::vector<size_t>>& shells,
                 const std::vector<std::vector<size_t>>& components,
                 const std::vector<std::vector<size_t>>& cells)
{
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jPts;
    for (const auto& iso_pt : mesh_pts) {
        jPts.push_back(json(iso_pt));
    }
    //
    json jFaces;
    for (const auto& iso_face : mesh_faces) {
        jFaces.push_back(json(iso_face.vert_indices));
    }
    //
    json jPatches;
    for (const auto& patch : patches) {
        jPatches.push_back(json(patch));
    }
    //
    json jEdges;
    for (const auto& edge : edges) {
        jEdges.push_back({edge.v1, edge.v2});
    }
    //
    json jChains;
    for (const auto& chain : chains) {
        jChains.push_back(json(chain));
    }
    //
    json jCorners;
    for (size_t i = 0; i < non_manifold_edges_of_vert.size(); i++) {
        if (non_manifold_edges_of_vert[i].size() > 2 || non_manifold_edges_of_vert[i].size() == 1) {
            jCorners.push_back(i);
        }
    }
    //
//    json jHalfPatchsList;
//    for (const auto& i : half_patch_list) {
//        json jHalfPatchs;
//        for (const auto& j : i) {
//            jHalfPatchs.push_back(json(j));
//        }
//        jHalfPatchsList.push_back(jHalfPatchs);
//    }
    //
    json jHalfPatchPairList;
    for (const auto& i : half_patch_pair_list) {
        json jHalfPatchPairs;
        for (const auto& j : i) {
            jHalfPatchPairs.push_back(json(j));
        }
        jHalfPatchPairList.push_back(jHalfPatchPairs);
    }
    //
    json jShells;
    for (const auto& shell : shells) {
        jShells.push_back(json(shell));
    }
    //
    json jComponents;
    for (const auto& component : components) {
        jComponents.push_back(json(component));
    }
    //
    json jArrCells;
    for (const auto& arrangement_cell : cells) {
        jArrCells.push_back(json(arrangement_cell));
    }
    //
    json jOut;
    jOut.push_back(jPts);
    jOut.push_back(jFaces);
    jOut.push_back(jPatches);
    jOut.push_back(jEdges);
    jOut.push_back(jChains);
    jOut.push_back(jCorners);
//    jOut.push_back(jHalfPatchsList);
    jOut.push_back(jHalfPatchPairList);
    jOut.push_back(jShells);
    jOut.push_back(jComponents);
    jOut.push_back(jArrCells);
    fout << jOut << std::endl;
    fout.close();
    return true;
}

bool save_result_msh(const std::string& filename,
                     const std::vector<std::array<double, 3>>& mesh_pts,
                     const std::vector<PolygonFace>& mesh_faces,
                     const std::vector<std::vector<size_t>>& patches,
                     const std::vector<Edge>& edges,
                     const std::vector<std::vector<size_t>>& chains,
                     const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
                     const std::vector<std::vector<size_t>>& shells,
                     const std::vector<std::vector<size_t>>& components,
                     const std::vector<std::vector<size_t>>& cells)
{
    constexpr size_t INVALID = std::numeric_limits<size_t>::max();
    std::vector<size_t> vertex_map(mesh_pts.size(), INVALID);
    wmtk::MshData msh, msh2, msh3;

    // Save chains.
    auto extract_chain = [&](size_t i) {
        std::fill(vertex_map.begin(), vertex_map.end(), INVALID);
        const auto& chain = chains[i];
        size_t num_vertices = 0;
        for (auto eid : chain) {
            const auto& e = edges[eid];
            if (vertex_map[e.v1] == INVALID) {
                vertex_map[e.v1] = num_vertices;
                num_vertices++;
            }
            if (vertex_map[e.v2] == INVALID) {
                vertex_map[e.v2] = num_vertices;
                num_vertices++;
            }
        }

        std::vector<std::array<double, 3>> vertices(num_vertices);
        for (size_t j = 0; j < vertex_map.size(); j++) {
            if (vertex_map[j] == INVALID) continue;
            vertices[vertex_map[j]] = mesh_pts[j];
        }

        msh.add_edge_vertices(num_vertices, [&](size_t j) { return vertices[j]; });
        msh.add_edges(chain.size(), [&](size_t j) -> std::array<size_t, 2> {
            const auto eid = chain[j];
            const auto& e = edges[eid];
            return {vertex_map[e.v1], vertex_map[e.v2]};
        });
    };
    for (size_t i = 0; i < chains.size(); i++) {
        extract_chain(i);
    }
    msh.save(filename + "_chains.msh");

    // Save patches
    auto extract_patch = [&](size_t i) -> std::tuple<std::vector<std::array<double, 3>>,
            std::vector<std::array<size_t, 3>>,
            std::vector<size_t>> {
        std::fill(vertex_map.begin(), vertex_map.end(), INVALID);
        std::vector<std::array<double, 3>> vertices;
        std::vector<std::array<size_t, 3>> triangles;
        std::vector<size_t> polygon_ids;

        size_t num_vertices = 0;
        auto add_vertex = [&](size_t vi) {
            if (vertex_map[vi] == INVALID) {
                vertex_map[vi] = num_vertices;
                num_vertices++;
            }
        };

        const auto& patch = patches[i];
        triangles.reserve(patch.size() * 2);
        polygon_ids.reserve(patch.size() * 2);

        for (const auto poly_id : patch) {
            const auto& f = mesh_faces[poly_id];
            const size_t s = f.vert_indices.size();
            add_vertex(f.vert_indices[0]);
            add_vertex(f.vert_indices[1]);
            for (size_t j = 2; j < s; j++) {
                triangles.push_back({f.vert_indices[0], f.vert_indices[j - 1], f.vert_indices[j]});
                polygon_ids.push_back(poly_id);
                add_vertex(f.vert_indices[j]);
            }
        }

        vertices.resize(num_vertices);
        for (size_t j = 0; j < mesh_pts.size(); j++) {
            if (vertex_map[j] == INVALID) continue;
            const auto& p = mesh_pts[j];
            vertices[vertex_map[j]] = p;
        }

        for (auto& t : triangles) {
            t[0] = vertex_map[t[0]];
            t[1] = vertex_map[t[1]];
            t[2] = vertex_map[t[2]];
        }
        return {vertices, triangles, polygon_ids};
    };

    std::vector<size_t> patch_ids;
    std::vector<size_t> polygon_ids;
    for (size_t i = 0; i < patches.size(); i++) {
        auto r = extract_patch(i);
        const auto& vertices = std::get<0>(r);
        const auto& triangles = std::get<1>(r);
        const auto& local_polygon_ids = std::get<2>(r);
        msh2.add_face_vertices(vertices.size(), [&](size_t j) { return vertices[j]; });
        msh2.add_faces(triangles.size(), [&](size_t j) { return triangles[j]; });
        patch_ids.insert(patch_ids.end(), triangles.size(), i);
        polygon_ids.insert(polygon_ids.end(), local_polygon_ids.begin(), local_polygon_ids.end());
    }

    if (patches.size() > 0) {
        msh2.add_face_attribute<1>("patch_id", [&](size_t i) { return patch_ids[i]; });
        msh2.add_face_attribute<1>("polygon_id", [&](size_t i) { return polygon_ids[i]; });
    }
    msh2.save(filename + "_patches.msh");

    auto extract_cell = [&](size_t i) {
        const auto& cell = cells[i];

        std::fill(vertex_map.begin(), vertex_map.end(), INVALID);
        std::vector<std::array<double, 3>> vertices;
        std::vector<std::array<size_t, 3>> triangles;

        size_t num_vertices = 0;
        auto add_vertex = [&](size_t vi) {
            if (vertex_map[vi] == INVALID) {
                vertex_map[vi] = num_vertices;
                num_vertices++;
            }
        };

        for (auto shell_id : cell) {
            const auto& shell = shells[shell_id];
            for (auto half_patch_id : shell) {
                size_t patch_id = half_patch_id / 2;
                const auto& patch = patches[patch_id];
                for (auto poly_id : patch) {
                    const auto& f = mesh_faces[poly_id];
                    const size_t s = f.vert_indices.size();
                    add_vertex(f.vert_indices[0]);
                    add_vertex(f.vert_indices[1]);
                    for (size_t j = 2; j < s; j++) {
                        triangles.push_back(
                                {f.vert_indices[0], f.vert_indices[j - 1], f.vert_indices[j]});
                        add_vertex(f.vert_indices[j]);
                    }
                }
            }
        }

        vertices.resize(num_vertices);
        for (size_t j = 0; j < mesh_pts.size(); j++) {
            if (vertex_map[j] == INVALID) continue;
            const auto& p = mesh_pts[j];
            vertices[vertex_map[j]] = p;
        }

        for (auto& t : triangles) {
            t[0] = vertex_map[t[0]];
            t[1] = vertex_map[t[1]];
            t[2] = vertex_map[t[2]];
        }

        msh3.add_face_vertices(num_vertices, [&](size_t j) { return vertices[j]; });
        msh3.add_faces(triangles.size(), [&](size_t j) { return triangles[j]; });
        return triangles.size();
    };

    std::vector<size_t> cell_ids;
    for (size_t i = 0; i < cells.size(); i++) {
        size_t num_faces = extract_cell(i);
        cell_ids.insert(cell_ids.end(), num_faces, i);
    }
    msh3.add_face_attribute<1>("cell_id", [&](size_t i) { return cell_ids[i]; });
    msh3.save(filename + "_cells.msh");
    return true;
}


bool save_timings(const std::string& filename,
                  const std::vector<std::string>& timing_labels,
                  const std::vector<double>& timings)
{
    // assert timing_labels.size() == timings.size()
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jOut;
    for (size_t i = 0; i < timings.size(); ++i) {
        jOut[timing_labels[i]] = timings[i];
    }
    //
    fout << jOut << std::endl;
    fout.close();
    return true;
}

bool save_statistics(const std::string& filename,
                     const std::vector<std::string>& stats_labels,
                     const std::vector<size_t>& stats)
{
    // assert stats_labels.size() == stats.size()
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jOut;
    for (size_t i = 0; i < stats.size(); ++i) {
        jOut[stats_labels[i]] = stats[i];
    }
    //
    fout << jOut << std::endl;
    fout.close();
    return true;
}
