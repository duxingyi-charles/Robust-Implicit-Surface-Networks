//
// Created by Charles Du on 7/20/22.
//

#include "io.h"
#include "msh_io.h"

#include <fstream>
#include <iostream>
#include <ghc/filesystem.hpp>

void process_path(const ghc::filesystem::path& base, ghc::filesystem::path& p)
{
    if (p.is_relative()) {
        p = ghc::filesystem::absolute(base / p);
    }
}

bool parse_config_file(const std::string& filename,
                       std::string& tet_mesh_file,
                       std::string& func_file,
                       std::string& output_dir,
                       bool& use_lookup,
                       bool& use_2func_lookup,
                       bool& use_topo_ray_shooting,
                       size_t& tet_mesh_resolution)
{
    using json = nlohmann::json;
    namespace fs = ghc::filesystem;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "configure file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();

    fs::path config_file(filename);
    fs::path config_path = config_file.parent_path();

    if (data.contains("tetMeshFile")) {
        fs::path tet_file(data["tetMeshFile"]);
        process_path(config_path, tet_file);
        tet_mesh_file = tet_file.string();
        tet_mesh_resolution = 0;
    } else {
        assert(data.contains("resolution"));
        tet_mesh_file = "";
        tet_mesh_resolution = data["resolution"];
    }

    fs::path function_file(data["funcFile"]);
    process_path(config_path, function_file);
    func_file = function_file.string();

    fs::path out_dir(data["outputDir"]);
    process_path(config_path, out_dir);
    output_dir = out_dir.string();

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
                          bool& use_topo_ray_shooting,
                          size_t& tet_mesh_resolution)
{
    using json = nlohmann::json;
    namespace fs = ghc::filesystem;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "configure file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();

    fs::path config_file(filename);
    fs::path config_path = config_file.parent_path();

    if (data.contains("tetMeshFile")) {
        fs::path tet_file(data["tetMeshFile"]);
        process_path(config_path, tet_file);
        tet_mesh_file = tet_file.string();
        tet_mesh_resolution = 0;
    } else {
        assert(data.contains("resolution"));
        tet_mesh_file = "";
        tet_mesh_resolution = data["resolution"];
    }

    fs::path mat_file(data["materialFile"]);
    process_path(config_path, mat_file);
    material_file = mat_file.string();

    fs::path out_dir(data["outputDir"]);
    process_path(config_path, out_dir);
    output_dir = out_dir.string();

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

bool generate_tet_mesh(size_t resolution,
                   std::vector<std::array<double, 3>> &pts,
                   std::vector<std::array<size_t, 4>> &tets)
{
    if (resolution == 0) return false;
    const size_t N = resolution + 1;
    pts.resize(N * N * N);
    for (size_t i=0; i<N; i++) {
        double x = (double)i / (double)(N-1) * 2 - 1;
        for (size_t j=0; j<N; j++) {
            double y = (double)j / (double)(N-1) * 2 - 1;
            for (size_t k=0; k<N; k++) {
                double z = (double)k / (double)(N-1) * 2 - 1;

                size_t idx = i*N*N + j*N + k;
                pts[idx] = {{x, y, z}};
            }
        }
    }

    tets.resize(resolution * resolution * resolution * 5);
    for (size_t i=0; i<resolution; i++) {
        for (size_t j=0; j<resolution; j++) {
            for (size_t k=0; k<resolution; k++) {
                size_t idx = (i * resolution * resolution + j * resolution + k) * 5;
                size_t v0 = i*N*N + j*N + k;
                size_t v1 = (i+1)*N*N + j*N + k;
                size_t v2 = (i+1)*N*N + (j+1)*N + k;
                size_t v3 = i*N*N + (j+1)*N + k;
                size_t v4 = i*N*N + j*N + k + 1;
                size_t v5 = (i+1)*N*N + j*N + k + 1;
                size_t v6 = (i+1)*N*N + (j+1)*N + k + 1;
                size_t v7 = i*N*N + (j+1)*N + k + 1;

                if ((i + j + k) % 2 == 0) {
                    tets[idx] = {{v4, v6, v1, v3}};
                    tets[idx + 1] = {{v6, v3, v4, v7}};
                    tets[idx + 2] = {{v1, v3, v0, v4}};
                    tets[idx + 3] = {{v3, v1, v2, v6}};
                    tets[idx + 4] = {{v4, v1, v6, v5}};
                } else {
                    tets[idx] = {{v7, v0, v2, v5}};
                    tets[idx + 1] = {{v2, v3, v0, v7}};
                    tets[idx + 2] = {{v5, v7, v0, v4}};
                    tets[idx + 3] = {{v7, v2, v6, v5}};
                    tets[idx + 4] = {{v0, v1, v2, v5}};
                }
            }
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
                 const std::vector<std::vector<size_t>>& shells,
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
    json jShells;
    for (const auto& shell : shells) {
        jShells.push_back(json(shell));
    }
    //
    json jCells;
    for (const auto& cell : cells) {
        jCells.push_back(json(cell));
    }
    //
    json jOut;
    jOut["points"] = jPts;
    jOut["faces"] = jFaces;
    jOut["patches"] = jPatches;
    jOut["edges"] = jEdges;
    jOut["chains"] = jChains;
    jOut["corners"] = jCorners;
    jOut["shells"] = jShells;
    jOut["cells"] = jCells;
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
