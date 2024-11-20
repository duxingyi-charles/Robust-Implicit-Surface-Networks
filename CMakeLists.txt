cmake_minimum_required(VERSION 3.14)

project(robust_implicit_networks LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(cli11)
include(nlohmann-json)
include(Eigen3)
include(mshio)
include(filesystem)
include(implicit_functions)
include(absl)
include(simplicial_arrangement)


# project library
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp ${CMAKE_CURRENT_LIST_DIR}/src/*.h)

add_library(robust_implicit_networks STATIC ${SRC_FILES})
target_include_directories(robust_implicit_networks PUBLIC
     ${CMAKE_CURRENT_LIST_DIR}/src)
target_compile_features(robust_implicit_networks PRIVATE cxx_std_20)
if (WIN32)
    target_compile_definitions(robust_implicit_networks PRIVATE NOMINMAX)
endif()
target_link_libraries(robust_implicit_networks PUBLIC
    simplicial_arrangement::simplicial_arrangement
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    absl::flat_hash_map
    absl::flat_hash_set
    mshio::mshio ghc::filesystem)


add_executable(impl_arrangement "app/implicit_arrangement.cpp")
target_compile_features(impl_arrangement PRIVATE cxx_std_20)
target_link_libraries(impl_arrangement simplicial_arrangement::simplicial_arrangement
    CLI11::CLI11 Eigen3::Eigen
    robust_implicit_networks implicit_functions::implicit_functions)

add_executable(csg "app/csg.cpp")
target_compile_features(csg PRIVATE cxx_std_20)
target_link_libraries(csg simplicial_arrangement::simplicial_arrangement
    CLI11::CLI11 Eigen3::Eigen
    robust_implicit_networks implicit_functions::implicit_functions)

 add_executable(material_interface "app/material_interface.cpp")
 target_compile_features(material_interface PRIVATE cxx_std_20)
 target_link_libraries(material_interface simplicial_arrangement::simplicial_arrangement
     CLI11::CLI11 Eigen3::Eigen
     robust_implicit_networks implicit_functions::implicit_functions)

option(IMPLICIT_NETWORK_UNIT_TESTS "Build implicit network unit tests" Off)
if (IMPLICIT_NETWORK_UNIT_TESTS)
    include(CTest)
    enable_testing()
    include(Catch2)

    file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
    add_executable(implicit_network_tests ${TEST_FILES})
    target_link_libraries(implicit_network_tests PRIVATE
             Catch2::Catch2 simplicial_arrangement::simplicial_arrangement
             Eigen3::Eigen robust_implicit_networks implicit_functions::implicit_functions)
    target_compile_features(implicit_network_tests PRIVATE cxx_std_20)
    target_compile_definitions(implicit_network_tests PRIVATE
        -DTEST_FILE="${CMAKE_CURRENT_LIST_DIR}/examples/tests"
    )
    catch_discover_tests(implicit_network_tests)
endif()
