cmake_minimum_required(VERSION 3.14)

project(robust_implicit_networks LANGUAGES CXX C)

include(FetchContent)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
# include(absl)
# include(implicit_predicates)
include(simplicial_arrangement)
include(cli11)
# include(nlohmann-json)



# implicit function utility
# add_subdirectory(implicit_functions)


# project library
# file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

# add_library(robust_implicit_networks STATIC ${SRC_FILES})
# target_include_directories(robust_implicit_networks PUBLIC
#     ${CMAKE_CURRENT_LIST_DIR}/src)
# target_link_libraries(robust_implicit_networks
#     PUBLIC absl::numeric absl::flat_hash_map nlohmann_json::nlohmann_json)
# target_compile_features(robust_implicit_networks PRIVATE cxx_std_17)



add_executable(impl_arrangement "app/implicit_arrangement.cpp")
target_compile_features(impl_arrangement PRIVATE cxx_std_17)
target_link_libraries(impl_arrangement simplicial_arrangement::simplicial_arrangement)
	# robust_implicit_networks implicit_functions CLI11::CLI11)
# if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
#     if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
#         target_compile_options(impl_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
#         target_link_options(impl_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
#     endif()
# endif()

# add_executable(material_interface "app/material_interface.cpp")
# target_compile_features(material_interface PRIVATE cxx_std_17)
# target_link_libraries(material_interface simplicial_arrangement::simplicial_arrangement implicit_functions CLI11::CLI11)
# if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
#     if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
#         target_compile_options(material_interface PRIVATE -fprofile-instr-generate -fcoverage-mapping)
#         target_link_options(material_interface PRIVATE -fprofile-instr-generate -fcoverage-mapping)
#     endif()
# endif()

