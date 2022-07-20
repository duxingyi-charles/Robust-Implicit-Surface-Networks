
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(nlohmann-json)
include(Eigen3)

add_library(implicit_functions STATIC "implicit_functions.h" "implicit_functions.cpp")
target_link_libraries(implicit_functions PRIVATE nlohmann_json::nlohmann_json Eigen3::Eigen)
target_compile_features(implicit_functions PRIVATE cxx_std_17)

target_include_directories(implicit_functions PUBLIC .)