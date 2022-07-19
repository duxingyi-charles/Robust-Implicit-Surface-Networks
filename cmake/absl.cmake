if (TARGET absl::base)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    absl
    GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
    GIT_TAG 20211102.0
    )

set(OLD_CXX_STANDARD ${CMAKE_CXX_STANDARD})
set(OLD_CXX_FLAGS ${CMAKE_CXX_FLAGS})
if (NOT MSVC)
    set(CMAKE_CXX_FLAGS "-Wno-unknown-warning-option")
endif()
set(CMAKE_CXX_STANDARD 17)
set(ABSL_PROPAGATE_CXX_STD On CACHE BOOL "")
FetchContent_MakeAvailable(absl)
set(CMAKE_CXX_FLAGS ${OLD_CXX_FLAGS})
set(CMAKE_CXX_STANDARD ${OLD_CXX_STANDARD})
