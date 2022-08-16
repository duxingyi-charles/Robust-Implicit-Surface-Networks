if(TARGET ghc::filesystem)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    filesystem
    GIT_REPOSITORY https://github.com/gulrak/filesystem.git
    GIT_TAG v1.5.12
)
FetchContent_MakeAvailable(filesystem)
add_library(ghc::filesystem ALIAS ghc_filesystem)
