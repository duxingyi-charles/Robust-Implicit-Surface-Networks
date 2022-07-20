if (TARGET simplicial_arrangement::simplicial_arrangement)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    simplicial_arrangement
    GIT_REPOSITORY git@github.com:qnzhou/simplicial_arrangement.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(simplicial_arrangement)
