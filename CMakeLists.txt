cmake_minimum_required(VERSION 3.14)

project(oiff
    VERSION 1.0.0
    DESCRIPTION "Optimizing the independent filter with FDR control"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(oiff INTERFACE)

target_include_directories(oiff INTERFACE include/)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()
