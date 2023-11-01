cmake_minimum_required(VERSION 3.23)
project(3_SSTree)

set(SFML_STATIC_LIBRARIES TRUE)
set(SFML_DIR D:/SFML/lib/cmake/SFML)
find_package(SFML COMPONENTS system window graphics audio network REQUIRED)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(3_SSTree
        indexing.cpp
        params.h
        Point.h
        SStree.cpp
        SStree.h)

target_include_directories(3_SSTree PRIVATE ${SFML_INCLUDE_DIR})
target_link_libraries(3_SSTree sfml-graphics sfml-window sfml-system)