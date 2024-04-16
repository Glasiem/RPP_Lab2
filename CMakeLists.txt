cmake_minimum_required(VERSION 3.27)
project(RPP_Lab2)

set(CMAKE_CXX_STANDARD 17)
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
    message("MPI found")
else()
    message("MPI Not found")
endif ()

add_executable(RPP_Lab2 main.cpp)
target_link_libraries(RPP_Lab2 PUBLIC MPI::MPI_CXX)
