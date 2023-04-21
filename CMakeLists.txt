cmake_minimum_required(VERSION 3.23)
project(CoilCalculations)
set(CMAKE_CXX_STANDARD 14)

find_package(Matplot++ REQUIRED)
find_package(GSL REQUIRED)
include_directories("${GSL_INCLUDE_DIRS}")

add_executable(CoilCalculations data_visualization.cpp)

target_link_libraries(CoilCalculations "/usr/local/Cellar/cuba/4.2.2/lib/libcuba.a")
target_link_libraries(CoilCalculations Matplot++::matplot)
target_link_libraries(CoilCalculations ${GSL_LIBRARIES})