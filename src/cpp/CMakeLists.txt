cmake_minimum_required(VERSION 3.14)
project(sgp4_propagator_plugin CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# dnwrnr/sgp4 library
file(GLOB SGP4_SOURCES deps/sgp4/libsgp4/*.cc)
add_library(sgp4 STATIC ${SGP4_SOURCES})
target_include_directories(sgp4 PUBLIC deps/sgp4/libsgp4)

# Plugin library
add_library(sgp4_prop STATIC
    src/gp_parser.cpp
    src/propagator.cpp
)
target_include_directories(sgp4_prop PUBLIC include deps/sgp4/libsgp4)
target_link_libraries(sgp4_prop sgp4)

# Tests
add_executable(test_sgp4_prop tests/test_propagator.cpp)
target_link_libraries(test_sgp4_prop sgp4_prop)
target_include_directories(test_sgp4_prop PRIVATE include)
