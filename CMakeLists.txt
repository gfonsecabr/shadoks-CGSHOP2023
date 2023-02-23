# cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=~/Documents/src/socg23challenge/cgal/ -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio221/ -DCMAKE_MODULE_PATH=.  CMakeLists.txt

cmake_minimum_required(VERSION 3.1...3.23)
project( solver )

add_compile_options(-fsplit-stack)
add_link_options(-fsplit-stack)
add_compile_options(-mpopcnt)
add_link_options(-mpopcnt)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )

find_package(CPLEX)
include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})

macro (add_executable _name)
    # invoke built-in add_executable
    _add_executable(${ARGV})
    if (TARGET ${_name})
        set_property(TARGET ${_name} PROPERTY CXX_STANDARD 20)
        target_link_libraries( ${_name} PRIVATE ${CPLEX_LIBRARIES})
        target_link_libraries( ${_name} PRIVATE dl)
    endif()
endmacro()

create_single_source_cgal_program( "solver.cpp" )

