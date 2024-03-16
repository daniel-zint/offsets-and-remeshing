cmake_minimum_required(VERSION 3.1...3.15)
project(voxeloffset LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

include(CMakeDependentOption)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Find CGAL and CGAL Qt5
find_package(CGAL REQUIRED)
# Find Eigen
find_package(Eigen3 REQUIRED) #(requires 3.1.0 or greater)
# glog
find_package(glog CONFIG REQUIRED)
# OpenMP
find_package(OpenMP)
# CLI
find_package(CLI11)

include(CGAL_Eigen3_support)

if(MSVC)
    # enable parallel build in Visual Studio
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

add_subdirectory(src)
