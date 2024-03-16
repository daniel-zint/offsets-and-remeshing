FILE(GLOB PROJ_FILES *.cpp *.h)

add_subdirectory(voxeloffset)
add_subdirectory(voxeloffset_main)
add_subdirectory(selfinter_repair)

cmake_dependent_option(BUILD_CLI "Build command line interface for voxeloffset" ON "CLI11_FOUND" OFF)
if(BUILD_CLI)
	add_subdirectory(voxeloffset_cli)
	add_subdirectory(mesh_cleanup)
endif()