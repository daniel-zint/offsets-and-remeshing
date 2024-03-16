FILE(GLOB MAIN_FILES *.cpp *.h)

add_executable(Voxeloffset_main 
  ${MAIN_FILES}
)

target_link_libraries(Voxeloffset_main PRIVATE 
  Voxeloffset::Voxeloffset
)

target_compile_definitions(Voxeloffset_main PRIVATE SOURCE_DIR="${CMAKE_SOURCE_DIR}")