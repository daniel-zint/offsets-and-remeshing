FILE(GLOB MAIN_FILES *.cpp *.h)

add_executable(Voxeloffset_cli 
  ${MAIN_FILES}
)

target_link_libraries(Voxeloffset_cli PRIVATE 
  Voxeloffset::Voxeloffset
  CLI11::CLI11
)

target_compile_definitions(Voxeloffset_cli PRIVATE SOURCE_DIR="${CMAKE_SOURCE_DIR}")