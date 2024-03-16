FILE(GLOB VOX_FILES *.cpp *.h)

add_library(Voxeloffset
  ${VOX_FILES}
)

target_link_libraries(Voxeloffset PUBLIC 
  # CGAL
  CGAL::CGAL
  # Eigen
  CGAL::Eigen3_support
  # glog
  glog::glog
  # OpenMP
  OpenMP::OpenMP_CXX
)

target_include_directories(Voxeloffset PUBLIC ..)

add_library(Voxeloffset::Voxeloffset ALIAS Voxeloffset)