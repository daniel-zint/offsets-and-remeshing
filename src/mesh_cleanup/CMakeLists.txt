FILE(GLOB MAIN_FILES *.cpp *.h)

add_executable(Mesh_cleanup 
  ${MAIN_FILES}
)

target_link_libraries(Mesh_cleanup PRIVATE 
    # CGAL
  CGAL::CGAL
  # Eigen
  CGAL::Eigen3_support
  # glog
  glog::glog
  # CLI11
  CLI11::CLI11
)