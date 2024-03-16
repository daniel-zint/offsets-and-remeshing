FILE(GLOB MAIN_FILES *.cpp *.h)

add_executable(Selfinter_repair 
  ${MAIN_FILES}
)

target_link_libraries(Selfinter_repair PRIVATE 
    # CGAL
  CGAL::CGAL
  # Eigen
  CGAL::Eigen3_support
  # glog
  glog::glog
  # CLI11
  CLI11::CLI11
)