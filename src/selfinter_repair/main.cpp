#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <fstream>
#include <algorithm>
#include <boost/foreach.hpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor Facet_handle;
typedef boost::graph_traits<Mesh>::vertex_descriptor Vertex_handle;
typedef boost::graph_traits<Mesh>::halfedge_descriptor Halfedge_handle;

namespace PMP = CGAL::Polygon_mesh_processing;
// COMPUTES THE NUMBER OF SELF INTERSECTIONS
int main(int argc, char* argv[])
{
  const char* filename = argv[1];
  const char* outputname = (argc > 2) ? argv[2] : "out.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }
  PMP::experimental::remove_self_intersections(mesh, CGAL::parameters::preserve_genus(false));
  CGAL::IO::write_OFF(outputname, mesh );
  return 0;
}
