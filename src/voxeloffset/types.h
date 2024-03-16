#pragma once

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#define TAG CGAL::Parallel_if_available_tag

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Triangle_3 Triangle;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

// Metamesh (=final output)
typedef CGAL::Surface_mesh<Point> SurfaceMesh;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor Vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor Halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor Edge_descriptor;
typedef SurfaceMesh::Face_index Face_descriptor;

typedef CGAL::Plane_3<Kernel> Plane_3;
namespace PMP = CGAL::Polygon_mesh_processing;

constexpr size_t INVALID_UNIQUE_ID = std::numeric_limits<size_t>::max();

constexpr FT FT_MAX = std::numeric_limits<FT>::max();

typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh> AABB_primitive;
typedef CGAL::AABB_traits<Kernel, AABB_primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> AABB_tree;

typedef CGAL::Side_of_triangle_mesh<SurfaceMesh, CGAL::GetGeomTraits<SurfaceMesh>::type> Sotm;
