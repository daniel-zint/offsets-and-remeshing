#pragma once

#include "Remeshing.h"
#include "types.h"

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <vector>

namespace Quality_measurements {

    inline bool does_self_intersect( const SurfaceMesh& m ) { return PMP::does_self_intersect( m ); }

    inline double hausdorff( const SurfaceMesh& offset, const SurfaceMesh& primal_mesh ) {
        return PMP::approximate_Hausdorff_distance<CGAL::Parallel_if_available_tag, SurfaceMesh>( offset, primal_mesh );
    }

    inline void normal_deviation( SurfaceMesh& offset, const SurfaceMesh& primal_mesh, const FT& iso_value ) {
        AABB_tree tree_primal( faces( primal_mesh ).first, faces( primal_mesh ).second, primal_mesh );
        AABB_tree tree_offset( faces( offset ).first, faces( offset ).second, offset );
        // std::vector<Point> sample_points;
        // PMP::sample_triangle_mesh( offset, std::back_inserter( sample_points ),
        //                            CGAL::parameters::do_sample_edges( false ).do_sample_vertices( false ) );

        auto [f_samples, f_samples_created] = offset.add_property_map<Face_descriptor, std::vector<Point>>( "f:samples_qt" );
        LOG_ASSERT( f_samples_created );
        for( Face_descriptor f: faces( offset ) ) {
            // std::array<Face_descriptor, 1> arr { f };
            // CGAL::Face_filtered_graph ffg( offset, arr );
            // PMP::sample_triangle_mesh(
            //     ffg, std::back_inserter( f_samples[f] ),
            //     CGAL::parameters::number_of_points_on_faces( 10 ).do_sample_edges( false ).do_sample_vertices( false ) );

            const Halfedge_descriptor h = offset.halfedge( f );
            const Vertex_descriptor v0  = offset.target( h );
            const Vertex_descriptor v1  = offset.target( offset.next( h ) );
            const Vertex_descriptor v2  = offset.source( h );
            const Point p0              = offset.point( v0 );
            const Point p1              = offset.point( v1 );
            const Point p2              = offset.point( v2 );

            const Vector p_vec0 = p0 - CGAL::ORIGIN;
            const Vector p_vec1 = p1 - CGAL::ORIGIN;
            const Vector p_vec2 = p2 - CGAL::ORIGIN;

            const Point s0 = CGAL::ORIGIN + ( 0.8 * p_vec0 + 0.1 * p_vec1 + 0.1 * p_vec2 );
            const Point s1 = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.8 * p_vec1 + 0.1 * p_vec2 );
            const Point s2 = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.1 * p_vec1 + 0.8 * p_vec2 );
            const Point s3 = CGAL::centroid( p0, p1, p2 );
            f_samples[f]   = { s0, s1, s2, s3 };
        }

        auto [f_normal, f_normal_created] = offset.add_property_map<Face_descriptor, Vector>( "f:normal_qt" );
        LOG_ASSERT( f_normal_created );
        PMP::compute_face_normals( offset, f_normal );

        // auto [f_nd, f_nd_created] = offset.add_property_map<Face_descriptor, double>( "f:nd" );
        // LOG_ASSERT( f_nd_created );
        std::vector<double> normal_deviations( offset.number_of_faces(), 0 );
        for( const Face_descriptor& f: offset.faces() ) {
            const auto& samples      = f_samples[f];
            const Vector& n_triangle = f_normal[f];
            for( const Point& p: samples ) {
                Point p_proj = tree_primal.closest_point( p );

                AABB_tree::Point_and_primitive_id papid_primal = tree_primal.closest_point_and_primitive( p );
                Face_descriptor triangle_primal                = papid_primal.second;

                Vector n_offset;
                if( iso_value == 0 ) {
                    n_offset = PMP::compute_face_normal( triangle_primal, primal_mesh );
                } else {
                    const auto [idx, type] = Distance::closest_triangle_primitive( p, triangle_primal, primal_mesh );

                    switch( type ) {
                    case 0:    // vertex
                    {
                        const Vertex_descriptor v( idx );
                        const auto s     = Distance::get_primitive_offset( primal_mesh, v, iso_value );
                        const Point proj = s.project( p );
                        n_offset         = normalize( proj - s.origin );
                        break;
                    }
                    case 1:    // edge
                    {
                        Edge_descriptor e( idx );
                        const auto c          = Distance::get_primitive_offset( primal_mesh, e, iso_value );
                        const Point proj      = c.project( p );
                        const auto proj_plane = c.projection_plane( p );
                        const Vector normal( proj_plane[0], proj_plane[1], proj_plane[2] );
                        n_offset = normal;
                        break;
                    }
                    case 2:    // triangle
                    {
                        const auto plane = Distance::get_primitive_offset( primal_mesh, triangle_primal, iso_value, p );
                        const Point proj = plane.project( p );
                        n_offset         = plane.normal;
                        break;
                    }
                    default:
                        LOG( FATAL ) << "Unknown id type: " << type;
                        return;
                    }
                }
                const double nd      = CGAL::approximate_angle( n_triangle, n_offset );
                normal_deviations[f] = CGAL::max( normal_deviations[f], nd );
            }

            // if( normal_deviations[f] > 45 ) {
            //     std::cout << "Face " << f << " has a normal deviation of " << normal_deviations[f] << std::endl;
            //     std::cout << "  pos = " << samples[3] << std::endl;
            //     for( const auto& v: offset.vertices_around_face( offset.halfedge( f ) ) ) {
            //         std::cout << "  v " << v << std::endl;
            //     }
            // }
        }

        std::sort( normal_deviations.begin(), normal_deviations.end(), std::greater<>() );
        for( size_t i = 0; i < 10; ++i ) {
            std::cout << normal_deviations[i] << std::endl;
        }

        offset.remove_property_map( f_samples );
        offset.remove_property_map( f_normal );
    }

    inline std::vector<FT> distance_to_real_offset( const SurfaceMesh& offset, const SurfaceMesh& primal_mesh, const FT& iso_value ) {
        // sample offset
        std::vector<Point> sample_points;
        PMP::sample_triangle_mesh( offset, std::back_inserter( sample_points ) );

        // get distances
        AABB_tree tree( faces( primal_mesh ).first, faces( primal_mesh ).second, primal_mesh );
        // AABB_tree::Point_and_primitive_id papid = tree.closest_point_and_primitive( p );
        std::vector<FT> distances;
        distances.reserve( sample_points.size() );

        for( const Point& p: sample_points ) {
            const Point pp = tree.closest_point( p );
            const FT d     = CGAL::sqrt( ( p - pp ).squared_length() );
            distances.push_back( d - iso_value );
        }

        return distances;
    }

    inline void analyze_distance( const SurfaceMesh& offset, const SurfaceMesh& primal_mesh, const FT& iso_value, std::vector<FT> sample_distances ) {
        std::sort( sample_distances.begin(), sample_distances.end() );

        for( auto& d: sample_distances ) {
            d = 100 * d / iso_value;
        }

        std::cout << "Distance Information (in %)" << std::endl;
        std::cout << "  min = " << sample_distances[0] << std::endl;
        std::cout << "  max = " << sample_distances[sample_distances.size() - 1] << std::endl;
        std::cout << "  median = " << sample_distances[( sample_distances.size() - 1 ) / 2] << std::endl;
        double sum  = std::accumulate( sample_distances.begin(), sample_distances.end(), 0.0 );
        double mean = sum / sample_distances.size();
        std::cout << "  mean = " << mean << std::endl;

        std::vector<double> diff( sample_distances.size() );
        std::transform( sample_distances.begin(), sample_distances.end(), diff.begin(), [mean]( double x ) { return x - mean; } );
        double sq_sum = std::inner_product( diff.begin(), diff.end(), diff.begin(), 0.0 );
        double stdev  = std::sqrt( sq_sum / sample_distances.size() );
        std::cout << "  standard deviation = " << stdev << std::endl;
    }

}    // namespace Quality_measurements