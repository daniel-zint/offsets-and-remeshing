#pragma once
/**
 * @file   Distance.h
 * @brief  Functions for computing the distance of a point to triangles, edges, or vertices.
 *
 * There are also functions that do not just return the distance but also the the type and index of the closest primitive, e.g., if the distance
 * between a point and a triangle is computed, the nearest point on the triangle might be located on the edge of the triangle.
 *
 * @author Daniel Zint
 * @date   December 2022
 */

#include "types.h"
#include "utilities.h"

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <array>
#include <glog/logging.h>

namespace Distance {

    /** Squared distance between two points. */
    inline FT point_to_point_squared( const Point& p1, const Point& p2 ) { return ( p2 - p1 ).squared_length(); }

    /** Distance between two points. */
    inline FT point_to_point( const Point& p1, const Point& p2 ) { return CGAL::approximate_sqrt( point_to_point_squared( p1, p2 ) ); }

    /** Squared distance between the point p and the line through (a,b). */
    inline FT point_to_line_squared( const Point& p, const Point& a, const Point& b ) {
        return CGAL::cross_product( b - a, a - p ).squared_length() / ( b - a ).squared_length();
    }

    /** Distance between the point p and the line through (a,b). */
    inline FT point_to_line( const Point& p, const Point& a, const Point& b ) {
        auto dist_squared = point_to_line_squared( p, a, b );
        return CGAL::approximate_sqrt( dist_squared );
    }

    inline Point project_to_line( const Point& p, const Point& a, const Point& b ) {
        // get dot product of e1, e2
        const Vector e1 = b - a;
        const Vector e2 = p - a;
        const FT e1e2   = CGAL::scalar_product( e1, e2 );
        return a + e1 * e1e2 / e1.squared_length();
    }

    /** Squared distance between the point p and the plane through the points (a,b,c). */
    inline FT point_to_plane_squared( const Point& p, const Point& a, const Point& b, const Point& c ) {
        return CGAL::squared_distance( p, Triangle( a, b, c ) );
    }

    /** Distance between the point p and the plane through the points (a,b,c). */
    inline FT point_to_plane( const Point& p, const Point& a, const Point& b, const Point& c ) {
        FT squared_dist = point_to_plane_squared( p, a, b, c );
        FT dist         = CGAL::approximate_sqrt( squared_dist );
        return dist;
    }

    inline Point project_to_plane( const Point& p, const Point& a, const Point& b, const Point& c ) {
        const Vector n = normalize( CGAL::cross_product( b - a, c - a ) );
        const FT dist  = CGAL::scalar_product( n, p - a );
        return p - n * dist;
    }

    /**
     * Find the closest primitive of an edge to a query point. This can be either the edge itself or one of its incident vertices.
     *
     * @param p
     * @param h
     * @param mesh
     * @return {edge index, 1}
     */
    inline std::tuple<size_t, int> closest_edge_primitive( const Point& p, const Halfedge_descriptor& h, const SurfaceMesh& mesh ) {
        const Vertex_descriptor& va = mesh.target( h );
        const Vertex_descriptor& vb = mesh.target( mesh.opposite( h ) );
        const Point& a              = mesh.point( va );
        const Point& b              = mesh.point( vb );

        // cases when point is not within the cylinder of the edge
        if( CGAL::scalar_product( p - a, b - a ) < 0 ) {
            // point maps to a
            return { va.idx(), 0 };
        }
        if( CGAL::scalar_product( p - b, a - b ) < 0 ) {
            // point maps to b
            return { vb.idx(), 0 };
        }

        return { mesh.edge( h ).idx(), 1 };
    }

    /**
     * Find the closest primitive of a triangle to a query point. This can be either the triangle itself or one of its incident edges or vertices.
     *
     * @param p query point
     * @param triangle face handle of the triangle
     * @param mesh the mesh the triangle belongs to
     * @return {index of primitive, type of primitive}. The types are: 0=vertex, 1=edge, 2=face
     */
    inline std::tuple<size_t, int> closest_triangle_primitive( const Point& p, const Face_descriptor& triangle, const SurfaceMesh& mesh ) {
        Halfedge_descriptor ha = mesh.halfedge( triangle );
        Halfedge_descriptor hb = mesh.next( ha );
        Halfedge_descriptor hc = mesh.next( hb );
        Vertex_descriptor va   = mesh.target( hb );
        Vertex_descriptor vb   = mesh.target( hc );
        Vertex_descriptor vc   = mesh.target( ha );
        Point a                = mesh.point( va );
        Point b                = mesh.point( vb );
        Point c                = mesh.point( vc );

        const Vector area_vec = CGAL::cross_product( b - a, c - a );

        Vector normal = area_vec;
        if( CGAL::scalar_product( p - a, normal ) < 0 ) {
            normal = -normal;
            std::swap( hb, hc );
            std::swap( vb, vc );
            ha = mesh.opposite( ha );
            hb = mesh.opposite( hb );
            hc = mesh.opposite( hc );
            b  = mesh.point( vb );
            c  = mesh.point( vc );
        }

        // orientation tests
        int slab_ha = on_right_slab_of_triangle_edge( p, normal, b, c );
        int slab_hb = on_right_slab_of_triangle_edge( p, normal, c, a );
        int slab_hc = on_right_slab_of_triangle_edge( p, normal, a, b );

        if( slab_ha == -1 && slab_hb == -1 && slab_hc == -1 ) {
            return { triangle.idx(), 2 };
        }

        if( slab_ha == 0 ) {
            return { mesh.edge( ha ).idx(), 1 };
        }
        if( slab_hb == 0 ) {
            return { mesh.edge( hb ).idx(), 1 };
        }
        if( slab_hc == 0 ) {
            return { mesh.edge( hc ).idx(), 1 };
        }

        // TODO shorten that down by combining if conditions
        std::vector<bool> v_is_candidate( 3, false );
        bool a_is_candidate = false;
        bool b_is_candidate = false;
        bool c_is_candidate = false;
        if( slab_ha == 1 ) {
            b_is_candidate = true;
        }
        if( slab_ha == 2 ) {
            c_is_candidate = true;
        }
        if( slab_hb == 1 ) {
            c_is_candidate = true;
        }
        if( slab_hb == 2 ) {
            a_is_candidate = true;
        }
        if( slab_hc == 1 ) {
            a_is_candidate = true;
        }
        if( slab_hc == 2 ) {
            b_is_candidate = true;
        }

        if( !a_is_candidate && !b_is_candidate && !c_is_candidate ) {
            return { triangle.idx(), 2 };
        }

        Vertex_descriptor v_closest;
        FT d = FT_MAX;
        if( a_is_candidate ) {
            const FT dist = ( p - mesh.point( va ) ).squared_length();
            if( dist < d ) {
                v_closest = va;
                d         = dist;
            }
        }
        if( b_is_candidate ) {
            const FT dist = ( p - mesh.point( vb ) ).squared_length();
            if( dist < d ) {
                v_closest = vb;
                d         = dist;
            }
        }
        if( c_is_candidate ) {
            const FT dist = ( p - mesh.point( vc ) ).squared_length();
            if( dist < d ) {
                v_closest = vc;
                d         = dist;
            }
        }
        return { v_closest.idx(), 0 };

        // if( !vertex_candidates.empty() ) {
        //     Vertex_descriptor v_closest;
        //     FT d = FT_MAX;
        //     for( const auto& v: vertex_candidates ) {
        //         const Point& pv = mesh.point( v );
        //         const FT dist   = ( p - pv ).squared_length();
        //         if( dist < d ) {
        //             v_closest = v;
        //             d         = dist;
        //         }
        //     }
        //
        //     return { v_closest.idx(), 0 };
        // }
        //
        // return { triangle.idx(), 2 };
    }

    /**
     * Call `point_to_point` for an array of query points.
     *
     * @param p_array
     * @param p
     * @return
     */
    template<int N>
    inline std::array<FT, N> points_to_point( const std::array<Point, N>& p_array, const Point& p ) {
        std::array<FT, N> d;
        std::transform( p_array.begin(), p_array.end(), d.begin(), [&p]( const auto& e ) { return point_to_point( p, e ); } );
        return d;
    }

    /** Squared distance between a point p and an edge (a,b) */
    inline FT point_to_edge_squared( const Point& p, const Point& a, const Point& b ) {
        // cases when point is not within the cylinder of the edge
        if( CGAL::scalar_product( p - a, b - a ) < 0 ) {
            return point_to_point_squared( p, a );
        }
        if( CGAL::scalar_product( p - b, a - b ) < 0 ) {
            return point_to_point_squared( p, b );
        }

        // point - line distance
        return CGAL::cross_product( b - a, a - p ).squared_length() / ( b - a ).squared_length();
    }

    /** Distance between a point p and an edge (a,b) */
    inline FT point_to_edge( const Point& p, const Point& a, const Point& b ) {
        auto dist_squared = point_to_edge_squared( p, a, b );
        return CGAL::approximate_sqrt( dist_squared );
    }

    /**
     * Squared distance between a point p and an edge.
     *
     * @param p point
     * @param h halfedge
     * @param mesh mesh the halfedge belongs to
     * @return {index of primitive, type of primitive}. The types are: 0=vertex, 1=edge
     */
    inline std::tuple<FT, size_t, int> point_to_edge_squared( const Point& p, const Halfedge_descriptor& h, const SurfaceMesh& mesh ) {
        FT dist_squared;
        auto [idx, type] = closest_edge_primitive( p, h, mesh );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            dist_squared = point_to_point_squared( p, mesh.point( v ) );
            break;
        }
        case 1:    // edge
        {
            const Vertex_descriptor& va = mesh.target( h );
            const Vertex_descriptor& vb = mesh.source( h );
            dist_squared                = point_to_line_squared( p, mesh.point( va ), mesh.point( vb ) );
            break;
        }
        default:
            return { std::numeric_limits<FT>::max(), std::numeric_limits<size_t>::max(), -1 };
        }

        return { dist_squared, idx, type };
    }

    /**
     * Distance between a point p and an edge.
     *
     * @param p point
     * @param h halfedge
     * @param mesh mesh the halfedge belongs to
     * @return {index of primitive, type of primitive}. The types are: 0=vertex, 1=edge
     */
    inline std::tuple<FT, size_t, int> point_to_edge( const Point& p, const Halfedge_descriptor& h, const SurfaceMesh& mesh ) {
        FT dist;
        auto [idx, type] = closest_edge_primitive( p, h, mesh );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            dist = point_to_point( p, mesh.point( v ) );
            break;
        }
        case 1:    // edge
        {
            const Vertex_descriptor& va = mesh.target( h );
            const Vertex_descriptor& vb = mesh.source( h );
            dist                        = point_to_line( p, mesh.point( va ), mesh.point( vb ) );
            break;
        }
        default:
            return { std::numeric_limits<FT>::max(), std::numeric_limits<size_t>::max(), -1 };
        }

        return { dist, idx, type };
    }

    inline Point project_to_edge( const Point& p, const Halfedge_descriptor& h, const SurfaceMesh& mesh ) {
        auto [idx, type] = closest_edge_primitive( p, h, mesh );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            return mesh.point( v );
        }
        case 1:    // edge
        {
            const Vertex_descriptor& va = mesh.target( h );
            const Vertex_descriptor& vb = mesh.source( h );
            return project_to_line( p, mesh.point( va ), mesh.point( vb ) );
        }
        default:
            return Point();
        }
    }

    /**
     * Call `point_to_edge` for an array of query points.
     *
     * @param p_array
     * @param h
     * @param mesh
     * @param is_inside_slab
     * @return
     */
    template<int N>
    inline std::array<FT, N> points_to_edge( const std::array<Point, N>& p_array, const Halfedge_descriptor& h, const SurfaceMesh& mesh,
                                             bool& is_inside_slab ) {
        std::array<FT, N> d;
        is_inside_slab = false;
        std::transform( p_array.begin(), p_array.end(), d.begin(), [&h, &mesh, &is_inside_slab]( const auto& e ) {
            auto [d, i, t] = point_to_edge( e, h, mesh );
            if( t == 1 ) {
                is_inside_slab = true;
            }
            return d;
        } );
        return d;
    }

    /**
     * Squared distance between a point p and a triangle (a,b,c).
     *
     * @param p
     * @param a
     * @param b
     * @param c
     * @return
     */
    inline FT point_to_triangle_squared( const Point& p, const Point& a, const Point& b, const Point& c ) {
        Plane_3 plane( a, b, c );
        Point p_proj = plane.projection( p );

        // barycentric coordinates
        const Vector area_vec  = CGAL::cross_product( b - a, c - a );
        const Vector alpha_vec = CGAL::cross_product( b - p, c - p );
        const Vector beta_vec  = CGAL::cross_product( p - a, c - a );
        const Vector gamma_vec = CGAL::cross_product( b - a, p - a );

        // if area_vec and ..._vec show in the same direction (dot-product > 0), the corresponding barycentric coordinate is positive
        const bool is_alpha = CGAL::scalar_product( alpha_vec, area_vec ) > 0;
        const bool is_beta  = CGAL::scalar_product( beta_vec, area_vec ) > 0;
        const bool is_gamma = CGAL::scalar_product( gamma_vec, area_vec ) > 0;

        if( is_alpha && !is_beta && !is_gamma ) {
            return point_to_point_squared( p, a );
        }
        if( !is_alpha && is_beta && !is_gamma ) {
            return point_to_point_squared( p, b );
        }
        if( !is_alpha && !is_beta && is_gamma ) {
            return point_to_point_squared( p, c );
        }

        if( !is_alpha && is_beta && is_gamma ) {
            return point_to_edge_squared( p, b, c );
        }
        if( is_alpha && !is_beta && is_gamma ) {
            return point_to_edge_squared( p, a, c );
        }
        if( is_alpha && is_beta && !is_gamma ) {
            return point_to_edge_squared( p, a, b );
        }

        LOG_ASSERT( is_alpha && is_beta && is_gamma );

        // point is inside the prism of the triangle
        return ( p_proj - p ).squared_length();
    }

    /**
     * Distance between a point p and a triangle (a,b,c).
     *
     * @param p
     * @param a
     * @param b
     * @param c
     * @return
     */
    inline FT point_to_triangle( const Point& p, const Point& a, const Point& b, const Point& c ) {
        Plane_3 plane( a, b, c );
        Point p_proj = plane.projection( p );

        // barycentric coordinates
        const Vector area_vec  = CGAL::cross_product( b - a, c - a );
        const Vector alpha_vec = CGAL::cross_product( b - p, c - p );
        const Vector beta_vec  = CGAL::cross_product( p - a, c - a );
        const Vector gamma_vec = CGAL::cross_product( b - a, p - a );

        // if area_vec and ..._vec show in the same direction (dot-product > 0), the corresponding barycentric coordinate is positive
        const bool is_alpha = CGAL::scalar_product( alpha_vec, area_vec ) > 0;
        const bool is_beta  = CGAL::scalar_product( beta_vec, area_vec ) > 0;
        const bool is_gamma = CGAL::scalar_product( gamma_vec, area_vec ) > 0;

        if( is_alpha && !is_beta && !is_gamma ) {
            return point_to_point( p, a );
        }
        if( !is_alpha && is_beta && !is_gamma ) {
            return point_to_point( p, b );
        }
        if( !is_alpha && !is_beta && is_gamma ) {
            return point_to_point( p, c );
        }

        if( !is_alpha && is_beta && is_gamma ) {
            return point_to_edge( p, b, c );
        }
        if( is_alpha && !is_beta && is_gamma ) {
            return point_to_edge( p, a, c );
        }
        if( is_alpha && is_beta && !is_gamma ) {
            return point_to_edge( p, a, b );
        }

        LOG_ASSERT( is_alpha && is_beta && is_gamma );

        // point is inside the prism of the triangle
        return CGAL::approximate_sqrt( ( p_proj - p ).squared_length() );
    }

    /**
     * Squared distance between a point p and a triangle of a surface mesh.
     *
     * @param p
     * @param triangle
     * @param mesh
     * @return {distance, index of primitive, type of primitive}. The types are: 0=vertex, 1=edge, 2=face
     */
    inline std::tuple<FT, size_t, int> point_to_triangle_squared( const Point& p, const Face_descriptor& triangle, const SurfaceMesh& mesh ) {
        FT dist_squared;
        auto [idx, type] = closest_triangle_primitive( p, triangle, mesh );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            dist_squared = point_to_point_squared( p, mesh.point( v ) );
            break;
        }
        case 1:    // edge
        {
            const Edge_descriptor e( idx );
            const Halfedge_descriptor h = mesh.halfedge( e, 0 );
            const Vertex_descriptor& va = mesh.target( h );
            const Vertex_descriptor& vb = mesh.source( h );
            dist_squared                = point_to_line_squared( p, mesh.point( va ), mesh.point( vb ) );
            break;
        }
        case 2:    // face
        {
            const Halfedge_descriptor& h = mesh.halfedge( triangle );
            const Vertex_descriptor& va  = mesh.source( h );
            const Vertex_descriptor& vb  = mesh.target( h );
            const Vertex_descriptor& vc  = mesh.target( mesh.next( h ) );
            dist_squared                 = point_to_plane_squared( p, mesh.point( va ), mesh.point( vb ), mesh.point( vc ) );
            break;
        }
        default:
            return { std::numeric_limits<FT>::max(), std::numeric_limits<size_t>::max(), -1 };
        }

        return { dist_squared, idx, type };
    }

    /**
     * Distance between a point p and a triangle of a surface mesh.
     *
     * @param p
     * @param triangle
     * @param mesh
     * @return {distance, index of primitive, type of primitive}. The types are: 0=vertex, 1=edge, 2=face
     */
    inline std::tuple<FT, size_t, int> point_to_triangle( const Point& p, const Face_descriptor& triangle, const SurfaceMesh& mesh ) {
        FT dist;
        auto [idx, type] = closest_triangle_primitive( p, triangle, mesh );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            dist = point_to_point( p, mesh.point( v ) );
            break;
        }
        case 1:    // edge
        {
            const Edge_descriptor e( idx );
            const Halfedge_descriptor h = mesh.halfedge( e, 0 );
            const Vertex_descriptor& va = mesh.target( h );
            const Vertex_descriptor& vb = mesh.source( h );
            dist                        = point_to_line( p, mesh.point( va ), mesh.point( vb ) );
            break;
        }
        case 2:    // face
        {
            const Halfedge_descriptor& h = mesh.halfedge( triangle );
            const Vertex_descriptor& va  = mesh.source( h );
            const Vertex_descriptor& vb  = mesh.target( h );
            const Vertex_descriptor& vc  = mesh.target( mesh.next( h ) );
            dist                         = point_to_plane( p, mesh.point( va ), mesh.point( vb ), mesh.point( vc ) );
            break;
        }
        default:
            return { std::numeric_limits<FT>::max(), std::numeric_limits<size_t>::max(), -1 };
        }

        return { dist, idx, type };
    }

    inline Point project_to_triangle( const Point& p, const Face_descriptor& triangle, const SurfaceMesh& mesh ) {
        auto [idx, type] = closest_triangle_primitive( p, triangle, mesh );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            return mesh.point( v );
        }
        case 1:    // edge
        {
            const Edge_descriptor e( idx );
            const Halfedge_descriptor h = mesh.halfedge( e, 0 );
            const Vertex_descriptor& va = mesh.target( h );
            const Vertex_descriptor& vb = mesh.source( h );
            return project_to_line( p, mesh.point( va ), mesh.point( vb ) );
        }
        case 2:    // face
        {
            const Halfedge_descriptor& h = mesh.halfedge( triangle );
            const Vertex_descriptor& va  = mesh.source( h );
            const Vertex_descriptor& vb  = mesh.target( h );
            const Vertex_descriptor& vc  = mesh.target( mesh.next( h ) );
            return project_to_plane( p, mesh.point( va ), mesh.point( vb ), mesh.point( vc ) );
        }
        default:
            return Point();
        }
    }

    /**
     * Call `point_to_triangle` for an array of query points.
     *
     * @param p_array
     * @param triangle
     * @param mesh
     * @param is_inside_slab
     * @return
     */
    template<int N>
    inline std::array<FT, N> points_to_triangle( const std::array<Point, N>& p_array, const Face_descriptor& triangle, const SurfaceMesh& mesh,
                                                 bool& is_inside_slab ) {
        std::array<FT, N> d;
        std::transform( p_array.begin(), p_array.end(), d.begin(), [&triangle, &mesh, &is_inside_slab]( const auto& e ) {
            auto [d, i, t] = point_to_triangle( e, triangle, mesh );
            if( t == 2 ) {
                is_inside_slab = true;
            }
            return d;
        } );
        return d;
    }

    /**
     * Call `point_to_triangle_squared` for an array of query points.
     *
     * @param p_array
     * @param triangle
     * @param mesh
     * @param is_inside_slab
     * @return
     */
    template<int N>
    inline std::array<FT, N> points_to_triangle_squared( const std::array<Point, N>& p_array, const Face_descriptor& triangle,
                                                         const SurfaceMesh& mesh, bool& is_inside_slab ) {
        std::array<FT, N> d;
        std::transform( p_array.begin(), p_array.end(), d.begin(), [&triangle, &mesh, &is_inside_slab]( const auto& e ) {
            auto [d, i, t] = point_to_triangle_squared( e, triangle, mesh );
            if( t == 2 ) {
                is_inside_slab = true;
            }
            return d;
        } );
        return d;
    }

    inline Point project_to_primitive( const Point& p, const Face_descriptor& f, const SurfaceMesh& mesh ) {
        const auto [idx, type] = closest_triangle_primitive( p, f, mesh );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            return mesh.point( v );
        }
        case 1:    // edge
        {
            Edge_descriptor e( idx );
            return project_to_edge( p, mesh.halfedge( e, 0 ), mesh );
        }
        case 2:    // face
        {
            return project_to_triangle( p, f, mesh );
        }
        default:
            return Point();
        }
    }

    /**
     * Compute the normal vector pointing from vertex v to point p.
     *
     * @param v
     * @param m
     * @param p
     * @return
     */
    inline Vector get_vertex_normal( const Vertex_descriptor& v, const SurfaceMesh& m, const Point& p ) {
        const Point v_point = m.point( v );
        return normalize( p - v_point );
    }

    /**
     * Compute the normal vector that is orthogonal to edge e and points towards point p.
     *
     * @param e
     * @param m
     * @param p
     * @return
     */
    inline Vector get_edge_normal( const Edge_descriptor& e, const SurfaceMesh& m, const Point& p ) {
        const Halfedge_descriptor h = m.halfedge( e );
        const Vertex_descriptor v0  = m.source( h );
        const Vertex_descriptor v1  = m.target( h );
        const Point p0              = m.point( v0 );
        const Point p1              = m.point( v1 );

        if( CGAL::scalar_product( p - p0, p1 - p0 ) < 0 ) {
            return normalize( p - p0 );
        }
        if( CGAL::scalar_product( p - p1, p0 - p1 ) < 0 ) {
            return normalize( p - p1 );
        }

        Vector orth = CGAL::cross_product( p - p0, p1 - p0 );
        Vector n    = CGAL::cross_product( p1 - p0, orth );

        return normalize( n );
    }

    /**
     * Compute the normal vector from the projection point of p onto face f pointing to point p.
     *
     * @param f
     * @param m
     * @param p
     * @return
     */
    inline Vector get_face_normal( const Face_descriptor& f, const SurfaceMesh& m, const Point& p ) {
        auto [idx, type] = Distance::closest_triangle_primitive( p, f, m );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            return get_vertex_normal( v, m, p );
        }
        case 1:    // edge
        {
            const Edge_descriptor e( idx );
            return get_edge_normal( e, m, p );
        }
        case 2:    // face
        {
            return PMP::compute_face_normal( f, m );
        }
        default:
            LOG( ERROR ) << "Unknown primitive type";
            break;
        }
        return { 0, 0, 0 };
    }

    /**
     * Compute the normal vector pointing towards p from its projection point on the primitive with the given unique id.
     *
     * @param id
     * @param m
     * @param p
     * @return
     */
    inline Vector get_prim_normal( const size_t& id, const SurfaceMesh& m, const Point& p ) {
        const size_t type = identify_unique_id_type( id, m );
        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( id );
            return get_vertex_normal( v, m, p );
        }
        case 1:    // edge
        {
            const Edge_descriptor e( id - m.number_of_vertices() );
            return get_edge_normal( e, m, p );
        }
        case 2:    // face
        {
            const Face_descriptor f( id - m.number_of_vertices() - m.number_of_edges() );
            return get_face_normal( f, m, p );
        }
        default:
            LOG( ERROR ) << "Unknown primitive type";
            break;
        }
        return { 0, 0, 0 };
    }

    struct Sphere {
        Point origin;
        FT radius;

        Point project( const Point& p ) const {
            LOG_ASSERT( ( p - origin ).squared_length() != 0 );
            const Vector proj_vec = normalize( p - origin );
            return origin + radius * proj_vec;
        }

        std::array<FT, 4> projection_plane( const Point& p ) const {
            const Vector n     = normalize( p - origin );
            const Point p_proj = origin + radius * n;
            const FT d         = -CGAL::scalar_product( n, p_proj - CGAL::ORIGIN );
            return { n.x(), n.y(), n.z(), d };
        }

        FT signed_distance( const Point& x ) const { return CGAL::approximate_sqrt( ( x - origin ).squared_length() ) - radius; }
    };

    struct Cylinder {
        Point origin;
        Vector direction;
        FT radius;

        Point project( const Point& p ) const {
            const Vector po = p - origin;
            const Vector xo = direction * CGAL::scalar_product( po, direction );
            const Point x   = origin + xo;
            LOG_ASSERT( ( p - x ).squared_length() != 0 );
            const Vector xp = normalize( p - x );
            return x + radius * xp;
        }

        std::array<FT, 4> projection_plane( const Point& p ) const {
            const Vector po    = p - origin;
            const Vector xo    = direction * CGAL::scalar_product( po, direction );
            const Point x      = origin + xo;
            const Vector n     = normalize( p - x );
            const Point p_proj = x + radius * n;
            const FT d         = -CGAL::scalar_product( n, p_proj - CGAL::ORIGIN );
            return { n.x(), n.y(), n.z(), d };
        }

        FT signed_distance( const Point& x ) const {
            return CGAL::approximate_sqrt( CGAL::cross_product( direction, ( x - origin ) ).squared_length() ) - radius;
        }
    };

    struct Plane {
        Point origin;
        Vector normal;

        Point project( const Point& p ) const {
            const Vector po = p - origin;
            const Vector xp = -normal * CGAL::scalar_product( po, normal );
            return p + xp;
        }

        std::array<FT, 4> projection_plane( const Point& p ) const {
            const FT d = -CGAL::scalar_product( normal, origin - CGAL::ORIGIN );
            return { normal.x(), normal.y(), normal.z(), d };
        }

        FT signed_distance( const Point& x ) const { return CGAL::scalar_product( ( x - origin ), normal ); }
    };

    inline Sphere get_primitive_offset( const SurfaceMesh& m, const Vertex_descriptor& v, const FT& dist ) {
        const Point o = m.point( v );
        return { o, dist };
    }

    inline Cylinder get_primitive_offset( const SurfaceMesh& m, const Edge_descriptor& e, const FT& dist ) {
        const Halfedge_descriptor h = m.halfedge( e );
        const Point p0              = m.point( m.source( h ) );
        const Point p1              = m.point( m.target( h ) );
        Vector direction            = ( p1 - p0 );
        direction /= CGAL::approximate_sqrt( direction.squared_length() );
        return { p0, direction, dist };
    }

    inline Plane get_primitive_offset( const SurfaceMesh& m, const Face_descriptor& face, const FT& dist, const Point& p ) {
        const Halfedge_descriptor& ha = m.halfedge( face );
        const Halfedge_descriptor& hb = m.next( ha );
        const Halfedge_descriptor& hc = m.next( hb );
        const Vertex_descriptor& va   = m.target( hb );
        const Vertex_descriptor& vb   = m.target( hc );
        const Vertex_descriptor& vc   = m.target( ha );
        const Point& a                = m.point( va );
        const Point& b                = m.point( vb );
        const Point& c                = m.point( vc );

        const Vector area_vec = CGAL::cross_product( b - a, c - a );
        Vector normal         = area_vec / CGAL::approximate_sqrt( area_vec.squared_length() );

        if( CGAL::scalar_product( p - a, normal ) < 0 ) {
            normal = -normal;
        }

        const Point o = a + normal * dist;

        return { o, normal };
    }

}    // namespace Distance