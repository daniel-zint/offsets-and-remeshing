#pragma once
/**
 * @file   utilities.h
 * @brief  Utility functions.
 *
 * @author Daniel Zint
 * @date   December 2022
 */

#include "types.h"

#include <CGAL/Distance_3/Point_3_Triangle_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <filesystem>
#include <glog/logging.h>

/**
 * Check if an edge is strictly convex using `CGAL::orientation`.
 *
 * @param h halfedge
 * @param m mesh
 * @return true if orientation is `CGAL::NEGATIVE`
 */
inline bool is_strictly_convex( const Halfedge_descriptor& h, const SurfaceMesh& m ) {
    if( m.is_border( h ) )
        return false;

    // get three points of mate facet
    const Point& a = m.point( m.target( h ) );
    const Point& b = m.point( m.target( m.next( h ) ) );
    const Point& c = m.point( m.source( h ) );

    // opposite point on neighboring face
    const Point& d   = m.point( m.target( m.next( m.opposite( h ) ) ) );
    auto orientation = CGAL::orientation( a, b, c, d );
    if( orientation == CGAL::NEGATIVE ) {
        return true;
    } else {
        return false;
    }
}

/**
 * Normalize vector v.
 *
 * @param v vector to normalize
 * @return normalized vector v. If v is of length 0, the `CGAL::NULL_VECTOR` is returned.
 */
inline Vector normalize( const Vector& v ) {
    // TODO implement Quake's fast inverse square root
    const FT sqnorm = v.squared_length();
    if( sqnorm != 0.0 )
        return v / std::sqrt( sqnorm );
    else
        return CGAL::NULL_VECTOR;
}

/**
 * Compute the centroid of a bbox.
 *
 * @param bbox
 * @return
 */
inline Point centroid_of_bbox( const CGAL::Bbox_3& bbox ) {
    const Point p0( bbox.xmin(), bbox.ymin(), bbox.zmin() );
    const Point p1( bbox.xmax(), bbox.ymax(), bbox.zmax() );
    return CGAL::midpoint( p0, p1 );
}

/**
 * Returns true iff pt is on the negative side of the plane defined by (ep0, ep1) and normal.
 *
 * @param p
 * @param normal
 * @param ep0
 * @param ep1
 * @return
 */
inline bool on_left_of_triangle_edge( const Point& p, const Vector& normal, const Point& ep0, const Point& ep1 ) {
    return CGAL::internal::on_left_of_triangle_edge( p, normal, ep0, ep1, Kernel() );
}

/**
 * Check if the point is within the slab defined by the edge (ep0,ep1) and the normal. The slab is bounded by three planes:
 * 1. The plane that contains the edge and the normal as tangent vectors.
 * 2. The plane through ep0 with the normal vector (ep1-ep0)
 * 3. The plane through ep1 with the normal vector (ep0-ep1)
 * @param pt
 * @param normal
 * @param ep0
 * @param ep1
 * @return -1: left side of triangle, 0: edge, 1: ep0, 2: ep1
 */
inline int on_right_slab_of_triangle_edge( const Point& pt, const Vector& normal, const Point& ep0, const Point& ep1 ) {
    // In case the normal is 0, the result of on_left_of_triangle_edge is garbage.
    // A 0 normal appears at degenerated triangles. In that case, the nearest point is never on the triangle but either on the edge or a vertex.
    if( normal.squared_length() > 0 && CGAL::internal::on_left_of_triangle_edge( pt, normal, ep0, ep1, Kernel() ) ) {
        return -1;
    }
    if( ( ep1 - ep0 ).squared_length() == 0 ) {
        // edge is degenerated, just return any point
        return 1;
    }
    if( CGAL::scalar_product( ep1 - ep0, pt - ep0 ) < 0 ) {
        return 1;
    }
    if( CGAL::scalar_product( ep0 - ep1, pt - ep1 ) < 0 ) {
        return 2;
    }

    return 0;
}

/**
 * Print the id of a primitive and the ids of its incident vertices. This method is purely for debugging.
 *
 * @param id unique id
 * @param m mesh
 */
inline void identify_unique_id( const size_t& id, const SurfaceMesh& m ) {
    std::cout << "Unique id " << id << std::endl;

    if( id < m.number_of_vertices() ) {
        // primitive is a vertex
        std::cout << "  Vertex" << std::endl;
    } else if( id < m.number_of_vertices() + m.number_of_edges() ) {
        // primitive is an edge
        Edge_descriptor e( id - m.number_of_vertices() );
        Halfedge_descriptor h = m.halfedge( e );
        Vertex_descriptor v0  = m.source( h );
        Vertex_descriptor v1  = m.target( h );
        std::cout << "  Edge (v0,v1) = (" << v0 << "," << v1 << ")" << std::endl;
    } else {
        // primitive is a face
        Face_descriptor f( id - m.number_of_vertices() - m.number_of_edges() );
        Halfedge_descriptor h = m.halfedge( f );
        Vertex_descriptor v0  = m.source( h );
        Vertex_descriptor v1  = m.target( h );
        Vertex_descriptor v2  = m.target( m.next( h ) );
        std::cout << "  Face (v0,v1,v2) = (" << v0 << "," << v1 << "," << v2 << ")" << std::endl;
        std::cout << "  face id = " << f.idx() << std::endl;
    }
}

/**
 * __DEPRECATED__ Compute the primitive type (vertex, edge, face) of an element.
 * 0: vertex
 * 1: edge
 * 2: face
 * @param id unique id
 * @param m mesh
 * @return 0: vertex, 1: edge, 2: face, MAX: invalid unique id
 */
inline size_t identify_unique_id_type( const size_t& id, const SurfaceMesh& m ) {
    LOG( FATAL ) << "This function is deprecated!";

    if( id < m.number_of_vertices() ) {
        // primitive is a vertex
        return 0;
    } else if( id < m.number_of_vertices() + m.number_of_edges() ) {
        // primitive is an edge
        return 1;
    } else if( id < m.number_of_vertices() + m.number_of_edges() + m.number_of_faces() ) {
        // primitive is a face
        return 2;
    } else {
        // primitive does not exist
        return std::numeric_limits<size_t>::max();
    }
}

/**
 * __DEPRECATED__ Check if two primitives are incident.
 *
 * @param id1 unique id of the first primitive
 * @param id2 unique id of the second primitive
 * @param m mesh
 * @return true, if the two primitives are incident
 */
inline bool are_ids_incident( const size_t& id1, const size_t& id2, const SurfaceMesh& m ) {
    LOG( FATAL ) << "This function is deprecated!";
    auto [v_unique_id, v_unique_id_found] = m.property_map<Vertex_descriptor, size_t>( "v:unique_id" );
    LOG_ASSERT( v_unique_id_found );
    auto [e_unique_id, e_unique_id_found] = m.property_map<Edge_descriptor, size_t>( "e:unique_id" );
    LOG_ASSERT( e_unique_id_found );
    auto [f_unique_id, f_unique_id_found] = m.property_map<Face_descriptor, size_t>( "f:unique_id" );
    LOG_ASSERT( f_unique_id_found );

    const size_t id1_type = identify_unique_id_type( id1, m );
    switch( id1_type ) {
    case 0: {
        const Vertex_descriptor v( id1 );
        // check incident edges and faces
        for( const auto& h: m.halfedges_around_target( m.halfedge( v ) ) ) {
            const auto& e    = m.edge( h );
            const auto& e_id = e_unique_id[e];
            if( e_id == id2 ) {
                return true;
            }

            if( !m.is_border( e ) ) {
                const auto& f    = m.face( h );
                const auto& f_id = f_unique_id[f];
                if( f_id == id2 ) {
                    return true;
                }
            }
        }
        break;
    }
    case 1: {
        const Edge_descriptor e( id1 - m.number_of_vertices() );
        const auto& h0 = m.halfedge( e, 0 );
        const auto& h1 = m.halfedge( e, 1 );

        // check incident vertices
        const auto& v0 = m.source( h0 );
        const auto& v1 = m.target( h0 );
        if( v_unique_id[v0] == id2 ) {
            return true;
        }
        if( v_unique_id[v1] == id2 ) {
            return true;
        }

        // check incident faces
        if( !m.is_border( h0 ) && f_unique_id[m.face( h0 )] == id2 ) {
            return true;
        }
        if( !m.is_border( h1 ) && f_unique_id[m.face( h1 )] == id2 ) {
            return true;
        }
        break;
    }
    case 2: {
        const Face_descriptor f( id1 - m.number_of_vertices() - m.number_of_edges() );

        // check incident edges and vertices
        for( const auto& h: m.halfedges_around_face( m.halfedge( f ) ) ) {
            const auto& e = m.edge( h );
            if( e_unique_id[e] == id2 ) {
                return true;
            }
            const auto& v = m.source( h );
            if( v_unique_id[v] == id2 ) {
                return true;
            }
        }

        break;
    }
    default: {
        LOG( FATAL ) << "ID type unknown";
        break;
    }
    }

    return false;
}

/**
 * __DEPRECATED__ Get all neighbors (incident and adjacent) of prim that are not complex or behind a complex prim.
 *
 * @param prim
 * @param m
 * @return
 */
inline std::set<size_t> get_non_complex_neighbors( const size_t& prim, const SurfaceMesh& m ) {
    LOG( FATAL ) << "This function is deprecated!";
    auto [v_is_complex, v_is_complex_found] = m.property_map<Vertex_descriptor, bool>( "v:is_complex" );
    LOG_ASSERT( v_is_complex_found );
    auto [e_is_complex, e_is_complex_found] = m.property_map<Edge_descriptor, bool>( "e:is_complex" );
    LOG_ASSERT( e_is_complex_found );

    auto [v_unique_id, v_unique_id_found] = m.property_map<Vertex_descriptor, size_t>( "v:unique_id" );
    LOG_ASSERT( v_unique_id_found );
    auto [e_unique_id, e_unique_id_found] = m.property_map<Edge_descriptor, size_t>( "e:unique_id" );
    LOG_ASSERT( e_unique_id_found );
    auto [f_unique_id, f_unique_id_found] = m.property_map<Face_descriptor, size_t>( "f:unique_id" );
    LOG_ASSERT( f_unique_id_found );

    std::set<size_t> prim_neighbors;
    const size_t prim_type = identify_unique_id_type( prim, m );
    switch( prim_type ) {
    case 0: {
        const Vertex_descriptor v( prim );
        if( v_is_complex[v] ) {
            break;
        }

        for( const auto& h: m.halfedges_around_target( m.halfedge( v ) ) ) {
            prim_neighbors.insert( e_unique_id[m.edge( h )] );
            prim_neighbors.insert( f_unique_id[m.face( h )] );
        }
        break;
    }
    case 1: {
        const Edge_descriptor e( prim - m.number_of_vertices() );
        if( e_is_complex[e] ) {
            break;
        }

        const auto& h0 = m.halfedge( e, 0 );
        const auto& h1 = m.halfedge( e, 1 );
        const auto& v0 = m.source( h0 );
        const auto& v1 = m.target( h0 );
        if( !v_is_complex[v0] ) {
            prim_neighbors.insert( v_unique_id[v0] );
        }
        if( !v_is_complex[v1] ) {
            prim_neighbors.insert( v_unique_id[v1] );
        }
        prim_neighbors.insert( f_unique_id[m.face( h0 )] );
        prim_neighbors.insert( f_unique_id[m.face( h1 )] );
        break;
    }
    case 2: {
        const Face_descriptor f( prim - m.number_of_vertices() - m.number_of_edges() );
        for( const auto& h: m.halfedges_around_face( m.halfedge( f ) ) ) {
            const auto& e = m.edge( h );
            if( e_is_complex[e] ) {
                continue;
            }
            prim_neighbors.insert( e_unique_id[e] );
            prim_neighbors.insert( f_unique_id[m.face( m.opposite( h ) )] );
            const auto& v = m.target( h );
            if( !v_is_complex[v] ) {
                prim_neighbors.insert( v_unique_id[v] );
            }
        }

        break;
    }
    default:
        LOG( FATAL ) << "ID type unknown";
        break;
    }

    return prim_neighbors;
}

/**
 * Get all incident faces of prim.
 *
 * @param prim
 * @param m
 * @return
 */
inline std::set<Face_descriptor> get_incident_faces( const size_t& prim, const SurfaceMesh& m ) {
    std::set<Face_descriptor> prim_neighbors;
    const size_t prim_type = identify_unique_id_type( prim, m );
    switch( prim_type ) {
    case 0: {
        const Vertex_descriptor v( prim );
        for( const auto& f: m.faces_around_face( m.halfedge( v ) ) ) {
            prim_neighbors.insert( f );
        }
        break;
    }
    case 1: {
        const Edge_descriptor e( prim - m.number_of_vertices() );
        prim_neighbors.insert( m.face( m.halfedge( e, 0 ) ) );
        prim_neighbors.insert( m.face( m.halfedge( e, 1 ) ) );
        break;
    }
    case 2: {
        const Face_descriptor f( prim - m.number_of_vertices() - m.number_of_edges() );
        prim_neighbors.insert( f );
        break;
    }
    default:
        LOG( FATAL ) << "ID type unknown";
        break;
    }

    return prim_neighbors;
}

/**
 * Check if the primitive is marked as complex. Only vertices or edges can be complex. For faces, the result is always false.
 *
 * @param prim
 * @param m
 * @return
 */
inline bool is_prim_complex( const size_t& prim, const SurfaceMesh& m ) {
    auto [v_is_complex, v_is_complex_found] = m.property_map<Vertex_descriptor, bool>( "v:is_complex" );
    LOG_ASSERT( v_is_complex_found );
    auto [e_is_complex, e_is_complex_found] = m.property_map<Edge_descriptor, bool>( "e:is_complex" );
    LOG_ASSERT( e_is_complex_found );

    const size_t prim_type = identify_unique_id_type( prim, m );
    switch( prim_type ) {
    case 0: {
        const Vertex_descriptor v( prim );
        return v_is_complex[v];
    }
    case 1: {
        const Edge_descriptor e( prim - m.number_of_vertices() );
        return e_is_complex[e];
    }
    case 2: {
        const Face_descriptor f( prim - m.number_of_vertices() - m.number_of_edges() );
        return false;
    }
    default:
        LOG( FATAL ) << "ID type unknown";
        return true;
    }
}

/**
 * Check if the entries in vector v1 are all contained in vector v2.
 *
 * @param v1
 * @param v2
 * @return
 */
inline bool is_v1_subset_of_v2( const std::vector<size_t>& v1, const std::vector<size_t>& v2 ) {
    if( v1.size() > v2.size() ) {
        return false;
    }
    std::vector<size_t> intersection;
    std::set_intersection( v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter( intersection ) );
    if( intersection.size() == v1.size() ) {
        return true;
    } else {
        return false;
    }
}

/**
 * Find the element with the most occurences in an unordered vector and return this element and the number of occurences.
 *
 * @param v
 * @return {most frequent element, number of occurences}
 */
inline std::tuple<std::size_t, std::size_t> maxFreq( const std::vector<Face_descriptor>& v ) {
    std::unordered_map<std::size_t, std::size_t> hash;
    for( int i = 0; i < v.size(); i++ ) {
        hash[v[i].idx()]++;
    }

    // find the max frequency
    int max_count = 0, res = -1;
    for( const auto& i: hash ) {
        if( max_count < i.second ) {
            res       = i.first;
            max_count = i.second;
        }
    }

    return std::make_tuple( hash[res], max_count );
}

/**
 * splits the halfedge `h` and afterwards:
 * - calls `split_face()` with `prev(h, g)` and `next(h, g)`, if `h` is not a border halfedge
 * - calls `split_face()` with `opposite(h, g)` and `next(next(h, g), g)`, if `oppiste(h, g)`
 *   is not a border halfedge
 *
 * @tparam Graph must be a `MutableFaceGraph`
 *
 * @returns the halfedge `hnew` pointing to the inserted vertex in the edge split.
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor split_edge_and_incident_faces( typename boost::graph_traits<Graph>::halfedge_descriptor h,
                                                                                        Graph& g ) {
    typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

    halfedge_descriptor res = CGAL::Euler::split_edge( h, g );

    if( !is_border( res, g ) )
        CGAL::Euler::split_face( res, next( h, g ), g );

    halfedge_descriptor opp_h = opposite( h, g );
    if( !is_border( opp_h, g ) )
        CGAL::Euler::split_face( opp_h, next( next( opp_h, g ), g ), g );

    return res;
}

/**
 * Mean ration metric of a triangle projected onto the plane that is defined through the normal.
 *
 * @param triangle
 * @param normal
 * @return
 */
inline float meanRatioMetric( const std::array<Point, 3> triangle, const Vector& normal ) {
    const Vector e12 = triangle[0] - triangle[1];
    const Vector e23 = triangle[1] - triangle[2];
    const Vector e13 = triangle[0] - triangle[2];

    const FT l            = e12.squared_length() + e23.squared_length() + e13.squared_length();
    const Vector area_vec = CGAL::cross_product( e12, e23 );
    FT area               = CGAL::scalar_product( area_vec, normal );

    return 2.f * std::sqrt( 3.f ) * area / l;
}

/**
 * Mean ratio metric multiplied by the normalized deviation of the triangle normal and the given normal.
 * The normal deviation is made heavy by beeing multiplied three times.
 *
 * @param triangle array of three points describing the triangle
 * @param normal the *correct* normal
 * @return
 */
inline float meanRatioMetric_heavy_normal( const std::array<Point, 3> triangle, const Vector& normal ) {
    const Vector e12 = triangle[0] - triangle[1];
    const Vector e23 = triangle[1] - triangle[2];
    const Vector e13 = triangle[0] - triangle[2];

    const FT l            = e12.squared_length() + e23.squared_length() + e13.squared_length();
    const Vector area_vec = CGAL::cross_product( e12, e23 );
    const FT area         = CGAL::sqrt( area_vec.squared_length() );
    const FT angle        = CGAL::approximate_angle( area_vec, normal );

    const FT angle_normalized = ( ( 90.0 - angle ) / 90.0 );

    return ( 2.f * std::sqrt( 3.f ) * area / l ) * angle_normalized * angle_normalized * angle_normalized;
}

inline void print_mesh( const SurfaceMesh& mesh, const std::filesystem::path& filename ) {
    CGAL::IO::write_OFF( filename.string(), mesh, CGAL::parameters::stream_precision( 15 ) );
}

inline void normalize_mesh( SurfaceMesh& m ) {
    CGAL::Bbox_3 aabb = PMP::bbox( m );

    const FT d = CGAL::max( aabb.x_span(), CGAL::max( aabb.y_span(), aabb.z_span() ) );

    for( const auto& v: m.vertices() ) {
        Point& p   = m.point( v );
        const FT x = ( p.x() - aabb.xmin() ) / d;
        const FT y = ( p.y() - aabb.ymin() ) / d;
        const FT z = ( p.z() - aabb.zmin() ) / d;
        p          = Point( x, y, z );
    }
}