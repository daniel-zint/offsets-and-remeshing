#ifndef CGAL_DUAL_CONTOURING_3_H
#define CGAL_DUAL_CONTOURING_3_H

#include "Distance.h"
#include "Octree_wrapper.h"
#include "types.h"

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/centroid.h>
#include <Eigen/SVD>
#include <array>
#include <map>
#include <vector>

namespace CGAL {

    namespace Dual_contouring_3 {
        namespace Positioning {
            template<bool use_bbox = false>
            class QEM_SVD_hermite {
                std::function<bool( const Point& )> is_inside;
                std::unique_ptr<Sotm> sotm;

              public:
                QEM_SVD_hermite( const SurfaceMesh& primal_mesh ) {
                    if( CGAL::is_closed( primal_mesh ) ) {
                        sotm      = std::make_unique<Sotm>( primal_mesh );
                        is_inside = [this]( const Point& p ) { return ( *( this->sotm ) )( p ) != CGAL::ON_UNBOUNDED_SIDE; };
                    } else {
                        is_inside = []( const Point& p ) { return false; };
                    }
                }

                /// <summary>
                /// Compute vertex position for Dual Contouring
                /// </summary>
                /// <returns> true, if there is a point in the cell</returns>
                int position( const Octree_wrapper& octree, SurfaceMesh& offset_mesh, const SurfaceMesh& primal_mesh, const FT iso_value,
                              const Vertex_descriptor& v ) const {
                    auto [v_primitives, v_primitives_found] =
                        offset_mesh.property_map<Vertex_descriptor, std::vector<Face_descriptor>>( "v:primitives" );
                    LOG_ASSERT( v_primitives_found );
                    auto [v_voxel_handle, v_voxel_handle_found] = offset_mesh.property_map<Vertex_descriptor, size_t>( "v:voxel_handle" );
                    LOG_ASSERT( v_voxel_handle );
                    auto [v_bbox, v_bbox_found] = offset_mesh.property_map<Vertex_descriptor, CGAL::Bbox_3>( "v:bbox" );
                    LOG_ASSERT( v_bbox_found );

                    const FT max_step_size = std::min( v_bbox[v].x_span(), std::min( v_bbox[v].y_span(), v_bbox[v].z_span() ) );

                    const auto& prims = v_primitives[v];

                    std::array<FT, Tables::N_VERTICES> s = octree.voxel_values( v_voxel_handle[v] );

                    std::array<bool, Tables::N_VERTICES> b;
                    std::transform( s.begin(), s.end(), b.begin(), [iso_value]( const auto& e ) { return e <= iso_value; } );

                    unsigned int cubeindex = 0;
                    // set bit if corresponding corner is below iso
                    for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                        cubeindex |= b[i] << i;
                    }

                    if( cubeindex == 0 || cubeindex == 255 ) {
                        return 0;
                    }

                    std::array<Point, Tables::N_VERTICES> p = octree.voxel_vertex_positions( v_voxel_handle[v] );
                    std::array<Vector, Tables::N_VERTICES> pos;
                    std::transform( p.begin(), p.end(), pos.begin(), []( const auto& e ) { return e - CGAL::ORIGIN; } );
                    std::array<FT, Tables::N_VERTICES> sign;
                    std::transform( p.begin(), p.end(), sign.begin(), [this]( const Point& e ) { return is_inside( e ) ? -1 : 1; } );

                    // compute edge intersections
                    std::vector<Point> edge_intersections;
                    std::vector<Vector> edge_intersection_normals;

                    for( int i = 0; i < Tables::N_EDGES; ++i ) {
                        const auto& v0 = Tables::edge_to_vertex[i][0];
                        const auto& v1 = Tables::edge_to_vertex[i][1];

                        if( b[v0] != b[v1] ) {    // e0
                            const FT u   = ( s[v0] - iso_value ) / ( s[v0] - s[v1] );
                            Point p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[v0] + u * pos[v1] );

                            // only use prims that intersect the edge for interpolation
                            std::vector<Face_descriptor> intersecting_prims;
                            for( const auto& prim: prims ) {
                                const Face_descriptor f( prim );
                                // check if prim intersect this edge
                                // const auto [dist0, id0, type0] = Distance::point_to_triangle( p[v0], f, primal_mesh );
                                // const auto [dist1, id1, type1] = Distance::point_to_triangle( p[v1], f, primal_mesh );
                                const Point proj0 = Distance::project_to_triangle( p[v0], f, primal_mesh );
                                const Point proj1 = Distance::project_to_triangle( p[v1], f, primal_mesh );
                                const FT dist0    = CGAL::sqrt( ( p[v0] - proj0 ).squared_length() );
                                const FT dist1    = CGAL::sqrt( ( p[v1] - proj1 ).squared_length() );
                                const bool b0     = ( sign[v0] * dist0 ) <= iso_value;
                                const bool b1     = ( sign[v1] * dist1 ) <= iso_value;
                                if( b0 != b1 ) {
                                    intersecting_prims.push_back( f );
                                    continue;
                                }
                                // TODO project to line
                                // check if projection onto line is inside the radius
                                const Point e_proj0 = Distance::project_to_line( proj0, p[v0], p[v1] );
                                if( ( proj0 - e_proj0 ).squared_length() < iso_value * iso_value ) {
                                    intersecting_prims.push_back( f );
                                    continue;
                                }
                                const Point e_proj1 = Distance::project_to_line( proj1, p[v0], p[v1] );
                                if( ( proj1 - e_proj1 ).squared_length() < iso_value * iso_value ) {
                                    intersecting_prims.push_back( f );
                                    continue;
                                }
                            }

                            if( intersecting_prims.empty() ) {
                                continue;
                            }

                            auto get_distance_to_intersecting_prims = [&intersecting_prims, &primal_mesh]( const Point& p ) {
                                FT dist_squared_min = FT_MAX;
                                for( const auto& prim: intersecting_prims ) {
                                    const auto [dist_squared, id, type] = Distance::point_to_triangle_squared( p, prim, primal_mesh );
                                    dist_squared_min                    = CGAL::min( dist_squared_min, dist_squared );
                                }
                                return CGAL::sqrt( dist_squared_min );
                            };

                            // find better approximation by performing a bisection search with only intersecting prims
                            Point p0  = p[v0];
                            Point p1  = p[v1];
                            FT s0     = s[v0];
                            FT s1     = s[v1];
                            FT s_lerp = get_distance_to_intersecting_prims( p_lerp );

                            for( int i = 0; i < 4; ++i ) {
                                if( CGAL::abs( s_lerp * s_lerp - iso_value * iso_value ) < max_step_size * 1e-8 ) {
                                    break;
                                }
                                const Point p_mid = CGAL::midpoint( p0, p1 );
                                const FT s_mid    = get_distance_to_intersecting_prims( p_mid );

                                const bool b0    = s0 <= iso_value;
                                const bool b1    = s1 <= iso_value;
                                const bool b_mid = s_mid <= iso_value;
                                if( b0 == b_mid ) {
                                    s0 = s_mid;
                                    p0 = p_mid;
                                } else if( b1 == b_mid ) {
                                    s1 = s_mid;
                                    p1 = p_mid;
                                } else {
                                    LOG( ERROR ) << "Could not identify mid point in QEM_SVD_hermite";
                                    break;
                                }

                                const FT u = ( s0 - iso_value ) / ( s0 - s1 );
                                p_lerp     = p0 + u * ( p1 - p0 );
                                s_lerp     = get_distance_to_intersecting_prims( p_lerp );
                            }

                            edge_intersections.push_back( p_lerp );

                            Vector normal( 0, 0, 0 );
                            FT min_dist_squared     = FT_MAX;
                            auto& prims_for_normals = intersecting_prims.empty() ? prims : intersecting_prims;
                            for( const auto& prim: prims_for_normals ) {
                                const Point p_lerp_proj = Distance::project_to_primitive( p_lerp, prim, primal_mesh );
                                const FT squared_dist   = ( p_lerp - p_lerp_proj ).squared_length();
                                // const auto [dist, id, type] = Distance::point_to_primitive( p_lerp, prim, primal_mesh );
                                if( squared_dist < min_dist_squared ) {
                                    min_dist_squared = squared_dist;
                                    // normal           = Distance::get_prim_normal( prim, primal_mesh, p_lerp );
                                    normal = normalize( p_lerp - p_lerp_proj );
                                }
                            }
                            edge_intersection_normals.push_back( normal );
                        }
                    }

                    // LOG_ASSERT( !edge_intersections.empty() ) << " prims size = " << prims.size();
                    if( edge_intersections.empty() ) {
                        LOG( WARNING ) << "No edge intersections for " << v;
                        // there is supposed to be a vertex but there were no cell intersections
                        // vertex cannot be placed using qem
                        return 1;
                    }

                    // MC Polygon Center of Mass
                    Point point = CGAL::centroid( edge_intersections.begin(), edge_intersections.end() );

                    // SVD QEM
                    bool qem_success = qem( point, edge_intersections, edge_intersection_normals );

                    bool use_projection = !qem_success;
                    // bbox
                    if constexpr( use_bbox ) {
                        CGAL::Bbox_3 bbox = p[0].bbox() + p[7].bbox();
                        if( point.x() < bbox.xmin() || point.x() > bbox.xmax() ) {
                            use_projection = true;
                        } else if( point.y() < bbox.ymin() || point.y() > bbox.ymax() ) {
                            use_projection = true;
                        } else if( point.z() < bbox.zmin() || point.z() > bbox.zmax() ) {
                            use_projection = true;
                        }
                    }

                    if( use_projection ) {
                        point = CGAL::midpoint( p[0], p[7] );    // set point to voxel center

                        // in case of clipping, project cell center onto the offset of the nearest prim
                        Point p_nearest;
                        FT min_dist_squared = FT_MAX;
                        for( const auto& prim: prims ) {
                            const Point p_proj    = Distance::project_to_primitive( point, prim, primal_mesh );
                            const FT squared_dist = ( point - p_proj ).squared_length();
                            if( squared_dist < min_dist_squared ) {
                                min_dist_squared = squared_dist;
                                p_nearest        = p_proj;
                            }
                        }

                        Vector n = normalize( point - p_nearest );
                        if( is_inside( point ) == CGAL::ON_BOUNDED_SIDE ) {
                            n *= -1;
                        }
                        point = p_nearest + iso_value * n;

                        if constexpr( use_bbox ) {
                            CGAL::Bbox_3 bbox = p[0].bbox() + p[7].bbox();
                            FT x              = std::min<FT>( std::max<FT>( point.x(), bbox.xmin() ), bbox.xmax() );
                            FT y              = std::min<FT>( std::max<FT>( point.y(), bbox.ymin() ), bbox.ymax() );
                            FT z              = std::min<FT>( std::max<FT>( point.z(), bbox.zmin() ), bbox.zmax() );
                            point             = Point( x, y, z );
                        }
                    }

                    offset_mesh.point( v ) = point;

                    if( use_projection ) {
                        // there is a vertex but it could not be placed inside the cell with qem
                        return 1;
                    } else {
                        // this vertex could be placed inside the cell with qem
                        return 2;
                    }
                }

                int position( const Octree_wrapper& octree, const SurfaceMesh& primal_mesh, const FT iso_value,
                              const Octree_wrapper::Voxel_handle& vh, Point& point ) const {
                    const auto& prims = octree.voxel_prims( vh );

                    if( prims.empty() ) {
                        return false;
                    }

                    std::array<FT, Tables::N_VERTICES> s = octree.voxel_values( vh );

                    for( const auto& e: s ) {
                        if( e == std::numeric_limits<FT>::max() ) {
                            return false;
                        }
                    }

                    std::array<bool, Tables::N_VERTICES> b;
                    std::transform( s.begin(), s.end(), b.begin(), [iso_value]( const auto& e ) { return e <= iso_value; } );

                    unsigned int cubeindex = 0;
                    // set bit if corresponding corner is below iso
                    for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                        cubeindex |= b[i] << i;
                    }

                    if( cubeindex == 0 || cubeindex == 255 ) {
                        // no intersection of the cell with the surface visible
                        return 0;
                    }

                    std::array<Point, Tables::N_VERTICES> p = octree.voxel_vertex_positions( vh );
                    std::array<Vector, Tables::N_VERTICES> pos;
                    std::transform( p.begin(), p.end(), pos.begin(), []( const auto& e ) { return e - CGAL::ORIGIN; } );
                    std::array<FT, Tables::N_VERTICES> sign;
                    std::transform( p.begin(), p.end(), sign.begin(), [this]( const Point& e ) { return is_inside( e ) ? -1 : 1; } );

                    const FT max_step_size = ( p[0] - p[7] ).squared_length();

                    // compute edge intersections
                    std::vector<Point> edge_intersections;
                    std::vector<Vector> edge_intersection_normals;

                    for( int i = 0; i < Tables::N_EDGES; ++i ) {
                        const auto& v0 = Tables::edge_to_vertex[i][0];
                        const auto& v1 = Tables::edge_to_vertex[i][1];

                        if( b[v0] != b[v1] ) {    // e0
                            const FT u   = ( s[v0] - iso_value ) / ( s[v0] - s[v1] );
                            Point p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[v0] + u * pos[v1] );

                            // only use prims that intersect the edge for interpolation
                            std::vector<Face_descriptor> intersecting_prims;
                            for( const auto& prim: prims ) {
                                const Face_descriptor f( prim );
                                // check if prim intersect this edge
                                // const auto [dist0, id0, type0] = Distance::point_to_triangle( p[v0], f, primal_mesh );
                                // const auto [dist1, id1, type1] = Distance::point_to_triangle( p[v1], f, primal_mesh );
                                const Point proj0 = Distance::project_to_triangle( p[v0], f, primal_mesh );
                                const Point proj1 = Distance::project_to_triangle( p[v1], f, primal_mesh );
                                const FT dist0    = CGAL::sqrt( ( p[v0] - proj0 ).squared_length() );
                                const FT dist1    = CGAL::sqrt( ( p[v1] - proj1 ).squared_length() );
                                const bool b0     = ( sign[v0] * dist0 ) <= iso_value;
                                const bool b1     = ( sign[v1] * dist1 ) <= iso_value;
                                if( b0 != b1 ) {
                                    intersecting_prims.push_back( f );
                                    continue;
                                }
                                // TODO project to line
                                // check if projection onto line is inside the radius
                                const Point e_proj0 = Distance::project_to_line( proj0, p[v0], p[v1] );
                                if( ( proj0 - e_proj0 ).squared_length() < iso_value * iso_value ) {
                                    intersecting_prims.push_back( f );
                                    continue;
                                }
                                const Point e_proj1 = Distance::project_to_line( proj1, p[v0], p[v1] );
                                if( ( proj1 - e_proj1 ).squared_length() < iso_value * iso_value ) {
                                    intersecting_prims.push_back( f );
                                    continue;
                                }
                            }

                            if( intersecting_prims.empty() ) {
                                continue;
                            }

                            auto get_distance_to_intersecting_prims = [&intersecting_prims, &primal_mesh]( const Point& p ) {
                                FT dist_squared_min = FT_MAX;
                                for( const auto& prim: intersecting_prims ) {
                                    const auto [dist_squared, id, type] = Distance::point_to_triangle_squared( p, prim, primal_mesh );
                                    dist_squared_min                    = CGAL::min( dist_squared_min, dist_squared );
                                }
                                return CGAL::sqrt( dist_squared_min );
                            };

                            // find better approximation by performing a bisection search with only intersecting prims
                            Point p0  = p[v0];
                            Point p1  = p[v1];
                            FT s0     = s[v0];
                            FT s1     = s[v1];
                            FT s_lerp = get_distance_to_intersecting_prims( p_lerp );

                            for( int i = 0; i < 4; ++i ) {
                                if( CGAL::abs( s_lerp * s_lerp - iso_value * iso_value ) < max_step_size * 1e-8 ) {
                                    break;
                                }
                                const Point p_mid = CGAL::midpoint( p0, p1 );
                                const FT s_mid    = get_distance_to_intersecting_prims( p_mid );

                                const bool b0    = s0 <= iso_value;
                                const bool b1    = s1 <= iso_value;
                                const bool b_mid = s_mid <= iso_value;
                                if( b0 == b_mid ) {
                                    s0 = s_mid;
                                    p0 = p_mid;
                                } else if( b1 == b_mid ) {
                                    s1 = s_mid;
                                    p1 = p_mid;
                                } else {
                                    LOG( ERROR ) << "Could not identify mid point in QEM_SVD_hermite";
                                    break;
                                }

                                const FT u = ( s0 - iso_value ) / ( s0 - s1 );
                                p_lerp     = p0 + u * ( p1 - p0 );
                                s_lerp     = get_distance_to_intersecting_prims( p_lerp );
                            }

                            edge_intersections.push_back( p_lerp );

                            Vector normal( 0, 0, 0 );
                            FT min_dist_squared     = FT_MAX;
                            auto& prims_for_normals = intersecting_prims.empty() ? prims : intersecting_prims;
                            for( const auto& prim: prims_for_normals ) {
                                const Point p_lerp_proj = Distance::project_to_primitive( p_lerp, prim, primal_mesh );
                                const FT squared_dist   = ( p_lerp - p_lerp_proj ).squared_length();
                                // const auto [dist, id, type] = Distance::point_to_primitive( p_lerp, prim, primal_mesh );
                                if( squared_dist < min_dist_squared ) {
                                    min_dist_squared = squared_dist;
                                    // normal           = Distance::get_prim_normal( prim, primal_mesh, p_lerp );
                                    normal = normalize( p_lerp - p_lerp_proj );
                                }
                            }
                            edge_intersection_normals.push_back( normal );
                        }
                    }

                    // LOG_ASSERT( !edge_intersections.empty() ) << " prims size = " << prims.size();
                    if( edge_intersections.empty() ) {
                        LOG( WARNING ) << "No edge intersections for cell " << vh;
                        // there is supposed to be a vertex but there were no cell intersections
                        // vertex cannot be placed using qem
                        return 1;
                    }

                    // MC Polygon Center of Mass
                    point = CGAL::centroid( edge_intersections.begin(), edge_intersections.end() );

                    // check if normals point in a similar direction
                    FT min_scalar_prod = FT_MAX;
                    for( size_t i = 0; i < edge_intersection_normals.size(); ++i ) {
                        for( size_t j = 0; j < edge_intersection_normals.size(); ++j ) {
                            FT prod         = CGAL::scalar_product( edge_intersection_normals[i], edge_intersection_normals[j] );
                            min_scalar_prod = CGAL::min( prod, min_scalar_prod );
                        }
                    }

                    bool use_projection = false;
                    if( min_scalar_prod > -0.5 ) {
                        // SVD QEM
                        qem( point, edge_intersections, edge_intersection_normals );
                    } else {
                        use_projection = true;
                    }

                    // bbox
                    if constexpr( use_bbox ) {
                        CGAL::Bbox_3 bbox = p[0].bbox() + p[7].bbox();
                        if( point.x() < bbox.xmin() || point.x() > bbox.xmax() ) {
                            use_projection = true;
                        } else if( point.y() < bbox.ymin() || point.y() > bbox.ymax() ) {
                            use_projection = true;
                        } else if( point.z() < bbox.zmin() || point.z() > bbox.zmax() ) {
                            use_projection = true;
                        }
                    }

                    if( use_projection ) {
                        point = CGAL::midpoint( p[0], p[7] );    // set point to voxel center

                        // in case of clipping, project cell center onto the offset of the nearest prim
                        Point p_nearest;
                        FT min_dist_squared = FT_MAX;
                        for( const auto& prim: prims ) {
                            const Point p_proj    = Distance::project_to_primitive( point, prim, primal_mesh );
                            const FT squared_dist = ( point - p_proj ).squared_length();
                            if( squared_dist < min_dist_squared ) {
                                min_dist_squared = squared_dist;
                                p_nearest        = p_proj;
                            }
                        }

                        Vector n = normalize( point - p_nearest );
                        if( is_inside( point ) == CGAL::ON_BOUNDED_SIDE ) {
                            n *= -1;
                        }
                        point = p_nearest + iso_value * n;

                        if constexpr( use_bbox ) {
                            CGAL::Bbox_3 bbox = p[0].bbox() + p[7].bbox();
                            FT x              = std::min<FT>( std::max<FT>( point.x(), bbox.xmin() ), bbox.xmax() );
                            FT y              = std::min<FT>( std::max<FT>( point.y(), bbox.ymin() ), bbox.ymax() );
                            FT z              = std::min<FT>( std::max<FT>( point.z(), bbox.zmin() ), bbox.zmax() );
                            point             = Point( x, y, z );
                        }
                    }

                    if( min_scalar_prod > -0.5 ) {
                        // this vertex could be placed with qem
                        return 2;
                    } else {
                        // there is a vertex but it could not be placed with qem
                        return 1;
                    }
                }

              private:
                bool qem( Point& point, const std::vector<Point>& edge_intersections, const std::vector<Vector>& edge_intersection_normals ) const {
                    Eigen::Matrix3d A;
                    A.setZero();
                    Eigen::Vector3d b;
                    b.setZero();
                    for( int i = 0; i < edge_intersections.size(); ++i ) {
                        Eigen::Vector3d n_k = { edge_intersection_normals[i].x(), edge_intersection_normals[i].y(),
                                                edge_intersection_normals[i].z() };
                        Eigen::Vector3d p_k = { edge_intersections[i].x(), edge_intersections[i].y(), edge_intersections[i].z() };
                        double d_k          = n_k.transpose() * p_k;

                        Eigen::Matrix3d A_k = n_k * n_k.transpose();
                        Eigen::Vector3d b_k = d_k * n_k;
                        A += A_k;
                        b += b_k;
                    }

                    Eigen::JacobiSVD<Eigen::Matrix3d> svd( A, Eigen::ComputeFullU | Eigen::ComputeFullV );
                    // set threshold as in Peter Lindstrom's paper, "Out-of-Core
                    // Simplification of Large Polygonal Models"
                    svd.setThreshold( 1e-3 );

                    // Init x hat
                    Eigen::Vector3d x_hat;
                    x_hat << point.x(), point.y(), point.z();

                    // Lindstrom formula for QEM new position for singular matrices
                    Eigen::Vector3d v_svd = x_hat + svd.solve( b - A * x_hat );

                    // compute residual
                    const double residual = ( b - A * v_svd ).squaredNorm();

                    if( residual > 1e-12 ) {
                        return false;
                    }

                    point = Point( v_svd[0], v_svd[1], v_svd[2] );
                    return true;
                }
            };

            class Voxel_center {
              public:
                /// <summary>
                /// Compute vertex position for Dual Contouring
                /// </summary>
                /// <typeparam name="Domain_"></typeparam>
                /// <param name="domain"></param>
                /// <param name="iso_value"></param>
                /// <param name="i"></param>
                /// <param name="j"></param>
                /// <param name="k"></param>
                /// <returns> true, if there is a point in the cell</returns>
                bool position( const Octree_wrapper& octree, const FT iso_value, const Octree_wrapper::Voxel_handle& vh, Point& point ) const {
                    const auto& prims = octree.voxel_prims( vh );

                    if( prims.empty() ) {
                        return false;
                    }

                    std::array<FT, Tables::N_VERTICES> s = octree.voxel_values( vh );

                    for( const auto& e: s ) {
                        if( e == std::numeric_limits<FT>::max() ) {
                            return false;
                        }
                    }

                    std::array<bool, Tables::N_VERTICES> b;
                    std::transform( s.begin(), s.end(), b.begin(), [iso_value]( const auto& e ) { return e <= iso_value; } );

                    unsigned int cubeindex = 0;
                    // set bit if corresponding corner is below iso
                    for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                        cubeindex |= b[i] << i;
                    }

                    if( cubeindex == 0 || cubeindex == 255 ) {
                        return false;
                    }

                    std::array<Point, Tables::N_VERTICES> p = octree.voxel_vertex_positions( vh );
                    std::array<Vector, Tables::N_VERTICES> pos;
                    std::transform( p.begin(), p.end(), pos.begin(), []( const auto& e ) { return e - CGAL::ORIGIN; } );

                    point = CGAL::ORIGIN + ( pos[0] + 0.5 * ( pos[7] - pos[0] ) );    // set point to voxel center

                    return true;
                }
            };

            class MC_polygon_center {
              public:
                /// <summary>
                /// Compute vertex position for Dual Contouring
                /// </summary>
                /// <typeparam name="Domain_"></typeparam>
                /// <param name="domain"></param>
                /// <param name="iso_value"></param>
                /// <param name="i"></param>
                /// <param name="j"></param>
                /// <param name="k"></param>
                /// <returns> true, if there is a point in the cell</returns>
                bool position( const Octree_wrapper& octree, const FT iso_value, const Octree_wrapper::Voxel_handle& vh, Point& point ) const {
                    const auto& prims = octree.voxel_prims( vh );

                    if( prims.empty() ) {
                        return false;
                    }

                    std::array<FT, Tables::N_VERTICES> s = octree.voxel_values( vh );

                    for( const auto& e: s ) {
                        if( e == std::numeric_limits<FT>::max() ) {
                            return false;
                        }
                    }

                    std::array<bool, Tables::N_VERTICES> b;
                    std::transform( s.begin(), s.end(), b.begin(), [iso_value]( const auto& e ) { return e <= iso_value; } );

                    unsigned int cubeindex = 0;
                    // set bit if corresponding corner is below iso
                    for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                        cubeindex |= b[i] << i;
                    }

                    if( cubeindex == 0 || cubeindex == 255 ) {
                        return false;
                    }

                    std::array<Point, Tables::N_VERTICES> p = octree.voxel_vertex_positions( vh );
                    std::array<Vector, Tables::N_VERTICES> pos;
                    std::transform( p.begin(), p.end(), pos.begin(), []( const auto& e ) { return e - CGAL::ORIGIN; } );

                    point = CGAL::midpoint( p[0], p[7] );    // set point to voxel center

                    // compute edge intersections
                    std::vector<Point> edge_intersections;

                    for( int i = 0; i < Tables::N_EDGES; ++i ) {
                        const auto& v0 = Tables::edge_to_vertex[i][0];
                        const auto& v1 = Tables::edge_to_vertex[i][1];

                        if( b[v0] != b[v1] ) {    // e0
                            const FT u         = ( s[v0] - iso_value ) / ( s[v0] - s[v1] );
                            const Point p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[v0] + u * pos[v1] );
                            edge_intersections.push_back( p_lerp );
                        }
                    }

                    point = CGAL::centroid( edge_intersections.begin(), edge_intersections.end() );

                    return true;
                }
            };

        }    // namespace Positioning

        const size_t INVALID_ID = -1;

        struct Halfedge {
            size_t id       = INVALID_ID;
            size_t face     = INVALID_ID;
            size_t from     = INVALID_ID;    // from vertex
            size_t to       = INVALID_ID;    // to vertex
            size_t next     = INVALID_ID;    // next halfedge id
            size_t opposite = INVALID_ID;    // opposite halfedge id
        };

        struct Face {
            size_t id       = INVALID_ID;    // face id
            size_t halfedge = INVALID_ID;    // halfedge id
            Vector normal;
        };

        struct Vertex {
            size_t id       = INVALID_ID;    // vertex id
            size_t halfedge = INVALID_ID;    // ingoing halfedge id
        };

        inline SurfaceMesh resolve_non_manifoldness_on_polygon_soup( std::vector<Point>& points, std::vector<std::vector<size_t>>& polygons ) {
            // implement halfedge data structure

            std::vector<Vertex> vertices;
            vertices.reserve( points.size() );
            for( const auto& p: points ) {
                vertices.push_back( { vertices.size() } );
            }

            std::vector<Face> faces;
            faces.reserve( polygons.size() );
            std::vector<Halfedge> halfedges;
            halfedges.reserve( polygons.size() * 4 );    // assume input are quads

            std::multimap<std::pair<size_t, size_t>, size_t> map_vertices_to_halfedge_id;

            // generate halfedge data-structure (opposite halfedges are added later)
            for( const auto& pol: polygons ) {
                size_t h0_id = halfedges.size();
                for( int i = 0; i < pol.size(); ++i ) {
                    size_t v0 = pol[i];
                    size_t v1 = pol[( i + 1 ) % pol.size()];
                    halfedges.push_back( { h0_id + i, faces.size(), v0, v1, h0_id + ( ( i + 1 ) % pol.size() ) } );
                    map_vertices_to_halfedge_id.insert( { { v0, v1 }, h0_id + i } );
                    if( vertices[v1].halfedge == INVALID_ID ) {
                        vertices[v1].halfedge = h0_id + i;
                    }
                }

                // compute normal
                const Point& p0 = points[pol[0]];
                const Point& p1 = points[pol[1]];
                const Point& p2 = points[pol[2]];
                const Vector n  = normalize( CGAL::cross_product( ( p1 - p0 ), ( p2 - p0 ) ) );

                faces.push_back( { faces.size(), h0_id, n } );
            }

            // set opposite halfedges
            for( auto& he: halfedges ) {
                if( he.opposite != INVALID_ID ) {
                    continue;
                }

                const size_t v0 = he.from;
                const size_t v1 = he.to;
                if( map_vertices_to_halfedge_id.count( { v1, v0 } ) == 1 ) {
                    he.opposite = map_vertices_to_halfedge_id.find( { v1, v0 } )->second;
                } else {
                    // resolve ambiguity /  non-manifold edges
                    const Point p0           = points[v0];
                    const Point p1           = points[v1];
                    const Vector edge_vector = normalize( p1 - p0 );
                    const Vector n1          = faces[he.face].normal;

                    std::multimap<FT, size_t> candidates;

                    auto it = map_vertices_to_halfedge_id.find( { v1, v0 } );
                    while( it->first == std::pair { v1, v0 } ) {
                        const Halfedge& he_opp = halfedges[it->second];
                        const Vector n2        = faces[he_opp.face].normal;
                        // atan2((Va x Vb) . Vn, Va . Vb) <-- from
                        // https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
                        FT angle = std::atan2( CGAL::scalar_product( CGAL::cross_product( n1, n2 ), edge_vector ), CGAL::scalar_product( n1, n2 ) );
                        candidates.insert( { -angle, he_opp.id } );
                        ++it;
                    }

                    // size_t he_opposite = candidates.begin()->second;
                    size_t he_opposite = INVALID_ID;
                    for( const auto [angle, cand]: candidates ) {
                        if( halfedges[cand].opposite == INVALID_ID ) {
                            he_opposite = cand;
                            break;
                        }
                    }
                    LOG_ASSERT( he_opposite != INVALID_ID );
                    // take the one with the smallest angle
                    he.opposite                     = he_opposite;
                    halfedges[he_opposite].opposite = he.id;
                }
            }

            // sanity check, all halfedges have an opposite now
            for( const auto& he: halfedges ) {
                LOG_ASSERT( he.opposite != INVALID_ID );
            }

            // split non-manifold edges
            size_t halfedges_size_with_nmes = halfedges.size();
            for( size_t he = 0; he < halfedges_size_with_nmes; ++he ) {
                const size_t v0 = halfedges[he].from;
                const size_t v1 = halfedges[he].to;
                if( map_vertices_to_halfedge_id.count( { v0, v1 } ) < 2 ) {
                    continue;
                }

                const size_t he_opp = halfedges[he].opposite;

                const size_t f0 = halfedges[he].face;
                const size_t f1 = halfedges[he_opp].face;

                const size_t h_new     = halfedges.size();
                const size_t h_new_opp = halfedges.size() + 1;

                const Point& p0    = points[v0];
                const Point& p1    = points[v1];
                const Point p_mid  = CGAL::midpoint( p0, p1 );
                const size_t v_mid = vertices.size();
                vertices.push_back( { v_mid, h_new } );
                points.push_back( p_mid );

                // insert new halfedges
                halfedges.push_back( { h_new, f0, v0, v_mid, he, h_new_opp } );
                halfedges.push_back( { h_new_opp, f1, v_mid, v0, halfedges[he_opp].next, h_new } );
                // modify old halfedges
                halfedges[he].from     = v_mid;
                halfedges[he_opp].next = h_new_opp;
                halfedges[he_opp].to   = v_mid;
                size_t he_it           = he;
                do {
                    if( halfedges[he_it].next == he ) {
                        halfedges[he_it].next = h_new;
                        break;
                    }
                    he_it = halfedges[he_it].next;
                } while( he_it != he );
                // modify vertex
                vertices[v0].halfedge = h_new_opp;
            }

            // detect and duplicate non-manifold vertices
            std::map<size_t, std::set<size_t>> vertex_faces;
            for( const auto& f: faces ) {
                const size_t he_init = f.halfedge;
                size_t he            = he_init;
                do {
                    vertex_faces[halfedges[he].to].insert( f.id );
                    he = halfedges[he].next;
                } while( he != he_init );
            }

            size_t vertices_size_with_nmvs = vertices.size();
            for( size_t v_id = 0; v_id < vertices_size_with_nmvs; ++v_id ) {
                auto& incident_face_ids = vertex_faces[v_id];

                while( !incident_face_ids.empty() ) {
                    const size_t he_init = vertices[v_id].halfedge;
                    size_t he            = he_init;

                    do {
                        size_t f_id = halfedges[he].face;

                        // make sure f_id exists, then remove it from incident face list
                        auto it = incident_face_ids.find( f_id );
                        LOG_ASSERT( it != incident_face_ids.end() );
                        incident_face_ids.erase( it );

                        // go to next face
                        he = halfedges[he].next;
                        he = halfedges[he].opposite;
                    } while( he != he_init );

                    if( !incident_face_ids.empty() ) {
                        // this is a non-manifold vertex --> duplicate vertex and replace its id in this umbrella
                        const size_t duplicate = vertices.size();
                        vertices.push_back( vertices[v_id] );
                        vertices[duplicate].id = duplicate;
                        points.push_back( points[v_id] );

                        // change ids in umbrella halfedges
                        do {
                            halfedges[he].to   = duplicate;
                            he                 = halfedges[he].next;
                            halfedges[he].from = duplicate;
                            he                 = halfedges[he].opposite;
                        } while( he != he_init );

                        // fix halfedge for v_id
                        const size_t he_init = faces[*( incident_face_ids.begin() )].halfedge;
                        size_t he            = he_init;
                        do {
                            if( halfedges[he].to == v_id ) {
                                vertices[v_id].halfedge = he;
                                break;
                            }
                            he = halfedges[he].next;
                        } while( he != he_init );
                    }
                }
            }

            // convert into SurfaceMesh
            SurfaceMesh m;
            for( const auto& p: points ) {
                m.add_vertex( p );
            }

            // add faces
            for( const auto& f: faces ) {
                std::vector<Vertex_descriptor> f_vertices;
                size_t he_init = f.halfedge;
                size_t he      = he_init;
                do {
                    f_vertices.push_back( Vertex_descriptor( halfedges[he].to ) );
                    he = halfedges[he].next;
                } while( he != he_init );

                if( m.add_face( f_vertices ) == SurfaceMesh::null_face() ) {
                    LOG( WARNING ) << "Could not add face";
                }
            }

            size_t n_remaining_nmvs = PMP::duplicate_non_manifold_vertices( m );
            if( n_remaining_nmvs > 0 ) {
                LOG( INFO ) << "PMP removed another " << n_remaining_nmvs << " non-manifold vertices";
            }

            return m;
        }

        inline std::set<size_t> find_non_manifoldness_on_polygon_soup( const std::vector<Point>& points,
                                                                       const std::vector<std::vector<size_t>>& polygons ) {
            std::set<size_t> nmvs;

            std::vector<Vertex> vertices;
            vertices.reserve( points.size() );
            for( const auto& p: points ) {
                vertices.push_back( { vertices.size() } );
            }

            std::vector<Face> faces;
            faces.reserve( polygons.size() );
            std::vector<Halfedge> halfedges;
            halfedges.reserve( polygons.size() * 4 );    // assume input are quads

            std::multimap<std::pair<size_t, size_t>, size_t> map_vertices_to_halfedge_id;

            // generate halfedge data-structure (opposite halfedges are added later)
            for( const auto& pol: polygons ) {
                size_t h0_id = halfedges.size();
                for( int i = 0; i < pol.size(); ++i ) {
                    size_t v0 = pol[i];
                    size_t v1 = pol[( i + 1 ) % pol.size()];
                    halfedges.push_back( { h0_id + i, faces.size(), v0, v1, h0_id + ( ( i + 1 ) % pol.size() ) } );
                    map_vertices_to_halfedge_id.insert( { { v0, v1 }, h0_id + i } );
                    if( vertices[v1].halfedge == INVALID_ID ) {
                        vertices[v1].halfedge = h0_id + i;
                    }
                }

                // compute normal
                const Point& p0 = points[pol[0]];
                const Point& p1 = points[pol[1]];
                const Point& p2 = points[pol[2]];
                const Vector n  = normalize( CGAL::cross_product( ( p1 - p0 ), ( p2 - p0 ) ) );

                faces.push_back( { faces.size(), h0_id, n } );
            }

            // set opposite halfedges
            for( auto& he: halfedges ) {
                if( he.opposite != INVALID_ID ) {
                    continue;
                }

                const size_t v0 = he.from;
                const size_t v1 = he.to;
                if( map_vertices_to_halfedge_id.count( { v1, v0 } ) == 1 ) {
                    he.opposite = map_vertices_to_halfedge_id.find( { v1, v0 } )->second;
                } else {
                    nmvs.insert( v0 );
                    nmvs.insert( v1 );

                    // resolve ambiguity /  non-manifold edges
                    const Point p0           = points[v0];
                    const Point p1           = points[v1];
                    const Vector edge_vector = normalize( p1 - p0 );
                    const Vector n1          = faces[he.face].normal;

                    std::multimap<FT, size_t> candidates;

                    auto it = map_vertices_to_halfedge_id.find( { v1, v0 } );
                    while( it->first == std::pair { v1, v0 } ) {
                        const Halfedge& he_opp = halfedges[it->second];
                        const Vector n2        = faces[he_opp.face].normal;
                        // atan2((Va x Vb) . Vn, Va . Vb) <-- from
                        // https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
                        FT angle = std::atan2( CGAL::scalar_product( CGAL::cross_product( n1, n2 ), edge_vector ), CGAL::scalar_product( n1, n2 ) );
                        candidates.insert( { -angle, he_opp.id } );
                        ++it;
                    }

                    // size_t he_opposite = candidates.begin()->second;
                    size_t he_opposite = INVALID_ID;
                    for( const auto [angle, cand]: candidates ) {
                        if( halfedges[cand].opposite == INVALID_ID ) {
                            he_opposite = cand;
                            break;
                        }
                    }
                    LOG_ASSERT( he_opposite != INVALID_ID );

                    // take the one with the smallest angle
                    he.opposite                     = he_opposite;
                    halfedges[he_opposite].opposite = he.id;
                }
            }

            // sanity check, all halfedges have an opposite now
            for( const auto& he: halfedges ) {
                LOG_ASSERT( he.opposite != INVALID_ID );
            }

            // detect and duplicate non-manifold vertices
            std::map<size_t, std::set<size_t>> vertex_faces;
            for( const auto& f: faces ) {
                const size_t he_init = f.halfedge;
                size_t he            = he_init;
                do {
                    vertex_faces[halfedges[he].to].insert( f.id );
                    he = halfedges[he].next;
                } while( he != he_init );
            }

            for( size_t v_id = 0; v_id < vertices.size(); ++v_id ) {
                auto& incident_face_ids = vertex_faces[v_id];

                const size_t he_init = vertices[v_id].halfedge;
                size_t he            = he_init;

                do {
                    size_t f_id = halfedges[he].face;
                    // make sure f_id exists, then remove it from incident face list
                    auto it = incident_face_ids.find( f_id );
                    LOG_ASSERT( it != incident_face_ids.end() ) << "\nv = " << v_id;
                    incident_face_ids.erase( it );

                    // go to next face
                    he = halfedges[he].next;
                    he = halfedges[he].opposite;
                } while( he != he_init );

                // if( v_id == 2045 ) {
                //     LOG( WARNING ) << "DEBUG " << v_id << ", p = " << points[v_id];
                //     LOG( WARNING ) << "is manifold: " << incident_face_ids.empty();
                // }

                if( !incident_face_ids.empty() ) {
                    // this is a non-manifold vertex
                    nmvs.insert( v_id );
                }
            }

            return nmvs;
        }

    }    // namespace Dual_contouring_3

    struct Visitor : public CGAL::Polygon_mesh_processing::Default_orientation_visitor {
        void non_manifold_edge( std::size_t id1, std::size_t id2, std::size_t nb_poly ) {
            std::cout << "The edge " << id1 << ", " << id2 << " is not manifold: " << nb_poly << " incident polygons." << std::endl;
        }
        void non_manifold_vertex( std::size_t id, std::size_t nb_cycles ) {
            std::cout << "The vertex " << id << " is not manifold: " << nb_cycles << " connected components of vertices in the link." << std::endl;
        }
        void duplicated_vertex( std::size_t v1, std::size_t v2 ) {
            std::cout << "The vertex " << v1 << " has been duplicated, its new id is " << v2 << "." << std::endl;
        }
        void vertex_id_in_polygon_replaced( std::size_t p_id, std::size_t i1, std::size_t i2 ) {
            std::cout << "In the polygon " << p_id << ", the index " << i1 << " has been replaced by " << i2 << "." << std::endl;
        }
        void polygon_orientation_reversed( std::size_t p_id ) { std::cout << "The polygon " << p_id << " has been reversed." << std::endl; }
    };

    template<class Positioning>
    void make_quad_mesh_using_dual_contouring( const Octree_wrapper& octree, const FT iso_value, SurfaceMesh& m, const Positioning& positioning ) {
        // compute dc-vertices
        std::map<size_t, size_t> map_voxel_to_point_id;
        std::map<size_t, Point> map_voxel_to_point;
        size_t points_counter = 0;
        std::map<Octree_wrapper::Edge_handle, std::array<size_t, 4>> quads;

        const std::size_t n_voxel = octree.n_voxels();

        // save all points
        for( int i = 0; i < n_voxel; ++i ) {
            const auto& vh = octree.voxels( i );
            Point p;
            if( positioning.position( octree, iso_value, vh, p ) ) {
                map_voxel_to_point[vh]    = p;
                map_voxel_to_point_id[vh] = points_counter++;
            }
        }

        const std::size_t n_edges = octree.n_edges();

        // save all quads
        for( int i = 0; i < n_edges; ++i ) {
            const auto& e        = octree.edges( i );
            const auto& [s0, s1] = octree.edge_values( e );

            if( s0 == std::numeric_limits<FT>::max() || s1 == std::numeric_limits<FT>::max() ) {
                continue;
            }

            if( s0 <= iso_value && s1 > iso_value ) {
                const auto voxels = octree.voxels_incident_to_edge( e );
                quads[e]          = voxels;
            } else if( s1 <= iso_value && s0 > iso_value ) {
                const auto voxels = octree.voxels_incident_to_edge( e );
                quads[e]          = { voxels[0], voxels[3], voxels[2], voxels[1] };
            }
        }

        // add points to empty voxels
        size_t n_empty_voxels = 0;
        std::map<Octree_wrapper::Edge_handle, std::array<size_t, 4>> quads_buf;
        for( const auto& [id, q]: quads ) {
            bool ignore_quad = false;

            // ignore quads that only have empty voxels
            bool has_non_empty_voxels = false;
            for( const auto& cell_id: q ) {
                if( !octree.voxel_prims( cell_id ).empty() ) {
                    has_non_empty_voxels = true;
                    break;
                }
            }

            if( !has_non_empty_voxels ) {
                ignore_quad = true;
            } else {
                for( const auto& cell_id: q ) {
                    if( map_voxel_to_point_id.count( cell_id ) == 0 ) {
                        // LOG( WARNING ) << "Voxel " << cell_id << " does not have a point. Add point.";
                        ++n_empty_voxels;
                        // voxels_without_prims.push_back( cell_id );
                        const auto points              = octree.voxel_vertex_positions( cell_id );
                        const Point_3 p                = points[0] + 0.5 * ( points[7] - points[0] );    // set point to voxel center;
                        map_voxel_to_point[cell_id]    = p;
                        map_voxel_to_point_id[cell_id] = points_counter++;
                    }
                }
            }

            if( !ignore_quad ) {
                quads_buf[id] = q;
            }
        }
        quads = quads_buf;
        if( n_empty_voxels > 0 ) {
            LOG( WARNING ) << n_empty_voxels << " cells were empty. A point was added to them.";
        }

        std::vector<Point> mesh_points;
        mesh_points.reserve( map_voxel_to_point.size() );
        std::map<size_t, size_t> map_cell_id_to_mesh_point;
        std::map<size_t, size_t> map_point_id_to_cell_id;
        std::vector<std::vector<size_t>> mesh_polygons;
        mesh_polygons.reserve( quads.size() );

        for( const auto& [id, p]: map_voxel_to_point ) {
            map_cell_id_to_mesh_point[id]               = mesh_points.size();
            map_point_id_to_cell_id[mesh_points.size()] = id;
            mesh_points.push_back( p );
        }
        for( const auto& [id, q]: quads ) {
            std::vector<size_t> vertices;
            std::set<size_t> cell_set;    // make sure that vertices only appear once
            for( const auto& cell_id: q ) {
                if( map_voxel_to_point_id.count( cell_id ) == 0 ) {
                    LOG( WARNING ) << "Voxel " << cell_id << " does not have a point. Ignoring quad.";
                    continue;
                }
                if( cell_set.find( cell_id ) == cell_set.end() ) {
                    vertices.push_back( map_cell_id_to_mesh_point[cell_id] );
                    cell_set.insert( cell_id );
                }
            }
            mesh_polygons.push_back( vertices );
        }

        // sanity check - mesh_polygons only hold vertices that do exist
        {
            std::vector<size_t> point_used( mesh_points.size(), false );
            for( const auto& poly: mesh_polygons ) {
                for( const auto& v: poly ) {
                    if( v >= mesh_points.size() ) {
                        LOG( ERROR ) << "Vertex v=" << v << " must not exist. Sanity check failed!";
                    }
                    point_used[v] = true;
                }
            }

            for( int i = 0; i < point_used.size(); ++i ) {
                if( !point_used[i] ) {
                    LOG( ERROR ) << "Vertex v=" << i << "in cell " << map_point_id_to_cell_id[i] << " was never used. Sanity check failed!";
                }
            }
        }

        //////////////////////////////////
        //// find non-manifold vertices //
        // std::vector<std::vector<size_t>> vertex_polygons( map_voxel_to_point.size() );
        // for( size_t q_it = 0; q_it != mesh_polygons.size(); ++q_it ) {
        //     for( const auto& v: mesh_polygons[q_it] ) {
        //         vertex_polygons[v].push_back( q_it );
        //     }
        // }
        //
        // auto polygon_next_vertex_id = []( const size_t& vid, const std::vector<size_t>& quad ) {
        //    for( size_t i = 0; i < quad.size(); ++i ) {
        //        if( quad[i] == vid ) {
        //            return quad[( i + 1 ) % quad.size()];
        //        }
        //    }
        //    LOG( WARNING ) << "vertex id not found in quad";
        //    return std::numeric_limits<size_t>::max();
        //};
        //
        // auto polygon_prev_vertex_id = []( const size_t& vid, const std::vector<size_t>& quad ) {
        // for( size_t i = 0; i < quad.size(); ++i ) {
        //    if( quad[( i + 1 ) % quad.size()] == vid ) {
        //        return quad[i];
        //    }
        //}
        // LOG( WARNING ) << "vertex id not found in quad";
        // return std::numeric_limits<size_t>::max();
        //};
        //
        // auto polygon_has_he = []( const size_t& source, const size_t& target, const std::vector<size_t>& quad ) {
        // for( size_t i = 0; i < quad.size(); ++i ) {
        //    if( quad[i] == source && quad[( i + 1 ) % quad.size()] == target ) {
        //        return true;
        //    }
        //}
        // return false;
        //}
        //;
        //
        // size_t n_duplicate_vertices = 0;
        // for( size_t point_id = 0; point_id != vertex_polygons.size(); ++point_id ) {
        //     std::vector<size_t> polygon_ids = vertex_polygons[point_id];
        //     std::vector<size_t> umbrella_ids;    // oriented umbrella around vertex
        //     umbrella_ids.push_back( polygon_ids[0] );
        //     polygon_ids.erase( polygon_ids.begin() );
        //     // find the rest of the umbrella
        //     while( true ) {
        //         // ingoing halfedge is (vid,vo)
        //         const auto& curr_quad = mesh_polygons[umbrella_ids[umbrella_ids.size() - 1]];
        //         size_t vo             = polygon_prev_vertex_id( point_id, curr_quad );
        //
        //         // search for next polygon in umbrella
        //         size_t next_polygon_id = std::numeric_limits<size_t>::max();
        //         for( size_t i = 0; i < umbrella_ids.size() - 1; ++i ) {
        //             const auto& next_quad = mesh_polygons[umbrella_ids[i]];
        //             if( polygon_has_he( point_id, vo, next_quad ) ) {
        //                 next_polygon_id = i;
        //                 break;
        //             }
        //         }
        //
        //         // umbrella is closed --> break;
        //         if( next_polygon_id != std::numeric_limits<size_t>::max() ) {
        //             if( !polygon_ids.empty() ) {
        //                 // if next_polygon_id is non-zero, it means that some ids actually belong to the other fan
        //                 while( next_polygon_id != 0 ) {
        //                     polygon_ids.push_back( umbrella_ids[0] );
        //                     umbrella_ids.erase( umbrella_ids.begin() );
        //                     next_polygon_id--;
        //                 }
        //
        //                 // vertex is non-manifold --> duplicate it for all remaining polygon_ids
        //                 map_point_id_to_cell_id[mesh_points.size()] = map_point_id_to_cell_id[point_id];
        //                 mesh_points.push_back( mesh_points[point_id] );
        //                 // LOG( INFO ) << "Duplicate vertex to resolve non-manifoldness";
        //                 ++n_duplicate_vertices;
        //
        //                 for( const auto& poly_id: polygon_ids ) {
        //                     for( auto& v: mesh_polygons[poly_id] ) {
        //                         if( v == point_id ) {
        //                             v = mesh_points.size() - 1;
        //                         }
        //                     }
        //                 }
        //             }
        //             break;
        //         }
        //
        //         // find next quad
        //         for( size_t i = 0; i < polygon_ids.size(); ++i ) {
        //             const auto& next_quad = mesh_polygons[polygon_ids[i]];
        //             if( polygon_has_he( point_id, vo, next_quad ) ) {
        //                 next_polygon_id = i;
        //                 break;
        //             }
        //         }
        //
        //         // could not find next polygon (this happens for boundary vertices)
        //         if( next_polygon_id == std::numeric_limits<size_t>::max() ) {
        //             // LOG( WARNING ) << "Could not close umbrella for vertex " << point_id;
        //             break;
        //         }
        //
        //         umbrella_ids.push_back( polygon_ids[next_polygon_id] );
        //         polygon_ids.erase( polygon_ids.begin() + next_polygon_id );
        //     }
        // }
        // if( n_duplicate_vertices > 0 ) {
        //     LOG( INFO ) << n_duplicate_vertices << " vertices were dupicated to resolve non-manifoldness";
        // }

        //// create surface mesh with topology information
        // std::vector<std::vector<Vertex_descriptor>> surface_mesh_vertices( mesh_points.size() );
        // auto [v_voxel_handle, v_voxel_handle_created] = m.add_property_map<Vertex_descriptor, size_t>( "v:voxel_handle" );
        // LOG_ASSERT( v_voxel_handle_created );
        //
        // for( const auto& p: mesh_points ) {
        //     const Vertex_descriptor v = m.add_vertex( p );
        //     v_voxel_handle[v]         = map_point_id_to_cell_id[v.idx()];
        //     surface_mesh_vertices[v.idx()].push_back( v );
        // }
        // size_t n_null_faces = 0;
        // for( const auto& q: mesh_polygons ) {
        //     std::vector<Vertex_descriptor> vertex_ids;
        //     for( const auto& v_id: q ) {
        //         vertex_ids.push_back( surface_mesh_vertices[v_id][0] );
        //     }
        //
        //     if( m.add_face( vertex_ids ) == SurfaceMesh::null_face() ) {
        //         LOG( WARNING ) << "Could not add face (" << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << ") in dual contouring";
        //         ++n_null_faces;
        //     }
        // }
        // if( n_null_faces > 0 ) {
        //     LOG( WARNING ) << n_null_faces << " faces could not be added to the mesh";
        // }

        m = Dual_contouring_3::resolve_non_manifoldness_on_polygon_soup( mesh_points, mesh_polygons );

        std::map<Point, size_t> map_point_to_cell_id;
        for( const auto& [point_id, cell_id]: map_point_id_to_cell_id ) {
            map_point_to_cell_id[mesh_points[point_id]] = cell_id;
        }

        // Visitor visitor;
        // PMP::orient_polygon_soup( mesh_points, mesh_polygons, CGAL::parameters::visitor( visitor ) );
        // PMP::polygon_soup_to_polygon_mesh( mesh_points, mesh_polygons, m );
        auto [v_voxel_handle, v_voxel_handle_created] = m.add_property_map<Vertex_descriptor, size_t>( "v:voxel_handle" );
        LOG_ASSERT( v_voxel_handle_created );
        for( const auto& v: m.vertices() ) {
            const auto& p = m.point( v );
            if( map_point_to_cell_id.count( p ) == 1 ) {
                v_voxel_handle[v] = map_point_to_cell_id[p];
            } else {
                v_voxel_handle[v] = -1;
            }
        }
    }

    inline std::vector<std::size_t> find_non_manifold_cells( const Octree_wrapper& octree, const SurfaceMesh& primal_mesh, const FT& iso_value ) {
        std::vector<std::size_t> nm_cells;

        //// compute dc-vertices
        // std::map<Octree_wrapper::Edge_handle, std::array<size_t, 4>> quads;
        //
        // const std::size_t n_edges = octree.n_edges();
        //
        //// save all quads
        // std::set<std::size_t> active_cells;
        // for( int i = 0; i < n_edges; ++i ) {
        //     const auto& e        = octree.edges( i );
        //     const auto& [s0, s1] = octree.edge_values( e );
        //
        //     if( s0 == std::numeric_limits<FT>::max() || s1 == std::numeric_limits<FT>::max() ) {
        //         continue;
        //     }
        //
        //     if( s0 <= iso_value && s1 > iso_value ) {
        //         const auto voxels = octree.voxels_incident_to_edge( e );
        //         quads[e]          = voxels;
        //         active_cells.insert( voxels.begin(), voxels.end() );
        //     } else if( s1 <= iso_value && s0 > iso_value ) {
        //         const auto voxels = octree.voxels_incident_to_edge( e );
        //         quads[e]          = { voxels[0], voxels[3], voxels[2], voxels[1] };
        //         active_cells.insert( voxels.begin(), voxels.end() );
        //     }
        // }
        //
        //// add all cells to the list that do not intersect the iso value
        // Dual_contouring_3::Positioning::Voxel_center positioning;
        // for( const std::size_t& cell: active_cells ) {
        //     Point p;
        //     if( !positioning.position( octree, iso_value, cell, p ) ) {
        //         nm_cells.push_back( cell );
        //     }
        // }
        // if( nm_cells.size() > 0 ) {
        //     LOG( INFO ) << "Add " << nm_cells.size() << " empty cells for subdivision";
        // }
        //
        // std::vector<std::vector<size_t>> mesh_polygons;
        // mesh_polygons.reserve( quads.size() );
        //
        // for( const auto& [id, q]: quads ) {
        //     std::vector<size_t> vertices;
        //     std::set<size_t> cell_set;    // make sure that vertices only appear once
        //     for( const auto& cell_id: q ) {
        //         if( cell_set.find( cell_id ) == cell_set.end() ) {
        //             vertices.push_back( cell_id );
        //             cell_set.insert( cell_id );
        //         }
        //     }
        //     mesh_polygons.push_back( vertices );
        // }

        //////////////////////////////////
        //// find non-manifold vertices //
        // std::map<size_t, std::vector<size_t>> vertex_polygons;
        // for( size_t q_it = 0; q_it != mesh_polygons.size(); ++q_it ) {
        //     for( const auto& v: mesh_polygons[q_it] ) {
        //         vertex_polygons[v].push_back( q_it );
        //     }
        // }
        //
        // auto polygon_next_vertex_id = []( const size_t& vid, const std::vector<size_t>& quad ) {
        //     for( size_t i = 0; i < quad.size(); ++i ) {
        //         if( quad[i] == vid ) {
        //             return quad[( i + 1 ) % quad.size()];
        //         }
        //     }
        //     LOG( WARNING ) << "vertex id not found in quad";
        //     return std::numeric_limits<size_t>::max();
        // };
        //
        // auto polygon_prev_vertex_id = []( const size_t& vid, const std::vector<size_t>& quad ) {
        //     for( size_t i = 0; i < quad.size(); ++i ) {
        //         if( quad[( i + 1 ) % quad.size()] == vid ) {
        //             return quad[i];
        //         }
        //     }
        //     LOG( WARNING ) << "vertex id not found in quad";
        //     return std::numeric_limits<size_t>::max();
        // };
        //
        // auto polygon_has_he = []( const size_t& source, const size_t& target, const std::vector<size_t>& quad ) {
        //     for( size_t i = 0; i < quad.size(); ++i ) {
        //         if( quad[i] == source && quad[( i + 1 ) % quad.size()] == target ) {
        //             return true;
        //         }
        //     }
        //     return false;
        // };
        //
        // for( const auto& [point_id, polygons]: vertex_polygons ) {
        //     std::vector<size_t> polygon_ids = polygons;
        //     std::vector<size_t> umbrella_ids;    // oriented umbrella around vertex
        //     umbrella_ids.push_back( polygon_ids[0] );
        //     polygon_ids.erase( polygon_ids.begin() );
        //     // find the rest of the umbrella
        //     while( true ) {
        //         // ingoing halfedge is (vid,vo)
        //         const auto& curr_quad = mesh_polygons[umbrella_ids[umbrella_ids.size() - 1]];
        //         size_t vo             = polygon_prev_vertex_id( point_id, curr_quad );
        //
        //         // search for next polygon in umbrella
        //         size_t next_polygon_id = std::numeric_limits<size_t>::max();
        //         for( size_t i = 0; i < umbrella_ids.size() - 1; ++i ) {
        //             const auto& next_quad = mesh_polygons[umbrella_ids[i]];
        //             if( polygon_has_he( point_id, vo, next_quad ) ) {
        //                 next_polygon_id = i;
        //                 break;
        //             }
        //         }
        //
        //         // umbrella is closed --> break;
        //         if( next_polygon_id != std::numeric_limits<size_t>::max() ) {
        //             if( !polygon_ids.empty() ) {
        //                 // if next_polygon_id is non-zero, it means that some ids actually belong to the other fan
        //                 while( next_polygon_id != 0 ) {
        //                     polygon_ids.push_back( umbrella_ids[0] );
        //                     umbrella_ids.erase( umbrella_ids.begin() );
        //                     next_polygon_id--;
        //                 }
        //
        //                 // vertex is non-manifold --> duplicate it for all remaining polygon_ids
        //                 nm_cells.push_back( point_id );
        //                 // LOG( INFO ) << "Non-manifold vertex in cell " << point_id;
        //             }
        //             break;
        //         }
        //
        //         // find next quad
        //         for( size_t i = 0; i < polygon_ids.size(); ++i ) {
        //             const auto& next_quad = mesh_polygons[polygon_ids[i]];
        //             if( polygon_has_he( point_id, vo, next_quad ) ) {
        //                 next_polygon_id = i;
        //                 break;
        //             }
        //         }
        //
        //         // could not find next polygon
        //         if( next_polygon_id == std::numeric_limits<size_t>::max() ) {
        //             LOG( WARNING ) << "Could not close umbrella for cell " << point_id;
        //             break;
        //         }
        //
        //         umbrella_ids.push_back( polygon_ids[next_polygon_id] );
        //         polygon_ids.erase( polygon_ids.begin() + next_polygon_id );
        //     }
        // }

        Dual_contouring_3::Positioning::QEM_SVD_hermite<true> positioning( primal_mesh );

        // compute dc-vertices
        std::map<size_t, size_t> map_voxel_to_point_id;
        std::map<size_t, Point> map_voxel_to_point;
        size_t points_counter = 0;
        std::map<Octree_wrapper::Edge_handle, std::array<size_t, 4>> quads;

        const std::size_t n_voxel = octree.n_voxels();

        // save all points
        for( int i = 0; i < n_voxel; ++i ) {
            const auto& vh = octree.voxels( i );
            Point p;
            int position_result = positioning.position( octree, primal_mesh, iso_value, vh, p );
            if( position_result > 0 ) {
                map_voxel_to_point[vh]    = p;
                map_voxel_to_point_id[vh] = points_counter++;
                if( position_result == 1 ) {
                    nm_cells.push_back( vh );
                }
            }
        }

        const std::size_t n_edges = octree.n_edges();

        // save all quads
        for( int i = 0; i < n_edges; ++i ) {
            const auto& e        = octree.edges( i );
            const auto& [s0, s1] = octree.edge_values( e );

            if( s0 == std::numeric_limits<FT>::max() || s1 == std::numeric_limits<FT>::max() ) {
                continue;
            }

            if( s0 <= iso_value && s1 > iso_value ) {
                const auto voxels = octree.voxels_incident_to_edge( e );
                quads[e]          = voxels;
            } else if( s1 <= iso_value && s0 > iso_value ) {
                const auto voxels = octree.voxels_incident_to_edge( e );
                quads[e]          = { voxels[0], voxels[3], voxels[2], voxels[1] };
            }
        }

        // add points to empty voxels
        size_t n_empty_voxels = 0;
        std::map<Octree_wrapper::Edge_handle, std::array<size_t, 4>> quads_buf;
        for( const auto& [id, q]: quads ) {
            bool ignore_quad = false;

            // ignore quads that only have empty voxels
            bool has_non_empty_voxels = false;
            for( const auto& cell_id: q ) {
                if( !octree.voxel_prims( cell_id ).empty() ) {
                    has_non_empty_voxels = true;
                    break;
                }
            }

            if( !has_non_empty_voxels ) {
                ignore_quad = true;
            } else {
                for( const auto& cell_id: q ) {
                    if( map_voxel_to_point_id.count( cell_id ) == 0 ) {
                        // LOG( WARNING ) << "Voxel " << cell_id << " does not have a point. Add point.";
                        ++n_empty_voxels;
                        // voxels_without_prims.push_back( cell_id );
                        const auto points              = octree.voxel_vertex_positions( cell_id );
                        const Point_3 p                = points[0] + 0.5 * ( points[7] - points[0] );    // set point to voxel center;
                        map_voxel_to_point[cell_id]    = p;
                        map_voxel_to_point_id[cell_id] = points_counter++;
                        nm_cells.push_back( cell_id );
                    }
                }
            }

            if( !ignore_quad ) {
                quads_buf[id] = q;
            }
        }
        quads = quads_buf;
        // if( n_empty_voxels > 0 ) {
        //     LOG( WARNING ) << n_empty_voxels << " cells were empty. A point was added to them.";
        // }

        std::vector<Point> mesh_points;
        mesh_points.reserve( map_voxel_to_point.size() );
        std::map<size_t, size_t> map_cell_id_to_mesh_point;
        std::map<size_t, size_t> map_point_id_to_cell_id;
        std::vector<std::vector<size_t>> mesh_polygons;
        mesh_polygons.reserve( quads.size() );

        for( const auto& [id, p]: map_voxel_to_point ) {
            map_cell_id_to_mesh_point[id]               = mesh_points.size();
            map_point_id_to_cell_id[mesh_points.size()] = id;
            mesh_points.push_back( p );
        }
        for( const auto& [id, q]: quads ) {
            std::vector<size_t> vertices;
            std::set<size_t> cell_set;    // make sure that vertices only appear once
            for( const auto& cell_id: q ) {
                if( map_voxel_to_point_id.count( cell_id ) == 0 ) {
                    LOG( WARNING ) << "Voxel " << cell_id << " does not have a point. Ignoring quad.";
                    continue;
                }
                if( cell_set.find( cell_id ) == cell_set.end() ) {
                    vertices.push_back( map_cell_id_to_mesh_point[cell_id] );
                    cell_set.insert( cell_id );
                }
            }
            mesh_polygons.push_back( vertices );
        }

        // sanity check - mesh_polygons only hold vertices that do exist
        {
            std::vector<size_t> point_used( mesh_points.size(), false );
            for( const auto& poly: mesh_polygons ) {
                for( const auto& v: poly ) {
                    if( v >= mesh_points.size() ) {
                        LOG( ERROR ) << "Vertex v=" << v << " must not exist. Sanity check failed!";
                    }
                    point_used[v] = true;
                }
            }

            for( int i = 0; i < point_used.size(); ++i ) {
                if( !point_used[i] ) {
                    LOG( ERROR ) << "Vertex v=" << i << "in cell " << map_point_id_to_cell_id[i] << " was never used. Sanity check failed!";
                }
            }
        }

        std::set<size_t> nmvs = Dual_contouring_3::find_non_manifoldness_on_polygon_soup( mesh_points, mesh_polygons );

        for( const size_t v_id: nmvs ) {
            nm_cells.push_back( map_point_id_to_cell_id[v_id] );
        }

        std::sort( nm_cells.begin(), nm_cells.end() );
        nm_cells.erase( std::unique( nm_cells.begin(), nm_cells.end() ), nm_cells.end() );

        return nm_cells;
    }

}    // namespace CGAL

#endif    // CGAL_DUAL_CONTOURING_3_H
