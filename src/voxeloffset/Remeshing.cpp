#include "Remeshing.h"

#include "Quality_measurements.h"

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

// get the normal deviation only in the given direction
FT oriented_normal_deviation( const Vector& direction, const std::vector<Vector> normals ) {
    FT min_angle = FT_MAX;
    FT max_angle = -FT_MAX;
    for( const Vector& n: normals ) {
        // compute the angle between the normal and the plane described by 'direction'
        const FT dot   = CGAL::scalar_product( n, direction );
        const FT angle = std::asin( dot ) * 180.0 / CGAL_PI;
        min_angle      = CGAL::min( angle, min_angle );
        max_angle      = CGAL::max( angle, max_angle );
    }
    return max_angle - min_angle;
}

bool print_mesh( const std::filesystem::path& filename, SurfaceMesh& m ) {
    bool is_ok = CGAL::IO::write_OFF( filename.string(), m );
    return is_ok;
}

inline bool is_flip_topologically_allowed( const Edge_descriptor& e, const SurfaceMesh& m ) {
    Halfedge_descriptor h = halfedge( e, m );
    return !halfedge( target( next( h, m ), m ), target( next( opposite( h, m ), m ), m ), m ).second;
}

void Remeshing::remesh( const size_t n_iterations ) {
    // PMP::keep_largest_connected_components( mesh_, 1 );
    // mesh_.collect_garbage();

    auto [f_normals, f_normals_added] = mesh_.add_property_map<Face_descriptor, Vector>( "f:normal" );
    LOG_ASSERT( f_normals_added );
    PMP::compute_face_normals( mesh_, f_normals );

    auto [f_samples, f_samples_added] = mesh_.add_property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
    LOG_ASSERT( f_samples_added );

    auto [f_last_operation, f_last_operation_added] = mesh_.add_property_map<Face_descriptor, size_t>( "f:last_operation" );
    LOG_ASSERT( f_last_operation_added );

    for( n_remeshing_iterations_ = 0; n_remeshing_iterations_ < n_iterations; ++n_remeshing_iterations_ ) {
        LOG( INFO ) << "Remeshing " << ( n_remeshing_iterations_ + 1 ) << " / " << n_iterations;
        if( std::filesystem::exists( output_path_ ) ) {
            print_mesh( output_path_ / ( "pre_remesh_" + std::to_string( n_remeshing_iterations_ ) + ".off" ), mesh_ );
        }
        size_t n_splits = split();
        LOG( INFO ) << "# splits = " << n_splits;
        if( std::filesystem::exists( output_path_ ) ) {
            print_mesh( output_path_ / ( "split_" + std::to_string( n_remeshing_iterations_ ) + ".off" ), mesh_ );
        }
        // size_t n_collapse = collapse();
        size_t n_collapse = collapse_halfedge();
        LOG( INFO ) << "# collapses = " << n_collapse;
        // size_t n_collapse_low = collapse_low_quality_triangles();
        // LOG( INFO ) << "# collapses low = " << n_collapse_low;
        if( std::filesystem::exists( output_path_ ) ) {
            mesh_.collect_garbage();
            print_mesh( output_path_ / ( "collapsed_" + std::to_string( n_remeshing_iterations_ ) + ".off" ), mesh_ );
        }
        size_t n_flips = flip_delaunay();
        LOG( INFO ) << "# flips = " << n_flips;
        if( std::filesystem::exists( output_path_ ) ) {
            print_mesh( output_path_ / ( "flipped_" + std::to_string( n_remeshing_iterations_ ) + ".off" ), mesh_ );
        }
        smooth();
        if( std::filesystem::exists( output_path_ ) ) {
            print_mesh( output_path_ / ( "smoothed_" + std::to_string( n_remeshing_iterations_ ) + ".off" ), mesh_ );
        }
        if( iso_value_ > 0 ) {
            size_t n_collapse_untangle = collapse_halfedge_untangle();
            LOG( INFO ) << "# collapses untangling = " << n_collapse_untangle;
            if( std::filesystem::exists( output_path_ ) ) {
                mesh_.collect_garbage();
                print_mesh( output_path_ / ( "collapsed_untangling_" + std::to_string( n_remeshing_iterations_ ) + ".off" ), mesh_ );
            }
            smooth_convex_sharp_edges( 100, 10 );
            if( std::filesystem::exists( output_path_ ) ) {
                print_mesh( output_path_ / ( "smoothed_convex_" + std::to_string( n_remeshing_iterations_ ) + ".off" ), mesh_ );
            }
        }

        bool was_operation_performed = false;
        for( const auto& f: mesh_.faces() ) {
            if( f_last_operation[f] == n_remeshing_iterations_ ) {
                was_operation_performed = true;
                break;
            }
        }
        if( !was_operation_performed ) {
            LOG( INFO ) << "no more operations were performed in iteration " + std::to_string( n_remeshing_iterations_ + 1 ) + " --> break";
            break;
        }
    }

    if( iso_value_ > 0 ) {
        smooth_convex_sharp_edges( 100, 10 );
    }

    PMP::keep_large_connected_components( mesh_, 5 );

    mesh_.collect_garbage();
    // print_mesh( "../../../remeshed.off", mesh_ );

    mesh_.remove_property_map( f_normals );
    mesh_.remove_property_map( f_samples );
    mesh_.remove_property_map( f_last_operation );
}

size_t Remeshing::split() {
    auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
    LOG_ASSERT( f_samples_found );

    // LOG( INFO ) << "split all";
    typedef boost::bimap<boost::bimaps::set_of<Halfedge_descriptor>, boost::bimaps::multiset_of<FT, std::greater<FT>>> Boost_bimap;
    typedef typename Boost_bimap::value_type long_edge;

    // collect edges with large normal deviation in incident face
    Boost_bimap long_edges;
    if( iso_value_ > 0 ) {
        for( const auto& f: mesh_.faces() ) {
            if( !is_operation_required( mesh_.halfedge( f ), Remeshing_operation_types::split ) ) {
                continue;
            }
            const FT nd = normal_deviation_sampled( f );
            if( nd < max_normal_deviation_ ) {
                continue;
            }

            if( nd > 150 ) {
                // do not split faces that are tangled
                continue;
            }

            // add the longest edge
            Edge_descriptor e;
            FT l = 0;
            for( const auto& h: mesh_.halfedges_around_face( mesh_.halfedge( f ) ) ) {
                const Point p0 = mesh_.point( mesh_.source( h ) );
                const Point p1 = mesh_.point( mesh_.target( h ) );
                const FT lh    = ( p1 - p0 ).squared_length();
                if( lh > l ) {
                    e = mesh_.edge( h );
                    l = lh;
                }
            }

            // const auto& sample_normals  = f_samples[f].normals;
            // const Halfedge_descriptor h = mesh_.halfedge( f );
            // const FT nd0                = CGAL::approximate_angle( sample_normals[0], sample_normals[1] );
            // const FT nd1                = CGAL::approximate_angle( sample_normals[1], sample_normals[2] );
            // const FT nd2                = CGAL::approximate_angle( sample_normals[2], sample_normals[0] );
            //
            // if( nd0 > nd1 && nd0 > nd2 ) {
            //     const Halfedge_descriptor he = h;
            //     e                            = mesh_.edge( he );
            //     const Point p0               = mesh_.point( mesh_.source( he ) );
            //     const Point p1               = mesh_.point( mesh_.target( he ) );
            //     l                            = ( p1 - p0 ).squared_length();
            // } else if( nd1 > nd0 && nd1 > nd2 ) {
            //     const Halfedge_descriptor he = mesh_.next( h );
            //     e                            = mesh_.edge( he );
            //     const Point p0               = mesh_.point( mesh_.source( he ) );
            //     const Point p1               = mesh_.point( mesh_.target( he ) );
            //     l                            = ( p1 - p0 ).squared_length();
            // } else {
            //     const Halfedge_descriptor he = mesh_.prev( h );
            //     e                            = mesh_.edge( he );
            //     const Point p0               = mesh_.point( mesh_.source( he ) );
            //     const Point p1               = mesh_.point( mesh_.target( he ) );
            //     l                            = ( p1 - p0 ).squared_length();
            // }

            // Do not add it, if it is smaller than lMin_
            if( l < lMin_ * lMin_ ) {
                continue;
            }

            long_edges.insert( long_edge( halfedge( e, mesh_ ), nd ) );
        }
    } else {
        for( const auto& e: mesh_.edges() ) {
            const Halfedge_descriptor h = mesh_.halfedge( e );
            const Point p0              = mesh_.point( mesh_.source( h ) );
            const Point p1              = mesh_.point( mesh_.target( h ) );
            const FT l                  = ( p1 - p0 ).squared_length();

            if( l > lMax_ * lMax_ ) {
                long_edges.insert( long_edge( halfedge( e, mesh_ ), l ) );
            }
        }
    }

    // split long edges
    size_t nb_splits = 0;
    while( !long_edges.empty() ) {
        // the edge with longest length
        typename Boost_bimap::right_map::iterator eit = long_edges.right.begin();
        const Halfedge_descriptor he                  = eit->second;
        const FT normal_deviation                     = eit->first;
        long_edges.right.erase( eit );

        const Vertex_descriptor& v0 = mesh_.source( he );
        const Vertex_descriptor& v1 = mesh_.target( he );
        const Point& p0             = mesh_.point( v0 );
        const Point& p1             = mesh_.point( v1 );
        const Point p_mid           = CGAL::midpoint( p0, p1 );

        Point p_mid_proj;
        {
            std::vector<Point> p_proj;
            std::vector<Vector> n_proj;
            std::vector<FT> w;

            const Face_descriptor& f0 = mesh_.face( he );
            const Face_descriptor& f1 = mesh_.face( mesh_.opposite( he ) );

            std::vector<Face_descriptor> fv0 { f0 };
            const FT area0 = PMP::area( fv0, mesh_ );
            std::vector<Face_descriptor> fv1 { f1 };
            const FT area1 = PMP::area( fv1, mesh_ );

            // only use the two sample points nearest to the split point
            // FT dist_s0_mid = FT_MAX;
            // Point s0;
            // Vector sn0;
            // for( int i = 0; i < f_samples[f0].points.size(); ++i ) {
            //    const Point& sp = f_samples[f0].points[i];
            //    const FT d      = ( sp - p_mid ).squared_length();
            //    if( d < FT_MAX ) {
            //        dist_s0_mid = d;
            //        s0          = sp;
            //        sn0         = f_samples[f0].normals[i];
            //    }
            //}
            //
            // FT dist_s1_mid = FT_MAX;
            // Point s1;
            // Vector sn1;
            // for( int i = 0; i < f_samples[f1].points.size(); ++i ) {
            //    const Point& sp = f_samples[f1].points[i];
            //    const FT d      = ( sp - p_mid ).squared_length();
            //    if( d < FT_MAX ) {
            //        dist_s1_mid = d;
            //        s1          = sp;
            //        sn1         = f_samples[f1].normals[i];
            //    }
            //}
            //
            // p_proj.push_back( s0 );
            // n_proj.push_back( sn0 );
            // w.push_back( 1 );
            //
            // p_proj.push_back( s1 );
            // n_proj.push_back( sn1 );
            // w.push_back( 1 );

            for( const auto& sp: f_samples[f0].points ) {
                p_proj.push_back( sp );
            }
            for( const auto& sn: f_samples[f0].normals ) {
                n_proj.push_back( sn );
                w.push_back( area0 );
            }
            for( const auto& sp: f_samples[f1].points ) {
                p_proj.push_back( sp );
            }
            for( const auto& sn: f_samples[f1].normals ) {
                n_proj.push_back( sn );
                w.push_back( area1 );
            }
            p_mid_proj = qem_weighted( p_mid, p_proj, n_proj, w );
        }

        {
            // check for quality after split
            FT current_length = ( p1 - p0 ).squared_length();
            if( ( p_mid_proj - p0 ).squared_length() > current_length || ( p_mid_proj - p1 ).squared_length() > current_length ) {
                // new vertex position is far away --> only move it in the direction by "current_length"
                Vector dir = normalize( p_mid_proj - p_mid );
                p_mid_proj = p_mid + CGAL::sqrt( current_length ) * dir;
            }
        }

        // split edge
        Halfedge_descriptor hnew     = CGAL::Euler::split_edge( he, mesh_ );
        Vertex_descriptor vNew       = mesh_.target( hnew );
        Halfedge_descriptor hnew_opp = opposite( hnew, mesh_ );
        ++nb_splits;

        mesh_.point( vNew ) = p_mid_proj;

        // insert new edges to keep triangular faces
        if( !mesh_.is_border( hnew ) ) {
            CGAL::Euler::split_face( hnew, next( next( hnew, mesh_ ), mesh_ ), mesh_ );
        }
        if( !mesh_.is_border( hnew_opp ) ) {
            CGAL::Euler::split_face( prev( hnew_opp, mesh_ ), next( hnew_opp, mesh_ ), mesh_ );
        }

        // update normals and samples
        for( const auto& f: mesh_.faces_around_target( mesh_.halfedge( vNew ) ) ) {
            update_triangle_normal( f );
            update_triangle_sampling( f );
            update_last_operation( f );
        }
    }

    return nb_splits;
}

size_t Remeshing::collapse() {
    // LOG( INFO ) << "collapse all";

    auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
    LOG_ASSERT( f_samples_found );

    typedef boost::bimap<boost::bimaps::set_of<Halfedge_descriptor>, boost::bimaps::multiset_of<double, std::less<double>>> Boost_bimap;
    typedef typename Boost_bimap::value_type short_edge;

    auto [f_normal_deviation, f_normal_deviation_created] = mesh_.add_property_map<Face_descriptor, FT>( "f:normal_deviation" );
    LOG_ASSERT( f_normal_deviation_created );
    for( const auto& f: mesh_.faces() ) {
        f_normal_deviation[f] = normal_deviation_sampled( f );
    }

    Boost_bimap short_edges;
    for( const auto& e: mesh_.edges() ) {
        if( !is_operation_required( mesh_.halfedge( e ), Remeshing_operation_types::collapse ) ) {
            continue;
        }
        if( mesh_.is_border( e ) ) {
            continue;
        }
        // add edges that have a small normal deviation and are below lMin_
        const FT sqlen = PMP::squared_edge_length( e, mesh_ );
        // if( sqlen > lMin_ * lMin_ ) {
        //    continue;
        //}
        const Halfedge_descriptor h = mesh_.halfedge( e );
        const FT fnd1               = f_normal_deviation[mesh_.face( h )];
        const FT fnd2               = f_normal_deviation[mesh_.face( mesh_.opposite( h ) )];
        const FT fnd                = std::max( fnd1, fnd2 );
        if( fnd > min_normal_deviation_ && sqlen > lMin_ * lMin_ ) {
            continue;
        }
        short_edges.insert( short_edge( mesh_.halfedge( e ), sqlen ) );
    }

    mesh_.remove_property_map( f_normal_deviation );

    size_t nb_collapses = 0;
    while( !short_edges.empty() ) {
        // the edge with shortest length
        typename Boost_bimap::right_map::iterator eit = short_edges.right.begin();
        const Halfedge_descriptor he                  = eit->second;
        const Halfedge_descriptor he_opp              = mesh_.opposite( he );
        short_edges.right.erase( eit );

        const Edge_descriptor e     = edge( he, mesh_ );
        const Edge_descriptor e_opp = edge( he_opp, mesh_ );

        if( !CGAL::Euler::does_satisfy_link_condition( e, mesh_ ) ) {
            continue;
        }

        const Edge_descriptor e_collapse     = e;
        const Halfedge_descriptor h_collapse = mesh_.halfedge( e_collapse );
        const Vertex_descriptor v0           = source( h_collapse, mesh_ );
        const Vertex_descriptor v1           = target( h_collapse, mesh_ );

        // do not collapse if edge becomes too long
        bool collapse_ok = true;
        for( const auto& h: mesh_.halfedges_around_target( mesh_.halfedge( v0 ) ) ) {
            const Vertex_descriptor v = mesh_.source( h );
            if( ( mesh_.point( v ) - mesh_.point( v1 ) ).squared_length() > lMax_ * lMax_ ) {
                collapse_ok = false;
                break;
            }
        }
        if( !collapse_ok ) {
            continue;
        }

        Point p_optimal;
        {
            // compute optimal point position with qem
            std::set<Face_descriptor> neighs;
            for( const auto& f: mesh_.faces_around_target( h_collapse ) ) {
                neighs.insert( f );
            }
            for( const auto& f: mesh_.faces_around_target( mesh_.opposite( h_collapse ) ) ) {
                neighs.insert( f );
            }
            std::vector<Point> p_proj;
            std::vector<Vector> n_proj;
            std::vector<FT> w;
            for( const auto& f: neighs ) {
                std::vector<Face_descriptor> fv { f };
                const FT area  = PMP::area( fv, mesh_ );
                const auto& fs = f_samples[f];
                for( const auto& sp: fs.points ) {
                    p_proj.push_back( sp );
                }
                for( const auto& sn: fs.normals ) {
                    n_proj.push_back( sn );
                    w.push_back( area );
                }
            }
            p_optimal = qem_weighted( CGAL::midpoint( mesh_.point( v0 ), mesh_.point( v1 ) ), p_proj, n_proj, w );
        }

        // do not collapse if elements go below min_quality_
        // const Point p_current = mesh_.point( v0 );
        // mesh_.point( v0 )     = mesh_.point( v1 );
        const Point p0_current = mesh_.point( v0 );
        const Point p1_current = mesh_.point( v1 );
        mesh_.point( v0 )      = p_optimal;
        mesh_.point( v1 )      = p_optimal;
        for( const auto& h: mesh_.halfedges_around_target( mesh_.halfedge( v0 ) ) ) {
            if( mesh_.source( h ) == v1 || mesh_.target( mesh_.next( h ) ) == v1 ) {
                continue;
            }
            const FT normal_deviation = face_normal_deviation( mesh_.face( h ) );
            if( normal_deviation > max_normal_deviation_ ) {
                collapse_ok = false;
                break;
            }
        }
        for( const auto& h: mesh_.halfedges_around_target( mesh_.halfedge( v1 ) ) ) {
            if( mesh_.source( h ) == v0 || mesh_.target( mesh_.next( h ) ) == v0 ) {
                continue;
            }
            const FT normal_deviation = face_normal_deviation( mesh_.face( h ) );
            if( normal_deviation > max_normal_deviation_ ) {
                collapse_ok = false;
                break;
            }
        }
        // mesh_.point( v0 ) = p_current;
        if( !collapse_ok ) {
            mesh_.point( v0 ) = p0_current;
            mesh_.point( v1 ) = p1_current;
            continue;
        }

        //"collapse v0 into v1 along e"
        // remove edges incident to va and vb, because their lengths will change
        for( Halfedge_descriptor ha: halfedges_around_target( v0, mesh_ ) ) {
            short_edges.left.erase( ha );
            short_edges.left.erase( opposite( ha, mesh_ ) );
        }
        for( Halfedge_descriptor hb: halfedges_around_target( v1, mesh_ ) ) {
            short_edges.left.erase( hb );
            short_edges.left.erase( opposite( hb, mesh_ ) );
        }

        Vertex_descriptor vkept = CGAL::Euler::collapse_edge( e_collapse, mesh_ );

        ++nb_collapses;

        // update normals and samples
        for( Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
            update_triangle_normal( mesh_.face( ht ) );
            update_triangle_sampling( mesh_.face( ht ) );
            update_last_operation( mesh_.face( ht ) );
        }

        // insert new/remaining short edges
        for( Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
            if( mesh_.is_border( mesh_.edge( ht ) ) ) {
                continue;
            }
            double sqlen = PMP::squared_edge_length( mesh_.edge( ht ), mesh_ );
            // if( sqlen > lMin_ * lMin_ ) {
            //    continue;
            //}

            // const FT fnd1 = face_normal_deviation( mesh_.face( ht ) );
            // const FT fnd2 = face_normal_deviation( mesh_.face( mesh_.opposite( ht ) ) );
            const FT fnd1 = normal_deviation_sampled( mesh_.face( ht ) );
            const FT fnd2 = normal_deviation_sampled( mesh_.face( mesh_.opposite( ht ) ) );
            const FT fnd  = std::max( fnd1, fnd2 );
            if( fnd > min_normal_deviation_ && sqlen > lMin_ * lMin_ ) {
                continue;
            }
            short_edges.insert( short_edge( ht, sqlen ) );
        }
    }

    return nb_collapses;
}

size_t Remeshing::collapse_halfedge() {
    // LOG( INFO ) << "collapse all";

    auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
    LOG_ASSERT( f_samples_found );

    typedef boost::bimap<boost::bimaps::set_of<Halfedge_descriptor>, boost::bimaps::multiset_of<double, std::less<double>>> Boost_bimap;
    typedef typename Boost_bimap::value_type short_edge;

    Boost_bimap short_edges;
    if( iso_value_ > 0 ) {
        for( const auto& e: mesh_.edges() ) {
            if( !is_operation_required( mesh_.halfedge( e ), Remeshing_operation_types::collapse ) ) {
                continue;
            }
            if( mesh_.is_border( e ) ) {
                continue;
            }
            // add halfedges that have a small normal deviation and are shorter than lMax_
            const FT sqlen = PMP::squared_edge_length( e, mesh_ );
            if( sqlen > lMax_ * lMax_ ) {
                continue;
            }
            const Halfedge_descriptor h = mesh_.halfedge( e );
            const Vertex_descriptor v0  = mesh_.source( h );
            const Vertex_descriptor v1  = mesh_.target( h );
            const Point& p0             = mesh_.point( v0 );
            const Point& p1             = mesh_.point( v1 );
            const Vector ev             = normalize( p1 - p0 );

            std::vector<Vector> normals0;
            for( const Face_descriptor f: mesh_.faces_around_target( mesh_.halfedge( v0 ) ) ) {
                normals0.insert( normals0.end(), f_samples[f].normals.begin(), f_samples[f].normals.end() );
            }
            FT angle_diff0 = oriented_normal_deviation( ev, normals0 );

            std::vector<Vector> normals1;
            for( const Face_descriptor f: mesh_.faces_around_target( mesh_.halfedge( v1 ) ) ) {
                normals1.insert( normals1.end(), f_samples[f].normals.begin(), f_samples[f].normals.end() );
            }
            FT angle_diff1 = oriented_normal_deviation( ev, normals1 );

            if( angle_diff0 <= angle_diff1 && angle_diff0 < max_normal_deviation_ ) {
                short_edges.insert( short_edge( h, angle_diff0 ) );
            }

            if( angle_diff1 < angle_diff0 && angle_diff1 < max_normal_deviation_ ) {
                short_edges.insert( short_edge( mesh_.opposite( h ), angle_diff1 ) );
            }
        }
    } else {
        for( const auto& e: mesh_.edges() ) {
            const Halfedge_descriptor h = mesh_.halfedge( e );
            const Point p0              = mesh_.point( mesh_.source( h ) );
            const Point p1              = mesh_.point( mesh_.target( h ) );
            const FT l                  = ( p1 - p0 ).squared_length();

            if( l < lMin_ * lMin_ ) {
                short_edges.insert( short_edge( halfedge( e, mesh_ ), l ) );
            }
        }
    }

    size_t nb_collapses = 0;
    while( !short_edges.empty() ) {
        // the edge with shortest length
        typename Boost_bimap::right_map::iterator eit = short_edges.right.begin();
        const Halfedge_descriptor he                  = eit->second;
        const Halfedge_descriptor he_opp              = mesh_.opposite( he );
        short_edges.right.erase( eit );

        const Edge_descriptor e     = edge( he, mesh_ );
        const Edge_descriptor e_opp = edge( he_opp, mesh_ );

        if( !CGAL::Euler::does_satisfy_link_condition( e, mesh_ ) ) {
            continue;
        }

        const Edge_descriptor e_collapse     = e;
        const Halfedge_descriptor h_collapse = he;
        const Vertex_descriptor v0           = source( h_collapse, mesh_ );
        const Vertex_descriptor v1           = target( h_collapse, mesh_ );

        // do not collapse if edge becomes too long
        bool collapse_ok = true;
        for( const auto& h: mesh_.halfedges_around_target( mesh_.halfedge( v0 ) ) ) {
            const Vertex_descriptor v = mesh_.source( h );
            if( ( mesh_.point( v ) - mesh_.point( v1 ) ).squared_length() > lMax_ * lMax_ ) {
                collapse_ok = false;
                break;
            }
        }
        if( !collapse_ok ) {
            continue;
        }

        // mesh_.point( v0 ) = mesh_.point( v1 );

        //"collapse v0 into v1 along e"
        // remove edges incident to va and vb, because their lengths will change
        for( Halfedge_descriptor ha: halfedges_around_target( v0, mesh_ ) ) {
            short_edges.left.erase( ha );
            short_edges.left.erase( opposite( ha, mesh_ ) );
        }
        for( Halfedge_descriptor hb: halfedges_around_target( v1, mesh_ ) ) {
            short_edges.left.erase( hb );
            short_edges.left.erase( opposite( hb, mesh_ ) );
        }

        Vertex_descriptor vkept = CGAL::Euler::collapse_edge( e_collapse, mesh_ );

        ++nb_collapses;

        // update normals and samples
        for( Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
            update_triangle_normal( mesh_.face( ht ) );
            update_triangle_sampling( mesh_.face( ht ) );
            update_last_operation( mesh_.face( ht ) );
        }

        // insert new/remaining short edges
        if( iso_value_ > 0 ) {
            for( const Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
                const Edge_descriptor e = mesh_.edge( ht );
                if( mesh_.is_border( e ) ) {
                    continue;
                }
                // add halfedges that have a small normal deviation and are shorter than lMax_
                const FT sqlen = PMP::squared_edge_length( e, mesh_ );
                if( sqlen > lMax_ * lMax_ ) {
                    continue;
                }
                const Halfedge_descriptor h = mesh_.halfedge( e );
                const Vertex_descriptor v0  = mesh_.source( h );
                const Vertex_descriptor v1  = mesh_.target( h );
                const Point& p0             = mesh_.point( v0 );
                const Point& p1             = mesh_.point( v1 );
                const Vector ev             = normalize( p1 - p0 );

                FT min_angle0 = FT_MAX;
                FT max_angle0 = -FT_MAX;
                for( const Face_descriptor f: mesh_.faces_around_target( mesh_.halfedge( v0 ) ) ) {
                    for( const Vector& n: f_samples[f].normals ) {
                        const FT dot   = CGAL::scalar_product( n, ev );
                        const FT angle = std::asin( dot ) * 180.0 / CGAL_PI;
                        min_angle0     = CGAL::min( angle, min_angle0 );
                        max_angle0     = CGAL::max( angle, max_angle0 );
                    }
                }
                const FT angle_diff0 = max_angle0 - min_angle0;

                FT min_angle1 = FT_MAX;
                FT max_angle1 = -FT_MAX;
                for( const Face_descriptor f: mesh_.faces_around_target( mesh_.halfedge( v1 ) ) ) {
                    for( const Vector& n: f_samples[f].normals ) {
                        const FT dot   = CGAL::scalar_product( n, ev );
                        const FT angle = std::asin( dot ) * 180.0 / CGAL_PI;
                        min_angle1     = CGAL::min( angle, min_angle1 );
                        max_angle1     = CGAL::max( angle, max_angle1 );
                    }
                }
                const FT angle_diff1 = max_angle1 - min_angle1;

                if( angle_diff0 <= angle_diff1 && angle_diff0 < max_normal_deviation_ ) {
                    short_edges.insert( short_edge( h, angle_diff0 ) );
                }

                if( angle_diff1 < angle_diff0 && angle_diff1 < max_normal_deviation_ ) {
                    short_edges.insert( short_edge( mesh_.opposite( h ), angle_diff1 ) );
                }
            }
        } else {
            for( const Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
                const Point p0               = mesh_.point( mesh_.source( ht ) );
                const Point p1               = mesh_.point( mesh_.target( ht ) );
                const FT l                   = ( p1 - p0 ).squared_length();

                if( l < lMin_ * lMin_ ) {
                    short_edges.insert( short_edge( halfedge( e, mesh_ ), l ) );
                }
            }
        }
    }

    return nb_collapses;
}

size_t Remeshing::collapse_halfedge_untangle() {
    auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
    LOG_ASSERT( f_samples_found );

    typedef boost::bimap<boost::bimaps::set_of<Halfedge_descriptor>, boost::bimaps::multiset_of<size_t, std::less<size_t>>> Boost_bimap;
    typedef typename Boost_bimap::value_type short_edge;

    auto e_is_sharp = [this]( const Halfedge_descriptor h0, const Halfedge_descriptor h1 ) {
        // f0 (v0,v1,v2)
        // f1 (v3,v4,v5)
        const Vertex_descriptor v0 = mesh_.source( h0 );
        const Vertex_descriptor v1 = mesh_.target( h0 );
        const Vertex_descriptor v2 = mesh_.target( mesh_.next( h0 ) );
        const Vertex_descriptor v3 = mesh_.source( h1 );
        const Vertex_descriptor v4 = mesh_.target( h1 );
        const Vertex_descriptor v5 = mesh_.target( mesh_.next( h1 ) );
        const Point& p0            = mesh_.point( v0 );
        const Point& p1            = mesh_.point( v1 );
        const Point& p2            = mesh_.point( v2 );
        const Point& p3            = mesh_.point( v3 );
        const Point& p4            = mesh_.point( v4 );
        const Point& p5            = mesh_.point( v5 );

        LOG_ASSERT( p0 == p4 );
        LOG_ASSERT( p1 == p3 );

        const Vector n0  = normalize( CGAL::cross_product( p1 - p0, p2 - p0 ) );
        const Vector n1  = normalize( CGAL::cross_product( p4 - p3, p5 - p3 ) );
        const Vector nxn = CGAL::cross_product( n0, n1 );
        if( CGAL::scalar_product( nxn, p1 - p0 ) < 0 ) {
            // edge is concave --> not sharp
            return false;
        }
        FT angle = CGAL::approximate_angle( n0, n1 );
        if( angle < 150 ) {
            return false;
        }
        return true;
    };

    auto untangling_potential = [this, &e_is_sharp]( const Halfedge_descriptor h ) {
        const Vertex_descriptor v0 = mesh_.source( h );

        size_t n_tangled_edges_before = 0;
        for( const Halfedge_descriptor& hh: mesh_.halfedges_around_target( mesh_.halfedge( v0 ) ) ) {
            // ingoing halfedges
            if( e_is_sharp( hh, mesh_.opposite( hh ) ) ) {
                ++n_tangled_edges_before;
            }
            // prev
            const Halfedge_descriptor hh_prev = mesh_.prev( hh );
            if( e_is_sharp( hh_prev, mesh_.opposite( hh_prev ) ) ) {
                ++n_tangled_edges_before;
            }
        }

        if( n_tangled_edges_before == 0 ) {
            return n_tangled_edges_before;
        }

        size_t n_tangled_edges_after = 0;
        const Point p0_buf           = mesh_.point( v0 );
        mesh_.point( v0 )            = mesh_.point( mesh_.target( h ) );
        const Face_descriptor f0     = mesh_.face( h );
        const Face_descriptor f1     = mesh_.face( mesh_.opposite( h ) );
        for( const Halfedge_descriptor& hh: mesh_.halfedges_around_target( mesh_.halfedge( v0 ) ) ) {
            if( mesh_.face( hh ) == f1 ) {
                // that edge is getting collapsed
                continue;
            }
            if( mesh_.face( hh ) == f0 ) {
                Halfedge_descriptor hh_corrected = mesh_.opposite( mesh_.prev( hh ) );
                if( e_is_sharp( hh_corrected, mesh_.opposite( hh ) ) ) {
                    ++n_tangled_edges_after;
                }
            } else if( mesh_.face( mesh_.opposite( hh ) ) == f1 ) {
                Halfedge_descriptor hh_opp_corrected = mesh_.opposite( mesh_.next( mesh_.opposite( hh ) ) );
                if( e_is_sharp( hh, hh_opp_corrected ) ) {
                    ++n_tangled_edges_after;
                }
            } else {
                // ingoing halfedges
                if( e_is_sharp( hh, mesh_.opposite( hh ) ) ) {
                    ++n_tangled_edges_after;
                }
                // prev
                const Halfedge_descriptor hh_prev = mesh_.prev( hh );
                if( e_is_sharp( hh_prev, mesh_.opposite( hh_prev ) ) ) {
                    ++n_tangled_edges_after;
                }
            }
        }
        mesh_.point( v0 ) = p0_buf;

        return n_tangled_edges_before - n_tangled_edges_after;
    };

    Boost_bimap short_edges;
    for( const auto& e: mesh_.edges() ) {
        if( !is_operation_required( mesh_.halfedge( e ), Remeshing_operation_types::collapse ) ) {
            continue;
        }
        if( mesh_.is_border( e ) ) {
            continue;
        }

        const Halfedge_descriptor h0 = mesh_.halfedge( e );
        const Halfedge_descriptor h1 = mesh_.opposite( h0 );

        // find halfedgeedge collapses that removes sharpness
        size_t h0_potential = untangling_potential( h0 );
        size_t h1_potential = untangling_potential( h1 );

        if( h0_potential > 0 ) {
            short_edges.insert( short_edge( h0, h0_potential ) );
        }
        if( h1_potential > 0 ) {
            short_edges.insert( short_edge( h1, h1_potential ) );
        }
    }

    size_t nb_collapses = 0;
    while( !short_edges.empty() ) {
        // the edge with shortest length
        typename Boost_bimap::right_map::iterator eit = short_edges.right.begin();
        const Halfedge_descriptor he                  = eit->second;
        const Halfedge_descriptor he_opp              = mesh_.opposite( he );
        short_edges.right.erase( eit );

        const Edge_descriptor e     = edge( he, mesh_ );
        const Edge_descriptor e_opp = edge( he_opp, mesh_ );

        if( !CGAL::Euler::does_satisfy_link_condition( e, mesh_ ) ) {
            continue;
        }

        const Edge_descriptor e_collapse     = e;
        const Halfedge_descriptor h_collapse = he;
        const Vertex_descriptor v0           = source( h_collapse, mesh_ );
        const Vertex_descriptor v1           = target( h_collapse, mesh_ );

        // do not collapse if edge becomes too long
        bool collapse_ok = true;
        for( const auto& h: mesh_.halfedges_around_target( mesh_.halfedge( v0 ) ) ) {
            const Vertex_descriptor v = mesh_.source( h );
            if( ( mesh_.point( v ) - mesh_.point( v1 ) ).squared_length() > lMax_ * lMax_ ) {
                collapse_ok = false;
                break;
            }
        }
        if( !collapse_ok ) {
            continue;
        }

        // mesh_.point( v0 ) = mesh_.point( v1 );

        //"collapse v0 into v1 along e"
        // remove edges incident to va and vb, because their lengths will change
        for( Halfedge_descriptor ha: halfedges_around_target( v0, mesh_ ) ) {
            short_edges.left.erase( ha );
            short_edges.left.erase( opposite( ha, mesh_ ) );
        }
        for( Halfedge_descriptor hb: halfedges_around_target( v1, mesh_ ) ) {
            short_edges.left.erase( hb );
            short_edges.left.erase( opposite( hb, mesh_ ) );
        }

        Vertex_descriptor vkept = CGAL::Euler::collapse_edge( e_collapse, mesh_ );

        ++nb_collapses;

        // update normals and samples
        for( const Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
            update_triangle_normal( mesh_.face( ht ) );
            update_triangle_sampling( mesh_.face( ht ) );
            update_last_operation( mesh_.face( ht ) );
        }

        // insert new/remaining short edges
        for( Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
            const Edge_descriptor e = mesh_.edge( ht );
            if( mesh_.is_border( e ) ) {
                continue;
            }

            const Halfedge_descriptor h0 = mesh_.halfedge( e );
            const Halfedge_descriptor h1 = mesh_.opposite( h0 );

            // find halfedgeedge collapses that removes sharpness
            size_t h0_potential = untangling_potential( h0 );
            size_t h1_potential = untangling_potential( h1 );

            if( h0_potential > 0 ) {
                short_edges.insert( short_edge( h0, h0_potential ) );
            }
            if( h1_potential > 0 ) {
                short_edges.insert( short_edge( h1, h1_potential ) );
            }
        }
    }

    return nb_collapses;
}

size_t Remeshing::collapse_low_quality_triangles() {
    // LOG( INFO ) << "collapse low quality triangles";

    typedef boost::bimap<boost::bimaps::set_of<Face_descriptor>, boost::bimaps::multiset_of<double, std::less<double>>> Boost_bimap;
    typedef typename Boost_bimap::value_type trash_triangle;

    // Boost_bimap short_edges;
    Boost_bimap trash_triangles;
    for( const Face_descriptor& f: mesh_.faces() ) {
        // const Halfedge_descriptor h = mesh_.halfedge( f );
        // const Vertex_descriptor v0  = mesh_.source( h );
        // const Vertex_descriptor v1  = mesh_.target( h );
        // const Vertex_descriptor v2  = mesh_.target( mesh_.next( h ) );
        // const Point p0              = mesh_.point( v0 );
        // const Point p1              = mesh_.point( v1 );
        // const Point p2              = mesh_.point( v2 );
        //
        // const Vector normal = normalize( CGAL::cross_product( p1 - p0, p2 - p0 ) );
        const FT q = mean_ratio_metric_with_sampling( f );
        if( q < min_quality_ ) {
            trash_triangles.insert( trash_triangle( f, q ) );
        }
    }

    size_t nb_collapses = 0;
    while( !trash_triangles.empty() ) {
        // the edge with shortest length
        typename Boost_bimap::right_map::iterator eit = trash_triangles.right.begin();
        const Face_descriptor f                       = eit->second;
        trash_triangles.right.erase( eit );

        // sort edges by length
        std::multimap<FT, Edge_descriptor> face_edges;
        for( const auto& h: mesh_.halfedges_around_face( mesh_.halfedge( f ) ) ) {
            const Vertex_descriptor v0 = mesh_.source( h );
            const Vertex_descriptor v1 = mesh_.target( h );
            const Point p0             = mesh_.point( v0 );
            const Point p1             = mesh_.point( v1 );
            const FT l_squared         = ( p1 - p0 ).squared_length();
            face_edges.insert( { l_squared, mesh_.edge( h ) } );
        }

        Edge_descriptor e_collapse;
        for( const auto& [l, e]: face_edges ) {
            if( CGAL::Euler::does_satisfy_link_condition( e, mesh_ ) ) {
                e_collapse = e;
                break;
            }
        }

        if( !e_collapse.is_valid() ) {
            continue;
        }

        const Halfedge_descriptor h_collapse = mesh_.halfedge( e_collapse );
        Vertex_descriptor v0                 = source( h_collapse, mesh_ );
        Vertex_descriptor v1                 = target( h_collapse, mesh_ );

        //"collapse v0 into v1 along e"
        // remove edges incident to va and vb, because their lengths will change
        for( Halfedge_descriptor ha: halfedges_around_target( v0, mesh_ ) ) {
            trash_triangles.left.erase( mesh_.face( ha ) );
        }
        for( Halfedge_descriptor hb: halfedges_around_target( v1, mesh_ ) ) {
            trash_triangles.left.erase( mesh_.face( hb ) );
        }

        Vertex_descriptor vkept = CGAL::Euler::collapse_edge( e_collapse, mesh_ );

        ++nb_collapses;

        // update normals and samples
        for( Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
            update_triangle_normal( mesh_.face( ht ) );
            update_triangle_sampling( mesh_.face( ht ) );
            update_last_operation( mesh_.face( ht ) );
        }

        // insert new/remaining short edges
        for( Halfedge_descriptor ht: halfedges_around_target( vkept, mesh_ ) ) {
            const FT q = mean_ratio_metric_with_sampling( mesh_.face( ht ) );
            if( q < min_quality_ ) {
                trash_triangles.insert( trash_triangle( mesh_.face( ht ), q ) );
            }
        }
    }

    return nb_collapses;
}

size_t Remeshing::flip_delaunay() {
    // LOG( INFO ) << "flip all to Delaunay";
    size_t n_flips     = 0;
    size_t n_flips_old = n_flips;

    // current orientation:
    //    2 - 1
    //    | / |
    //    0 - 3

    auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
    LOG_ASSERT( f_samples_found );
    auto [f_normals, f_normals_found] = mesh_.property_map<Face_descriptor, Vector>( "f:normal" );
    LOG_ASSERT( f_normals_found );

    auto [f_normal_samples, f_normal_created] = mesh_.add_property_map<Face_descriptor, Vector>( "f:normal_samples" );
    LOG_ASSERT( f_normal_created );

    for( const auto& f: mesh_.faces() ) {
        // f_normal[f] = compute_f_normal( f );
        f_normal_samples[f] = f_samples[f].normals[3];
    }

    typedef boost::bimap<boost::bimaps::set_of<Halfedge_descriptor>, boost::bimaps::multiset_of<double, std::greater<double>>> Boost_bimap;
    typedef typename Boost_bimap::value_type flip_edge;
    Boost_bimap flip_edges;

    for( const auto& e: mesh_.edges() ) {
        if( !is_operation_required( mesh_.halfedge( e ), Remeshing_operation_types::flip ) ) {
            continue;
        }
        if( mesh_.is_border( e ) ) {
            continue;
        }
        const auto& h = mesh_.halfedge( e );

        const auto& v0 = mesh_.source( h );
        const auto& v1 = mesh_.target( h );
        const auto& p0 = mesh_.point( v0 );
        const auto& p1 = mesh_.point( v1 );
        if( !is_flip_topologically_allowed( e, mesh_ ) ) {
            continue;
        }

        const auto& v2 = mesh_.target( mesh_.next( h ) );
        const auto& v3 = mesh_.target( mesh_.next( mesh_.opposite( h ) ) );
        const auto& p2 = mesh_.point( v2 );
        const auto& p3 = mesh_.point( v3 );

        const Face_descriptor f_012 = mesh_.face( h );
        const Face_descriptor f_031 = mesh_.face( mesh_.opposite( h ) );

        std::vector<Vector> normals;
        normals.insert( normals.end(), f_samples[f_012].normals.begin(), f_samples[f_012].normals.end() );
        normals.insert( normals.end(), f_samples[f_031].normals.begin(), f_samples[f_031].normals.end() );
        const Vector e01_direction = normalize( p1 - p0 );
        FT e01_nd                  = oriented_normal_deviation( e01_direction, normals );

        const Vector e23_direction = normalize( p3 - p2 );
        FT e23_nd                  = oriented_normal_deviation( e23_direction, normals );

        if( e23_nd > max_normal_deviation_ && e01_nd < max_normal_deviation_ ) {
            continue;
        }

        // check min quality
        const FT q012    = meanRatioMetric_heavy_normal( { p0, p1, p2 }, f_normal_samples[f_012] );
        const FT q031    = meanRatioMetric_heavy_normal( { p0, p3, p1 }, f_normal_samples[f_031] );
        const FT qMinOld = std::min( q012, q031 );

        const Point p_center_032 = CGAL::centroid( p0, p3, p2 );
        const Point p_center_312 = CGAL::centroid( p3, p1, p2 );
        FT qMinNew;
        if( f_samples[f_012].has_one_primitive && f_samples[f_031].has_one_primitive &&
            f_samples[f_012].primitive_id == f_samples[f_031].primitive_id ) {
            // directly use this primitive, do not use the AABB_tree
            const auto [p_center_032_proj, n_032, t_032] = get_offset_projection_and_normal( p_center_032, f_samples[f_012].primitive_id );
            const auto [p_center_312_proj, n_312, t_312] = get_offset_projection_and_normal( p_center_312, f_samples[f_012].primitive_id );
            const FT q032                                = meanRatioMetric_heavy_normal( { p0, p3, p2 }, n_032 );
            const FT q312                                = meanRatioMetric_heavy_normal( { p3, p1, p2 }, n_312 );
            qMinNew                                      = std::min( q032, q312 );
        } else {
            const auto [p_center_032_proj, n_032] = get_offset_projection_and_normal( p_center_032 );
            const auto [p_center_312_proj, n_312] = get_offset_projection_and_normal( p_center_312 );
            const FT q032                         = meanRatioMetric_heavy_normal( { p0, p3, p2 }, n_032 );
            const FT q312                         = meanRatioMetric_heavy_normal( { p3, p1, p2 }, n_312 );
            qMinNew                               = std::min( q032, q312 );
        }

        if( qMinNew <= qMinOld ) {
            continue;
        }

        flip_edges.insert( flip_edge( mesh_.halfedge( e ), qMinNew - qMinOld ) );
    }

    while( !flip_edges.empty() ) {
        // the edge with shortest length
        typename Boost_bimap::right_map::iterator eit = flip_edges.right.begin();
        const Halfedge_descriptor h                   = eit->second;
        const Halfedge_descriptor h_opp               = mesh_.opposite( h );
        flip_edges.right.erase( eit );

        const Edge_descriptor e = edge( h, mesh_ );

        if( !is_flip_topologically_allowed( e, mesh_ ) ) {
            continue;
        }

        const Face_descriptor f0 = mesh_.face( h );
        const Face_descriptor f1 = mesh_.face( mesh_.opposite( h ) );

        if( normal_deviation_sampled( f0 ) > 179 || normal_deviation_sampled( f1 ) > 179 ) {
            continue;
        }

        const Vertex_descriptor& v0 = source( h, mesh_ );
        const Vertex_descriptor& v1 = target( h, mesh_ );
        const Vertex_descriptor& v2 = mesh_.target( mesh_.next( h ) );
        const Vertex_descriptor& v3 = mesh_.target( mesh_.next( mesh_.opposite( h ) ) );

        std::vector<Halfedge_descriptor> surrounding_halfedges;
        surrounding_halfedges.push_back( mesh_.next( h ) );
        surrounding_halfedges.push_back( mesh_.prev( h ) );
        surrounding_halfedges.push_back( mesh_.next( mesh_.opposite( h ) ) );
        surrounding_halfedges.push_back( mesh_.prev( mesh_.opposite( h ) ) );
        const auto& p0 = mesh_.point( v0 );
        const auto& p1 = mesh_.point( v1 );
        const auto& p2 = mesh_.point( v2 );
        const auto& p3 = mesh_.point( v3 );

        //"flip edge from v0v1 to v2v3"
        // remove edges around flipped edge
        for( const auto& hh: surrounding_halfedges ) {
            flip_edges.left.erase( hh );
            flip_edges.left.erase( mesh_.opposite( hh ) );
        }

        CGAL::Euler::flip_edge( h, mesh_ );

        Face_descriptor f2 = mesh_.face( h );
        Face_descriptor f3 = mesh_.face( mesh_.opposite( h ) );

        update_triangle_normal( f2 );
        update_triangle_sampling( f2 );
        update_triangle_normal( f3 );
        update_triangle_sampling( f3 );

        ++n_flips;

        update_last_operation( f2 );
        update_last_operation( f3 );

        // const Point p_center_032              = CGAL::centroid( p0, p3, p2 );
        // const auto [p_center_032_proj, n_032] = get_offset_projection_and_normal( p_center_032 );
        // const Point p_center_312              = CGAL::centroid( p3, p1, p2 );
        // const auto [p_center_312_proj, n_312] = get_offset_projection_and_normal( p_center_312 );
        //
        // f_normal[f1] = n_312;
        // f_normal[f2] = n_032;
        f_normal_samples[f2] = f_samples[f2].normals[3];
        f_normal_samples[f3] = f_samples[f3].normals[3];

        // insert new/remaining short edges
        for( const auto& hh: surrounding_halfedges ) {
            const auto& v0 = mesh_.source( hh );
            const auto& v1 = mesh_.target( hh );
            const auto& p0 = mesh_.point( v0 );
            const auto& p1 = mesh_.point( v1 );
            if( !is_flip_topologically_allowed( mesh_.edge( hh ), mesh_ ) ) {
                continue;
            }

            const auto& v2 = mesh_.target( mesh_.next( hh ) );
            const auto& v3 = mesh_.target( mesh_.next( mesh_.opposite( hh ) ) );
            const auto& p2 = mesh_.point( v2 );
            const auto& p3 = mesh_.point( v3 );

            const Face_descriptor f_012 = mesh_.face( hh );
            const Face_descriptor f_031 = mesh_.face( mesh_.opposite( hh ) );

            std::vector<Vector> normals;
            normals.insert( normals.end(), f_samples[f_012].normals.begin(), f_samples[f_012].normals.end() );
            normals.insert( normals.end(), f_samples[f_031].normals.begin(), f_samples[f_031].normals.end() );
            const Vector e01_direction = normalize( p1 - p0 );
            FT e01_nd                  = oriented_normal_deviation( e01_direction, normals );

            const Vector e23_direction = normalize( p3 - p2 );
            FT e23_nd                  = oriented_normal_deviation( e23_direction, normals );

            if( e23_nd > max_normal_deviation_ && e01_nd < max_normal_deviation_ ) {
                continue;
            }

            // check min quality
            const FT q012    = meanRatioMetric_heavy_normal( { p0, p1, p2 }, f_normal_samples[f_012] );
            const FT q031    = meanRatioMetric_heavy_normal( { p0, p3, p1 }, f_normal_samples[f_031] );
            const FT qMinOld = std::min( q012, q031 );

            const Point p_center_032 = CGAL::centroid( p0, p3, p2 );
            const Point p_center_312 = CGAL::centroid( p3, p1, p2 );
            FT qMinNew;
            if( f_samples[f_012].has_one_primitive && f_samples[f_031].has_one_primitive &&
                f_samples[f_012].primitive_id == f_samples[f_031].primitive_id ) {
                // directly use this primitive, do not use the AABB_tree
                const auto [p_center_032_proj, n_032, t_032] = get_offset_projection_and_normal( p_center_032, f_samples[f_012].primitive_id );
                const auto [p_center_312_proj, n_312, t_312] = get_offset_projection_and_normal( p_center_312, f_samples[f_012].primitive_id );
                const FT q032                                = meanRatioMetric_heavy_normal( { p0, p3, p2 }, n_032 );
                const FT q312                                = meanRatioMetric_heavy_normal( { p3, p1, p2 }, n_312 );
                qMinNew                                      = std::min( q032, q312 );
            } else {
                const auto [p_center_032_proj, n_032] = get_offset_projection_and_normal( p_center_032 );
                const auto [p_center_312_proj, n_312] = get_offset_projection_and_normal( p_center_312 );
                const FT q032                         = meanRatioMetric_heavy_normal( { p0, p3, p2 }, n_032 );
                const FT q312                         = meanRatioMetric_heavy_normal( { p3, p1, p2 }, n_312 );
                qMinNew                               = std::min( q032, q312 );
            }

            if( qMinNew <= qMinOld ) {
                continue;
            }

            flip_edges.insert( flip_edge( hh, qMinNew - qMinOld ) );
        }
    }

    mesh_.remove_property_map( f_normal_samples );

    return n_flips;
}

void Remeshing::smooth() {
    LOG( INFO ) << "smooth all";

    auto [f_normals, f_normals_found] = mesh_.property_map<Face_descriptor, Vector>( "f:normal" );
    LOG_ASSERT( f_normals_found );
    auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
    LOG_ASSERT( f_samples_found );

    const int n_iter = 5;
    for( int iter = 0; iter < n_iter; ++iter ) {
        LOG( INFO ) << "smooth iteration " << iter + 1 << "/" << n_iter;
        for( const auto& v: mesh_.vertices() ) {
            if( !is_operation_required( mesh_.halfedge( v ), Remeshing_operation_types::smooth ) ) {
                continue;
            }
            const Point p_current = mesh_.point( v );

            std::vector<Vertex_descriptor> neighs;
            neighs.reserve( 6 );
            FT squared_edge_lenth = 0;
            for( const auto& vv: mesh_.vertices_around_target( mesh_.halfedge( v ) ) ) {
                neighs.push_back( vv );
                const Point p_neigh = mesh_.point( vv );
                squared_edge_lenth += ( p_neigh - p_current ).squared_length();
            }
            squared_edge_lenth /= neighs.size();

            std::vector<Point> points;
            std::transform( neighs.begin(), neighs.end(), std::back_inserter( points ), [this]( Vertex_descriptor v ) { return mesh_.point( v ); } );

            const Point p_smooth = CGAL::centroid( points.begin(), points.end() );
            // mesh_.point( v )     = p_smooth;

            const Point p_proj = project_vertex_back_by_sampling( v, p_smooth );

            if( ( p_proj - p_current ).squared_length() > 0.001 * squared_edge_lenth ) {
                // store current triangle normals and samples
                double max_nd_buf = 0;
                std::map<Face_descriptor, Vector> normal_buf;
                std::map<Face_descriptor, Triangle_sampling> samples_buf;
                for( const auto& f: mesh_.faces_around_target( mesh_.halfedge( v ) ) ) {
                    normal_buf[f]  = f_normals[f];
                    samples_buf[f] = f_samples[f];
                    double nd      = normal_deviation_sampled( f );
                    max_nd_buf     = CGAL::max( max_nd_buf, nd );
                }

                mesh_.point( v ) = p_proj;

                double max_nd_new = 0;
                for( const auto& f: mesh_.faces_around_target( mesh_.halfedge( v ) ) ) {
                    update_triangle_normal( f );
                    update_triangle_sampling( f );
                    double nd  = normal_deviation_sampled( f );
                    max_nd_new = CGAL::max( max_nd_new, nd );
                }

                // check if normal deviation got worse
                // if( max_nd_new > max_nd_buf && max_nd_new > max_normal_deviation_ ) {
                if( max_nd_new > max_normal_deviation_ && max_nd_buf < max_normal_deviation_ ) {
                    // reverse changes
                    for( const auto& f: mesh_.faces_around_target( mesh_.halfedge( v ) ) ) {
                        f_normals[f] = normal_buf[f];
                        f_samples[f] = samples_buf[f];
                    }
                    mesh_.point( v ) = p_current;
                    continue;
                }

                for( const auto& f: mesh_.faces_around_target( mesh_.halfedge( v ) ) ) {
                    update_last_operation( f );
                }
            }
        }
        if( std::filesystem::exists( output_path_ ) ) {
            print_mesh( output_path_ / ( "smoothed_" + std::to_string( n_remeshing_iterations_ ) + "_" + std::to_string( iter ) + ".off" ), mesh_ );
        }
    }
}

void Remeshing::smooth_convex_sharp_edges( const FT& e_sharpness_threshold, const size_t& n_iterations ) {
    LOG( INFO ) << "smooth convex sharp edges";

    auto [e_is_sharp, e_is_sharp_created] = mesh_.add_property_map<Edge_descriptor, bool>( "e:is_sharp" );
    LOG_ASSERT( e_is_sharp_created );
    auto [v_smooth, v_smooth_created] = mesh_.add_property_map<Vertex_descriptor, bool>( "v:smooth" );
    LOG_ASSERT( v_smooth_created );

    for( int iter = 0; iter < n_iterations; ++iter ) {
        LOG( INFO ) << "smooth iteration " << iter + 1 << "/" << n_iterations;

        bool smoothed_something = false;

        for( const auto& v: mesh_.vertices() ) {
            v_smooth[v] = false;
        }

        for( const auto& e: mesh_.edges() ) {
            const Halfedge_descriptor h = mesh_.halfedge( e );
            const Vertex_descriptor v0  = mesh_.source( h );
            const Vertex_descriptor v1  = mesh_.target( h );
            const Vertex_descriptor v2  = mesh_.target( mesh_.next( h ) );
            const Vertex_descriptor v3  = mesh_.target( mesh_.next( mesh_.opposite( h ) ) );
            const Point p0              = mesh_.point( v0 );
            const Point p1              = mesh_.point( v1 );
            const Point p2              = mesh_.point( v2 );
            const Point p3              = mesh_.point( v3 );

            const Vector n0  = normalize( CGAL::cross_product( p1 - p0, p2 - p0 ) );
            const Vector n1  = normalize( CGAL::cross_product( p3 - p0, p1 - p0 ) );
            const Vector nxn = CGAL::cross_product( n0, n1 );
            if( CGAL::scalar_product( nxn, p1 - p0 ) < 0 ) {
                // edge is concave --> not sharp
                e_is_sharp[e] = false;
            } else {
                FT angle = CGAL::approximate_angle( n0, n1 );
                if( angle < e_sharpness_threshold ) {
                    e_is_sharp[e] = false;
                } else {
                    e_is_sharp[e] = true;
                }
            }

            if( e_is_sharp[e] ) {
                v_smooth[v0] = true;
                v_smooth[v1] = true;
                v_smooth[v2] = true;
                v_smooth[v3] = true;
            }
        }

        for( const auto& v: mesh_.vertices() ) {
            if( !v_smooth[v] ) {
                continue;
            }

            std::vector<Vertex_descriptor> neighs;
            neighs.reserve( 6 );
            for( const auto& vv: mesh_.vertices_around_target( mesh_.halfedge( v ) ) ) {
                neighs.push_back( vv );
            }
            std::vector<Point> points;
            std::transform( neighs.begin(), neighs.end(), std::back_inserter( points ), [this]( Vertex_descriptor v ) { return mesh_.point( v ); } );

            const Point p_new = CGAL::centroid( points.begin(), points.end() );
            mesh_.point( v )  = CGAL::midpoint( mesh_.point( v ), p_new );

            for( const auto& f: mesh_.faces_around_target( mesh_.halfedge( v ) ) ) {
                update_triangle_normal( f );
                update_triangle_sampling( f );
            }

            smoothed_something = true;
        }
        // print_dual_contouring_mesh_complex_vertices( "../../../smoothed.off", mesh_ );

        if( !smoothed_something ) {
            LOG( INFO ) << "no convex sharp edges to smooth --> break";
            break;
        }
    }

    mesh_.remove_property_map( e_is_sharp );
    mesh_.remove_property_map( v_smooth );
}