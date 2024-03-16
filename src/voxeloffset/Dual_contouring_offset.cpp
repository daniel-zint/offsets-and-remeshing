#include "Dual_contouring_offset.h"

#include "Distance.h"
#include "Dual_contouring_octree_3.h"
#include "Octree_split_predicates.h"
#include "Remeshing.h"
#include "utilities.h"

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/helpers.h>
/* public */

void Dual_contouring_offset::compute_offset( const SurfaceMesh& primal_mesh ) {
    if( std::filesystem::exists( output_path_ ) ) {
        octree_->print( ( output_path_ / "octree.off" ).string(), iso_value_ );
    }
    compute_mesh_from_octree( primal_mesh );

    if( std::filesystem::exists( output_path_ ) ) {
        print_mesh( mesh_, output_path_ / "dual_contouring_offset.off" );
        octree_->print( ( output_path_ / "octree.off" ).string(), iso_value_ );
    }

    CGAL::is_valid_polygon_mesh( mesh_, true );
    LOG( INFO ) << "is quad mesh: " << CGAL::is_quad_mesh( mesh_ );
    LOG( INFO ) << "is closed mesh: " << CGAL::is_closed( mesh_ );

    auto [v_voxel_handle, v_voxel_handle_found] = mesh_.property_map<Vertex_descriptor, size_t>( "v:voxel_handle" );
    LOG_ASSERT( v_voxel_handle );
    auto [v_primitives, v_primitives_found] = mesh_.property_map<Vertex_descriptor, std::vector<Face_descriptor>>( "v:primitives" );
    LOG_ASSERT( v_primitives_found );

    // add primitives of incident voxels to vertices that come from empty voxels
    for( const auto& v: mesh_.vertices() ) {
        if( v_primitives[v].size() == 0 ) {
            std::set<Face_descriptor> prims_set;
            for( const auto& vv: mesh_.vertices_around_target( mesh_.halfedge( v ) ) ) {
                prims_set.insert( v_primitives[vv].begin(), v_primitives[vv].end() );
            }
            std::vector<Face_descriptor> prims( prims_set.begin(), prims_set.end() );
            v_primitives[v] = prims;
        }
    }

    auto [v_bbox, v_bbox_created] = mesh_.add_property_map<Vertex_descriptor, CGAL::Bbox_3>( "v:bbox" );
    LOG_ASSERT( v_bbox_created );
    for( const auto& v: mesh_.vertices() ) {
        if( v_voxel_handle[v] == -1 ) {
            continue;
        }
        auto p    = octree_->node_points( octree_->get_node( v_voxel_handle[v] ) );
        v_bbox[v] = p[0].bbox() + p[1].bbox() + p[2].bbox() + p[3].bbox() + p[4].bbox() + p[5].bbox() + p[6].bbox() + p[7].bbox();
    }

    // triangulate mesh
    LOG( INFO ) << "Triangulate";
    triangulate();
    CGAL::is_valid_polygon_mesh( mesh_, true );
    if( CGAL::is_closed( mesh_ ) ) {
        LOG( INFO ) << "Mesh is closed";
    } else {
        LOG( WARNING ) << "Mesh is not closed. Perform hole filling.";
        for( const auto& h: mesh_.halfedges() ) {
            if( !mesh_.is_border( h ) ) {
                continue;
            }
            std::vector<Face_descriptor> patch_facets;
            PMP::triangulate_hole( mesh_, h, std::back_inserter( patch_facets ) );
        }
        if( !CGAL::is_closed( mesh_ ) ) {
            LOG( ERROR ) << "Could not fill holes!";
        }
    }
    if( std::filesystem::exists( output_path_ ) ) {
        print_mesh( mesh_, output_path_ / "triangulated.off" );
    }

    reposition_vertices( primal_mesh );

    if( std::filesystem::exists( output_path_ ) ) {
        LOG( INFO ) << "Print rep mesh";
        print_mesh( mesh_, output_path_ / "triangulated_rep.off" );
    }
}

void Dual_contouring_offset::print_distance_metrics( const SurfaceMesh& primal_mesh ) {
    auto [v_primitives, v_primitives_found] = mesh_.property_map<Vertex_descriptor, std::vector<Face_descriptor>>( "v:primitives" );
    LOG_ASSERT( v_primitives_found );

    FT max_dist = 0;
    FT avg_dist = 0;
    for( const auto& v: mesh_.vertices() ) {
        const auto& p     = mesh_.point( v );
        const auto& prims = v_primitives[v];

        FT vertex_max_dist = 0;
        for( const auto& prim: prims ) {
            auto [dist, idx, type] = Distance::point_to_triangle( p, prim, primal_mesh );
            vertex_max_dist        = std::max( dist, vertex_max_dist );
        }
        max_dist = std::max( max_dist, vertex_max_dist );
        avg_dist += vertex_max_dist;
    }
    avg_dist /= mesh_.number_of_vertices();

    LOG( INFO ) << "Iso-value = " << iso_value_;
    // LOG( INFO ) << "Grid spacing = " << grid_spacing;
    LOG( INFO ) << "AVG dist = " << std::abs( avg_dist - iso_value_ );
    LOG( INFO ) << "MAX dist = " << std::abs( max_dist - iso_value_ );
    LOG( INFO ) << "relative AVG dist = " << std::abs( ( avg_dist - iso_value_ ) / iso_value_ );
    LOG( INFO ) << "relative MAX dist = " << std::abs( ( max_dist - iso_value_ ) / iso_value_ );

    FT hausdorff = PMP::approximate_Hausdorff_distance<TAG>( mesh_, primal_mesh, CGAL::parameters::number_of_points_per_area_unit( 4000 ) );
    LOG( INFO ) << "Hausdorff distance offset -> mesh = " << hausdorff;
    LOG( INFO ) << "relative Hausdorff = " << ( hausdorff - iso_value_ ) / iso_value_;    // TODO replace iso_value by iso_value_

    // LOG( INFO ) << "Hausdorff distance mesh -> offset = "
    //            << PMP::approximate_Hausdorff_distance<TAG>( primal_mesh, dual_contouring_mesh_,
    //                                                         CGAL::parameters::number_of_points_per_area_unit( 4000 ) );
}

/* private */

void Dual_contouring_offset::init_octree( const SurfaceMesh& primal_mesh ) {
    LOG( INFO ) << "Init octree";

    // compute bounding box
    CGAL::Bbox_3 aabb = PMP::bbox( primal_mesh );
    const FT aabb_min = CGAL::min( aabb.xmin(), CGAL::min( aabb.ymin(), aabb.zmin() ) ) - CGAL::abs( iso_value_ ) * 1.01;
    const FT aabb_max = CGAL::max( aabb.xmax(), CGAL::max( aabb.ymax(), aabb.zmax() ) ) + CGAL::abs( iso_value_ ) * 1.01;
    aabb              = CGAL::Bbox_3( aabb_min, aabb_min, aabb_min, aabb_max, aabb_max, aabb_max );

    octree_ = std::make_shared<Octree_wrapper>( aabb, max_depth_ );
    octree_->init_from_primal_mesh( primal_mesh, iso_value_, depth_, mid_depth_ );

    std::vector<std::size_t> nm_cells;
    for( int i = 0; i < max_depth_ - depth_; ++i ) {
        nm_cells.clear();
        nm_cells = CGAL::find_non_manifold_cells( *octree_, primal_mesh, iso_value_ );
        if( nm_cells.empty() ) {
            LOG( INFO ) << "No non-manifold vertices found";
            break;
        }
        LOG( INFO ) << "Splitting " << nm_cells.size() << " cells with non-manifold vertices";
        octree_->refine_cell_vector_to_convergence( primal_mesh, iso_value_, max_depth_, nm_cells );
    }
}

void Dual_contouring_offset::compute_mesh_from_octree( const SurfaceMesh& primal_mesh ) {
    LOG( INFO ) << "Compute DC mesh";
    // create DC mesh
    CGAL::Dual_contouring_3::Positioning::Voxel_center dc_positioning;

    CGAL::make_quad_mesh_using_dual_contouring( *octree_, iso_value_, mesh_, dc_positioning );

    auto [v_voxel_handle, v_voxel_handle_found] = mesh_.property_map<Vertex_descriptor, size_t>( "v:voxel_handle" );
    LOG_ASSERT( v_voxel_handle );

    auto [v_primitives, v_primitives_created] = mesh_.add_property_map<Vertex_descriptor, std::vector<Face_descriptor>>( "v:primitives" );
    LOG_ASSERT( v_primitives_created );

    for( const auto& v: mesh_.vertices() ) {
        if( v_voxel_handle[v] != -1 ) {
            v_primitives[v] = octree_->voxel_prims( v_voxel_handle[v] );
        }
    }

    // print_dual_contouring_mesh( "../../../resultOct.off" );

    // refine / compute point positions
}

void Dual_contouring_offset::reposition_vertices( const SurfaceMesh& primal_mesh ) {
    LOG( INFO ) << "Reposition vertices";

    auto [v_voxel_handle, v_voxel_handle_found] = mesh_.property_map<Vertex_descriptor, size_t>( "v:voxel_handle" );
    LOG_ASSERT( v_voxel_handle );

    auto [v_primitives, v_primitives_found] = mesh_.property_map<Vertex_descriptor, std::vector<Face_descriptor>>( "v:primitives" );
    LOG_ASSERT( v_primitives_found );

    CGAL::Dual_contouring_3::Positioning::QEM_SVD_hermite<true> dc_positioning( primal_mesh );

    // #pragma omp parallel for schedule( dynamic )
    for( int v_index = 0; v_index < mesh_.number_of_vertices(); ++v_index ) {
        Vertex_descriptor v( v_index );

        if( v_voxel_handle[v] == -1 || v_primitives[v].empty() ) {
            continue;
        }

        const bool is_converged = false;
        if( !is_converged ) {
            dc_positioning.position( *octree_, mesh_, primal_mesh, iso_value_, v );
        }
    }
    LOG( INFO ) << "Reposition vertices done";
}

void Dual_contouring_offset::triangulate() {
    auto& m                                       = mesh_;
    auto [e_is_from_quad, e_is_from_quad_created] = mesh_.add_property_map<Edge_descriptor, bool>( "e:is_from_quad" );
    LOG_ASSERT( e_is_from_quad_created );
    for( const auto& e: m.edges() ) {
        e_is_from_quad[e] = true;
    }

    if( !PMP::triangulate_faces( m ) ) {
        LOG( INFO ) << "Triangulate degenerated elements";
        // triangulate degenerated quads
        for( const auto& f: m.faces() ) {
            size_t degree = m.degree( f );
            auto h1       = m.halfedge( f );
            for( size_t i = 0; i + 3 < degree; ++i ) {
                auto h2 = m.next( m.next( h1 ) );
                CGAL::Euler::split_face( h1, m.prev( m.prev( h1 ) ), m );
                h1 = h2;
            }
        }
    }
    LOG_ASSERT( CGAL::is_triangle_mesh( m ) );
    auto [v_primitives, v_primitives_found] = m.property_map<Vertex_descriptor, std::vector<Face_descriptor>>( "v:primitives" );
    LOG_ASSERT( v_primitives_found );

    for( const auto& e: m.edges() ) {
        if( e_is_from_quad[e] ) {
            continue;
        }
        const auto& h = m.halfedge( e );

        const auto& f0 = m.face( h );
        const auto& f1 = m.face( m.opposite( h ) );

        // v2 --- v1   v2 --- v1
        //  |ta / |     | \td |
        //  | /tb |     |tc \ |
        // v0 --- v3   v0 --- v3

        std::array<Vertex_descriptor, 4> v;
        v[0] = m.source( h );
        v[1] = m.target( h );
        v[2] = m.target( m.next( h ) );
        v[3] = m.target( m.next( m.opposite( h ) ) );

        std::size_t ta_prim  = 0;
        std::size_t tb_prim  = 0;
        std::size_t tc_prim  = 0;
        std::size_t td_prim  = 0;
        std::size_t ta_count = 0;
        std::size_t tb_count = 0;
        std::size_t tc_count = 0;
        std::size_t td_count = 0;

        std::vector<Face_descriptor> ta_prims;
        ta_prims.insert( ta_prims.end(), v_primitives[v[0]].begin(), v_primitives[v[0]].end() );
        ta_prims.insert( ta_prims.end(), v_primitives[v[1]].begin(), v_primitives[v[1]].end() );
        ta_prims.insert( ta_prims.end(), v_primitives[v[2]].begin(), v_primitives[v[2]].end() );
        std::tie( ta_prim, ta_count ) = maxFreq( ta_prims );

        std::vector<Face_descriptor> tb_prims;
        tb_prims.insert( tb_prims.end(), v_primitives[v[0]].begin(), v_primitives[v[0]].end() );
        tb_prims.insert( tb_prims.end(), v_primitives[v[3]].begin(), v_primitives[v[3]].end() );
        tb_prims.insert( tb_prims.end(), v_primitives[v[1]].begin(), v_primitives[v[1]].end() );
        std::tie( tb_prim, tb_count ) = maxFreq( tb_prims );

        std::vector<Face_descriptor> tc_prims;
        tc_prims.insert( tc_prims.end(), v_primitives[v[0]].begin(), v_primitives[v[0]].end() );
        tc_prims.insert( tc_prims.end(), v_primitives[v[3]].begin(), v_primitives[v[3]].end() );
        tc_prims.insert( tc_prims.end(), v_primitives[v[2]].begin(), v_primitives[v[2]].end() );
        std::tie( tc_prim, tc_count ) = maxFreq( tc_prims );

        std::vector<Face_descriptor> td_prims;
        td_prims.insert( td_prims.end(), v_primitives[v[1]].begin(), v_primitives[v[1]].end() );
        td_prims.insert( td_prims.end(), v_primitives[v[2]].begin(), v_primitives[v[2]].end() );
        td_prims.insert( td_prims.end(), v_primitives[v[3]].begin(), v_primitives[v[3]].end() );
        std::tie( td_prim, td_count ) = maxFreq( td_prims );

        if( tc_count + td_count > ta_count + tb_count ) {
            CGAL::Euler::flip_edge( h, m );
        }
    }

    // make sure there are no degree 2 vertices
    for( const auto& v: m.vertices() ) {
        if( m.degree( v ) == 2 ) {
            // flip one edge towards this vertex
            Halfedge_descriptor h = m.halfedge( v );
            h                     = m.opposite( h );
            h                     = m.next( h );
            CGAL::Euler::flip_edge( h, m );
        }
    }

    m.remove_property_map( e_is_from_quad );
}