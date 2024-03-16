#pragma once
/**
 * @file   Octree_wrapper.h
 * @brief  This file contains a class that wraps the CGAL Octree and adds vertex and edge handles.
 *
 * @author Daniel Zint
 * @date   December 2022
 */

#include "Distance.h"
#include "Octree_tables.h"
#include "types.h"

// #include <CGAL/Min_sphere_annulus_d_traits_3.h>
// #include <CGAL/Min_sphere_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Octree.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/centroid.h>

/**
 * A wrapper for CGAL Octree with additional functionality for offset generation.
 *
 * This wrapper provides handles for voxels, edges, and vertices on the leaf level. It also contains functions specific to offset generation.
 *
 * Naming convention from "A parallel dual marching cubes approach to quad only surface reconstruction - Grosso & Zint":
 * ```
 *        ^ y
 *        |
 *       v2------e2------v3
 *       /|             /|
 *     e11|           e10|
 *     /  e3          /  e1
 *   v6------e6------v7  |
 *    |   |          |   |
 *    |  v0------e0--|---v1 --> x
 *    e7 /           e5 /
 *    | e8           | e9
 *    |/             |/
 *   v4------e4------v5
 *   /
 *  < z
 * ```
 */
class Octree_wrapper {
  public:
    typedef CGAL::Octree<Kernel, std::vector<Point>> Octree;

    typedef std::size_t Vertex_handle;
    typedef std::tuple<std::size_t, std::size_t> Edge_handle;
    typedef std::size_t Voxel_handle;

    typedef Octree::Node::Global_coordinates Uniform_coords;    // coordinates on max depth level

    typedef CGAL::Min_sphere_of_points_d_traits_3<Kernel, FT> Min_sphere_traits;
    typedef CGAL::Min_sphere_of_spheres_d<Min_sphere_traits> Min_sphere;

  private:
    std::size_t min_depth_ = 0;
    std::size_t max_depth_ = 0;

    FT offset_x_ = 0;
    FT offset_y_ = 0;
    FT offset_z_ = 0;

    CGAL::Bbox_3 bbox_;

    std::size_t dim_ = 1;

    FT hx_ = 0;
    FT hy_ = 0;
    FT hz_ = 0;

    std::vector<Point> pseudo_point_cloud_;

    Octree octree_;

    std::vector<Voxel_handle> leaf_voxels_;
    std::vector<Edge_handle> leaf_edges_;
    std::vector<Vertex_handle> leaf_vertices_;
    std::unordered_map<Vertex_handle, FT> vertex_values_;
    std::unordered_map<Vertex_handle, bool> vertex_is_inside_;

    std::unordered_map<Voxel_handle, std::vector<Face_descriptor>> voxel_prims_;
    std::unordered_map<Voxel_handle, bool> voxel_is_converged_;

  public:
    /**
     * Init the octree from a bbox and with a max depth.
     *
     * The Octree requires a point cloud. Therefore, the Octree wrapper constructs a pseudo point cloud from the bbox with only two points.
     *
     * @param bbox
     * @param max_depth Naximum depth of the octree. The octree **must not** exceed the maximum depth or the indexing of voxels and vertices
     * will be wrong.
     */
    Octree_wrapper( const CGAL::Bbox_3& bbox, const std::size_t& max_depth )
        : max_depth_( max_depth ), bbox_( bbox ), offset_x_( bbox.xmin() ), offset_y_( bbox.ymin() ),
          offset_z_( bbox.zmin() ), pseudo_point_cloud_ { { bbox.xmin(), bbox.ymin(), bbox.zmin() }, { bbox.xmax(), bbox.ymax(), bbox.zmax() } },
          octree_( pseudo_point_cloud_ ) {
        dim_ = std::size_t( 1 ) << max_depth_;
        hx_  = bbox_.x_span() / dim_;
        hy_  = bbox_.y_span() / dim_;
        hz_  = bbox_.z_span() / dim_;
    }

    /**
     * Initialize the octree from the primal mesh and the iso value.
     *
     * @tparam offset_side Use `CGAL::Side_of_triangle_mesh` to generate only the inside or outside of the offset. Possible values are 0 (= double
     * sided offset, default behavior), `CGAL::ON_BOUNDED_SIDE` (inside), `CGAL::ON_UNBOUNDED_SIDE` (outside).
     *
     * @param primal_mesh The mesh from which the offset is generated.
     * @param iso_value The offset value.
     * @param depth The maximum depth reached in the initialization. This must not be lower than the max_depth defined in the constructor.
     */
    void init_from_primal_mesh( const SurfaceMesh& primal_mesh, const FT& iso_value, const std::size_t& depth_min, const std::size_t& depth ) {
        LOG_ASSERT( octree_.root().is_leaf() );
        std::function<bool( const Point& )> is_inside;
        // Sotm sotm( primal_mesh );
        std::unique_ptr<Sotm> sotm;
        if( CGAL::is_closed( primal_mesh ) ) {
            sotm      = std::make_unique<Sotm>( primal_mesh );
            is_inside = [&sotm]( const Point& p ) { return ( *sotm )( p ) != CGAL::ON_UNBOUNDED_SIDE; };
        } else {
            is_inside = []( const Point& p ) { return false; };
        }

        update_vertex_is_inside( octree_.root(), primal_mesh, is_inside );

        // add all prims to root
        for( const auto& f: primal_mesh.faces() ) {
            voxel_prims_[0].push_back( f );
        }

        std::queue<Octree::Node> q;
        for( const Octree::Node& node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            q.push( node );
        }

        std::size_t iteration_counter = 0;
        while( !q.empty() ) {
            iteration_counter++;
            if( iteration_counter % 10000 == 0 ) {
                LOG( INFO ) << std::to_string( q.size() ) << " nodes in queue. Current node is on depth " << std::to_string( q.front().depth() );
            }
            // Get the next element
            Octree::Node current = q.front();
            q.pop();

            std::vector<Point> projection_points = get_triangle_projection_points( primal_mesh, current );
            if( projection_points.empty() ) {
                continue;
            }
            bool voxel_is_empty = ball_criterion( current, projection_points, iso_value, is_inside );
            if( voxel_is_empty ) {
                continue;
            }

            if( current.depth() == depth_min ) {
                continue;
            }

            const std::size_t current_lex = node_to_lex_index( current );

            const std::vector<Face_descriptor> current_prims = voxel_prims_[current_lex];

            // split
            octree_.split( current );

            // Process each of its children
            for( int child_id = 0; child_id < Octree::Degree::value; ++child_id ) {
                // update voxel prims
                const Uniform_coords uniform_coords = uniform_coordinates( current[child_id] );
                const std::size_t lex               = lex_index( uniform_coords[0], uniform_coords[1], uniform_coords[2], max_depth_ );
                voxel_prims_[lex]                   = current_prims;
                voxel_is_converged_[lex]            = true;    // assume convergence of voxel until disproven

                update_vertex_is_inside( current[child_id], primal_mesh, is_inside );
            }

            // Add children to queue
            for( int child_id = 0; child_id < Octree::Degree::value; ++child_id ) {
                q.push( current[child_id] );
            }
            // print( "../../../octree.off", iso_value );
        }

        LOG( INFO ) << "Min refinement done. Use ball and disk criterion from now on.";
        LOG_ASSERT( q.empty() );
        for( const Octree::Node& node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            q.push( node );
        }

        while( !q.empty() ) {
            iteration_counter++;
            if( iteration_counter % 10000 == 0 ) {
                LOG( INFO ) << std::to_string( q.size() ) << " nodes in queue. Current node is on depth " << std::to_string( q.front().depth() );
            }
            // Get the next element
            Octree::Node current = q.front();
            q.pop();

            std::vector<Point> projection_points = get_triangle_projection_points( primal_mesh, current );
            if( projection_points.empty() ) {
                continue;
            }

            const std::size_t current_lex = node_to_lex_index( current );
            // if( current_lex == 38678829568 ) {
            //     LOG( WARNING ) << "DEBUG " << current_lex;
            //     for( const auto& prim: voxel_prims_[current_lex] ) {
            //         identify_unique_id( prim, primal_mesh );
            //     }
            // }

            bool voxel_is_empty = ball_criterion( current, projection_points, iso_value, is_inside );
            if( voxel_is_empty ) {
                continue;
            }

            if( current.depth() == depth ) {
                continue;
            }

            // if( current_lex == 38678829568 ) {
            //     LOG( WARNING ) << "DEBUG " << current_lex;
            //     for( const auto& prim: voxel_prims_[current_lex] ) {
            //         identify_unique_id( prim, primal_mesh );
            //     }
            // }

            const bool prims_connected = are_prims_connected( primal_mesh, iso_value, current );
            // const bool prims_connected = false;
            if( prims_connected || do_primitive_offsets_overlap( current, projection_points, iso_value, primal_mesh ) ) {
                if( is_voxel_intersected_by_its_prims( primal_mesh, current, iso_value ) ) {
                    continue;
                }
            }

            const std::vector<Face_descriptor> current_prims = voxel_prims_[current_lex];

            // split
            octree_.split( current );

            // Process each of its children
            for( int child_id = 0; child_id < Octree::Degree::value; ++child_id ) {
                // update voxel prims
                const Uniform_coords uniform_coords = uniform_coordinates( current[child_id] );
                const std::size_t lex               = lex_index( uniform_coords[0], uniform_coords[1], uniform_coords[2], max_depth_ );
                voxel_prims_[lex]                   = current_prims;
                voxel_is_converged_[lex]            = true;    // assume convergence of voxel until disproven

                update_vertex_is_inside( current[child_id], primal_mesh, is_inside );
            }

            // Add children to queue
            for( int child_id = 0; child_id < Octree::Degree::value; ++child_id ) {
                q.push( current[child_id] );
            }
            // print( "../../../octree.off", iso_value );
        }

        // split triangles into primitives
        // split_triangles_into_primitives( primal_mesh, iso_value );
        // aggressive_voronoi_filtering( primal_mesh, iso_value );

        for( const Octree::Node& node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            update_vertex_values( node, primal_mesh, iso_value );
            const size_t vox = node_to_lex_index( node );
        }

        // unsplit nodes where all children contain the same information
        // unsplit_nodes( iso_value );

        // re-compute leaf vectors
        update_leaf_vectors();
    }

    /**
     * Refine a given vector of cells independent of any criteria besides the given depth.
     *
     * Refinement is not performed if a cell is on the given depth.
     *
     * @param primal_mesh
     * @param iso_value
     * @param depth
     * @param nm_cells Vector with cells that should be refined
     */
    void refine_cell_vector_to_convergence( const SurfaceMesh& primal_mesh, const FT& iso_value, const std::size_t& depth,
                                            const std::vector<std::size_t>& nm_cells ) {
        std::function<bool( const Point& )> is_inside;
        std::unique_ptr<Sotm> sotm;
        if( CGAL::is_closed( primal_mesh ) ) {
            sotm      = std::make_unique<Sotm>( primal_mesh );
            is_inside = [&sotm]( const Point& p ) { return ( *sotm )( p ) != CGAL::ON_UNBOUNDED_SIDE; };
        } else {
            is_inside = []( const Point& p ) { return false; };
        }

        // Initialize a queue of nodes that potentially need to be refined
        std::queue<Octree::Node> q;
        for( const auto& cell: nm_cells ) {
            q.push( get_node( cell ) );
        }

        std::size_t iteration_counter = 0;
        while( !q.empty() ) {
            iteration_counter++;
            if( iteration_counter % 10000 == 0 ) {
                LOG( INFO ) << std::to_string( q.size() ) << " nodes in queue. Current node is on depth " << std::to_string( q.front().depth() );
            }
            // Get the next element
            Octree::Node current = q.front();
            q.pop();

            std::vector<Point> projection_points = get_triangle_projection_points( primal_mesh, current );
            if( projection_points.empty() ) {
                continue;
            }
            bool voxel_is_empty = ball_criterion( current, projection_points, iso_value, is_inside );
            if( voxel_is_empty ) {
                continue;
            }

            if( current.depth() == depth ) {
                continue;
            }

            const std::size_t current_lex                    = node_to_lex_index( current );
            const std::vector<Face_descriptor> current_prims = voxel_prims_[current_lex];

            // split
            octree_.split( current );

            // Process each of its children
            for( int child_id = 0; child_id < Octree::Degree::value; ++child_id ) {
                // update voxel prims
                const Uniform_coords uniform_coords = uniform_coordinates( current[child_id] );
                const std::size_t lex               = lex_index( uniform_coords[0], uniform_coords[1], uniform_coords[2], max_depth_ );
                voxel_prims_[lex]                   = current_prims;
                voxel_is_converged_[lex]            = true;    // assume convergence of voxel until disproven

                update_vertex_is_inside( current[child_id], primal_mesh, is_inside );
            }
            // print( "../../../octree.off", iso_value );
        }

        // aggressive_voronoi_filtering( primal_mesh, iso_value );

        for( const Octree::Node& node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            update_vertex_values( node, primal_mesh, iso_value );
            const size_t vox = node_to_lex_index( node );
            // check if voxel is intersected by isosurface
            // if( !is_voxel_intersected_by_isosurface( vox, iso_value ) ) {
            //    voxel_prims_[vox].clear();
            //}
        }

        // re-compute leaf vectors
        update_leaf_vectors();
    }

    /**
     * Update vertex values of the given node.
     *
     * @param node
     * @param primal_mesh
     * @param iso_value
     */
    void update_vertex_values( const Octree::Node& node, const SurfaceMesh& primal_mesh, const FT& iso_value ) {
        const Uniform_coords uniform_coords = uniform_coordinates( node );
        const std::size_t lex               = node_to_lex_index( node );
        const std::size_t df                = depth_factor( node.depth() );
        std::array<Point, 8> points         = voxel_vertex_positions( lex );
        // update vertex values using the prims of the voxel
        for( const auto& prim: voxel_prims_[lex] ) {
            bool is_inside_slab;
            const auto dists_squared = Distance::points_to_triangle_squared<8>( points, prim, primal_mesh, is_inside_slab );

            for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                const std::size_t i_point   = uniform_coords[0] + df * Tables::local_vertex_position[i][0];
                const std::size_t j_point   = uniform_coords[1] + df * Tables::local_vertex_position[i][1];
                const std::size_t k_point   = uniform_coords[2] + df * Tables::local_vertex_position[i][2];
                const std::size_t lex_point = lex_index( i_point, j_point, k_point, max_depth_ );

                if( vertex_values_.count( lex_point ) == 0 ) {
                    // point does not exist in the list --> add it
                    const int dist_sign       = vertex_is_inside_[lex_point] ? -1 : 1;
                    const FT val              = CGAL::to_double( dists_squared[i] );
                    vertex_values_[lex_point] = dist_sign * CGAL::sqrt( val );
                } else {
                    if( dists_squared[i] < vertex_values_[lex_point] * vertex_values_[lex_point] ) {
                        // point exists but has a larger value --> replace value
                        const int dist_sign       = vertex_is_inside_[lex_point] ? -1 : 1;
                        const FT val              = CGAL::to_double( dists_squared[i] );
                        vertex_values_[lex_point] = dist_sign * CGAL::sqrt( val );
                    }
                }
            }
        }

        // give voxel vertices also values if they are empty
        if( voxel_prims_[lex].empty() ) {
            for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                const std::size_t i_point   = uniform_coords[0] + df * Tables::local_vertex_position[i][0];
                const std::size_t j_point   = uniform_coords[1] + df * Tables::local_vertex_position[i][1];
                const std::size_t k_point   = uniform_coords[2] + df * Tables::local_vertex_position[i][2];
                const std::size_t lex_point = lex_index( i_point, j_point, k_point, max_depth_ );

                if( vertex_values_.count( lex_point ) == 0 ) {
                    // point does not exist in the list --> add it
                    vertex_values_[lex_point] = FT_MAX;
                }
            }
        }
    }

    void update_vertex_is_inside( const Octree::Node& node, const SurfaceMesh& primal_mesh, const std::function<bool( const Point& )> is_inside ) {
        const Uniform_coords uniform_coords = uniform_coordinates( node );
        const std::size_t lex               = node_to_lex_index( node );
        const std::size_t df                = depth_factor( node.depth() );
        // std::array<Point, 8> points         = voxel_vertex_positions( lex );
        for( int i = 0; i < Tables::N_VERTICES; ++i ) {
            const std::size_t i_point   = uniform_coords[0] + df * Tables::local_vertex_position[i][0];
            const std::size_t j_point   = uniform_coords[1] + df * Tables::local_vertex_position[i][1];
            const std::size_t k_point   = uniform_coords[2] + df * Tables::local_vertex_position[i][2];
            const std::size_t lex_point = lex_index( i_point, j_point, k_point, max_depth_ );

            if( vertex_is_inside_.count( lex_point ) == 0 ) {
                Point p = point_structured( lex_point );
                if( is_inside( p ) ) {
                    vertex_is_inside_[lex_point] = true;
                } else {
                    vertex_is_inside_[lex_point] = false;
                }
            }
        }
    }

    /**
     * **DEPRECATED** Split node, update vertex values and perform Voronoi filtering on primitives.
     *
     * @param node
     * @param primal_mesh
     * @param iso_value
     */
    void split_node( Octree::Node& node, const SurfaceMesh& primal_mesh, const FT& iso_value ) {
        // CGAL::Side_of_triangle_mesh<SurfaceMesh, CGAL::GetGeomTraits<SurfaceMesh>::type> sotm( primal_mesh );

        const Uniform_coords current_uniform = uniform_coordinates( node );
        const std::size_t current_lex        = lex_index( current_uniform[0], current_uniform[1], current_uniform[2], max_depth_ );
        // std::cout << "Split node (" << node.global_coordinates()[0] << "," << node.global_coordinates()[1] << "," << node.global_coordinates()[2]
        //          << ") depth = " << std::to_string( node.depth() ) << std::endl;

        const std::vector<Face_descriptor> current_prims = voxel_prims_[current_lex];

        // split
        octree_.split( node );

        // Process each of its children
        for( int child_id = 0; child_id < Octree::Degree::value; ++child_id ) {
            // update voxel prims
            const Uniform_coords uniform_coords = uniform_coordinates( node[child_id] );
            const std::size_t lex               = lex_index( uniform_coords[0], uniform_coords[1], uniform_coords[2], max_depth_ );
            voxel_prims_[lex]                   = current_prims;
            voxel_is_converged_[lex]            = true;    // assume convergence of voxel until disproven

            const std::size_t df = depth_factor( node[child_id].depth() );

            std::array<Point, 8> points = voxel_vertex_positions( lex );

            // init vertex values
            for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                const std::size_t i_point   = uniform_coords[0] + df * Tables::local_vertex_position[i][0];
                const std::size_t j_point   = uniform_coords[1] + df * Tables::local_vertex_position[i][1];
                const std::size_t k_point   = uniform_coords[2] + df * Tables::local_vertex_position[i][2];
                const std::size_t lex_point = lex_index( i_point, j_point, k_point, max_depth_ );

                // if point does not exist in the list --> add it
                if( vertex_values_.count( lex_point ) == 0 ) {
                    vertex_values_[lex_point] = FT_MAX;
                }
            }

            // remove voxels that do not intersect the iso-value
            if( !voxel_prims_[lex].empty() ) {
                int pattern                    = 0;
                const std::array<FT, 8> values = voxel_values( lex );

                for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                    pattern += ( values[i] <= iso_value ) << i;
                }

                if( pattern == 0 || pattern == 255 ) {
                    voxel_prims_[lex].clear();
                }
            }

            // Voronoi filtering
            if( voxel_prims_[lex].size() > 1 ) {
                const auto& prims = voxel_prims_[lex];

                std::vector<size_t> prim_types( prims.size() );
                for( size_t idx = 0; idx < prims.size(); ++idx ) {
                    prim_types[idx] = identify_unique_id_type( prims[idx], primal_mesh );
                }

                // const std::array<Point, 8>& cell_positions = points;
                const std::array<Point, 8>& cell_positions = voxel_vertex_positions_inflated( lex );

                std::vector<bool> is_covered( prims.size(), false );
                bool covered_something = false;

                std::vector<std::array<FT, 8>> cell_dists( prims.size() );
                for( size_t idx = 0; idx < prims.size(); ++idx ) {
                    bool is_inside_slab;
                    cell_dists[idx] = Distance::points_to_triangle<8>( cell_positions, prims[idx], primal_mesh, is_inside_slab );

                    // Remove primitives that do not intersect voxel
                    if( !is_inside_slab ) {
                        is_covered[idx]   = true;
                        covered_something = true;
                    }
                    int pattern = 0;
                    for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                        FT value = cell_dists[idx][i];

                        bool s = value <= iso_value;
                        pattern += s << i;
                    }

                    if( pattern == 0 || pattern == 255 ) {
                        is_covered[idx]   = true;
                        covered_something = true;
                    }
                }

                for( size_t idx1 = 0; idx1 < prims.size(); ++idx1 ) {
                    const size_t& prim1      = prims[idx1];
                    const size_t& prim1_type = prim_types[idx1];
                    if( is_covered[idx1] ) {
                        continue;
                    }
                    const auto& cell_dist1 = cell_dists[idx1];
                    for( size_t idx2 = idx1 + 1; idx2 < prims.size(); ++idx2 ) {
                        const size_t& prim2      = prims[idx2];
                        const size_t& prim2_type = prim_types[idx2];
                        if( is_covered[idx2] ) {
                            continue;
                        }

                        bool may_1_cover_2 = true;    // an element of smaller or equal dimension can always cover the other
                        bool may_2_cover_1 = false;

                        if( are_ids_incident( prim1, prim2, primal_mesh ) ) {
                            if( prim1_type == 0 && prim2_type == 1 ) {    // vertex - edge
                                may_2_cover_1 = true;
                            }
                            if( prim1_type == 1 && prim2_type == 2 ) {    // edge - triangle
                                may_2_cover_1 = true;
                            }
                        }
                        if( prim1_type == 0 ) {
                            if( prim2_type == 0 ) {    // vertex - vertex
                                may_2_cover_1 = true;
                            }
                        } else if( prim1_type == 1 ) {
                            if( prim2_type == 1 ) {    // edge - edge
                                may_2_cover_1 = true;
                            }
                        } else if( prim1_type == 2 ) {    // triangle - triangle
                            may_2_cover_1 = true;
                        }

                        const auto& cell_dist2 = cell_dists[idx2];

                        // compare distances of id1 and id2 in all corners
                        int comp_pattern  = 0;
                        int pattern_1_l_2 = 0;
                        int pattern_2_l_1 = 0;
                        int pattern_eq    = 0;

                        for( size_t l = 0; l < 8; ++l ) {
                            comp_pattern += ( cell_dist1[l] <= cell_dist2[l] ) << l;

                            pattern_1_l_2 += ( cell_dist1[l] < cell_dist2[l] ) << l;
                            pattern_2_l_1 += ( cell_dist1[l] > cell_dist2[l] ) << l;
                            pattern_eq += ( cell_dist1[l] == cell_dist2[l] ) << l;
                        }

                        bool does_1_cover_2, does_2_cover_1;
                        if( prim_types[idx1] == prim_types[idx2] ) {
                            does_1_cover_2 = ( pattern_1_l_2 | pattern_eq ) == 255;
                            does_2_cover_1 = ( pattern_2_l_1 | pattern_eq ) == 255;
                        } else {
                            does_1_cover_2 = ( pattern_1_l_2 | pattern_eq ) == 255;
                            does_2_cover_1 = ( pattern_2_l_1 ) == 255;
                        }

                        if( does_1_cover_2 && may_1_cover_2 ) {
                            is_covered[idx2]  = true;
                            covered_something = true;
                        } else if( does_2_cover_1 && may_2_cover_1 ) {
                            is_covered[idx1]  = true;
                            covered_something = true;
                        }
                    }
                }

                if( covered_something ) {
                    std::vector<Face_descriptor> ids_new;
                    ids_new.reserve( prims.size() );
                    for( size_t i = 0; i < prims.size(); ++i ) {
                        if( !is_covered[i] ) {
                            ids_new.push_back( prims[i] );
                        }
                    }
                    voxel_prims_[lex] = ids_new;
                }
            }
        }
    }

    /** Unsplit nodes where all children contain the same primitives. */
    void unsplit_nodes( const FT& iso_value ) {
        // TODO change this to a queue
        std::set<Octree::Node> nodes_to_unsplit;
        do {
            nodes_to_unsplit.clear();
            for( const Octree::Node& node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
                if( node.depth() == min_depth_ ) {
                    continue;
                }
                const auto [p, r] = node_center_and_radius( node );

                // only refine until the voxel diameter is smaller than the offset
                if( 2 * r < iso_value ) {
                    continue;
                }

                Octree::Node parent = node.parent();
                // collect all child prims and check if any child is not a leaf node
                std::set<Face_descriptor> prims;
                bool can_unsplit = true;
                for( int i = 0; i < 8; ++i ) {
                    const Octree::Node& child = parent[i];
                    if( !child.is_leaf() ) {
                        can_unsplit = false;
                        break;
                    }

                    const size_t child_lex  = node_to_lex_index( child );
                    const auto& child_prims = voxel_prims_[child_lex];

                    if( !child_prims.empty() && !voxel_is_converged_[child_lex] ) {
                        can_unsplit = false;
                        break;
                    }

                    std::vector<Face_descriptor> intersection;
                    std::set_intersection( prims.begin(), prims.end(), child_prims.begin(), child_prims.end(), std::back_inserter( intersection ) );
                    if( intersection.size() == child_prims.size() ) {
                        // child prims are already in prims --> do nothing
                    } else if( intersection.size() == prims.size() ) {
                        prims.insert( child_prims.begin(), child_prims.end() );
                    } else {
                        can_unsplit = false;
                    }
                }
                if( !can_unsplit ) {
                    continue;
                }

                // clean up child prims
                for( int i = 0; i < 8; ++i ) {
                    voxel_prims_.erase( node_to_lex_index( parent[i] ) );
                }

                voxel_prims_[node_to_lex_index( parent )] = std::vector<Face_descriptor>( prims.begin(), prims.end() );
                // unsplit node
                nodes_to_unsplit.insert( parent );
            }
            LOG( INFO ) << "nodes to unsplit: " << nodes_to_unsplit.size();
            for( auto n: nodes_to_unsplit ) {
                n.unsplit();
            }
        } while( !nodes_to_unsplit.empty() );
    }

    /**
     * Update leaf vectors.
     * This is required for vertex and edge handles to work properly and must be called whenever the octree was modified.
     */
    void update_leaf_vectors() {
        // re-compute leaf vectors
        std::unordered_set<Voxel_handle> leaf_voxels_set;
        std::set<Edge_handle> leaf_edges_set;
        for( const Octree::Node& node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            const auto& coords_uniform = uniform_coordinates( node );
            // write all leaf nodes in a set
            const std::size_t vox_lex = lex_index( coords_uniform[0], coords_uniform[1], coords_uniform[2], max_depth_ );
            // if( voxel_prims_[vox_lex].empty() ) {
            //    continue;
            //}
            leaf_voxels_set.insert( vox_lex );

            // write all leaf edges in a set
            const auto& coords_global = node.global_coordinates();
            const auto& depth         = node.depth();
            const auto& df            = depth_factor( node.depth() );
            for( const auto& edge_voxels: Tables::edge_to_voxel_neighbor ) {
                bool are_all_voxels_leafs = true;
                for( const auto& node_ijk: edge_voxels ) {
                    const std::size_t x = coords_uniform[0] + df * node_ijk[0];
                    const std::size_t y = coords_uniform[1] + df * node_ijk[1];
                    const std::size_t z = coords_uniform[2] + df * node_ijk[2];
                    // check for overflow / ignore edges on boundary
                    if( x >= dim_ || y >= dim_ || z >= dim_ ) {
                        are_all_voxels_leafs = false;
                        break;
                    }

                    const Octree::Node n = get_node( x, y, z );
                    if( n.depth() > depth ) {
                        are_all_voxels_leafs = false;
                        break;
                    }
                }
                if( are_all_voxels_leafs ) {
                    // add to leaf edge set
                    std::size_t e_gl = e_glIndex( edge_voxels[0][3], coords_global[0], coords_global[1], coords_global[2], depth );
                    leaf_edges_set.insert( { e_gl, depth } );
                }
            }
        }

        leaf_voxels_ = std::vector<Voxel_handle>( leaf_voxels_set.begin(), leaf_voxels_set.end() );
        std::sort( leaf_voxels_.begin(), leaf_voxels_.end() );
        leaf_edges_ = std::vector<Edge_handle>( leaf_edges_set.begin(), leaf_edges_set.end() );
    }

    /** perform aggressive voronoi filtering on max depth level */
    void aggressive_voronoi_filtering( const SurfaceMesh& primal_mesh, const FT& iso_value, const std::size_t& depth ) {
        for( const Octree::Node& node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            if( node.depth() != depth ) {
                continue;
            }
            const std::size_t lex = node_to_lex_index( node );
            if( voxel_is_converged_[lex] ) {
                continue;
            }
            std::array<Point, 8> points = voxel_vertex_positions_inflated( lex );
            // Voronoi filtering
            if( voxel_prims_[lex].size() > 1 ) {
                const auto& prims = voxel_prims_[lex];

                std::vector<size_t> prim_types( prims.size() );
                for( size_t idx = 0; idx < prims.size(); ++idx ) {
                    prim_types[idx] = identify_unique_id_type( prims[idx], primal_mesh );
                }

                const std::array<Point, 8>& cell_positions = points;

                std::vector<bool> is_covered( prims.size(), false );
                bool covered_something = false;

                std::vector<std::array<FT, 8>> cell_dists( prims.size() );
                for( size_t idx = 0; idx < prims.size(); ++idx ) {
                    bool is_inside_slab;
                    cell_dists[idx] = Distance::points_to_triangle<8>( cell_positions, prims[idx], primal_mesh, is_inside_slab );

                    // Remove primitives that do not intersect voxel
                    if( !is_inside_slab ) {
                        is_covered[idx]   = true;
                        covered_something = true;
                    }
                    int pattern = 0;
                    for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                        FT value = cell_dists[idx][i];

                        bool s = value <= iso_value;
                        pattern += s << i;
                    }

                    if( pattern == 0 || pattern == 255 ) {
                        is_covered[idx]   = true;
                        covered_something = true;
                    }
                }

                for( size_t idx1 = 0; idx1 < prims.size(); ++idx1 ) {
                    const size_t& prim1      = prims[idx1];
                    const size_t& prim1_type = prim_types[idx1];
                    if( is_covered[idx1] ) {
                        continue;
                    }
                    const auto& cell_dist1 = cell_dists[idx1];
                    for( size_t idx2 = idx1 + 1; idx2 < prims.size(); ++idx2 ) {
                        const size_t& prim2      = prims[idx2];
                        const size_t& prim2_type = prim_types[idx2];
                        if( is_covered[idx2] ) {
                            continue;
                        }

                        const auto& cell_dist2 = cell_dists[idx2];

                        // compare distances of id1 and id2 in all corners
                        int pattern_1_l_2 = 0;
                        int pattern_2_l_1 = 0;
                        int pattern_eq    = 0;

                        for( size_t l = 0; l < 8; ++l ) {
                            pattern_1_l_2 += ( cell_dist1[l] < cell_dist2[l] ) << l;
                            pattern_2_l_1 += ( cell_dist1[l] > cell_dist2[l] ) << l;
                            pattern_eq += ( cell_dist1[l] == cell_dist2[l] ) << l;
                        }

                        bool does_1_cover_2, does_2_cover_1;
                        if( prim_types[idx1] == prim_types[idx2] ) {
                            does_1_cover_2 = ( pattern_1_l_2 | pattern_eq ) == 255;
                            does_2_cover_1 = ( pattern_2_l_1 | pattern_eq ) == 255;
                        } else {
                            does_1_cover_2 = ( pattern_1_l_2 | pattern_eq ) == 255;
                            does_2_cover_1 = ( pattern_2_l_1 ) == 255;
                        }

                        if( does_1_cover_2 ) {
                            is_covered[idx2]  = true;
                            covered_something = true;
                        } else if( does_2_cover_1 ) {
                            is_covered[idx1]  = true;
                            covered_something = true;
                        }
                    }
                }

                if( covered_something ) {
                    std::vector<Face_descriptor> ids_new;
                    ids_new.reserve( prims.size() );
                    for( size_t i = 0; i < prims.size(); ++i ) {
                        if( !is_covered[i] ) {
                            ids_new.push_back( prims[i] );
                        }
                    }
                    voxel_prims_[lex] = ids_new;

                    // recompute vertex position
                    Point point;
                    // const bool is_point_on_offset = position_point_in_voxel( primal_mesh, iso_value, lex, point );
                    // voxel_is_converged_[lex]      = is_point_on_offset;
                }
            }
        }
    }

    /** Aggressive Voronoi filtering on all leaf voxels. */
    void aggressive_voronoi_filtering( const SurfaceMesh& primal_mesh, const FT& iso_value ) {
        for( const Octree::Node& node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            const std::size_t lex       = node_to_lex_index( node );
            std::array<Point, 8> points = voxel_vertex_positions_inflated( lex );
            const auto& prims           = voxel_prims_[lex];
            if( prims.empty() ) {
                continue;
            }

            // Voronoi filtering
            std::vector<size_t> prim_types( prims.size() );
            for( size_t idx = 0; idx < prims.size(); ++idx ) {
                prim_types[idx] = identify_unique_id_type( prims[idx], primal_mesh );
            }

            const std::array<Point, 8>& cell_positions = points;

            std::vector<bool> is_covered( prims.size(), false );
            bool covered_something      = false;
            bool voxel_is_fully_covered = false;

            std::vector<std::array<FT, 8>> cell_dists( prims.size() );
            for( size_t idx = 0; idx < prims.size(); ++idx ) {
                bool is_inside_slab;
                cell_dists[idx] = Distance::points_to_triangle<8>( cell_positions, prims[idx], primal_mesh, is_inside_slab );

                // Remove primitives that do not intersect voxel
                if( !is_inside_slab ) {
                    is_covered[idx]   = true;
                    covered_something = true;
                }
                int pattern = 0;
                for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                    FT value = cell_dists[idx][i];

                    bool s = value <= iso_value;
                    pattern += s << i;
                }

                if( pattern == 0 ) {
                    is_covered[idx]   = true;
                    covered_something = true;
                } else if( pattern == 255 ) {
                    voxel_is_fully_covered = true;
                    break;
                }
            }

            if( voxel_is_fully_covered ) {
                voxel_prims_[lex].clear();
                voxel_is_converged_[lex] = true;
                continue;
            }

            for( size_t idx1 = 0; idx1 < prims.size(); ++idx1 ) {
                const size_t& prim1      = prims[idx1];
                const size_t& prim1_type = prim_types[idx1];
                if( is_covered[idx1] ) {
                    continue;
                }
                const auto& cell_dist1 = cell_dists[idx1];
                for( size_t idx2 = idx1 + 1; idx2 < prims.size(); ++idx2 ) {
                    const size_t& prim2      = prims[idx2];
                    const size_t& prim2_type = prim_types[idx2];
                    if( is_covered[idx2] ) {
                        continue;
                    }

                    const auto& cell_dist2 = cell_dists[idx2];

                    // compare distances of id1 and id2 in all corners
                    int pattern_1_l_2 = 0;
                    int pattern_2_l_1 = 0;
                    int pattern_eq    = 0;

                    for( size_t l = 0; l < 8; ++l ) {
                        pattern_1_l_2 += ( cell_dist1[l] < cell_dist2[l] ) << l;
                        pattern_2_l_1 += ( cell_dist1[l] > cell_dist2[l] ) << l;
                        pattern_eq += ( cell_dist1[l] == cell_dist2[l] ) << l;
                    }

                    bool does_1_cover_2, does_2_cover_1;
                    if( prim_types[idx1] == prim_types[idx2] ) {
                        does_1_cover_2 = ( pattern_1_l_2 | pattern_eq ) == 255;
                        does_2_cover_1 = ( pattern_2_l_1 | pattern_eq ) == 255;
                    } else {
                        does_1_cover_2 = ( pattern_1_l_2 | pattern_eq ) == 255;
                        does_2_cover_1 = ( pattern_2_l_1 ) == 255;
                    }

                    if( does_1_cover_2 ) {
                        is_covered[idx2]  = true;
                        covered_something = true;
                    } else if( does_2_cover_1 ) {
                        is_covered[idx1]  = true;
                        covered_something = true;
                    }
                }
            }

            if( covered_something ) {
                std::vector<Face_descriptor> ids_new;
                ids_new.reserve( prims.size() );
                for( size_t i = 0; i < prims.size(); ++i ) {
                    if( !is_covered[i] ) {
                        ids_new.push_back( prims[i] );
                    }
                }
                voxel_prims_[lex] = ids_new;

                // recompute vertex position
                Point point;
                // const bool is_point_on_offset = position_point_in_voxel( primal_mesh, iso_value, lex, point );
                // voxel_is_converged_[lex]      = is_point_on_offset;
            }
        }
    }

    /**
     * Check if isovalue is visible in the voxel.
     *
     * @param iso_value
     * @param vox
     * @return
     */
    bool isovalue_is_in_voxel( const FT& iso_value, const Voxel_handle& vox ) {
        std::array<FT, Tables::N_VERTICES> s = voxel_values( vox );
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
        return true;
    }

    /**
     * Check if the voxel contains a complex offset by checking for local features and global intersections.
     *
     * @param primal_mesh
     * @param vox
     * @return true, if voxel is complex
     */
    bool is_voxel_complex( const SurfaceMesh& primal_mesh, const Voxel_handle& vox ) const {
        auto [v_is_complex, v_is_complex_found] = primal_mesh.property_map<Vertex_descriptor, bool>( "v:is_complex" );
        LOG_ASSERT( v_is_complex_found );
        auto [e_is_complex, e_is_complex_found] = primal_mesh.property_map<Edge_descriptor, bool>( "e:is_complex" );
        LOG_ASSERT( e_is_complex_found );

        const auto& prims = voxel_prims( vox );

        if( prims.empty() ) {
            return false;
        }

        // check if any of the prims is complex
        for( const auto& prim: prims ) {
            const size_t prim_type = identify_unique_id_type( prim, primal_mesh );
            switch( prim_type ) {
            case 0: {
                const Vertex_descriptor v( prim );
                if( v_is_complex[v] ) {
                    return true;
                }
                break;
            }
            case 1: {
                const Edge_descriptor e( prim - primal_mesh.number_of_vertices() );
                if( e_is_complex[e] ) {
                    return true;
                }
                break;
            }
            default:
                break;
            }
        }

        if( prims.size() == 1 ) {
            return false;
        }

        std::set<Face_descriptor> all_faces;
        for( const auto& prim: prims ) {
            // get neighbor faces
            std::set<Face_descriptor> prim_neighbors = get_incident_faces( prim, primal_mesh );
            // merge all neighbor faces
            all_faces.insert( prim_neighbors.begin(), prim_neighbors.end() );
        }

        // Check if all neighbors form a connected component
        std::vector<bool> visited( primal_mesh.number_of_faces(), false );
        std::queue<Face_descriptor> q;

        visited[all_faces.begin()->idx()] = true;
        q.push( *all_faces.begin() );

        std::vector<Face_descriptor> connected_faces;
        while( !q.empty() ) {
            const Face_descriptor f = q.front();
            q.pop();

            if( all_faces.find( f ) == all_faces.end() ) {
                continue;
            }
            connected_faces.push_back( f );

            // get adjacent faces
            for( const auto& h: primal_mesh.halfedges_around_face( primal_mesh.halfedge( f ) ) ) {
                const auto& e = primal_mesh.edge( h );
                if( e_is_complex[e] ) {
                    continue;
                }
                const auto& ff = primal_mesh.face( primal_mesh.opposite( h ) );
                if( visited[ff.idx()] ) {
                    continue;
                }
                visited[ff.idx()] = true;
                q.push( ff );
            }
        }

        if( connected_faces.size() != all_faces.size() ) {
            // LOG( INFO ) << "Global intersection in voxel " << vox;
            return true;
        } else {
            return false;
        }
    }

    /**
     * **DEPRECATED** Remove primitives that are not a face. This should be done for non-complex voxels to avoid unneccesary computation.
     *
     * @param primal_mesh
     * @param vox
     */
    void remove_non_face_prims( const SurfaceMesh& primal_mesh, const Voxel_handle& vox ) {
        const auto& prims = voxel_prims( vox );

        if( prims.empty() ) {
            return;
        }

        std::vector<Face_descriptor> prims_new;
        for( const auto& prim: prims ) {
            const size_t prim_type = identify_unique_id_type( prim, primal_mesh );
            if( prim_type == 2 ) {
                prims_new.push_back( prim );
            }
        }
        voxel_prims_[vox] = prims_new;
    }

    std::vector<Point> get_triangle_projection_points( const SurfaceMesh& primal_mesh, const Octree::Node& node ) {
        const std::size_t vox = node_to_lex_index( node );
        auto& prims           = voxel_prims( vox );
        const auto [p, r]     = node_center_and_radius( node );

        std::vector<Point> projection_points_vec;
        projection_points_vec.reserve( prims.size() );
        for( const auto& prim: prims ) {
            const auto proj = Distance::project_to_primitive( p, prim, primal_mesh );
            projection_points_vec.push_back( proj );
        }

        return projection_points_vec;
    }

    Point project_point_to_cell( const Point& p, const Octree::Node& node ) {
        const auto& pts = node_points( node );
        const FT x      = CGAL::max( pts[0].x(), CGAL::min( pts[7].x(), p.x() ) );
        const FT y      = CGAL::max( pts[0].y(), CGAL::min( pts[7].y(), p.y() ) );
        const FT z      = CGAL::max( pts[0].z(), CGAL::min( pts[7].z(), p.z() ) );
        return { x, y, z };
    }

    /**
     * Remove triangles from the primitive list if they are too far away.
     *
     * @param primal_mesh
     * @param iso_value
     * @param node
     * @return true if a triangle covers the cell. In that case, all triangles are removed from the list
     */
    bool ball_criterion( const SurfaceMesh& primal_mesh, const FT& iso_value, const Octree::Node& node ) {
        const std::size_t vox = node_to_lex_index( node );
        auto& prims           = voxel_prims( vox );
        const auto [p, r]     = node_center_and_radius( node );

        bool node_is_covered_by_offset = false;
        std::vector<Face_descriptor> prims_new;
        for( const auto& prim: prims ) {
            const auto [dist, id, prim_type] = Distance::point_to_triangle( p, prim, primal_mesh );
            if( dist + r < iso_value ) {
                // voxel is fully covered by prim
                node_is_covered_by_offset = true;
                prims_new.clear();
                break;
            } else if( dist - r > iso_value ) {
                // prim does not intersect the ball around the voxel
                continue;
            } else {
                // keep prim
                prims_new.push_back( prim );
            }
        }

        prims = prims_new;

        bool node_contains_single_triangle_offset = ( prims.size() == 1 ) && ( 2 * r < iso_value );

        return node_is_covered_by_offset || node_contains_single_triangle_offset;
    }

    bool ball_criterion( const Octree::Node& node, std::vector<Point>& projection_points, const FT& iso_value,
                         const std::function<bool( const Point& )> is_inside ) {
        const std::size_t vox = node_to_lex_index( node );
        auto& prims           = voxel_prims( vox );
        auto [p, r]           = node_center_and_radius( node );
        const FT sign         = is_inside( p ) ? -1 : 1;

        // if( p.x() > 27.5 && p.x() < 34.4 && p.y() > -20.6 && p.y() < -13.76 && p.z() > 6.88 && p.z() < 13.76 ) {
        //     LOG( WARNING ) << "DEBUG " << p;
        // }

        bool node_is_covered_by_offset = false;
        std::vector<Face_descriptor> prims_new;
        std::vector<Point> projs_new;
        for( size_t i = 0; i < prims.size(); ++i ) {
            const auto& prim = prims[i];
            const FT dist    = CGAL::sqrt( ( projection_points[i] - p ).squared_length() );

            // const FT center   = sign * dist;
            const FT center   = dist;
            const FT ball_max = center + r;
            const FT ball_min = center - r;

            if( CGAL::abs( iso_value ) <= ball_max && CGAL::abs( iso_value ) >= ball_min ) {
                // prim offset intersects the ball --> keep prim
                prims_new.push_back( prim );
                projs_new.push_back( projection_points[i] );
                continue;
            }

            if( dist + r - CGAL::abs( iso_value ) < 0 ) {
                // voxel is fully covered by prim
                node_is_covered_by_offset = true;
                prims_new.clear();
                projs_new.clear();
                break;
            }

            // prim does not intersect the ball and does not cover it --> remove prim
        }
        prims             = prims_new;
        projection_points = projs_new;

        if( !node_is_covered_by_offset ) {
            node_is_covered_by_offset = true;
            for( size_t i = 0; i < prims.size(); ++i ) {
                const auto& prim = prims[i];
                const FT dist    = CGAL::sqrt( ( projection_points[i] - p ).squared_length() );

                const FT center   = sign * dist;
                const FT ball_max = center + r;
                const FT ball_min = center - r;

                if( iso_value <= ball_max && iso_value >= ball_min ) {
                    // prim offset intersects the ball --> keep voxel
                    node_is_covered_by_offset = false;
                    break;
                }
                // prim does not intersect the ball
            }

            if( node_is_covered_by_offset ) {
                prims.clear();
                projection_points.clear();
            }
        }

        // bool node_contains_single_triangle_offset = ( prims.size() == 1 ) && ( 2 * r < CGAL::abs( iso_value ) );
        // return node_is_covered_by_offset || node_contains_single_triangle_offset;
        return node_is_covered_by_offset || prims.empty();
    }

    bool do_primitive_offsets_overlap( const Octree::Node& node, const std::vector<Point>& projection_points, const FT& iso_value,
                                       const SurfaceMesh& primal_mesh ) {
        // typedef CGAL::Min_sphere_annulus_d_traits_3<Kernel> Traits;
        // typedef CGAL::Min_sphere_d<Traits> Min_sphere;

        std::set<Point> p( projection_points.begin(), projection_points.end() );

        // Min_sphere ms( p.begin(), p.end() );

        // alternative...
        Min_sphere ms( p.begin(), p.end() );

        // if( ms.squared_radius() < iso_value * iso_value ) {
        if( ms.radius() < iso_value ) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Check if all primitives in this cell are connected.
     *
     * @param primal_mesh
     * @param iso_value
     * @param node
     * @return true if they are connected, false if they are not and their offsets are not intersecting (global intersection)
     */
    bool are_prims_connected( const SurfaceMesh& primal_mesh, const FT& iso_value, const Octree::Node& node ) const {
        const Voxel_handle vox = node_to_lex_index( node );
        const auto& prims      = voxel_prims( vox );

        if( prims.size() < 2 ) {
            return true;
        }

        // std::set<Face_descriptor> all_faces;
        // for( const auto& prim: prims ) {
        //     // get neighbor faces
        //     std::set<Face_descriptor> prim_neighbors = get_incident_faces( prim, primal_mesh );
        //     // merge all neighbor faces
        //     all_faces.insert( prim_neighbors.begin(), prim_neighbors.end() );
        // }
        //
        // CGAL::Face_filtered_graph<SurfaceMesh> ffg( primal_mesh, all_faces );
        // const size_t n_borders = PMP::number_of_borders( ffg );
        //
        //// Check if connected faces has more than one boundary
        // if( n_borders > 1 ) {
        //     return false;
        // } else {
        //     return true;
        // }

        std::map<Vertex_descriptor, size_t> vm;

        for( const auto& f: prims ) {
            for( const auto& v: primal_mesh.vertices_around_face( primal_mesh.halfedge( f ) ) ) {
                vm[v]++;
            }
        }

        for( const auto& [v, i]: vm ) {
            if( i == prims.size() ) {
                return true;
            }
        }

        return false;
    }

    std::size_t dim() const { return dim_; }
    FT hx() const { return hx_; }
    FT hy() const { return hy_; }
    FT hz() const { return hz_; }
    FT offset_x() const { return offset_x_; }
    FT offset_y() const { return offset_y_; }
    FT offset_z() const { return offset_z_; }
    std::size_t max_depth() const { return max_depth_; }
    auto& octree() { return octree_; }

    std::size_t n_edges() const { return leaf_edges_.size(); }
    std::size_t n_vertices() const { return leaf_vertices_.size(); }
    std::size_t n_voxels() const { return leaf_voxels_.size(); }

    const std::vector<Edge_handle>& leaf_edges() const { return leaf_edges_; }
    const std::vector<Vertex_handle>& leaf_vertices() const { return leaf_vertices_; }
    const std::vector<Vertex_handle>& leaf_voxels() const { return leaf_voxels_; }

    const Edge_handle& edges( const std::size_t& i ) const { return leaf_edges_[i]; }
    const Vertex_handle& vertices( const std::size_t& i ) const { return leaf_vertices_[i]; }
    const Voxel_handle& voxels( const std::size_t& i ) const { return leaf_voxels_[i]; }

    FT value( const Vertex_handle& v ) const { return vertex_values_.at( v ); }
    FT& value( const Vertex_handle& v ) { return vertex_values_[v]; }

    const std::vector<Face_descriptor>& voxel_prims( const Voxel_handle& v ) const { return voxel_prims_.at( v ); }
    std::vector<Face_descriptor>& voxel_prims( const Voxel_handle& v ) { return voxel_prims_[v]; }

    std::size_t depth_factor( const std::size_t& depth ) const { return std::size_t( 1 ) << ( max_depth_ - depth ); }

    Uniform_coords uniform_coordinates( const Octree::Node& node ) const {
        auto coords          = node.global_coordinates();
        const std::size_t df = depth_factor( node.depth() );
        for( int i = 0; i < Octree::Node::Dimension::value; ++i ) {
            coords[i] *= df;
        }

        return coords;
    }

    std::array<Point, 8> node_points( const Octree::Node& node ) const {
        auto coords          = node.global_coordinates();
        const std::size_t df = depth_factor( node.depth() );

        const std::size_t i0 = coords[0] * df;
        const std::size_t j0 = coords[1] * df;
        const std::size_t k0 = coords[2] * df;
        const std::size_t i1 = ( coords[0] + std::size_t( 1 ) ) * df;
        const std::size_t j1 = ( coords[1] + std::size_t( 1 ) ) * df;
        const std::size_t k1 = ( coords[2] + std::size_t( 1 ) ) * df;

        std::array<Point, 8> points;
        points[0] = point_structured( lex_index( i0, j0, k0, max_depth_ ) );
        points[7] = point_structured( lex_index( i1, j1, k1, max_depth_ ) );
        points[1] = { points[7].x(), points[0].y(), points[0].z() };
        points[2] = { points[0].x(), points[7].y(), points[0].z() };
        points[3] = { points[7].x(), points[7].y(), points[0].z() };
        points[4] = { points[0].x(), points[0].y(), points[7].z() };
        points[5] = { points[7].x(), points[0].y(), points[7].z() };
        points[6] = { points[0].x(), points[7].y(), points[7].z() };
        return points;
    }

    std::array<Point, 8> node_points_inflated( const Octree::Node& node ) const {
        auto coords          = node.global_coordinates();
        const std::size_t df = depth_factor( node.depth() );

        const std::size_t i0 = coords[0] * df;
        const std::size_t j0 = coords[1] * df;
        const std::size_t k0 = coords[2] * df;
        const std::size_t i1 = ( coords[0] + std::size_t( 1 ) ) * df;
        const std::size_t j1 = ( coords[1] + std::size_t( 1 ) ) * df;
        const std::size_t k1 = ( coords[2] + std::size_t( 1 ) ) * df;

        std::array<Point, 8> points;
        points[0] = point_structured( lex_index( i0, j0, k0, max_depth_ ) );
        points[0] = { points[0].x() - 0.01 * hx_, points[0].y() - 0.01 * hy_, points[0].z() - 0.01 * hz_ };
        points[7] = point_structured( lex_index( i1, j1, k1, max_depth_ ) );
        points[7] = { points[7].x() + 0.01 * hx_, points[7].y() + 0.01 * hy_, points[7].z() + 0.01 * hz_ };
        points[1] = { points[7].x(), points[0].y(), points[0].z() };
        points[2] = { points[0].x(), points[7].y(), points[0].z() };
        points[3] = { points[7].x(), points[7].y(), points[0].z() };
        points[4] = { points[0].x(), points[0].y(), points[7].z() };
        points[5] = { points[7].x(), points[0].y(), points[7].z() };
        points[6] = { points[0].x(), points[7].y(), points[7].z() };
        return points;
    }

    std::tuple<Point, FT> node_center_and_radius( const Octree::Node& node ) const {
        auto coords          = node.global_coordinates();
        const std::size_t df = depth_factor( node.depth() );

        const std::size_t i0 = coords[0] * df;
        const std::size_t j0 = coords[1] * df;
        const std::size_t k0 = coords[2] * df;
        const std::size_t i1 = ( coords[0] + std::size_t( 1 ) ) * df;
        const std::size_t j1 = ( coords[1] + std::size_t( 1 ) ) * df;
        const std::size_t k1 = ( coords[2] + std::size_t( 1 ) ) * df;

        const Point p0 = point_structured( lex_index( i0, j0, k0, max_depth_ ) );
        const Point p1 = point_structured( lex_index( i1, j1, k1, max_depth_ ) );

        const FT r = CGAL::approximate_sqrt( ( 0.5 * ( p1 - p0 ) ).squared_length() );

        return { CGAL::midpoint( p0, p1 ), r };
    }

    Point point_structured( const Uniform_coords& vertex_coordinates ) const {
        const FT x0 = offset_x_ + vertex_coordinates[0] * hx_;
        const FT y0 = offset_y_ + vertex_coordinates[1] * hy_;
        const FT z0 = offset_z_ + vertex_coordinates[2] * hz_;
        return { x0, y0, z0 };
    }

    Point point_structured( const Vertex_handle& v ) const {
        const auto [i, j, k] = ijk_index( v, max_depth_ );
        const FT x0          = offset_x_ + i * hx_;
        const FT y0          = offset_y_ + j * hy_;
        const FT z0          = offset_z_ + k * hz_;
        return { x0, y0, z0 };
    }

    Point position( const Vertex_handle& v ) const {
        const auto [i, j, k] = ijk_index( v, max_depth_ );
        const FT x0          = offset_x_ + i * hx_;
        const FT y0          = offset_y_ + j * hy_;
        const FT z0          = offset_z_ + k * hz_;
        return { x0, y0, z0 };
    }

    Uniform_coords vertex_uniform_coordinates( const Octree::Node& node, const Octree::Node::Local_coordinates local_coords ) const {
        const auto node_coords = node.global_coordinates();
        auto v_coords          = node_coords;
        for( int i = 0; i < Octree::Node::Dimension::value; ++i ) {
            v_coords[i] += std::size_t( local_coords[i] );
        }

        const auto df = depth_factor( node.depth() );
        for( int i = 0; i < Octree::Node::Dimension::value; ++i ) {
            v_coords[i] *= df;
        }

        return v_coords;
    }

    Octree::Node get_node( const std::size_t& i, const std::size_t& j, const std::size_t& k ) const {
        Octree::Node node    = octree_.root();
        const std::size_t& x = i;
        const std::size_t& y = j;
        const std::size_t& z = k;
        while( !node.is_leaf() ) {
            std::size_t dist_to_max = max_depth_ - node.depth() - 1;
            Octree::Node::Local_coordinates loc;
            if( x & ( std::size_t( 1 ) << dist_to_max ) ) {
                loc[0] = true;
            }
            if( y & ( std::size_t( 1 ) << dist_to_max ) ) {
                loc[1] = true;
            }
            if( z & ( std::size_t( 1 ) << dist_to_max ) ) {
                loc[2] = true;
            }
            node = node[loc.to_ulong()];
        }
        return node;
    }

    Octree::Node get_node( const std::size_t lex_index ) const {
        const auto [i, j, k] = ijk_index( lex_index, max_depth_ );
        return get_node( i, j, k );
    }

    std::size_t lex_index( const std::size_t& i, const std::size_t& j, const std::size_t& k, const std::size_t& depth ) const {
        std::size_t dim = ( std::size_t( 1 ) << depth ) + 1;
        return k * dim * dim + j * dim + i;
    }

    std::size_t i_index( const std::size_t& lex_index, const std::size_t& depth ) const {
        std::size_t dim = ( std::size_t( 1 ) << depth ) + 1;
        return lex_index % dim;
    }
    std::size_t j_index( const std::size_t& lex_index, const std::size_t& depth ) const {
        std::size_t dim = ( std::size_t( 1 ) << depth ) + 1;
        return ( ( lex_index / dim ) % dim );
    }
    std::size_t k_index( const std::size_t& lex_index, const std::size_t& depth ) const {
        std::size_t dim = ( std::size_t( 1 ) << depth ) + 1;
        return ( lex_index / ( dim * dim ) );
    }

    std::size_t node_to_lex_index( const Octree::Node& node ) const {
        const Uniform_coords uc = uniform_coordinates( node );
        const std::size_t lex   = lex_index( uc[0], uc[1], uc[2], max_depth_ );
        return lex;
    }

    std::tuple<std::size_t, std::size_t, std::size_t> ijk_index( const std::size_t& lex_index, const std::size_t& depth ) const {
        return std::make_tuple( i_index( lex_index, depth ), j_index( lex_index, depth ), k_index( lex_index, depth ) );
    }

    /// <summary>
    /// compute unique edge global index.
    /// </summary>
    /// <param name="e">local edge index</param>
    /// <param name="i_idx">i-index of cell</param>
    /// <param name="j_idx">j-index of cell</param>
    /// <param name="k_idx">k-index of cell</param>
    /// <returns></returns>
    std::size_t e_glIndex( const std::size_t& e, const std::size_t& i_idx, const std::size_t& j_idx, const std::size_t& k_idx,
                           const std::size_t& depth ) const {
        const unsigned long long gei_pattern_ = 670526590282893600ull;
        const size_t i                        = i_idx + (size_t)( ( gei_pattern_ >> 5 * e ) & 1 );            // global_edge_id[eg][0];
        const size_t j                        = j_idx + (size_t)( ( gei_pattern_ >> ( 5 * e + 1 ) ) & 1 );    // global_edge_id[eg][1];
        const size_t k                        = k_idx + (size_t)( ( gei_pattern_ >> ( 5 * e + 2 ) ) & 1 );    // global_edge_id[eg][2];
        const size_t offs                     = (size_t)( ( gei_pattern_ >> ( 5 * e + 3 ) ) & 3 );
        return ( 3 * lex_index( i, j, k, depth ) + offs );
    }

    std::array<FT, 8> voxel_values( const Voxel_handle& vox ) const {
        const auto [i, j, k] = ijk_index( vox, max_depth_ );
        Octree::Node node    = get_node( i, j, k );
        const auto& df       = depth_factor( node.depth() );

        std::array<Vertex_handle, 8> v;
        for( int v_id = 0; v_id < Tables::N_VERTICES; ++v_id ) {
            const auto& l  = Tables::local_vertex_position[v_id];
            const auto lex = lex_index( i + df * l[0], j + df * l[1], k + df * l[2], max_depth_ );
            v[v_id]        = lex;
        }

        std::array<FT, 8> s;
        std::transform( v.begin(), v.end(), s.begin(), [this]( const auto& e ) { return this->vertex_values_.at( e ); } );

        return s;
    }
    std::array<Point, 8> voxel_vertex_positions( const Voxel_handle& vox ) const {
        Octree::Node node = get_node( vox );
        return node_points( node );
    }
    std::array<Point, 8> voxel_vertex_positions_inflated( const Voxel_handle& vox ) const {
        Octree::Node node = get_node( vox );
        return node_points_inflated( node );
    }

    /// <summary>
    /// Get the values at the incident two vertices. Vertices are sorted in ascending order.
    /// </summary>
    /// <param name="e_id"></param>
    /// <returns></returns>
    std::array<FT, 2> edge_values( const Edge_handle& e_id ) const {
        const auto& [e_global_id, depth] = e_id;
        const auto df                    = depth_factor( depth );

        const size_t v0_lex_index = e_global_id / 3;
        const auto [i0, j0, k0]   = ijk_index( v0_lex_index, depth );

        // v1
        const std::size_t e_local_index = Tables::edge_store_index[e_global_id % 3];
        const auto& v1_local            = Tables::local_vertex_position[Tables::edge_to_vertex[e_local_index][1]];

        const std::size_t i1 = i0 + v1_local[0];
        const std::size_t j1 = j0 + v1_local[1];
        const std::size_t k1 = k0 + v1_local[2];

        const auto v0 = lex_index( df * i0, df * j0, df * k0, max_depth_ );
        const auto v1 = lex_index( df * i1, df * j1, df * k1, max_depth_ );

        return { value( v0 ), value( v1 ) };
    }

    /// <summary>
    /// Get the 4 voxels incident to an edge. If an edge has only three incident voxels, one will appear twice. The voxels are given with the uniform
    /// lexicographical index.
    /// </summary>
    /// <param name="e_id"></param>
    /// <returns></returns>
    std::array<std::size_t, 4> edge_voxels( const Edge_handle& e_id ) const {
        const auto& [e_global_id, depth] = e_id;
        const std::size_t e_local_index  = Tables::edge_store_index[e_global_id % 3];

        const auto df = depth_factor( depth );

        const size_t v0_lex_index = e_global_id / 3;
        auto [i, j, k]            = ijk_index( v0_lex_index, depth );
        i *= df;
        j *= df;
        k *= df;

        const auto& voxel_neighbors = Tables::edge_to_voxel_neighbor[e_local_index];
        Octree::Node n0             = get_node( i + voxel_neighbors[0][0], j + voxel_neighbors[0][1], k + voxel_neighbors[0][2] );
        Octree::Node n1             = get_node( i + voxel_neighbors[1][0], j + voxel_neighbors[1][1], k + voxel_neighbors[1][2] );
        Octree::Node n2             = get_node( i + voxel_neighbors[2][0], j + voxel_neighbors[2][1], k + voxel_neighbors[2][2] );
        Octree::Node n3             = get_node( i + voxel_neighbors[3][0], j + voxel_neighbors[3][1], k + voxel_neighbors[3][2] );

        const Uniform_coords n0_uniform_coords = uniform_coordinates( n0 );
        const Uniform_coords n1_uniform_coords = uniform_coordinates( n1 );
        const Uniform_coords n2_uniform_coords = uniform_coordinates( n2 );
        const Uniform_coords n3_uniform_coords = uniform_coordinates( n3 );

        std::size_t n0_lex = lex_index( n0_uniform_coords[0], n0_uniform_coords[1], n0_uniform_coords[2], max_depth_ );
        std::size_t n1_lex = lex_index( n1_uniform_coords[0], n1_uniform_coords[1], n1_uniform_coords[2], max_depth_ );
        std::size_t n2_lex = lex_index( n2_uniform_coords[0], n2_uniform_coords[1], n2_uniform_coords[2], max_depth_ );
        std::size_t n3_lex = lex_index( n3_uniform_coords[0], n3_uniform_coords[1], n3_uniform_coords[2], max_depth_ );

        return { n0_lex, n1_lex, n2_lex, n3_lex };

        // return { value( i0, j0, k0 ), value( i1, j1, k1 ) };
    }

    std::array<std::size_t, 4> voxels_incident_to_edge( const Edge_handle& e_id ) const { return edge_voxels( e_id ); }

    void print( const std::string& filename, const FT& iso_value ) const {
        SurfaceMesh m;
        for( Octree::Node node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            const auto uniform_coords = uniform_coordinates( node );
            const std::size_t lex     = lex_index( uniform_coords[0], uniform_coords[1], uniform_coords[2], max_depth_ );
            if( voxel_prims_.at( lex ).size() == 0 ) {
                continue;
            }

            const auto vals = voxel_values( lex );

            std::array<bool, 8> s;
            s[0]        = vals[0] > iso_value;
            s[1]        = vals[1] > iso_value;
            s[2]        = vals[2] > iso_value;
            s[3]        = vals[3] > iso_value;
            s[4]        = vals[4] > iso_value;
            s[5]        = vals[5] > iso_value;
            s[6]        = vals[6] > iso_value;
            s[7]        = vals[7] > iso_value;
            int pattern = 0;
            for( size_t it = 0; it < s.size(); ++it ) {
                pattern += static_cast<int>( s[it] ) << it;
            }
            if( pattern == 0 || pattern == 255 ) {
                continue;
            }

            auto coords            = node.global_coordinates();
            std::array<Point, 8> p = voxel_vertex_positions_inflated( lex );

            std::array<Vertex_descriptor, 8> v;
            std::transform( p.begin(), p.end(), v.begin(), [&m]( const Point& e ) { return m.add_vertex( e ); } );

            m.add_face( v[0], v[2], v[3], v[1] );
            m.add_face( v[0], v[4], v[6], v[2] );
            m.add_face( v[2], v[6], v[7], v[3] );
            m.add_face( v[1], v[3], v[7], v[5] );
            m.add_face( v[0], v[1], v[5], v[4] );
            m.add_face( v[7], v[6], v[4], v[5] );
        }

        CGAL::IO::write_OFF( filename, m, CGAL::parameters::stream_precision( 10 ) );
    }

  private:
    bool is_voxel_intersected_by_isosurface( const Voxel_handle& vox, const FT& iso_value ) const {
        const auto vals = voxel_values( vox );
        std::array<bool, 8> s;
        s[0]        = vals[0] > iso_value;
        s[1]        = vals[1] > iso_value;
        s[2]        = vals[2] > iso_value;
        s[3]        = vals[3] > iso_value;
        s[4]        = vals[4] > iso_value;
        s[5]        = vals[5] > iso_value;
        s[6]        = vals[6] > iso_value;
        s[7]        = vals[7] > iso_value;
        int pattern = 0;
        for( size_t it = 0; it < s.size(); ++it ) {
            pattern += static_cast<int>( s[it] ) << it;
        }
        if( pattern == 0 || pattern == 255 ) {
            return false;
        } else {
            return true;
        }
    }

    bool is_voxel_intersected_by_its_prims( const SurfaceMesh& primal_mesh, const Octree::Node& node, const FT& iso_value ) {
        const std::size_t vox = node_to_lex_index( node );
        const auto& prims     = voxel_prims_[vox];
        if( prims.empty() ) {
            return true;
        }

        std::vector<FT> vals( 8, FT_MAX );

        const auto points = voxel_vertex_positions( vox );
        for( const auto& prim: prims ) {
            bool is_inside_slab;
            const auto dists = Distance::points_to_triangle_squared<8>( points, prim, primal_mesh, is_inside_slab );
            for( int i = 0; i < Tables::N_VERTICES; ++i ) {
                vals[i] = std::min( vals[i], dists[i] );
            }
        }

        const Uniform_coords uniform_coords = uniform_coordinates( node );
        const std::size_t df                = depth_factor( node.depth() );
        for( int i = 0; i < Tables::N_VERTICES; ++i ) {
            const std::size_t i_point   = uniform_coords[0] + df * Tables::local_vertex_position[i][0];
            const std::size_t j_point   = uniform_coords[1] + df * Tables::local_vertex_position[i][1];
            const std::size_t k_point   = uniform_coords[2] + df * Tables::local_vertex_position[i][2];
            const std::size_t lex_point = lex_index( i_point, j_point, k_point, max_depth_ );
            if( vertex_is_inside_[lex_point] ) {
                vals[i] *= -1;
            }
        }

        const FT iso_value_sign = iso_value >= 0 ? 1 : -1;

        std::array<bool, 8> s;
        s[0]        = vals[0] > iso_value_sign * iso_value * iso_value;
        s[1]        = vals[1] > iso_value_sign * iso_value * iso_value;
        s[2]        = vals[2] > iso_value_sign * iso_value * iso_value;
        s[3]        = vals[3] > iso_value_sign * iso_value * iso_value;
        s[4]        = vals[4] > iso_value_sign * iso_value * iso_value;
        s[5]        = vals[5] > iso_value_sign * iso_value * iso_value;
        s[6]        = vals[6] > iso_value_sign * iso_value * iso_value;
        s[7]        = vals[7] > iso_value_sign * iso_value * iso_value;
        int pattern = 0;
        for( size_t it = 0; it < s.size(); ++it ) {
            pattern += static_cast<int>( s[it] ) << it;
        }
        if( pattern == 0 || pattern == 255 ) {
            return false;
        } else {
            return true;
        }
    }
};