#pragma once
/**
 * @file   Remeshing.h
 * @brief  Contains the Remeshing class which uses the normal deviation to remesh an offset.
 *
 * @author Daniel Zint
 * @date   January 2023
 */

#include "Distance.h"
#include "types.h"
#include "utilities.h"

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/helpers.h>
#include <Eigen/Eigenvalues>
#include <glog/logging.h>

enum Remeshing_operation_types { split = 0, collapse = 1, flip = 2, smooth = 3 };
const size_t n_remeshing_operation_types = 4;

/**
 * Remeshing for offset surfaces.
 *
 * If the given offset has a wrong topology, remeshing might fail.
 * The remesher requires an initial offset surface, the primal mesh, and the offset distance. The resolution of the remeshing is determined by four
 * variables, the minimum and maximum normal deviation, and the minimum and maximum edge length (which is currently computed automatically, based on
 * the normal deviation). Finally, there is a check for low quality (degenerated or flipped) triangles.
 *
 * If the offset distance goes towards zero, the edge lenth also goes towards zero and therefore, the result is not very reasonable. This behavior can
 * be adapted by changing the edge lengths.
 *
 */
class Remeshing {
    SurfaceMesh mesh_;                     // offset mesh
    SurfaceMesh primal_mesh_;              // primal mesh
    FT iso_value_            = -FT_MAX;    // offset distance
    FT min_quality_          = 0.05;       // minimum quality (mean ratio metric) every triangle must have.
    FT min_normal_deviation_ = 1;          // mimimum normal deviation
    FT max_normal_deviation_ = 3;          // maximum normal deviation
    FT lMin_                 = -FT_MAX;    // minimum edge length
    FT lMax_                 = -FT_MAX;    // maximum edge lenth
    AABB_tree tree_;                       // AABB tree for finding the closest primitive triangle

    size_t n_remeshing_iterations_ = 0;

    size_t n_sample_calls_ = 0;

    std::filesystem::path output_path_;

    struct Triangle_sampling {
        std::array<Point, 4> points;
        std::array<Vector, 4> normals;
        bool needs_update      = true;
        bool has_one_primitive = false;
        size_t primitive_id    = INVALID_UNIQUE_ID;
    };

  public:
    /**
     * Constructor for the Remeshing class.
     *
     * @param offset_mesh An initial offset mesh, generated for example with Dual Contouring
     * @param primal_mesh The primal mesh that is offsetted
     * @param offset_distance The signed offset distance
     */
    Remeshing( const SurfaceMesh& offset_mesh, const SurfaceMesh& primal_mesh, const FT& offset_distance )
        : mesh_( offset_mesh ), primal_mesh_( primal_mesh ), iso_value_( offset_distance ),
          tree_( faces( primal_mesh_ ).first, faces( primal_mesh_ ).second, primal_mesh_ ) {
        // LOG_ASSERT( iso_value_ > 0 );
        LOG_ASSERT( min_normal_deviation_ > 0 );
        LOG_ASSERT( max_normal_deviation_ > 0 );
        // compute lMin_ and lMax_ from max/min normal deviation
        lMin_ = 2 * iso_value_ * std::sin( min_normal_deviation_ * ( CGAL_PI / 180. ) );
        lMax_ = 2 * iso_value_ * std::sin( max_normal_deviation_ * ( CGAL_PI / 180. ) );
        // LOG_ASSERT( lMin_ > 0 );
        // LOG_ASSERT( lMax_ > 0 );
    }

    /**
     * Remesh the offset mesh.
     *
     * @param n_iterations Number of remeshing iterations
     */
    void remesh( const size_t n_iterations = 5 );

    /**
     * Set the maximum normal deviation.
     *
     * The minimum normal deviation and also l_min and l_max are adapted accordingly.
     *
     * @param n_iterations Number of remeshing iterations
     */
    void set_max_normal_deviation( const FT& max_normal_deviation ) {
        LOG_ASSERT( max_normal_deviation > 0 );
        max_normal_deviation_ = max_normal_deviation;
        min_normal_deviation_ = ( 4. / 9. ) * max_normal_deviation_;

        lMin_ = 2 * iso_value_ * std::sin( min_normal_deviation_ * ( CGAL_PI / 180. ) );
        lMax_ = 2 * iso_value_ * std::sin( max_normal_deviation_ * ( CGAL_PI / 180. ) );
    }

    const FT& min_normal_deviation() const { return min_normal_deviation_; }
    const FT& max_normal_deviation() const { return max_normal_deviation_; }

    FT& l_min() { return lMin_; }
    const FT& l_min() const { return lMin_; }
    FT& l_max() { return lMax_; }
    const FT& l_max() const { return lMax_; }

    auto& output_path() { return output_path_; }

    /** The offset mesh */
    const SurfaceMesh& mesh() const { return mesh_; }

    const size_t n_sample_calls() const { return n_sample_calls_; }

    /**
     * Smooth vertices in close proximity of convex sharp edges.
     * Convex sharp edges must not appear in offset meshes. However, if the initial offset has the wrong topology, these edges may appear. This method
     * removes tangling in such regions. If the initial offset is of correct topology, this method should not be required.
     *
     * It is only performed once after remeshing.
     */
    void smooth_convex_sharp_edges( const FT& e_sharpness_threshold, const size_t& n_iterations = 10 );

  private:
    /**
     * Split faces with large normal deviation at their longest edge.
     * The method does not re-insert faces after splitting. Otherwise, tangled regions would be refined indefinitely.
     */
    size_t split();
    /**
     * Collapse edges that have a small normal deviation at neighboring elements.
     * Edges are not collapsed, if an incident edge would become larger than lMax_, or if the collapse causes normal deviations larger than
     * max_normal_deviation_.
     * Adjacent edges are updated in the queue.
     */
    size_t collapse();

    size_t collapse_halfedge();

    size_t collapse_halfedge_untangle();

    /**
     * Remove low quality and inverted triangles by collapsing the shortest incident edge possible.
     * This method helps dealing with bad input meshes.
     */
    size_t collapse_low_quality_triangles();
    /**
     * Flip edges whenever this improves the mean ratio metric with heavy normal.
     * The metric is chosen such that in planar regions it defaults to flipping according to the Delaunay criterion. In curved regions however, the
     * normal deviation has a major impact. This behavior pulls edges into creases.
     */
    size_t flip_delaunay();
    /**
     * Performs the QEM back projection with the Laplacian smoothed vertex as initial position.
     */
    void smooth();

    /**
     * A simple sampling scheme for the triangle (p0,p1,p2).
     * Sampling consists of four points, one close to each corner and one in the center.
     *
     * @param p0
     * @param p1
     * @param p2
     * @return
     */
    Vector get_normal_by_sampling( const Point& p0, const Point& p1, const Point& p2 ) {
        // get normals at sample points
        const Vector p_vec0 = p0 - CGAL::ORIGIN;
        const Vector p_vec1 = p1 - CGAL::ORIGIN;
        const Vector p_vec2 = p2 - CGAL::ORIGIN;

        const Point s0              = CGAL::ORIGIN + ( 0.8 * p_vec0 + 0.1 * p_vec1 + 0.1 * p_vec2 );
        const Point s1              = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.8 * p_vec1 + 0.1 * p_vec2 );
        const Point s2              = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.1 * p_vec1 + 0.8 * p_vec2 );
        const Point s3              = CGAL::centroid( p0, p1, p2 );
        const auto [proj0, n_proj0] = get_offset_projection_and_normal( s0 );
        const auto [proj1, n_proj1] = get_offset_projection_and_normal( s1 );
        const auto [proj2, n_proj2] = get_offset_projection_and_normal( s2 );
        const auto [proj3, n_proj3] = get_offset_projection_and_normal( s3 );

        return normalize( n_proj0 + n_proj1 + n_proj2 + n_proj3 );
    }

    Vector get_normal_by_sampling( const Face_descriptor& f ) {
        const Halfedge_descriptor h = mesh_.halfedge( f );
        const Vertex_descriptor v0  = mesh_.target( h );
        const Vertex_descriptor v1  = mesh_.target( mesh_.next( h ) );
        const Vertex_descriptor v2  = mesh_.source( h );
        const Point p0              = mesh_.point( v0 );
        const Point p1              = mesh_.point( v1 );
        const Point p2              = mesh_.point( v2 );

        return get_normal_by_sampling( p0, p1, p2 );
    }

    /**
     * For any point in space, find the closest triangle and the projection point on its offset.
     *
     * After finding the closest triangle using the AABB tree, the closest triangle primitive is determined. Finally, the point is projected onto the
     * primitive offset and the offset normal is computed.
     *
     * TODO: The normal is always normalize(p - proj). The method could be reduced to get_offset_projection and another method that calls that one
     * could then also add the normal. That might have a positive impact on performance.
     *
     * @param p Query point
     * @return {projection point, offset normal}
     */
    std::tuple<Point, Vector> get_offset_projection_and_normal( const Point& p ) {
        ++n_sample_calls_;
        AABB_tree::Point_and_primitive_id papid = tree_.closest_point_and_primitive( p );
        Face_descriptor triangle                = papid.second;

        if( iso_value_ == 0 ) {
            const Vector normal = PMP::compute_face_normal( triangle, primal_mesh_ );

            // auto [f_unique_id, f_unique_id_found] = primal_mesh_.property_map<Face_descriptor, size_t>( "f:unique_id" );
            // LOG_ASSERT( f_unique_id_found );
            // const size_t id = f_unique_id[triangle];

            return { papid.first, normal };
        }

        const auto [idx, type] = Distance::closest_triangle_primitive( p, triangle, primal_mesh_ );

        switch( type ) {
        case 0:    // vertex
        {
            const Vertex_descriptor v( idx );
            const auto s        = Distance::get_primitive_offset( primal_mesh_, v, iso_value_ );
            const Point proj    = s.project( p );
            const Vector normal = normalize( proj - s.origin );

            return { proj, normal };
        }
        case 1:    // edge
        {
            Edge_descriptor e( idx );
            const auto c          = Distance::get_primitive_offset( primal_mesh_, e, iso_value_ );
            const Point proj      = c.project( p );
            const auto proj_plane = c.projection_plane( p );
            const Vector normal( proj_plane[0], proj_plane[1], proj_plane[2] );

            return { proj, normal };
        }
        case 2:    // triangle
        {
            const auto plane    = Distance::get_primitive_offset( primal_mesh_, triangle, iso_value_, p );
            const Point proj    = plane.project( p );
            const Vector normal = plane.normal;

            return { proj, normal };
        }
        default:
            LOG( FATAL ) << "Unknown id type: " << type;
            break;
        }
    }

    std::tuple<Point, Vector, int> get_offset_projection_and_normal( const Point& p, const size_t& unique_id ) {
        if( iso_value_ == 0 ) {
            ++n_sample_calls_;
            AABB_tree::Point_and_primitive_id papid = tree_.closest_point_and_primitive( p );
            Face_descriptor triangle                = papid.second;
            const Vector normal                     = PMP::compute_face_normal( triangle, primal_mesh_ );
            return { p, normal, 2 };
        }

        if( unique_id < primal_mesh_.number_of_vertices() ) {
            // primitive is a vertex
            Vertex_descriptor v( unique_id );
            const auto s        = Distance::get_primitive_offset( primal_mesh_, v, iso_value_ );
            const Point proj    = s.project( p );
            const Vector normal = normalize( proj - s.origin );
            return { proj, normal, 0 };
        } else if( unique_id < primal_mesh_.number_of_vertices() + primal_mesh_.number_of_edges() ) {
            // primitive is an edge
            Edge_descriptor e( unique_id - primal_mesh_.number_of_vertices() );
            const auto c          = Distance::get_primitive_offset( primal_mesh_, e, iso_value_ );
            const Point proj      = c.project( p );
            const auto proj_plane = c.projection_plane( p );
            const Vector normal( proj_plane[0], proj_plane[1], proj_plane[2] );
            return { proj, normal, 1 };
        } else {
            // primitive is a face
            Face_descriptor f( unique_id - primal_mesh_.number_of_vertices() - primal_mesh_.number_of_edges() );
            const auto plane    = Distance::get_primitive_offset( primal_mesh_, f, iso_value_, p );
            const Point proj    = plane.project( p );
            const Vector normal = plane.normal;
            return { proj, normal, 2 };
        }
    }

    Point qem( const Point& p, const std::vector<Point>& p_proj, const std::vector<Vector>& n_proj ) {
        Eigen::Matrix3d A;
        A.setZero();
        Eigen::Vector3d b;
        b.setZero();
        for( int i = 0; i < p_proj.size(); ++i ) {
            Eigen::Vector3d n_k = { n_proj[i].x(), n_proj[i].y(), n_proj[i].z() };
            Eigen::Vector3d p_k = { p_proj[i].x(), p_proj[i].y(), p_proj[i].z() };
            double d_k          = n_k.transpose() * p_k;

            Eigen::Matrix3d A_k = n_k * n_k.transpose();
            Eigen::Vector3d b_k = d_k * n_k;
            A += A_k;
            b += b_k;
        }

        Eigen::JacobiSVD<Eigen::Matrix3d> svd( A, Eigen::ComputeFullU | Eigen::ComputeFullV );
        svd.setThreshold( 1e-2 );

        // Init x hat
        Eigen::Vector3d x_hat;
        x_hat << p.x(), p.y(), p.z();

        // Lindstrom formula for QEM new position for singular matrices
        Eigen::Vector3d v_svd = x_hat + svd.solve( b - A * x_hat );

        const Point p_new = Point( v_svd[0], v_svd[1], v_svd[2] );

        return p_new;
    }

    Point qem_weighted( const Point& p, const std::vector<Point>& p_proj, const std::vector<Vector>& n_proj, const std::vector<FT>& w ) {
        Eigen::Matrix3d A;
        A.setZero();
        Eigen::Vector3d b;
        b.setZero();
        for( int i = 0; i < p_proj.size(); ++i ) {
            Eigen::Vector3d n_k = { n_proj[i].x(), n_proj[i].y(), n_proj[i].z() };
            Eigen::Vector3d p_k = { p_proj[i].x(), p_proj[i].y(), p_proj[i].z() };
            double d_k          = n_k.transpose() * p_k;

            Eigen::Matrix3d A_k = n_k * n_k.transpose();
            Eigen::Vector3d b_k = d_k * n_k;
            A += w[i] * A_k;
            b += w[i] * b_k;
        }

        Eigen::JacobiSVD<Eigen::Matrix3d> svd( A, Eigen::ComputeFullU | Eigen::ComputeFullV );
        svd.setThreshold( 1e-2 );

        // Init x hat
        Eigen::Vector3d x_hat;
        x_hat << p.x(), p.y(), p.z();

        // Lindstrom formula for QEM new position for singular matrices
        Eigen::Vector3d v_svd = x_hat + svd.solve( b - A * x_hat );

        const Point p_new = Point( v_svd[0], v_svd[1], v_svd[2] );

        return p_new;
    }

    /**
     * Reposition vertex **v** such that its umbrella approximates the offset optimally in the least squares sense.
     * The repositioning is performed with QEM and the SVD pseudo inverse, i.e. if planes are almost complanar, they are assumed as one plane. Planes
     * are constructed from sampling points in the umbrella.
     *
     * @param v vertex
     * @param p initial position that is used in case that eigenvalues in the SVD are smaller than the threshold
     * @return optimized position
     */
    Point project_vertex_back_by_sampling( const Vertex_descriptor& v, const Point& p ) {
        auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
        LOG_ASSERT( f_samples_found );

        bool has_invalid_face = false;

        const Point p0 = mesh_.point( v );
        // sample one-ring
        // std::vector<Point> samples;
        std::vector<Point> p_proj;
        std::vector<Vector> n_proj;
        std::vector<FT> w;
        for( const auto& h: mesh_.halfedges_around_target( mesh_.halfedge( v ) ) ) {
            std::vector<Face_descriptor> fv { mesh_.face( h ) };
            const FT area              = PMP::area( fv, mesh_ );
            const Vertex_descriptor v1 = mesh_.target( mesh_.next( h ) );
            const Vertex_descriptor v2 = mesh_.source( h );
            const Point p1             = mesh_.point( v1 );
            const Point p2             = mesh_.point( v2 );
            const Vector face_normal   = CGAL::cross_product( p1 - p0, p2 - p0 );

            const auto& fs = f_samples[mesh_.face( h )];
            // int n_valid_face_samples = 0;
            for( size_t i = 0; i < fs.points.size(); ++i ) {
                // do not use samples that are oriented completely different than this face
                if( CGAL::scalar_product( face_normal, fs.normals[i] ) <= 0 ) {
                    continue;
                }
                p_proj.push_back( fs.points[i] );
                n_proj.push_back( fs.normals[i] );
                w.push_back( area );
                //++n_valid_face_samples;
            }

            // if( n_valid_face_samples == 0 ) {
            //     has_invalid_face = true;
            //     break;
            // }
        }
        // if( !has_invalid_face ) {
        //     return qem_weighted( p, p_proj, n_proj, w );
        // } else {
        //     return p;
        // }
        if( p_proj.empty() ) {
            return p;
        } else {
            return qem_weighted( p, p_proj, n_proj, w );
        }
        // return qem( p, p_proj, n_proj );
    }

    /** Sample face and return the offset normals at this samples. */
    std::vector<Vector> get_face_normal_samples( const Point& p0, const Point& p1, const Point& p2 ) {
        const Vector p_vec0 = p0 - CGAL::ORIGIN;
        const Vector p_vec1 = p1 - CGAL::ORIGIN;
        const Vector p_vec2 = p2 - CGAL::ORIGIN;

        const Point s0              = CGAL::ORIGIN + ( 0.8 * p_vec0 + 0.1 * p_vec1 + 0.1 * p_vec2 );
        const Point s1              = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.8 * p_vec1 + 0.1 * p_vec2 );
        const Point s2              = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.1 * p_vec1 + 0.8 * p_vec2 );
        const Point s3              = CGAL::centroid( p0, p1, p2 );
        const auto [proj0, n_proj0] = get_offset_projection_and_normal( s0 );
        const auto [proj1, n_proj1] = get_offset_projection_and_normal( s1 );
        const auto [proj2, n_proj2] = get_offset_projection_and_normal( s2 );
        const auto [proj3, n_proj3] = get_offset_projection_and_normal( s3 );

        return { n_proj0, n_proj1, n_proj2, n_proj3 };
    }

    std::vector<Vector> get_face_normal_samples( const Face_descriptor& f ) {
        const Halfedge_descriptor h = mesh_.halfedge( f );
        const Vertex_descriptor v0  = mesh_.target( h );
        const Vertex_descriptor v1  = mesh_.target( mesh_.next( h ) );
        const Vertex_descriptor v2  = mesh_.source( h );
        const Point p0              = mesh_.point( v0 );
        const Point p1              = mesh_.point( v1 );
        const Point p2              = mesh_.point( v2 );
        return get_face_normal_samples( p0, p1, p2 );
    }

    /**
     * Compute the maximum deviation of the triangle normal to the offset normals in the triangle region.
     * Offset normals are computed by sampling the triangle.
     *
     * @param p0
     * @param p1
     * @param p2
     * @return
     */
    FT face_normal_deviation( const Point& p0, const Point& p1, const Point& p2 ) {
        // get face normal
        const Vector face_normal = normalize( CGAL::cross_product( p1 - p0, p2 - p0 ) );

        // get normals at sample points
        std::vector<Vector> normal_samples = get_face_normal_samples( p0, p1, p2 );
        std::vector<Vector> normal_samples_filtered;
        for( const auto& n: normal_samples ) {
            if( CGAL::scalar_product( face_normal, n ) > 0 ) {
                normal_samples_filtered.push_back( n );
            }
        }

        if( normal_samples_filtered.empty() ) {
            return 180;
        }

        // compute max angle between face normal and sampled normals
        double max_angle = 0;
        for( const auto& n: normal_samples_filtered ) {
            max_angle = CGAL::max( max_angle, CGAL::abs( CGAL::approximate_angle( face_normal, n ) ) );
        }

        return max_angle;
    }

    FT face_normal_deviation( const Face_descriptor& f ) {
        const Halfedge_descriptor h = mesh_.halfedge( f );
        const Vertex_descriptor v0  = mesh_.target( h );
        const Vertex_descriptor v1  = mesh_.target( mesh_.next( h ) );
        const Vertex_descriptor v2  = mesh_.source( h );
        const Point p0              = mesh_.point( v0 );
        const Point p1              = mesh_.point( v1 );
        const Point p2              = mesh_.point( v2 );
        return face_normal_deviation( p0, p1, p2 );
    }

    /**
     * Compute the normal deviation along the edge between the faces f(0,1,2) and f(0,3,2).
     * Sampling is performed in close proximity of the edge.
     *
     * @param p0
     * @param p1
     * @param p2
     * @param p3
     * @return
     */
    FT edge_normal_deviation( const Point& p0, const Point& p1, const Point& p2, const Point& p3 ) {
        const Vector face_normal_1 = CGAL::cross_product( p1 - p0, p2 - p0 );       // f(0,1,2)
        const Vector face_normal_2 = CGAL::cross_product( p3 - p0, p1 - p0 );       // f(0,3,2)
        const Vector edge_normal   = normalize( face_normal_1 + face_normal_2 );    // face area weighted (face normals were not normalized)
        // get samples
        const Vector p_vec0 = p0 - CGAL::ORIGIN;
        const Vector p_vec1 = p1 - CGAL::ORIGIN;
        const Vector p_vec2 = p2 - CGAL::ORIGIN;
        const Vector p_vec3 = p3 - CGAL::ORIGIN;

        const Point s_f1_0                  = CGAL::ORIGIN + ( 0.45 * p_vec0 + 0.45 * p_vec1 + 0.1 * p_vec2 );    // mid
        const Point s_f1_1                  = CGAL::ORIGIN + ( 0.8 * p_vec0 + 0.1 * p_vec1 + 0.1 * p_vec2 );      // close to p0
        const Point s_f1_2                  = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.8 * p_vec1 + 0.1 * p_vec2 );      // close to p1
        const auto [proj_f1_0, n_proj_f1_0] = get_offset_projection_and_normal( s_f1_0 );
        const auto [proj_f1_1, n_proj_f1_1] = get_offset_projection_and_normal( s_f1_1 );
        const auto [proj_f1_2, n_proj_f1_2] = get_offset_projection_and_normal( s_f1_2 );
        const Point s_f2_0                  = CGAL::ORIGIN + ( 0.45 * p_vec0 + 0.45 * p_vec1 + 0.1 * p_vec3 );    // mid
        const Point s_f2_1                  = CGAL::ORIGIN + ( 0.8 * p_vec0 + 0.1 * p_vec1 + 0.1 * p_vec3 );      // close to p0
        const Point s_f2_2                  = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.8 * p_vec1 + 0.1 * p_vec3 );      // close to p1
        const auto [proj_f2_0, n_proj_f2_0] = get_offset_projection_and_normal( s_f2_0 );
        const auto [proj_f2_1, n_proj_f2_1] = get_offset_projection_and_normal( s_f2_1 );
        const auto [proj_f2_2, n_proj_f2_2] = get_offset_projection_and_normal( s_f2_2 );
        std::vector<Vector> normal_samples { n_proj_f1_0, n_proj_f1_1, n_proj_f1_2, n_proj_f2_0, n_proj_f2_1, n_proj_f2_2 };

        // get normals at sample points
        std::vector<Vector> normal_samples_filtered;
        for( const auto& n: normal_samples ) {
            if( CGAL::scalar_product( edge_normal, n ) > 0 ) {
                normal_samples_filtered.push_back( n );
            }
        }

        if( normal_samples_filtered.empty() ) {
            return 180;
        }

        // compute max angle between face normal and sampled normals
        double max_angle = 0;
        for( const auto& n: normal_samples_filtered ) {
            max_angle = CGAL::max( max_angle, CGAL::abs( CGAL::approximate_angle( edge_normal, n ) ) );
        }

        return max_angle;
    }

    FT edge_normal_deviation( const Edge_descriptor& e ) {
        const Halfedge_descriptor h = mesh_.halfedge( e );
        const Vertex_descriptor v0  = mesh_.source( h );
        const Vertex_descriptor v1  = mesh_.target( h );
        const Vertex_descriptor v2  = mesh_.target( mesh_.next( h ) );
        const Vertex_descriptor v3  = mesh_.target( mesh_.next( mesh_.opposite( h ) ) );
        const Point p0              = mesh_.point( v0 );
        const Point p1              = mesh_.point( v1 );
        const Point p2              = mesh_.point( v2 );
        const Point p3              = mesh_.point( v3 );
        return edge_normal_deviation( p0, p1, p2, p3 );
    }

    FT max_face_normal_deviation_around_vertex( const Vertex_descriptor& v ) {
        FT nd = 0;
        for( const auto& f: mesh_.faces_around_target( mesh_.halfedge( v ) ) ) {
            nd = CGAL::max( nd, face_normal_deviation( f ) );
        }
        return nd;
    }

    float mean_ratio_metric_with_sampling( const Face_descriptor& f ) {
        const Halfedge_descriptor h = mesh_.halfedge( f );
        const Point p0              = mesh_.point( mesh_.source( h ) );
        const Point p1              = mesh_.point( mesh_.target( h ) );
        const Point p2              = mesh_.point( mesh_.target( mesh_.next( h ) ) );
        const Point p_center        = CGAL::centroid( p0, p1, p2 );
        const auto [p_proj, n]      = get_offset_projection_and_normal( p_center );
        const FT q                  = meanRatioMetric( { p0, p1, p2 }, n );
        return q;
    }

    float min_mean_ratio_metric_around_vertex_with_sampling( const Vertex_descriptor& v ) {
        const Point p0 = mesh_.point( v );
        // const auto [p_proj, n] = get_offset_projection_and_normal( p0 );
        FT min_quality = FT_MAX;
        for( const auto& h: mesh_.halfedges_around_target( mesh_.halfedge( v ) ) ) {
            if( mesh_.is_border( h ) ) {
                continue;
            }
            const Point p1         = mesh_.point( mesh_.target( mesh_.next( h ) ) );
            const Point p2         = mesh_.point( mesh_.target( mesh_.next( mesh_.next( h ) ) ) );
            const Point p_center   = CGAL::centroid( p0, p1, p2 );
            const auto [p_proj, n] = get_offset_projection_and_normal( p_center );
            const FT q             = meanRatioMetric( { p0, p1, p2 }, n );
            min_quality            = CGAL::min( min_quality, q );
        }

        return min_quality;
    }

    void update_triangle_sampling( const Face_descriptor& f ) {
        auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
        LOG_ASSERT( f_samples_found );
        auto& samples = f_samples[f];

        const Halfedge_descriptor h = mesh_.halfedge( f );
        const Vertex_descriptor v0  = mesh_.target( h );
        const Vertex_descriptor v1  = mesh_.target( mesh_.next( h ) );
        const Vertex_descriptor v2  = mesh_.source( h );
        const Point p0              = mesh_.point( v0 );
        const Point p1              = mesh_.point( v1 );
        const Point p2              = mesh_.point( v2 );

        const Vector p_vec0 = p0 - CGAL::ORIGIN;
        const Vector p_vec1 = p1 - CGAL::ORIGIN;
        const Vector p_vec2 = p2 - CGAL::ORIGIN;

        const Point s0 = CGAL::ORIGIN + ( 0.8 * p_vec0 + 0.1 * p_vec1 + 0.1 * p_vec2 );
        const Point s1 = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.8 * p_vec1 + 0.1 * p_vec2 );
        const Point s2 = CGAL::ORIGIN + ( 0.1 * p_vec0 + 0.1 * p_vec1 + 0.8 * p_vec2 );
        const Point s3 = CGAL::centroid( p0, p1, p2 );

        if( samples.has_one_primitive ) {
            // if this and all surrounding triangles have the same primitive, do not use AABB_tree
            const size_t& id = samples.primitive_id;
            std::set<Face_descriptor> neighs;
            for( const auto& v: mesh_.vertices_around_face( mesh_.halfedge( f ) ) ) {
                for( const auto& ff: mesh_.faces_around_target( mesh_.halfedge( v ) ) ) {
                    neighs.insert( ff );
                }
            }
            bool are_all_the_same = true;
            for( const auto& ff: neighs ) {
                if( !f_samples[ff].has_one_primitive || f_samples[ff].primitive_id != id ) {
                    are_all_the_same = false;
                    break;
                }
            }
            if( are_all_the_same ) {
                const auto [proj0, n_proj0, t0] = get_offset_projection_and_normal( s0, id );
                const auto [proj1, n_proj1, t1] = get_offset_projection_and_normal( s1, id );
                const auto [proj2, n_proj2, t2] = get_offset_projection_and_normal( s2, id );
                const auto [proj3, n_proj3, t3] = get_offset_projection_and_normal( s3, id );
                samples.points                  = { proj0, proj1, proj2, proj3 };
                samples.normals                 = { n_proj0, n_proj1, n_proj2, n_proj3 };
                samples.has_one_primitive       = true;
                samples.primitive_id            = id;
                samples.needs_update            = false;
                return;
            }
        }

        const auto [proj0, n_proj0] = get_offset_projection_and_normal( s0 );
        const auto [proj1, n_proj1] = get_offset_projection_and_normal( s1 );
        const auto [proj2, n_proj2] = get_offset_projection_and_normal( s2 );
        const auto [proj3, n_proj3] = get_offset_projection_and_normal( s3 );
        LOG_ASSERT( n_proj0.squared_length() > 0.5 && n_proj0.squared_length() < 1.5 );
        LOG_ASSERT( n_proj1.squared_length() > 0.5 && n_proj1.squared_length() < 1.5 );
        LOG_ASSERT( n_proj2.squared_length() > 0.5 && n_proj2.squared_length() < 1.5 );
        LOG_ASSERT( n_proj3.squared_length() > 0.5 && n_proj3.squared_length() < 1.5 );

        samples.points  = { proj0, proj1, proj2, proj3 };
        samples.normals = { n_proj0, n_proj1, n_proj2, n_proj3 };
        // samples.is_plane = { t0 == 2, t1 == 2, t2 == 2, t3 == 2 };
        //
        // if( id0 == id1 && id0 == id2 && id0 == id3 ) {
        //     samples.has_one_primitive = true;
        //     samples.primitive_id      = id0;
        // } else {
        //     samples.has_one_primitive = false;
        //     samples.primitive_id      = INVALID_UNIQUE_ID;
        // }

        samples.needs_update = false;
    }

    void update_triangle_normal( const Face_descriptor& f ) {
        auto [f_normals, f_normals_found] = mesh_.property_map<Face_descriptor, Vector>( "f:normal" );
        LOG_ASSERT( f_normals_found );
        f_normals[f] = PMP::compute_face_normal( f, mesh_ );
    }

    double normal_deviation_sampled( const Face_descriptor& f ) {
        auto [f_normals, f_normals_found] = mesh_.property_map<Face_descriptor, Vector>( "f:normal" );
        LOG_ASSERT( f_normals_found );

        auto [f_samples, f_samples_found] = mesh_.property_map<Face_descriptor, Triangle_sampling>( "f:samples" );
        LOG_ASSERT( f_samples_found );

        if( f_samples[f].needs_update ) {
            update_triangle_sampling( f );
        }
        // update_triangle_sampling( f );

        // get face normal
        const Vector& face_normal = f_normals[f];

        // get normals at sample points
        std::array<Vector, 4> normal_samples = f_samples[f].normals;
        std::vector<Vector> normal_samples_filtered;
        for( const auto& n: normal_samples ) {
            if( CGAL::scalar_product( face_normal, n ) > 0 ) {
                normal_samples_filtered.push_back( n );
            }
        }

        if( normal_samples_filtered.empty() ) {
            return 180;
        }

        // compute max angle between face normal and sampled normals
        double max_angle = 0;
        for( const auto& n: normal_samples_filtered ) {
            max_angle = CGAL::max( max_angle, CGAL::approximate_angle( face_normal, n ) );
        }

        return max_angle;
    }

    void update_last_operation( const Face_descriptor& f ) {
        auto [f_last_operation, f_last_operation_found] = mesh_.property_map<Face_descriptor, size_t>( "f:last_operation" );
        LOG_ASSERT( f_last_operation_found );

        f_last_operation[f] = n_remeshing_iterations_;
    }

    bool is_operation_required( const Halfedge_descriptor& h, const Remeshing_operation_types& t ) {
        auto [f_last_operation, f_last_operation_found] = mesh_.property_map<Face_descriptor, size_t>( "f:last_operation" );
        LOG_ASSERT( f_last_operation_found );

        std::vector<Face_descriptor> faces;
        switch( t ) {
        case Remeshing_operation_types::split:
            faces.push_back( mesh_.face( h ) );
        case Remeshing_operation_types::collapse:
            faces.push_back( mesh_.face( h ) );
            faces.push_back( mesh_.face( mesh_.opposite( h ) ) );
            break;
        case Remeshing_operation_types::flip:
            faces.push_back( mesh_.face( h ) );
            faces.push_back( mesh_.face( mesh_.opposite( h ) ) );
            break;
        case Remeshing_operation_types::smooth:
            for( const auto& f: mesh_.faces_around_target( h ) ) {
                faces.push_back( f );
            }
            break;
        default:
            LOG( ERROR ) << "Unknown operation type";
            break;
        }

        if( faces.empty() ) {
            return true;
        }

        for( const auto& f: faces ) {
            if( n_remeshing_iterations_ - f_last_operation[f] > 1 ) {
                // last OP is more than a full iteration ago
                continue;
            }
            // triangle was recently modified --> operation could be necessary
            return true;
        }
        return false;
    }
};