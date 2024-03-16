#include "types.h"

/**
 * collapses an edge in a graph that consists of triangles and other polygons.
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 * Let `h` be the halfedge of `e`, and let `v0` and `v1` be the source and target vertices of `h`.
 * Let `p_h` and `p_o_h` be respectively the edges of `prev(h,g)` and `prev(opposite(h, g), g)`.
 * Let `o_n_h` and `o_n_o_h` be respectively the edges of `opposite(next(h,g))` and `opposite(next(opposite(h, g), g))`.
 *
 * After the collapse of edge `e` the following holds:
 *   - The edge `e` is no longer in `g`.
 *   - The faces incident to edge `e` are no longer in `g`.
 *   - `v0` is no longer in `g`.
 *   - If `h` is not a border halfedge, `p_h` is no longer in `g` and is replaced by `o_n_h`.
 *   - If the opposite of `h` is not a border halfedge, `p_o_h` is no longer in `g` and is replaced by `o_n_o_h`.
 *   - The halfedges kept in `g` that had `v0` as target and source now have `v1` as target and source, respectively.
 *   - No other incidence information is changed in `g`.
 *
 * \returns vertex `v1`.
 * \pre g must be a triangulated graph
 * \pre `does_satisfy_link_condition(e,g) == true`.
 */
template<typename Graph>
typename boost::graph_traits<Graph>::vertex_descriptor collapse_hybrid_edge( typename boost::graph_traits<Graph>::halfedge_descriptor h, Graph& g ) {
    typedef boost::graph_traits<Graph> Traits;
    typedef typename Traits::vertex_descriptor vertex_descriptor;
    typedef typename Traits::halfedge_descriptor halfedge_descriptor;

    const halfedge_descriptor hop       = opposite( h, g );
    const halfedge_descriptor hnext     = next( hop, g );
    const halfedge_descriptor gnext     = next( h, g );
    const halfedge_descriptor hprev     = prev( hop, g );
    const halfedge_descriptor gprev     = prev( h, g );
    const halfedge_descriptor gprev_opp = opposite( gprev, g );
    const halfedge_descriptor hprev_opp = opposite( hprev, g );

    const bool h_face_is_triangle    = next( gnext, g ) == gprev;
    const bool hop_face_is_triangle  = next( hnext, g ) == hprev;
    const bool h_face_exists         = !is_border( h, g );
    const bool hop_face_exists       = !is_border( hop, g );
    const bool gprev_opp_face_exists = h_face_exists && !is_border( gprev_opp, g );
    const bool hprev_opp_face_exists = hop_face_exists && !is_border( hprev_opp, g );

    // CGAL_precondition( !h_face_exists || ( h_face_exists && ( degree( target( gprev_opp, g ), g ) > 2 ) ) );
    // CGAL_precondition( !hop_face_exists || ( hop_face_exists && ( degree( target( hprev_opp, g ), g ) > 2 ) ) );

    vertex_descriptor v1 = target( h, g );
    vertex_descriptor v0 = source( h, g );

    if( !CGAL::Euler::does_satisfy_link_condition( edge( h, g ), g ) ) {
        LOG( WARNING ) << "Edge does not satisfy link condition";
    }

    if( degree( v0, g ) == 2 ) {
        CGAL::Euler::remove_center_vertex( hop, g );
        return v1;
    }

    if( degree( v1, g ) == 2 ) {
        CGAL::Euler::remove_center_vertex( h, g );
        return v0;
    }

    bool lP_Erased = false;

    if( h_face_exists && h_face_is_triangle ) {
        CGAL_precondition( !is_border( opposite( gprev_opp, g ), g ) );    // v0-v1-t is a face of the mesh
        if( gprev_opp_face_exists ) {
            // CGAL_ECMS_TRACE(3, "Removing p-t E" << pt.idx() << " (V"
            //                << p.idx() << "->V" << target(pt, g).idx()
            //                << ") by joining top-left face" ) ;

            CGAL::Euler::join_face( gprev_opp, g );
        } else {
            // CGAL_ECMS_TRACE(3, "Removing p-t E" << pt.idx() << " (V" << p.idx()
            //                << "->V" << target(pt, g).idx() << ") by erasing top face" ) ;

            CGAL::Euler::remove_face( opposite( gprev_opp, g ), g );

            if( !hop_face_exists ) {
                // CGAL_ECMS_TRACE(3, "Bottom face doesn't exist so vertex v0 already removed" ) ;

                lP_Erased = true;
            }
        }
    }

    if( hop_face_exists && hop_face_is_triangle ) {
        CGAL_precondition( !is_border( opposite( hprev_opp, g ), g ) );    // v0-v1-b is a face of the mesh
        if( hprev_opp_face_exists ) {
            // CGAL_ECMS_TRACE(3, "Removing q-b E" << qb.idx() << " (V"
            //                << q.idx() << "->V" << target(qb, g).idx()
            //                << ") by joining bottom-right face" ) ;

            CGAL::Euler::join_face( hprev_opp, g );
        } else {
            // CGAL_ECMS_TRACE(3, "Removing q-b E" << qb.idx() << " (V"
            //                << q.idx() << "->V" << target(qb, g).idx()
            //                << ") by erasing bottom face" ) ;

            if( !h_face_exists ) {
                // CGAL_ECMS_TRACE(3, "Top face doesn't exist so vertex v1 already removed" ) ;
                lP_Erased = true;

                // v1 will be removed, swap v0 and v1
                CGAL::internal::swap_vertices( v0, v1, g );
            }

            CGAL::Euler::remove_face( opposite( hprev_opp, g ), g );
        }
    }

    if( !lP_Erased ) {
        // CGAL_ECMS_TRACE(3, "Removing vertex P by joining pQ" ) ;
        CGAL::Euler::join_vertex( h, g );
        lP_Erased = true;
    }

    CGAL_expensive_assertion( is_valid_polygon_mesh( g ) );

    return v1;
}