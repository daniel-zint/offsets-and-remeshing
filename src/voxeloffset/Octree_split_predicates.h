#pragma once

#include "Octree_wrapper.h"
#include "Dual_contouring_octree_3.h"

namespace Split_predicates {

    // Refine 1/8 of the octree. This split predicate is mostly for debugging purposes
    struct Refine_one_eighth {
        std::size_t min_depth_;
        std::size_t max_depth_;

        std::size_t octree_dim_;

        Octree_wrapper::Uniform_coords uniform_coordinates( const Octree_wrapper::Octree::Node& node ) const {
            auto coords                    = node.global_coordinates();
            const std::size_t depth_factor = std::size_t( 1 ) << ( max_depth_ - node.depth() );
            for( int i = 0; i < Octree_wrapper::Octree::Node::Dimension::value; ++i ) {
                coords[i] *= depth_factor;
            }

            return coords;
        }

        Refine_one_eighth( std::size_t min_depth, std::size_t max_depth ) : min_depth_( min_depth ), max_depth_( max_depth ) {
            octree_dim_ = std::size_t( 1 ) << max_depth_;
        }

        bool operator()( const Octree_wrapper::Octree::Node& n ) const {
            // n.depth()
            if( n.depth() < min_depth_ ) {
                return true;
            }
            if( n.depth() == max_depth_ ) {
                return false;
            }

            auto leaf_coords = uniform_coordinates( n );

            if( leaf_coords[0] >= octree_dim_ / 2 ) {
                return false;
            }
            if( leaf_coords[1] >= octree_dim_ / 2 ) {
                return false;
            }
            if( leaf_coords[2] >= octree_dim_ / 2 ) {
                return false;
            }
            return true;
        }
    };

    struct Split_uniformly {
        std::size_t depth_;

        Split_uniformly( const std::size_t& depth ) : depth_( depth ) {}

        bool operator()( const Octree_wrapper::Octree::Node& n ) const {
            if( n.depth() < depth_ ) {
                return true;
            }
            return false;
        }
    };
}    // namespace Split_predicates