#pragma once

#include "Octree_wrapper.h"
#include "types.h"

#include <filesystem>
#include <memory>

class Dual_contouring_offset {
    std::shared_ptr<Octree_wrapper> octree_;    // The underlying octree
    SurfaceMesh mesh_;                          // The dual contouring mesh

    std::size_t depth_     = 5;
    double iso_value_      = 0.1;      // TODO dummy values should be replaced by more generic values
    std::size_t mid_depth_ = 5 + 2;    // Mid depth level of the octree
    std::size_t max_depth_ = 5 + 3;    // Maximum depth level of the octree (only used on cells that would produce non-manifold vertices otherwise)

    std::filesystem::path output_path_;

  public:
    // constructors
    Dual_contouring_offset( const SurfaceMesh& primal_mesh, std::size_t depth, std::size_t mid_depth, std::size_t max_depth, double iso_value )
        : depth_ { depth }, mid_depth_ { mid_depth }, max_depth_ { max_depth }, iso_value_ { iso_value } {
        init_octree( primal_mesh );
    }

    /// <summary>
    /// Compute the Dual Contouring mesh and simplify it to the meta mesh using the primal element lists that are stored in the vertices.
    /// </summary>
    void compute_offset( const SurfaceMesh& primal_mesh );

    // access functions
    auto& mesh() { return mesh_; }
    auto& output_path() { return output_path_; }

    void print_distance_metrics( const SurfaceMesh& primal_mesh );

  private:
    void init_octree( const SurfaceMesh& primal_mesh );

    /// <summary>
    /// Perform Dual Contouring on the voxel grid. Vertices are placed in the center of the voxel. Non-manifoldness is resolved by duplicating
    /// vertices. Additionally, all vertices have the list of primal elements.
    /// </summary>
    void compute_mesh_from_octree( const SurfaceMesh& primal_mesh );

    void reposition_vertices( const SurfaceMesh& primal_mesh );

    /// <summary>
    /// Flip edges in the triangulated Dual Contouring mesh.
    /// </summary>
    void triangulate();
};
