#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Real_timer.h>
#include <CGAL/boost/graph/helpers.h>
#include <filesystem>
#include <fstream>
#include <glog/logging.h>
#include <voxeloffset/Dual_contouring_offset.h>
#include <voxeloffset/Quality_measurements.h>
#include <voxeloffset/Remeshing.h>
#include <voxeloffset/types.h>

std::filesystem::path SOURCE_PATH = SOURCE_DIR;
std::filesystem::path MESH_DIR = std::filesystem::canonical(SOURCE_PATH / "data");

int main(int argc, char* argv[]) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = true;
    FLAGS_colorlogtostderr = true;
    FLAGS_minloglevel = 0;    // Set log level: INFO = 0, WARNING = 1, ERROR = 2

    CGAL::Real_timer timer;
    timer.start();

    const size_t depth = 5;
    const size_t mid_depth = depth + 5;
    const size_t max_depth = depth + 7;
    const double iso_value = 0.5;
    // load some mesh
    SurfaceMesh primal_mesh;
    std::ifstream in(MESH_DIR / "cube.off");
    CGAL::IO::read_OFF(in, primal_mesh);
    LOG_ASSERT(in);

    print_mesh(primal_mesh, MESH_DIR / "primal.off");

    timer.stop();
    LOG(WARNING) << "Runtime init: " << timer.time() << " seconds";
    timer.reset();

    // Dual Contouring //
    std::cout << "===== Dual Contouring =====" << std::endl;
    timer.start();
    Dual_contouring_offset dual_contouring_offset(primal_mesh, depth, mid_depth, max_depth, iso_value);    // inner
    timer.stop();
    LOG(WARNING) << "Runtime Dual Contouring: " << timer.time() << " seconds";
    timer.reset();
    // offset
    dual_contouring_offset.output_path() = MESH_DIR;
    dual_contouring_offset.compute_offset(primal_mesh);

    SurfaceMesh offset_init = dual_contouring_offset.mesh();
    print_mesh(offset_init, MESH_DIR / "triangulated.off");

    if (iso_value < 0) {
        // reverse orientation if mesh is on the inside
        PMP::reverse_face_orientations(offset_init);
    }

    // Post Processing //
    std::cout << "===== Remeshing =====" << std::endl;
    timer.start();
    Remeshing remesher(offset_init, primal_mesh, CGAL::abs(iso_value));
    remesher.set_max_normal_deviation(7);
    remesher.l_max() *= 2;
    remesher.l_min() *= 0.1;
    remesher.remesh(5);
    timer.stop();
    LOG(WARNING) << "Runtime remeshing: " << timer.time() << " seconds";
    timer.reset();

    SurfaceMesh offset = remesher.mesh();
    print_mesh(offset, MESH_DIR / "postprocessed.off");

    LOG(INFO) << "Number of sample calls: " << remesher.n_sample_calls();

    // Quality Tests //
    std::cout << "===== Quality Tests =====" << std::endl;
    timer.start();
    LOG(INFO) << "Check for self intersections: ";
    if (Quality_measurements::does_self_intersect(offset)) {
        LOG(WARNING) << "Mesh has self intersections";
    }
    else {
        LOG(INFO) << "No self intersections";
    }

    LOG(INFO) << "Hausdorff distance: ";
    double hausdorff = Quality_measurements::hausdorff(offset, primal_mesh);
    LOG(INFO) << "Hausdorff = " << hausdorff << ". Offset value = " << CGAL::abs(iso_value);

    LOG(INFO) << "Normal deviation:";
    Quality_measurements::normal_deviation(offset, primal_mesh, CGAL::abs(iso_value));

    timer.stop();
    LOG(WARNING) << "Runtime quality tests: " << timer.time() << " seconds";
}