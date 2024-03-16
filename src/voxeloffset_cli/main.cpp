#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Real_timer.h>
#include <CLI/CLI.hpp>
#include <filesystem>
#include <fstream>
#include <glog/logging.h>
#include <iostream>
#include <voxeloffset/Dual_contouring_offset.h>
#include <voxeloffset/Quality_measurements.h>
#include <voxeloffset/Remeshing.h>
#include <voxeloffset/types.h>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = true;
    FLAGS_colorlogtostderr = true;
    FLAGS_minloglevel = 0;    // Set log level: INFO = 0, WARNING = 1, ERROR = 2

    CLI::App app{ "App description" };

    FT offset_radius = 0.1;
    fs::path primal_path, output_path, remesh_input, debug_output_folder, debug_output_folder_dc, debug_output_folder_remeshing;
    size_t n_remeshing_iterations = 20;
    size_t d0 = 0, d1 = 10, d2 = 12;
    bool run_dual_contouring = false;
    FT l_min = FT_MAX;
    FT l_max = -FT_MAX;
    FT l_min_factor = 0.1;
    FT l_max_factor = 2.0;
    FT max_normal_deviation = 7;
    bool use_normalize_mesh = false;
    bool perform_quality_tests = false;
    FT diag = FT_MAX;

    app.add_option("-j,--radius", offset_radius, "The offset radius/distance to the primal mesh");
    app.add_option("--diag", diag, "Set offset radius relative to the bbox diagonal");
    app.add_option("-p,--primal", primal_path, "primal mesh")->required();
    app.add_option("-o,--output", output_path, "output mesh")->required();

    // DC options
    app.add_flag("-d,--dualcontouring", run_dual_contouring, "Dual Contouring should be applied");
    app.add_option("--d0", d0, "minimum octree depth")->expected(0, 20);
    app.add_option("--d1", d1, "standard maximum octree depth")->expected(0, 20);
    app.add_option("--d2", d2, "octree depth for resolving non-manifold vertices")->expected(0, 20);

    // remeshing options
    app.add_option("-r,--remesh", n_remeshing_iterations, "number of remeshing iterations");
    app.add_option("--remeshinput", remesh_input, "remesh input mesh");
    app.add_option("--nd", max_normal_deviation, "maximum normal deviation")->expected(0, 45);
    app.add_option("--lminfactor", l_min_factor, "factor for minimal edge length");
    app.add_option("--lmaxfactor", l_max_factor, "factor for maximal edge length");
    app.add_option("--lmin", l_min, "minimal edge length (use if j=0, overrides l_min_factor)");
    app.add_option("--lmax", l_max, "maximal edge length (use if j=0, overrides l_max_factor)");

    // debug output
    app.add_option("--debugout", debug_output_folder, "folder for all intermediate results");
    app.add_option("--doutd", debug_output_folder_dc, "folder for all intermediate results in DC");
    app.add_option("--doutr", debug_output_folder_remeshing, "folder for all intermediate results in remeshing");

    // normalize input to [0,1] bbox
    app.add_flag("-n,--normalize", use_normalize_mesh, "Transform mesh to the unit cube");

    app.add_flag("-q", perform_quality_tests, "Perform quality tests");

    CLI11_PARSE(app, argc, argv);

    if (fs::exists(debug_output_folder)) {
        if (!fs::exists(debug_output_folder_dc)) {
            debug_output_folder_dc = debug_output_folder;
        }
        if (!fs::exists(debug_output_folder_remeshing)) {
            debug_output_folder_remeshing = debug_output_folder;
        }
    }

    // print settings
    std::cout << "Primal: " << primal_path << std::endl;
    std::cout << "Output: " << output_path << std::endl;
    std::cout << "Offset radius: " << offset_radius << std::endl;
    std::cout << "Normalize primal: " << use_normalize_mesh << std::endl;
    std::cout << "Run DC: " << run_dual_contouring << std::endl;
    if (run_dual_contouring) {
        std::cout << "  d0: " << d0 << std::endl;
        std::cout << "  d1: " << d1 << std::endl;
        std::cout << "  d2: " << d2 << std::endl;
    }
    std::cout << "Remeshing iterations: " << n_remeshing_iterations << std::endl;
    if (n_remeshing_iterations > 0) {
        if (!run_dual_contouring) {
            std::cout << " Remeshing input: " << remesh_input << std::endl;
        }
        std::cout << " Max normal deviation: " << max_normal_deviation << std::endl;
        std::cout << " l_min_factor: " << l_min_factor << std::endl;
        std::cout << " l_max_factor: " << l_max_factor << std::endl;
    }
    std::cout << "Test quality: " << perform_quality_tests << std::endl;
    if (fs::exists(debug_output_folder_dc)) {
        std::cout << "Debug output folder DC: " << debug_output_folder_dc << std::endl;
    }
    if (fs::exists(debug_output_folder_remeshing)) {
        std::cout << "Debug output folder remeshing: " << debug_output_folder_remeshing << std::endl;
    }

    if (!fs::exists(primal_path)) {
        LOG(ERROR) << "The primal mesh does not exist.\n"
            << "Input: " << fs::absolute(primal_path);
        return 0;
    }
    if (fs::exists(output_path)) {
        LOG(WARNING) << "The output mesh already exists. Abort.\n"
            << "Output: " << fs::absolute(output_path);
        // return 0;
    }

    SurfaceMesh primal_mesh;
    std::ifstream in(primal_path);
    CGAL::IO::read_OFF(in, primal_mesh);
    LOG_ASSERT(in);

    if (use_normalize_mesh) {
        normalize_mesh(primal_mesh);
    }

    if (diag < FT_MAX) {
        std::cout << "Re-compute offset radius from the bbox diagonal.";
        const CGAL::Bbox_3 bb = PMP::bbox(primal_mesh);
        const FT diag_length = CGAL::sqrt(bb.x_span() * bb.x_span() + bb.y_span() * bb.y_span() + bb.z_span() * bb.z_span());
        offset_radius = diag * diag_length;
        std::cout << "New Offset radius: " << offset_radius << std::endl;
    }

    // triangulate input
    if (!CGAL::is_triangle_mesh(primal_mesh) && !PMP::triangulate_faces(primal_mesh)) {
        LOG(FATAL) << "Could not triangulate input";
    }

    if (fs::exists(debug_output_folder_dc)) {
        print_mesh(primal_mesh, debug_output_folder_dc / "primal.off");
    }

    CGAL::Real_timer timer;
    CGAL::Real_timer timer_complete;
    timer_complete.start();

    SurfaceMesh offset_input_mesh;
    if (run_dual_contouring) {
        std::cout << "===== Octree Initialization =====" << std::endl;
        if (d2 == 0) {
            LOG(ERROR) << "Octree levels (--d0,--d1,--d2) are not specified.";
            return 0;
        }
        if (!(d0 <= d1 && d1 <= d2)) {
            LOG(ERROR) << "The three octree levels must have ascending values: d0 <= d1 <= d2.\n"
                << "d0 = " << d0 << "\n"
                << "d1 = " << d1 << "\n"
                << "d2 = " << d2 << "\n";
            return 0;
        }
        if (!remesh_input.string().empty()) {
            LOG(WARNING) << "Remesh input mesh will be ignored.\n"
                << "Remesh input: " << fs::absolute(remesh_input);
        }

        timer.start();
        Dual_contouring_offset dual_contouring_offset(primal_mesh, d0, d1, d2, offset_radius);
        timer.stop();
        std::cout << "===== Octree initialization took " << timer.time() << " seconds =====" << std::endl;
        timer.reset();

        std::cout << "===== Dual Contouring =====" << std::endl;
        if (fs::exists(debug_output_folder_dc)) {
            dual_contouring_offset.output_path() = debug_output_folder_dc;
        }
        timer.start();
        dual_contouring_offset.compute_offset(primal_mesh);
        timer.stop();
        std::cout << "===== Dual Contouring took " << timer.time() << " seconds =====" << std::endl;
        timer.reset();

        offset_input_mesh = dual_contouring_offset.mesh();

        if (offset_radius < 0) {
            // reverse orientation if mesh is on the inside
            PMP::reverse_face_orientations(offset_input_mesh);
        }

    }
    else {
        if (remesh_input.string().empty() || !fs::exists(remesh_input)) {
            LOG(ERROR) << "The remesh input mesh does not exist.\n"
                << "Remesh input: " << fs::absolute(remesh_input);
            return 0;
        }
        std::ifstream in(remesh_input);
        CGAL::IO::read_OFF(in, offset_input_mesh);
        LOG_ASSERT(in);
    }

    if (CGAL::is_closed(offset_input_mesh)) {
        LOG(INFO) << "Remeshing input is closed";
    }
    else {
        LOG(WARNING) << "Remeshing input contains holes";
    }

    // remeshing
    Remeshing remesher(offset_input_mesh, primal_mesh, CGAL::abs(offset_radius));
    if (fs::exists(debug_output_folder_remeshing)) {
        remesher.output_path() = debug_output_folder_remeshing;
    }
    if (n_remeshing_iterations > 0) {
        std::cout << "===== Remeshing =====" << std::endl;
        remesher.set_max_normal_deviation(max_normal_deviation);
        remesher.l_max() *= l_max_factor;
        remesher.l_min() *= l_min_factor;
        if (l_max != -FT_MAX) {
            remesher.l_max() = l_max;
        }
        if (l_min != FT_MAX) {
            remesher.l_min() = l_min;
        }
        timer.start();
        remesher.remesh(n_remeshing_iterations);
        timer.stop();
        std::cout << "===== Remeshing took " << timer.time() << " seconds =====" << std::endl;
        timer.reset();
    }

    SurfaceMesh offset = remesher.mesh();
    timer_complete.stop();
    print_mesh(offset, output_path);

    // quality tests
    if (perform_quality_tests) {
        std::cout << "===== Quality Tests =====" << std::endl;
        timer.start();
        LOG(INFO) << "Check for self intersections: ";
        if (Quality_measurements::does_self_intersect(offset)) {
            LOG(WARNING) << "Mesh has self intersections";
        }
        else {
            LOG(INFO) << "No self intersections";
        }

        LOG(INFO) << "Normal deviation:";
        Quality_measurements::normal_deviation(offset, primal_mesh, CGAL::abs(offset_radius));

        auto distances = Quality_measurements::distance_to_real_offset(offset, primal_mesh, offset_radius);
        Quality_measurements::analyze_distance(offset, primal_mesh, offset_radius, distances);

        timer.stop();
        std::cout << "===== Quality tests took " << timer.time() << " seconds =====" << std::endl;
    }

    std::cout << "===== Full offset generation took " << timer_complete.time() << " seconds for " << primal_mesh.number_of_faces()
        << " input triangles and " << offset.number_of_faces() << " output triangles =====" << std::endl;

    LOG(INFO) << "Program ended";
}