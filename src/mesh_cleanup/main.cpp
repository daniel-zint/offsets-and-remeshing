#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>
#include <CLI/CLI.hpp>
#include <filesystem>
#include <glog/logging.h>
#include <iostream>

#define TAG CGAL::Parallel_if_available_tag

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<size_t>> Polygon_range;

typedef CGAL::Surface_mesh<Point> SurfaceMesh;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor Vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor Halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor Edge_descriptor;
typedef SurfaceMesh::Face_index Face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = std::filesystem;

inline void normalize_mesh(SurfaceMesh& m) {
    CGAL::Bbox_3 aabb = PMP::bbox(m);

    const FT d = CGAL::max(aabb.x_span(), CGAL::max(aabb.y_span(), aabb.z_span()));

    for (const auto& v : m.vertices()) {
        Point& p = m.point(v);
        const FT x = (p.x() - aabb.xmin()) / d;
        const FT y = (p.y() - aabb.ymin()) / d;
        const FT z = (p.z() - aabb.zmin()) / d;
        p = Point(x, y, z);
    }
}

int main(int argc, char* argv[]) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = true;
    FLAGS_colorlogtostderr = true;
    FLAGS_minloglevel = 0;    // Set log level: INFO = 0, WARNING = 1, ERROR = 2

    CLI::App app{ "App description" };

    fs::path input_path, output_path;
    bool use_normalize_mesh = false;
    app.add_option("-i", input_path, "input mesh or folder")->required();
    app.add_option("-o", output_path, "output mesh or folder")->required();

    app.add_flag("-n,--normalize", use_normalize_mesh, "Transform mesh to the unit cube");

    CLI11_PARSE(app, argc, argv);

    if (!fs::exists(input_path)) {
        LOG(ERROR) << "The input path does not exist.\n"
            << "Input: " << fs::absolute(input_path);
        return 0;
    }

    LOG_IF(WARNING, fs::exists(output_path)) << "The output path already exists and will be overwritten.\n"
        << "Output: " << fs::absolute(output_path);

    std::vector<fs::path> input_files;
    if (fs::is_regular_file(input_path)) {
        input_files.push_back(input_path);
    }
    else if (fs::is_directory(input_path)) {
        fs::create_directories(output_path);
        for (const auto& f : fs::directory_iterator(input_path)) {
            if (fs::is_regular_file(f)) {
                input_files.push_back(f);
            }
        }
        LOG_IF(WARNING, input_files.empty()) << "No input files";
    }

    if (fs::is_directory(output_path)) {
        for (size_t i = 0; i < input_files.size(); ++i) {
            std::string subfolder = std::to_string(i / 10);
            fs::create_directories(output_path / subfolder);
        }
    }

    std::vector<fs::path> failed_files;
    for (size_t i = 0; i < input_files.size(); ++i) {
        const auto& f = input_files[i];

        if (fs::is_directory(output_path)) {
            std::string subfolder = std::to_string(i / 10);
            fs::path output_file = output_path / subfolder / (f.stem().string() + ".off");
            if (fs::exists(output_file)) {
                LOG(INFO) << i << "/" << input_files.size() << " Exists already: " << output_file;
                continue;
            }
        }

        // read file
        SurfaceMesh m;
        bool success = PMP::IO::read_polygon_mesh(f.string(), m);
        if (!success) {
            Point_range p;
            Polygon_range pol;
            success = CGAL::IO::read_polygon_soup(f.string(), p, pol);

            if (!success) {
                LOG(WARNING) << "Mesh could not be read or repaired\n"
                    << "Mesh: " << fs::absolute(f);
                failed_files.push_back(f);
                continue;
            }

            PMP::repair_polygon_soup(p, pol);
            PMP::orient_polygon_soup(p, pol);
            m.clear();
            PMP::polygon_soup_to_polygon_mesh(p, pol, m);
        }

        if (!CGAL::is_valid_polygon_mesh(m)) {
            LOG(WARNING) << "Mesh is not a valid polygon mesh\n"
                << "Mesh: " << fs::absolute(f);
            failed_files.push_back(f);
            continue;
        }

        if (use_normalize_mesh) {
            normalize_mesh(m);
        }
        if (!PMP::triangulate_faces(m)) {
            LOG(WARNING) << "Mesh could not be triangulated\n"
                << "Mesh: " << fs::absolute(f);
            failed_files.push_back(f);
            continue;
        }

        if (fs::is_directory(output_path)) {
            std::string subfolder = std::to_string(i / 10);
            // create an off file with the same name as the input
            fs::path output_file = output_path / subfolder / (f.stem().string() + ".off");
            LOG(INFO) << i << "/" << input_files.size() << " Write " << output_file;
            CGAL::IO::write_OFF(output_file.string(), m);
        }
        else {
            // assume output_path is the desired filename
            LOG(INFO) << i << "/" << input_files.size() << " Write " << output_path;
            CGAL::IO::write_OFF(output_path.string(), m);
        }
    }

    if (!failed_files.empty()) {
        LOG(WARNING) << "Some meshes could not be cleaned up:";
        for (const auto& f : failed_files) {
            std::cout << fs::absolute(f) << std::endl;
        }
    }
}