#include <getopt.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include "configuration/read_config.h"
#include "configuration/write_config.h"
#include "constant.h"
#include "helper.h"
#include "mesh/mesh.h"
#include "spatial/intersection.h"
#include "temporal/douglas_adi.h"
#include "temporal/lod.h"
#include "temporal/trapezoidal.h"
// User-defined surfaces
#include "surface/surface_cartesian.h"
// User-defined diffusion coefficient beta
#include "diffusion/beta.h"
// Users-defined functions
#include "equation/equation.h"

using namespace std;

int main(int argc, char *argv[]) {
  cout << "[INFO]: PROGRAM STARTS ... " << endl;

  auto t_start = chrono::high_resolution_clock::now();

  CubicDoub uh;
  string output_file, input_file;

  const char *short_opts = "i:o:h";
  static struct option long_opts[] = {
      {"input", required_argument, nullptr, 'i'},
      {"output", required_argument, nullptr, 'o'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, no_argument, nullptr, 0},
  };

  int opts;
  while ((opts = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
    switch (opts) {
      case 'i':
        input_file = static_cast<string>(optarg);
        cout << "[INFO]: reading configuration from file: " << static_cast<string>(optarg) << endl;
        break;
      case 'o':
        output_file = static_cast<string>(optarg);
        cout << "[INFO]: writing result to file: " << output_file << endl;
        break;
      case 'h':
        print_help(argv[0]);
        return EXIT_SUCCESS;
      default:
        print_help(argv[0]);
        return EXIT_FAILURE;
    }
  }

  if (input_file.empty() || output_file.empty()) {
    print_help(argv[0]);
    return EXIT_FAILURE;
  }

  // Read data from file
  ReadConfig config(input_file);

  // write configuration informations to result
  ofstream out_stream(output_file, ios::out | ios::app);
  if (!out_stream.good()) {
    LOG_FATAL("I/O failure when opening " + output_file);
  }
  write_configuration(out_stream, config);

  VecDoub domain = config.GetDomain();
  VecInt size = config.GetMesh();
  VecDoub time_info = config.GetTime();
  int diffusion_coef_type = config.GetDiffusionCoef();
  int equation_type = config.GetEquation();
  int spatial_accuracy = config.GetSpatialAccuracy();
  Surface_Type surface = config.GetSurface();
  Temporal_Method_Type temporal_method = config.GetTemporalMethod();
  Spatial_Method_Type spatial_method = config.GetSpatialmethod();

  unique_ptr<Beta> diffusion_coef_ptr = move(find_diffusion_coefficient(diffusion_coef_type));
  if (diffusion_coef_ptr.get() == nullptr) {
    LOG_FATAL("Diffusion coefficient pointer is NULL!");
  }
  unique_ptr<Surface_Cartesian> surface_ptr = move(find_surface(surface));
  if (surface_ptr.get() == nullptr) {
    LOG_FATAL("Surface coefficient pointer is NULL!");
  }

  // Mesh construction
  Mesh mesh(domain, size, *surface_ptr);
  // Spatial computation
  int index = (spatial_method == Spatial_Method_Type::MIB_V1) ? 1 : 2;
  Intersections inter(*surface_ptr, mesh, *diffusion_coef_ptr, spatial_accuracy, index, out_stream);

  // Temporal computation, ADI / LOD / Trapezoidal Splitting
  if (temporal_method == Temporal_Method_Type::ADI) {
    Call_ADI(equation_type, spatial_accuracy, *diffusion_coef_ptr, mesh, inter, time_info, uh,
             out_stream);
  } else if (temporal_method == Temporal_Method_Type::LOD_IE) {
    char index = 'I';
    Call_LOD(equation_type, *diffusion_coef_ptr, mesh, inter, time_info, uh, index, out_stream);
  } else if (temporal_method == Temporal_Method_Type::LOD_CN) {
    char index = 'C';
    Call_LOD(equation_type, *diffusion_coef_ptr, mesh, inter, time_info, uh, index, out_stream);
  } else if (temporal_method == Temporal_Method_Type::TS) {
    Call_TS(equation_type, *diffusion_coef_ptr, mesh, inter, time_info, uh, out_stream);
  } else {
    cout << "Currently, only LOD and Douglas-ADI are applied" << endl;
    exit(1);
  }

  auto t_end = chrono::high_resolution_clock::now();
  chrono::duration<double> period = t_end - t_start;

  out_stream << setprecision(1) << endl;
  out_stream << "CPU time cost: " << period.count() << " seconds" << endl
             << endl;

  Write_txt(inter, mesh, *diffusion_coef_ptr, uh, equation_type, time_info[1]);

  out_stream.close();

  cout << "[INFO]: PROGRAM FINISHED!" << endl;

  return EXIT_SUCCESS;
}
