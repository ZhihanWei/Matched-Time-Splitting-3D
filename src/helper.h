#pragma once

#include "data/data.h"

#define LOG_FATAL(msg)                                                                        \
  std::cerr << "[ERROR]: " << msg << " @ " << __FILE__ << " line: " << __LINE__ << std::endl; \
  exit(EXIT_FAILURE);

void print_help(char* binary_name) {
  cout << "Usage: " << binary_name << " -i input [-o output] \n \n";
  cout << "-i, --input [required] \n"
          "      input configuration file, sample can be found in directory \"example\" \n"
          "-o, --output [optional] \n"
          "      default: results.txt \n"
          "      output file \n \n";
}

void write_configuration(ofstream& out_stream, data& data) {
  VecDoub domain = data.GetDomain();
  VecDoub size = data.GetSize();
  VecInt time_info = data.GetTime();
  int beta_code = data.GetDiffusionCeff();
  int accuracy = data.GetSpatialAccuracy();
  string surface = data.GetSurface();
  string temporal_method = data.GetTemporalMethod();
  string spatial_method = data.GetSpatialmethod();
  int equation = data.GetEquation();

  out_stream << "------------ Configuration informations ------------" << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Mesh Size"
             << ": " << setw(5) << "NX = " << setiosflags(ios::left)
             << setw(5) << size[0] << setw(5)
             << " NY = " << setiosflags(ios::left) << setw(5) << size[1]
             << setw(5) << " NZ = " << setiosflags(ios::left) << setw(5)
             << size[2] << endl;
  out_stream << setprecision(1) << scientific;
  out_stream << setiosflags(ios::left) << setw(18) << "Time Step"
             << ": " << time_info[2] << endl;
  out_stream << fixed;
  out_stream << setiosflags(ios::left) << setw(18) << "Jump"
             << ": " << JP << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Surface"
             << ": " << surface << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Temporal Method"
             << ": " << temporal method << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Spatial Method"
             << ": L" << spatial_method << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Diffusion Coefficient Type "
             << ": " << beta_code << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Equation Type"
             << ": " << equation << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Order of Spatial Accuracy"
             << ": " << spatial_accuracy << endl
             << endl;
}
