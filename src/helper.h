#pragma once

#define LOG_FATAL(msg)                                                         \
  std::cerr << msg << " @ " << __FILE__ << " line: " << __LINE__ << std::endl; \
  exit(EXIT_FAILURE);
