cmake_minimum_required(VERSION 3.12)
project(99_99)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS " -mavx -mavx2 -O3 ")
add_executable(99_99
  parser/parser.hpp
  assembler.hpp
  cigar.hpp
  genetyper.hpp
  graph_wrapper.hpp
  haplotype.hpp
  intel_pairhmm.hpp
  intel_pairhmm.hpp
  math_utils.hpp
  pairhmm.hpp
  quality_utils.hpp
  read_clipper.hpp
  read_filter.hpp
  sam.hpp
  variant.hpp
  main.cpp)
set(BOOST_INCLUDE_DIRS "/usr/local/boost_1_73_0/include/")
set(BOOST_LIBRARY_DIRS "/usr/local/boost_1_73_0/lib")
include_directories(/usr/local/include)
include_directories(${BOOST_INCLUDE_DIRS})
target_link_libraries(99_99 hts)
target_link_libraries(99_99 omp)