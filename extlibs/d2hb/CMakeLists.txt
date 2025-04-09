# Set minimum version of CMake.
cmake_minimum_required(VERSION 3.11.4)

# Set project name and version
project(d2hb VERSION 1.0.0)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configure the option whether you want to use the systems own allocater
# WITH_D2HB_SYSTEM_ALLOCATOR=ON
# or our own written allocator strategy
# WITH_D2HB_SYSTEM_ALLOCATOR=OFF
option(WITH_D2HB_SYSTEM_ALLOCATOR "Use D2HB System Allocator" OFF)

add_library(d2hb STATIC)

set(public_headers
  include/d2hb/atomic_helpers.h
  include/d2hb/batcher.h
  include/d2hb/block.h
  include/d2hb/buckets.h
  include/d2hb/dynamic_doubly_hashed_blocks.h
  include/d2hb/dynamic_hashed_blocks.h
  include/d2hb/dynamic_hashtable.h
  include/d2hb/graph.h
  include/d2hb/hash_tools.h
  include/d2hb/integer_log2.h
  include/d2hb/submatrix.h
  include/d2hb/vec.h
)

target_sources(${PROJECT_NAME}
  PRIVATE
  src/buckets.cpp
  src/graph.cpp
  src/dynamic_hashtable.cpp
)

# ===================================
# Options check and set definitions
# ===================================
if(WITH_D2HB_SYSTEM_ALLOCATOR)
  MESSAGE(STATUS "Using D2HB System Allocator.")
  target_compile_definitions(d2hb PUBLIC -DD2HB_SYSTEM_ALLOCATOR)
endif()

set_property(TARGET d2hb PROPERTY POSITION_INDEPENDENT_CODE ON)

find_package(OpenMP)

target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

target_include_directories(${PROJECT_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS ${PROJECT_NAME}
  EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

install(EXPORT "${PROJECT_NAME}Targets"
  FILE "${PROJECT_NAME}Targets.cmake"
  NAMESPACE d2hb::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

add_library(${PROJECT_NAME}::${PROJECT_NAME} STATIC IMPORTED)

set_target_properties(${PROJECT_NAME}::${PROJECT_NAME} PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include"
)

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

# ===================================
# Make a Configuration Package
# ===================================
include(CMakePackageConfigHelpers)

export(EXPORT ${PROJECT_NAME}Targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Targets.cmake"
  NAMESPACE ${PROJECT_NAME}::
)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${version}"
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

set(version 1.0.0)

set_property(TARGET d2hb PROPERTY VERSION ${version})
set_property(TARGET d2hb PROPERTY SOVERSION 1)
set_property(TARGET d2hb PROPERTY
  INTERFACE_d2hb_MAJOR_VERSION 1)
set_property(TARGET d2hb APPEND PROPERTY
  COMPATIBLE_INTERFACE_STRING d2hb_MAJOR_VERSION
)

# ========================
# Test Target
# ========================
option(D2HB_TEST "Build test target." OFF)

if(D2HB_TEST)
  if(EXISTS "${PROJECT_SOURCE_DIR}/test/lib/Catch2")
    add_subdirectory(${PROJECT_SOURCE_DIR}/test/lib/Catch2)
  else()
    message(FATAL_ERROR
      "Missing Catch2 library in test/lib/Catch2 "
      "Please run `git submodule update --init` to fetch the submodule.")
  endif()

  if(EXISTS "${PROJECT_SOURCE_DIR}/test/graphs")
    add_definitions(-DD2HB_TEST_GRAPH_DIR="${PROJECT_SOURCE_DIR}/test/graphs")
  else()
    message(FATAL_ERROR
      "Could not find test graph directory.")
  endif()

  add_executable(d2hb_test
    test/atomic_helpers_test.cpp
    test/block_test.cpp
    test/graph_io.cpp
    test/graph_io.h
    test/hashtable_test.cpp
    test/hashtable_tooling_test.cpp
    test/hypersparse_matrix_test.cpp
    test/matrix_test.cpp
    test/submatrix_test.cpp
    test/vec_test.cpp
  )

  # needed for debugging
  # target_link_libraries(d2hb_test PRIVATE asan)
  target_link_libraries(d2hb_test PRIVATE d2hb Catch2::Catch2WithMain OpenMP::OpenMP_CXX stdc++fs)

  install(
    TARGETS d2hb_test
    RUNTIME DESTINATION bin
  )

  # Test setup. Use CTest + Catch2's CTest integration.
  include(CTest)
  include(${PROJECT_SOURCE_DIR}/test/lib/Catch2/extras/Catch.cmake)
  catch_discover_tests(d2hb_test)
endif()

option(D2HB_DEMO "Build demo program." OFF)

if(D2HB_DEMO)
  add_executable(demo
    demo/main.cpp
  )

  # needed for debugging
  # target_link_libraries(demo PRIVATE asan)
  target_link_libraries(demo PRIVATE d2hb OpenMP::OpenMP_CXX)

  install(
    TARGETS demo
    RUNTIME DESTINATION bin
  )
endif()
