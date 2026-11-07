cmake_minimum_required(VERSION 3.21)
project(potato_peeling_thesis VERSION 0.1.0 LANGUAGES CXX)

if (POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif ()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif ()

find_package(CGAL 6.0 REQUIRED)
find_package(Boost 1.74 CONFIG REQUIRED)

# ── C++ core ──────────────────────────────────────────────────────────────────
add_library(potato_peeling_geometry STATIC thesis_code/src/geometry.cpp)
target_include_directories(potato_peeling_geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/thesis_code/include)
target_compile_features(potato_peeling_geometry PUBLIC cxx_std_17)
target_link_libraries(potato_peeling_geometry PUBLIC CGAL::CGAL PRIVATE Boost::headers)

add_library(potato_peeling_chang_and_yap STATIC
        thesis_code/src/algorithms/chang_and_yap/chang_and_yap.cpp
        thesis_code/src/algorithms/chang_and_yap/memprobe.cpp
)
target_link_libraries(potato_peeling_chang_and_yap PUBLIC potato_peeling_geometry)

add_library(potato_peeling_aronov STATIC
        thesis_code/src/algorithms/aronov/mesh.cpp
        thesis_code/src/algorithms/aronov/predicates.cpp
        thesis_code/src/algorithms/aronov/sweep.cpp
        thesis_code/src/algorithms/aronov/aronov.cpp
)
target_link_libraries(potato_peeling_aronov PUBLIC potato_peeling_geometry)

add_library(potato_peeling_cabello STATIC
        thesis_code/src/algorithms/cabello/kernel.cpp
        thesis_code/src/algorithms/cabello/certify.cpp
        thesis_code/src/algorithms/cabello/sampling.cpp
        thesis_code/src/algorithms/cabello/visibility.cpp
        thesis_code/src/algorithms/cabello/convex_clique.cpp
        thesis_code/src/algorithms/cabello/a_bound.cpp
        thesis_code/src/algorithms/cabello/cabello.cpp
)
target_link_libraries(potato_peeling_cabello PUBLIC
        potato_peeling_geometry potato_peeling_aronov)

add_library(potato_peeling_refinement STATIC
        thesis_code/src/refinement/refinement.cpp
)
target_link_libraries(potato_peeling_refinement PUBLIC potato_peeling_geometry potato_peeling_cabello)

# ── Python module ─────────────────────────────────────────────────────────────
if (SKBUILD)
    find_package(Python COMPONENTS Interpreter Development.Module Development.SABIModule REQUIRED)
    execute_process(
            COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
    find_package(nanobind 3.0 CONFIG REQUIRED)

    nanobind_add_module(potato_peeling_python
            STABLE_ABI
            NB_STATIC
            thesis_code/bindings/module.cpp
            thesis_code/bindings/algorithms/chang_and_yap.cpp
            thesis_code/bindings/algorithms/aronov.cpp
            thesis_code/bindings/algorithms/cabello.cpp
            thesis_code/bindings/refinement.cpp
    )
    set_target_properties(potato_peeling_python PROPERTIES OUTPUT_NAME "potato_peeling")
    target_link_libraries(potato_peeling_python PRIVATE
            potato_peeling_chang_and_yap
            potato_peeling_aronov
            potato_peeling_cabello
            potato_peeling_refinement
    )
    install(TARGETS potato_peeling_python DESTINATION .)
endif ()
