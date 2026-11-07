#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>
#include "casters.hpp"

namespace nb = nanobind;

void bind_chang_and_yap(nb::module_ &m);
void bind_aronov(nb::module_ &m);
void bind_cabello(nb::module_ &m);
void bind_refinement(nb::module_ &m);

NB_MODULE(potato_peeling, m) {
    m.doc() = "potato_peeling: C++ convex skull algorithms via CGAL";

    // Submodules
    auto m_cy = m.def_submodule("chang_and_yap", "Chang & Yap (1986) algorithm");
    bind_chang_and_yap(m_cy);

    auto m_aronov = m.def_submodule("aronov", "Aronov et al (2011) meshed algorithm");
    bind_aronov(m_aronov);

    auto m_cabello
        = m.def_submodule("cabello", "Cabello et al (2017) randomized (1-eps)-approximation");
    bind_cabello(m_cabello);

    auto m_refinement = m.def_submodule("refinement", "Refinement of a convex polygon inside P.");
    bind_refinement(m_refinement);
}
