#ifndef PROBLEM_GENERATOR_H
#define PROBLEM_GENERATOR_H

#include "enums.h"
#include "global.h"

#include "utils/formatting.h"

#include "archetypes/particle_injector.h"
#include "archetypes/problem_generator.h"
#include "archetypes/traits.h"
#include "framework/domain/domain.h"
#include "framework/domain/metadomain.h"

#include <vector>

namespace user {
  using namespace ntt;

  template <SimEngine::type S, class M>
  struct PGen : public arch::ProblemGenerator<S, M> {

    // compatibility traits for the problem generator
    static constexpr auto engines {
      arch::traits::pgen::compatible_with<SimEngine::SRPIC, SimEngine::GRPIC>::value
    };
    static constexpr auto metrics {
      arch::traits::pgen::compatible_with<Metric::Minkowski,
                                          Metric::Spherical,
                                          Metric::QSpherical,
                                          Metric::Kerr_Schild,
                                          Metric::QKerr_Schild,
                                          Metric::Kerr_Schild_0>::value
    };
    static constexpr auto dimensions {
      arch::traits::pgen::compatible_with<Dim::_1D, Dim::_2D, Dim::_3D>::value
    };

    // for easy access to variables in the child class
    using arch::ProblemGenerator<S, M>::D;
    using arch::ProblemGenerator<S, M>::C;
    using arch::ProblemGenerator<S, M>::params;

    const Metadomain<S, M>& metadomain;

    inline PGen(const SimulationParams& p, const Metadomain<S, M>& metadomain)
      : arch::ProblemGenerator<S, M> { p }
      , metadomain { metadomain } {}

    inline void InitPrtls(Domain<S, M>& domain) {
      using vector_data_t      = std::vector<std::vector<real_t>>;
      const auto empty         = vector_data_t {};
      const auto sp1_positions = params.template get<vector_data_t>(
        "setup.sp1_positions",
        empty);
      const auto sp1_momenta = params.template get<vector_data_t>(
        "setup.sp1_momenta",
        empty);
      raise::ErrorIf(
        sp1_positions.size() != sp1_momenta.size(),
        "Particle positions and momenta size mismatch in problem generator.",
        HERE);

      std::map<std::string, std::vector<real_t>> sp1_data;
      for (const auto& pos : sp1_positions) {
        raise::ErrorIf(
          (pos.size() != static_cast<std::size_t>(D) and
           M::CoordType == Coord::Cart) or
            (pos.size() != 3u and M::CoordType != Coord::Cart),
          "Particle position dimension mismatch in problem generator.",
          HERE);
        for (auto d = 0u; d < pos.size(); ++d) {
          auto component_name = fmt::format("x%d", d + 1);
          if (M::CoordType != Coord::Cart and D != Dim::_3D and d == 2u) {
            component_name = "phi";
          }
          if (sp1_data.find(component_name) == sp1_data.end()) {
            sp1_data[component_name] = std::vector<real_t> {};
          }
          sp1_data[component_name].push_back(pos[d]);
        }
      }
      for (const auto& mom : sp1_momenta) {
        raise::ErrorIf(
          mom.size() != 3u,
          "Particle momenta dimension mismatch in problem generator.",
          HERE);
        for (auto d = 0u; d < mom.size(); ++d) {
          if (sp1_data.find(fmt::format("ux%d", d + 1)) == sp1_data.end()) {
            sp1_data[fmt::format("ux%d", d + 1)] = std::vector<real_t> {};
          }
          sp1_data[fmt::format("ux%d", d + 1)].push_back(mom[d]);
        }
      }

      arch::InjectGlobally<S, M>(metadomain, domain, (spidx_t)1, sp1_data);
    }
  };

} // namespace user

#endif
