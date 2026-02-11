
#ifndef PROBLEM_GENERATOR_H
#define PROBLEM_GENERATOR_H

#include "enums.h"
#include "global.h"

#include "arch/kokkos_aliases.h"
#include "utils/error.h"
#include "utils/numeric.h"

#include "archetypes/problem_generator.h"
#include "archetypes/traits.h"
#include "archetypes/utils.h"
#include "framework/domain/domain.h"
#include "framework/domain/metadomain.h"

/**
 *
 *      ^ x2              ^
 *      |        ^        |  B_background
 *      | B_wave |        |
 *      |
 *      -------------------------> x1
 *     /        /
 *    /        V E_wave
 *   /
 *  V  x3
 */

namespace user {
  using namespace ntt;

  template <Dimension D>
  struct InitFields {

    /*
      Sets up the initial background magnetic field for the simulation.

      @param B_background
    */
    InitFields(real_t B_background) : B_background { B_background } {}

    Inline auto bx2(const coord_t<D>&) const -> real_t {
      return B_background;
    }

  private:
    const real_t B_background;
  };

  template <SimEngine::type S, class M>
  struct PGen : public arch::ProblemGenerator<S, M> {

    // compatibility traits for the problem generator
    static constexpr auto engines =
      arch::traits::pgen::compatible_with<SimEngine::SRPIC>::value;
    static constexpr auto metrics =
      arch::traits::pgen::compatible_with<Metric::Minkowski>::value;
    static constexpr auto dimensions =
      arch::traits::pgen::compatible_with<Dim::_1D, Dim::_2D, Dim::_3D>::value;

    // for easy access to variables in the child class
    using arch::ProblemGenerator<S, M>::D;
    using arch::ProblemGenerator<S, M>::C;
    using arch::ProblemGenerator<S, M>::params;

    const real_t  temperature;
    const real_t  B_background, a0, omega;
    const real_t  t_transition, t_duration;
    InitFields<D> init_flds;

    inline PGen(const SimulationParams& p, const Metadomain<S, M>& global_domain)
      : arch::ProblemGenerator<S, M> { p }
      , temperature { p.template get<real_t>("setup.temperature",
                                             static_cast<real_t>(1e-4)) }
      , B_background { p.template get<real_t>("setup.B_background", ZERO) }
      , a0 { p.template get<real_t>("setup.a0", ONE) }
      , omega { p.template get<real_t>("setup.omega", ONE) }
      , t_transition { p.template get<real_t>("setup.t_transition", ZERO) }
      , t_duration { p.template get<real_t>("setup.t_duration", ZERO) }
      , init_flds { B_background } {}

    inline void InitPrtls(Domain<S, M>& domain) {
      const auto tot_density = ONE;
      arch::InjectUniformMaxwellian<S, M>(params,
                                          domain,
                                          tot_density,
                                          temperature,
                                          { 1, 2 });
    }

    /**
     * Sets up the driving field on the left boundary.
     *
     * @param bc_in Direction of the boundary (only be used for one side, so no need to check)
     * @param comp Electromagnetic component to set
     *
     * @note Because the wave is only set on the boundary, no coordinate dependency is needed.
     *
     * @note The fields are normalized to B0 (nominal magnetic field)
     *
     * @return Pair of (value to set, whether to set it or not)
     */
    auto FixFieldsConst(simtime_t time, const bc_in&, em comp) const
      -> std::pair<real_t, bool> {
      real_t phase     = time * omega;
      real_t amplitude = ZERO;
      if (time < t_transition) {
        amplitude = (time / t_transition);
      } else if (time < t_transition + t_duration) {
        amplitude = ONE;
      } else {
        amplitude = math::max(
          ONE - (static_cast<real_t>(time) - t_transition - t_duration) / t_transition,
          ZERO);
      }
      amplitude *= a0 * omega / params.template get<real_t>("scales.omegaB0");

      if (comp == em::ex3) {
        return { -amplitude * math::cos(phase), true };
      } else if (comp == em::bx2) {
        return { B_background + math::sin(phase) * amplitude, true };
      } else {
        return { ZERO, true };
      }
    }

    /**
     * Right boundary will simply match the fields to initial values smoothly
     */
    auto MatchFields(simtime_t) const -> InitFields<D> {
      return init_flds;
    }
  };

} // namespace user

#endif