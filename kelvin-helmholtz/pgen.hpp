#ifndef PROBLEM_GENERATOR_H
#define PROBLEM_GENERATOR_H

#include "enums.h"
#include "global.h"

#include "arch/traits.h"
#include "utils/numeric.h"

#include "archetypes/energy_dist.h"
#include "archetypes/particle_injector.h"
#include "archetypes/problem_generator.h"
#include "framework/domain/domain.h"
#include "framework/domain/metadomain.h"

namespace user {
  using namespace ntt;

  template <SimEngine::type S, class M>
  struct KHDistribution : public arch::EnergyDistribution<S, M> {
    KHDistribution(const M&              metric,
                   random_number_pool_t& pool,
                   real_t                temperature,
                   real_t                drift_vel,
                   in                    drift_dir,
                   real_t                width,
                   real_t                y1,
                   real_t                y2)
      : arch::EnergyDistribution<S, M> { metric }
      , pool { pool }
      , temperature { temperature }
      , drift_vel { drift_vel }
      , drift_dir { static_cast<unsigned short>(drift_dir) }
      , width { width }
      , y1 { y1 }
      , y2 { y2 } {}

    Inline void operator()(const coord_t<M::Dim>& x_Ph,
                           vec_t<Dim::_3D>&       v,
                           spidx_t = 0) const {
      arch::SampleFromMaxwellian<S, false>(v, pool, temperature);
      v[drift_dir] += drift_vel * (math::tanh((x_Ph[1] - y1) / width) -
                                   math::tanh((x_Ph[1] - y2) / width) - 1);
    }

  private:
    random_number_pool_t pool;

    const real_t         temperature, drift_vel;
    const unsigned short drift_dir;
    const real_t         width, y1, y2;
  };

  template <Dimension D>
  struct InitFields {
    InitFields(real_t bmag) : bmag { bmag } {}

    Inline auto bx1(const coord_t<D>&) const -> real_t {
      return bmag;
    }

  private:
    const real_t bmag;
  };

  template <SimEngine::type S, class M>
  struct PGen : public arch::ProblemGenerator<S, M> {
    static constexpr auto engines = traits::compatible_with<SimEngine::SRPIC>::value;
    static constexpr auto metrics = traits::compatible_with<Metric::Minkowski>::value;
    static constexpr auto dimensions = traits::compatible_with<Dim::_2D, Dim::_3D>::value;

    using arch::ProblemGenerator<S, M>::D;
    using arch::ProblemGenerator<S, M>::C;
    using arch::ProblemGenerator<S, M>::params;

    const real_t temperature, drift_vel, density;
    const in     drift_dir;
    const real_t width, y1, y2;

    InitFields<D> init_flds;

    inline PGen(const SimulationParams& p, const Metadomain<S, M>& global_domain)
      : arch::ProblemGenerator<S, M> { p }
      , temperature { p.template get<real_t>("setup.temperature", 0.01) }
      , drift_vel { p.template get<real_t>("setup.sonic_mach", 1.0) *
                    math::sqrt((5 / 3) * temperature) }
      , density { p.template get<real_t>("setup.density", 1.0) }
      , drift_dir { static_cast<in>(p.template get<int>("setup.drift_dir", 0)) }
      , width { p.template get<real_t>("setup.width", 0.01) }
      , y1 { p.template get<real_t>("setup.y1", -0.25) }
      , y2 { p.template get<real_t>("setup.y2", 0.25) }
      , init_flds { drift_vel * density *
                    p.template get<real_t>("setup.inv_alfven_mach", 0.0) /
                    math::sqrt(p.template get<real_t>("scales.sigma0")) } {}

    inline void InitPrtls(Domain<S, M>& domain) {
      const auto en_dist  = KHDistribution<S, M>(domain.mesh.metric,
                                                domain.random_pool,
                                                temperature,
                                                drift_vel,
                                                drift_dir,
                                                width,
                                                y1,
                                                y2);
      const auto injector = arch::UniformInjector<S, M, KHDistribution>(en_dist,
                                                                        { 1, 2 });

      arch::InjectUniform<S, M, decltype(injector)>(params,
                                                    domain,
                                                    injector,
                                                    density,
                                                    false);
    }
  };

} // namespace user

#endif
