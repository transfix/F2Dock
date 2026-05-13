#pragma once

#include "f3dock/flex/FlexKinematics.h"

#include <cstddef>
#include <string>
#include <vector>

namespace f3dock {
namespace flex {

// A single hinge axis with a discrete angle sweep. The sampler will
// enumerate one HingeSpec per angle in
//   { min_angle_radians + k * step_radians : k = 0..num_steps-1 }.
//
// `num_steps == 0` disables this hinge (it contributes no perturbations).
// `num_steps == 1` produces exactly one HingeSpec at `min_angle_radians`.
struct HingeAxisSpec {
  Point3 axis_point = {0.0, 0.0, 0.0};
  Point3 axis_direction = {0.0, 0.0, 1.0};
  double min_angle_radians = 0.0;
  double step_radians = 0.0;
  std::size_t num_steps = 0;
};

// A single shear plane with a discrete factor sweep. Analogous to
// HingeAxisSpec — enumerates one ShearSpec per factor in
//   { min_shear_factor + k * step : k = 0..num_steps-1 }.
struct ShearPlaneSpec {
  Point3 plane_normal = {0.0, 0.0, 1.0};
  Point3 shear_direction = {1.0, 0.0, 0.0};
  double min_shear_factor = 0.0;
  double step = 0.0;
  std::size_t num_steps = 0;
};

// Top-level configuration block read from the F2Dock parameter file
// (or set programmatically). Lives on PARAMS_IN when dockMode==kF3Dock.
//
// The receptor is sampled across the cross product of all hinge axes
// and shear planes. The first state in every sweep is the identity
// (angle = min_angle_radians, factor = min_shear_factor) — when those
// minimums are zero this means "no perturbation". Callers are expected
// to keep state 0 as the identity so F3Dock-mode top hits with all
// sweeps disabled match F2Dock-mode top hits.
struct FlexSamplingConfig {
  std::vector<HingeAxisSpec> hinges;
  std::vector<ShearPlaneSpec> shears;

  // True iff at least one hinge or shear has num_steps > 0. When this
  // is false, the sampler emits exactly one identity flex state.
  bool enabled() const;

  // Total number of flex states this configuration enumerates. Equal
  // to the product of (num_steps for each enabled hinge/shear), with
  // 1 when nothing is enabled. Bounded; callers should sanity-check
  // before turning on combinatorially large sweeps.
  std::size_t state_count() const;
};

// A concrete enumerated flex state. `state_id` is a stable index in
// [0, state_count()) that callers can attach to scoring records so
// downstream filters can group / re-rank by flex state. The two
// vectors describe the exact transforms to apply, in order:
// all hinges first, then all shears. The vectors are empty when the
// state is the identity / no perturbation.
struct FlexState {
  std::size_t state_id = 0;
  std::vector<HingeSpec> hinges;
  std::vector<ShearSpec> shears;
};

// Enumerates all FlexState combinations described by a config. The
// emitted order is deterministic (lexicographic over the per-knob
// step indices, hinges-major then shears-major).
//
// Returns false if `out` is null or the configuration is malformed
// (e.g. zero-length axis_direction on an enabled hinge). On success,
// `out` is overwritten with state_count() entries.
class FlexSampler {
public:
  static bool enumerate(const FlexSamplingConfig &config,
                        std::vector<FlexState> *out);

  // Parse a single `key value` pair from the F2Dock parameter file.
  // Recognised keys (case-insensitive, all values are space-separated
  // floats):
  //
  //   flexHingeAxis      ax ay az dx dy dz minAngleRad stepRad numSteps
  //   flexShearPlane     nx ny nz dx dy dz minFactor   step    numSteps
  //
  // Returns true if `key` was a flex-sampling key (whether parsing of
  // `value` succeeded or not — that is reported via `*ok`). Returns
  // false if the key was not for this subsystem, leaving `*ok`
  // untouched so the caller can fall through to other parsers.
  static bool parse_param(const char *key, const char *value,
                          FlexSamplingConfig *config, bool *ok);
};

} // namespace flex
} // namespace f3dock
