#pragma once

#include "f3dock/domain/DomainGraph.h"
#include "f3dock/domain/DomainSpec.h"

#include <cstddef>
#include <string>
#include <vector>

namespace f3dock {
namespace domain {

// Per-state override for a hinge joint's angle. The (parent_id,
// child_id) pair identifies the joint in the receptor's
// `ReceptorDomainConfig::joints` table. `angle_radians` replaces
// `JointSpec::initial_angle_radians` for this state only; all other
// joints remain at their rest pose.
struct HingeJointAngle {
  int parent_id = -1;
  int child_id = -1;
  double angle_radians = 0.0;
};

// Sweep description for a single hinge joint. The joint is identified
// by (parent_id, child_id); its angle is enumerated as
// min_angle_radians + k * step_radians for k in [0, num_steps). When
// `num_steps == 0` the joint contributes no states (i.e. is disabled).
// Mirrors `flex::HingeAxisSpec` semantics so callers can reason about
// the two sweeps the same way.
struct DomainHingeSweep {
  int parent_id = -1;
  int child_id = -1;
  double min_angle_radians = 0.0;
  double step_radians = 0.0;
  std::size_t num_steps = 0;
};

// Top-level domain-state sampling configuration. The state space is
// the cross product of per-joint hinge sweeps, in declaration order.
// The first emitted state always pins every joint at
// `min_angle_radians`; when those minimums equal the rest-pose
// `initial_angle_radians`, state 0 is the rest pose and F3Dock-mode
// top hits with no overrides match the no-sweep run bit-for-bit.
struct DomainSamplingConfig {
  std::vector<DomainHingeSweep> hinges;

  // True iff at least one hinge has num_steps > 0.
  bool enabled() const;

  // Total number of states this configuration enumerates. Equal to
  // the product of `num_steps` over enabled hinges, with 1 when
  // nothing is enabled. Callers should sanity-check before turning
  // on combinatorially large sweeps.
  std::size_t state_count() const;
};

// A concrete enumerated domain state. `state_id` is a stable index
// in [0, state_count()). `hinge_angles` lists the per-joint angle
// overrides for this state; joints absent from the list use their
// rest-pose `initial_angle_radians`.
struct DomainState {
  std::size_t state_id = 0;
  std::vector<HingeJointAngle> hinge_angles;
};

class DomainSampler {
public:
  // Enumerates all DomainState combinations described by `config`.
  // The emission order is lexicographic over per-hinge step indices
  // (first hinge varying slowest). Returns false on null `out`.
  static bool enumerate(const DomainSamplingConfig &config,
                        std::vector<DomainState> *out);

  // Parse a single `key value` pair. Recognised key (case-insensitive):
  //
  //   receptorJointHingeSweep parentId childId minAngleRad stepRad numSteps
  //
  // Returns true if `key` was a domain-sampling key (whether parsing
  // of `value` succeeded or not -- that is reported via `*ok`).
  // Returns false if the key was not for this subsystem.
  static bool parse_param(const char *key, const char *value,
                          DomainSamplingConfig *config, bool *ok);
};

// Build a `DomainGraph` from `config`'s rest pose, but replace each
// hinge joint's `initial_angle_radians` with the matching entry in
// `state.hinge_angles` (matched on (parent_id, child_id)). Joints not
// present in `state.hinge_angles` keep their rest-pose transform.
// Returns false (with `*error` populated) when an override refers to
// a joint that does not exist, or when the underlying validation in
// `ReceptorDomainConfig::build_graph` fails. On success, `*graph` is
// overwritten.
bool build_graph_with_state(const ReceptorDomainConfig &config,
                            const DomainState &state, DomainGraph *graph,
                            std::string *error);

} // namespace domain
} // namespace f3dock
