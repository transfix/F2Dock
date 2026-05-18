#pragma once

#include "f3dock/domain/DomainGraph.h"
#include "f3dock/domain/DomainSpec.h"

#include <cstddef>
#include <string>
#include <vector>

namespace f3dock {
namespace domain {

// Maps each receptor atom to the rigid sub-domain that owns it.
//
// Phase 4 task 2 turns the receptor's flat atom array into a labelled
// partition. Construction validates a `ReceptorDomainConfig` against
// the actual receptor size (`total_atoms`) and produces a dense
// atom-index -> domain-id lookup. Atoms that no `DomainSpec` claims
// retain id == -1 and are treated as belonging to the implicit root
// rigid body (matching legacy F2Dock behaviour for unconfigured
// receptors).
//
// `apply()` consumes a `DomainGraph` whose nodes carry the *current*
// per-joint pose (rest pose for now; per-state poses come in task 3)
// and rewrites baseline coordinates by composing each atom's owning
// world transform. This is the per-atom moral equivalent of
// `f3dock::flex::apply_flex_state` and is the seam the Phase 4 task 3
// scoring-loop refactor plugs into.
class DomainPartition {
public:
  // True iff `build()` has been called successfully with a non-empty
  // configuration. An "empty" partition is the legacy single-rigid-
  // body case; callers should skip `apply()` entirely in that case.
  bool enabled() const { return enabled_; }

  // Number of atoms the partition was built for.
  std::size_t size() const { return atom_to_domain_.size(); }

  // Validate `config` against `total_atoms` and build the dense
  // atom -> domain-id map. Validation rejects:
  //   * Any domain whose range falls outside [0, total_atoms).
  //   * Any domain with last < first.
  //   * Overlapping domain ranges.
  // Returns true on success. On failure, the partition is left
  // empty / disabled and (when non-null) `*error` carries the
  // diagnostic.
  bool build(const ReceptorDomainConfig &config, std::size_t total_atoms,
             std::string *error);

  // Returns the domain id that owns `atom_index`, or -1 when no
  // domain claims it (unclaimed atoms ride with the root rigid body).
  int domain_for_atom(std::size_t atom_index) const;

  // Apply per-domain world transforms (looked up in `graph`) to the
  // baseline atom arrays `x_in / y_in / z_in` (length `size()`),
  // writing the results into `x_out / y_out / z_out`. The output
  // pointers may alias the input pointers (per-atom independent
  // writes), so `apply(g, x, y, z, x, y, z)` is well-defined.
  // Unclaimed atoms are copied unchanged.
  //
  // Returns false (and populates `*error`) if `graph` is missing a
  // domain referenced by the partition.
  bool apply(const DomainGraph &graph, const double *x_in, const double *y_in,
             const double *z_in, double *x_out, double *y_out, double *z_out,
             std::string *error) const;

private:
  bool enabled_ = false;
  std::vector<int> atom_to_domain_; // -1 = unclaimed
};

} // namespace domain
} // namespace f3dock
