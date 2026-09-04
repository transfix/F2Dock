/*
  Copyright 2026 The University of Texas at Austin

        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of F2Dock.

  F2Dock is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  F2Dock is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/

#pragma once

#include "f3dock/domain/DomainGraph.h"

#include <cstddef>
#include <string>
#include <vector>

namespace f3dock {
namespace domain {

// Declares a single rigid sub-domain of the receptor. Atoms belonging
// to this domain are identified by a half-open index range
// [first_atom_index, last_atom_index) into the receptor atom array.
// Atom partitioning is consumed in Phase 4 task 2; for task 1 only the
// indices and the id need to round-trip through the parser.
struct DomainSpec {
  int id = -1;
  int first_atom_index = 0;
  int last_atom_index = 0; // exclusive
};

// Joint type between a parent domain and a child domain. `Fixed` is a
// rigid attachment with no extra degree of freedom; `Hinge` adds a
// single rotational DoF about the given axis. The hinge axis is
// interpreted in the parent's local frame; `initial_angle_radians`
// is the angle baked into `child_to_parent` and represents the
// "rest" configuration (often 0).
enum class JointType {
  Fixed = 0,
  Hinge = 1,
};

struct JointSpec {
  int parent_id = -1;
  int child_id = -1;
  JointType type = JointType::Fixed;
  Point3 axis_point = {0.0, 0.0, 0.0};
  Point3 axis_direction = {0.0, 0.0, 1.0};
  double initial_angle_radians = 0.0;
};

// Top-level configuration block, the domain-graph analogue of
// `FlexSamplingConfig`. Lives on `PARAMS_IN` when dockMode==kF3Dock.
//
// Validation rules enforced by `validate()` / `build_graph()`:
//   * Domain ids are unique.
//   * Every joint references existing domain ids.
//   * Exactly one domain has no parent (the root).
//   * `DomainGraph::setParent` rejects any cycle that would form.
//
// When `domains` is empty (the default) the receptor is treated as a
// single implicit rigid body, matching legacy F2Dock behaviour.
struct ReceptorDomainConfig {
  std::vector<DomainSpec> domains;
  std::vector<JointSpec> joints;

  // True iff at least one explicit domain has been declared. An
  // unconfigured receptor (no domains) is treated as a single rigid
  // body downstream.
  bool enabled() const;

  // Returns the root domain id (the only domain without an incoming
  // joint) on success, or -1 if the configuration is invalid or
  // disabled. When non-null, `*error` is populated with a one-line
  // diagnostic on failure.
  int find_root(std::string *error) const;

  // Validates the configuration and, on success, populates `*graph`
  // with one node per domain and one parent edge per joint. The
  // child-to-parent transform of each edge is the rest-pose joint
  // transform (the identity for `Fixed` joints; a rotation by
  // `initial_angle_radians` about the joint axis for `Hinge` joints).
  //
  // On failure, `*graph` is left untouched and `*error` (when non-
  // null) carries the diagnostic message.
  bool build_graph(DomainGraph *graph, std::string *error) const;
};

// Static-class style parser, parallel to FlexSampler. Recognised keys
// (case-insensitive, all numeric values space-separated):
//
//   receptorDomain         id firstAtom lastAtom
//   receptorJointFixed     parentId childId
//   receptorJointHinge     parentId childId ax ay az dx dy dz initAngleRad
//
// Returns true if `key` was a domain-spec key (whether parsing of
// `value` succeeded or not -- that is reported via `*ok`). Returns
// false if the key was not for this subsystem, leaving `*ok`
// untouched so the caller can fall through to other parsers.
class DomainSpecParser {
public:
  static bool parse_param(const char *key, const char *value,
                          ReceptorDomainConfig *config, bool *ok);
};

} // namespace domain
} // namespace f3dock
