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

#include "f3dock/domain/DomainSpec.h"

#include "f3dock/domain/DomainSampler.h"
#include "f3dock/util/StringUtil.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <sstream>
#include <string>

namespace f3dock {
namespace domain {

namespace {

// Parse `n` space-separated doubles from `value`. Returns true iff
// the string contains exactly `n` numeric tokens (any trailing
// whitespace is allowed).
bool parse_doubles(const char *value, double *out, std::size_t n) {
  if (value == nullptr) {
    return false;
  }
  const char *p = value;
  for (std::size_t i = 0; i < n; ++i) {
    while (*p == ' ' || *p == '\t') {
      ++p;
    }
    if (*p == '\0') {
      return false;
    }
    char *end = nullptr;
    out[i] = std::strtod(p, &end);
    if (end == p) {
      return false;
    }
    p = end;
  }
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') {
    ++p;
  }
  return *p == '\0';
}

// 3x3 rotation about an arbitrary axis through the origin (Rodrigues).
Matrix3 rotation_about_axis(const Point3 &axis, double angle_radians) {
  double nx = axis[0];
  double ny = axis[1];
  double nz = axis[2];
  const double len = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (len <= 0.0) {
    // Degenerate axis -> identity, leaving validation to the caller.
    Matrix3 i = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    return i;
  }
  nx /= len;
  ny /= len;
  nz /= len;
  const double c = std::cos(angle_radians);
  const double s = std::sin(angle_radians);
  const double C = 1.0 - c;
  Matrix3 m = {{{c + nx * nx * C, nx * ny * C - nz * s, nx * nz * C + ny * s},
                {ny * nx * C + nz * s, c + ny * ny * C, ny * nz * C - nx * s},
                {nz * nx * C - ny * s, nz * ny * C + nx * s, c + nz * nz * C}}};
  return m;
}

// Joint transform (child-local -> parent-local) corresponding to a
// rotation by `angle` about an axis line through `axis_point` with
// direction `axis_direction`.
RigidTransform hinge_transform(const Point3 &axis_point,
                               const Point3 &axis_direction, double angle) {
  RigidTransform t;
  t.rotation = rotation_about_axis(axis_direction, angle);
  // Translation = p - R p so that points on the axis are fixed.
  for (int i = 0; i < 3; ++i) {
    double Rp = t.rotation[i][0] * axis_point[0] +
                t.rotation[i][1] * axis_point[1] +
                t.rotation[i][2] * axis_point[2];
    t.translation[i] = axis_point[i] - Rp;
  }
  return t;
}

bool nonzero_axis(const Point3 &v) {
  return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) > 0.0;
}

} // namespace

bool ReceptorDomainConfig::enabled() const { return !domains.empty(); }

int ReceptorDomainConfig::find_root(std::string *error) const {
  if (!enabled()) {
    if (error != nullptr) {
      *error = "no domains declared";
    }
    return -1;
  }

  std::set<int> ids;
  for (const auto &d : domains) {
    if (!ids.insert(d.id).second) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "duplicate receptor domain id " << d.id;
        *error = os.str();
      }
      return -1;
    }
  }

  std::set<int> children;
  for (const auto &j : joints) {
    if (ids.count(j.parent_id) == 0 || ids.count(j.child_id) == 0) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "receptor joint references unknown domain (parent " << j.parent_id
           << " -> child " << j.child_id << ")";
        *error = os.str();
      }
      return -1;
    }
    if (!children.insert(j.child_id).second) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "receptor domain " << j.child_id
           << " has more than one parent joint";
        *error = os.str();
      }
      return -1;
    }
  }

  int root = -1;
  for (int id : ids) {
    if (children.count(id) == 0) {
      if (root != -1) {
        if (error != nullptr) {
          *error = "receptor domain graph has more than one root";
        }
        return -1;
      }
      root = id;
    }
  }

  if (root == -1) {
    if (error != nullptr) {
      *error = "receptor domain graph has no root (every domain has a parent)";
    }
    return -1;
  }

  return root;
}

bool ReceptorDomainConfig::build_graph(DomainGraph *graph,
                                       std::string *error) const {
  if (graph == nullptr) {
    if (error != nullptr) {
      *error = "build_graph called with null graph";
    }
    return false;
  }

  std::string local_err;
  if (find_root(&local_err) < 0) {
    if (error != nullptr) {
      *error = local_err;
    }
    return false;
  }

  DomainGraph tmp;
  for (const auto &d : domains) {
    if (!tmp.addDomain(d.id)) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "failed to add receptor domain " << d.id;
        *error = os.str();
      }
      return false;
    }
  }

  for (const auto &j : joints) {
    RigidTransform t;
    if (j.type == JointType::Hinge) {
      if (!nonzero_axis(j.axis_direction)) {
        if (error != nullptr) {
          std::ostringstream os;
          os << "receptor hinge joint (" << j.parent_id << " -> " << j.child_id
             << ") has zero-length axis direction";
          *error = os.str();
        }
        return false;
      }
      t = hinge_transform(j.axis_point, j.axis_direction,
                          j.initial_angle_radians);
    } else {
      t = RigidTransform::identity();
    }
    if (!tmp.setParent(j.child_id, j.parent_id, t)) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "failed to attach receptor domain " << j.child_id << " to parent "
           << j.parent_id << " (cycle or unknown id)";
        *error = os.str();
      }
      return false;
    }
  }

  *graph = std::move(tmp);
  return true;
}

bool build_graph_with_state(const ReceptorDomainConfig &config,
                            const DomainState &state, DomainGraph *graph,
                            std::string *error) {
  if (graph == nullptr) {
    if (error != nullptr) {
      *error = "build_graph_with_state called with null graph";
    }
    return false;
  }

  // Validate that every override references a real hinge joint in
  // the rest-pose configuration. We reject overrides that target
  // non-hinge joints or non-existent (parent, child) pairs so that
  // bad input files fail at config time rather than silently doing
  // nothing.
  for (const auto &ov : state.hinge_angles) {
    bool matched = false;
    for (const auto &j : config.joints) {
      if (j.parent_id == ov.parent_id && j.child_id == ov.child_id) {
        if (j.type != JointType::Hinge) {
          if (error != nullptr) {
            std::ostringstream os;
            os << "domain state " << state.state_id << " overrides joint ("
               << ov.parent_id << " -> " << ov.child_id
               << ") which is not a hinge";
            *error = os.str();
          }
          return false;
        }
        matched = true;
        break;
      }
    }
    if (!matched) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "domain state " << state.state_id << " overrides unknown joint ("
           << ov.parent_id << " -> " << ov.child_id << ")";
        *error = os.str();
      }
      return false;
    }
  }

  std::string local_err;
  if (config.find_root(&local_err) < 0) {
    if (error != nullptr) {
      *error = local_err;
    }
    return false;
  }

  DomainGraph tmp;
  for (const auto &d : config.domains) {
    if (!tmp.addDomain(d.id)) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "failed to add receptor domain " << d.id;
        *error = os.str();
      }
      return false;
    }
  }

  for (const auto &j : config.joints) {
    RigidTransform t;
    if (j.type == JointType::Hinge) {
      if (!nonzero_axis(j.axis_direction)) {
        if (error != nullptr) {
          std::ostringstream os;
          os << "receptor hinge joint (" << j.parent_id << " -> " << j.child_id
             << ") has zero-length axis direction";
          *error = os.str();
        }
        return false;
      }
      double angle = j.initial_angle_radians;
      for (const auto &ov : state.hinge_angles) {
        if (ov.parent_id == j.parent_id && ov.child_id == j.child_id) {
          angle = ov.angle_radians;
          break;
        }
      }
      t = hinge_transform(j.axis_point, j.axis_direction, angle);
    } else {
      t = RigidTransform::identity();
    }
    if (!tmp.setParent(j.child_id, j.parent_id, t)) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "failed to attach receptor domain " << j.child_id << " to parent "
           << j.parent_id << " (cycle or unknown id)";
        *error = os.str();
      }
      return false;
    }
  }

  *graph = std::move(tmp);
  return true;
}

bool DomainSpecParser::parse_param(const char *key, const char *value,
                                   ReceptorDomainConfig *config, bool *ok) {
  if (key == nullptr || config == nullptr || ok == nullptr) {
    return false;
  }

  const std::string_view key_sv(key);

  if (f3dock::util::iequals(key_sv, "receptorDomain")) {
    double v[3] = {0};
    if (!parse_doubles(value, v, 3)) {
      std::printf("Error: receptorDomain expects 3 integers: "
                  "id firstAtom lastAtom\n");
      *ok = false;
      return true;
    }
    DomainSpec d;
    d.id = static_cast<int>(v[0]);
    d.first_atom_index = static_cast<int>(v[1]);
    d.last_atom_index = static_cast<int>(v[2]);
    if (d.last_atom_index < d.first_atom_index) {
      std::printf(
          "Error: receptorDomain lastAtom (%d) must be >= firstAtom (%d)\n",
          d.last_atom_index, d.first_atom_index);
      *ok = false;
      return true;
    }
    config->domains.push_back(d);
    *ok = true;
    return true;
  }

  if (f3dock::util::iequals(key_sv, "receptorJointFixed")) {
    double v[2] = {0};
    if (!parse_doubles(value, v, 2)) {
      std::printf("Error: receptorJointFixed expects 2 integers: "
                  "parentId childId\n");
      *ok = false;
      return true;
    }
    JointSpec j;
    j.parent_id = static_cast<int>(v[0]);
    j.child_id = static_cast<int>(v[1]);
    j.type = JointType::Fixed;
    config->joints.push_back(j);
    *ok = true;
    return true;
  }

  if (f3dock::util::iequals(key_sv, "receptorJointHinge")) {
    double v[9] = {0};
    if (!parse_doubles(value, v, 9)) {
      std::printf("Error: receptorJointHinge expects 9 numbers: "
                  "parentId childId ax ay az dx dy dz initAngleRad\n");
      *ok = false;
      return true;
    }
    JointSpec j;
    j.parent_id = static_cast<int>(v[0]);
    j.child_id = static_cast<int>(v[1]);
    j.type = JointType::Hinge;
    j.axis_point = {v[2], v[3], v[4]};
    j.axis_direction = {v[5], v[6], v[7]};
    j.initial_angle_radians = v[8];
    config->joints.push_back(j);
    *ok = true;
    return true;
  }

  return false;
}

} // namespace domain
} // namespace f3dock
