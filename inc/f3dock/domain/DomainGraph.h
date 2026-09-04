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

#include <array>
#include <map>
#include <vector>

namespace f3dock {
namespace domain {

using Point3 = std::array<double, 3>;
using Matrix3 = std::array<std::array<double, 3>, 3>;

struct RigidTransform {
  Matrix3 rotation;
  Point3 translation;

  static RigidTransform identity();

  Point3 apply(const Point3 &point) const;
  RigidTransform compose(const RigidTransform &other) const;
};

class DomainGraph {
public:
  bool addDomain(int id);

  // Sets parent relationship and child-local transform (parent <- child).
  bool setParent(int child_id, int parent_id,
                 const RigidTransform &child_to_parent);

  bool worldTransform(int id, RigidTransform *out) const;

private:
  struct Node {
    int id = -1;
    int parent_id = -1;
    bool has_parent = false;
    RigidTransform child_to_parent = RigidTransform::identity();
  };

  bool hasCycleFrom(int id) const;

  std::map<int, Node> nodes_;
};

} // namespace domain
} // namespace f3dock
