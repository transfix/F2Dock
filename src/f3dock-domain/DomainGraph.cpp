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

#include "f3dock/domain/DomainGraph.h"

namespace {

using f3dock::domain::Matrix3;
using f3dock::domain::Point3;

Point3 add(const Point3 &a, const Point3 &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Point3 matVec(const Matrix3 &m, const Point3 &v) {
  return {
      m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
      m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
      m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
  };
}

Matrix3 matMul(const Matrix3 &a, const Matrix3 &b) {
  Matrix3 out = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    }
  }
  return out;
}

} // namespace

namespace f3dock {
namespace domain {

RigidTransform RigidTransform::identity() {
  return {{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
          {{0.0, 0.0, 0.0}}};
}

Point3 RigidTransform::apply(const Point3 &point) const {
  return add(matVec(rotation, point), translation);
}

RigidTransform RigidTransform::compose(const RigidTransform &other) const {
  // this.compose(other) means apply other first, then this.
  RigidTransform out;
  out.rotation = matMul(rotation, other.rotation);
  out.translation = add(matVec(rotation, other.translation), translation);
  return out;
}

bool DomainGraph::addDomain(int id) {
  if (nodes_.count(id) != 0) {
    return false;
  }

  Node node;
  node.id = id;
  nodes_[id] = node;
  return true;
}

bool DomainGraph::setParent(int child_id, int parent_id,
                            const RigidTransform &child_to_parent) {
  auto child_it = nodes_.find(child_id);
  auto parent_it = nodes_.find(parent_id);
  if (child_it == nodes_.end() || parent_it == nodes_.end() ||
      child_id == parent_id) {
    return false;
  }

  Node backup = child_it->second;
  child_it->second.parent_id = parent_id;
  child_it->second.has_parent = true;
  child_it->second.child_to_parent = child_to_parent;

  if (hasCycleFrom(child_id)) {
    child_it->second = backup;
    return false;
  }

  return true;
}

bool DomainGraph::worldTransform(int id, RigidTransform *out) const {
  if (out == nullptr) {
    return false;
  }

  auto it = nodes_.find(id);
  if (it == nodes_.end()) {
    return false;
  }

  RigidTransform acc = RigidTransform::identity();
  int current = id;
  while (true) {
    auto cur_it = nodes_.find(current);
    if (cur_it == nodes_.end()) {
      return false;
    }

    if (!cur_it->second.has_parent) {
      *out = acc;
      return true;
    }

    const RigidTransform &child_to_parent = cur_it->second.child_to_parent;
    acc = child_to_parent.compose(acc);
    current = cur_it->second.parent_id;
  }
}

bool DomainGraph::hasCycleFrom(int id) const {
  std::map<int, bool> visited;
  int current = id;
  while (true) {
    if (visited[current]) {
      return true;
    }
    visited[current] = true;

    auto it = nodes_.find(current);
    if (it == nodes_.end() || !it->second.has_parent) {
      return false;
    }
    current = it->second.parent_id;
  }
}

} // namespace domain
} // namespace f3dock
