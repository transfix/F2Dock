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

#include "f3dock/flex/FlexApply.h"

#include <vector>

namespace f3dock {
namespace flex {

bool apply_flex_state(const FlexState &state, double *x, double *y, double *z,
                      std::size_t n) {
  if (n == 0) {
    return true;
  }
  if (x == nullptr || y == nullptr || z == nullptr) {
    return false;
  }
  if (state.hinges.empty() && state.shears.empty()) {
    return true;
  }

  // Pack into the Point3 layout FlexKinematics consumes. Reused across
  // the hinge and shear sweeps as the running coordinate set.
  std::vector<Point3> pts(n);
  for (std::size_t i = 0; i < n; ++i) {
    pts[i] = {x[i], y[i], z[i]};
  }

  std::vector<Point3> tmp;
  tmp.reserve(n);

  for (const HingeSpec &h : state.hinges) {
    if (!FlexKinematics::applyHinge(pts, h, &tmp)) {
      return false;
    }
    pts.swap(tmp);
  }
  for (const ShearSpec &s : state.shears) {
    if (!FlexKinematics::applyShear(pts, s, &tmp)) {
      return false;
    }
    pts.swap(tmp);
  }

  for (std::size_t i = 0; i < n; ++i) {
    x[i] = pts[i][0];
    y[i] = pts[i][1];
    z[i] = pts[i][2];
  }
  return true;
}

} // namespace flex
} // namespace f3dock
