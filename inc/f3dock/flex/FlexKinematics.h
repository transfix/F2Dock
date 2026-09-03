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
#include <vector>

namespace f3dock {
namespace flex {

using Point3 = std::array<double, 3>;

struct HingeSpec {
  Point3 axis_point;
  Point3 axis_direction;
  double angle_radians;
};

struct ShearSpec {
  Point3 plane_normal;
  Point3 shear_direction;
  double shear_factor;
};

class FlexKinematics {
public:
  // Applies a rigid rotation around an arbitrary axis using Rodrigues formula.
  static bool applyHinge(const std::vector<Point3> &input,
                         const HingeSpec &hinge, std::vector<Point3> *output);

  // Applies planar shear: p' = p + dot(p, n) * factor * d.
  static bool applyShear(const std::vector<Point3> &input,
                         const ShearSpec &shear, std::vector<Point3> *output);
};

} // namespace flex
} // namespace f3dock
