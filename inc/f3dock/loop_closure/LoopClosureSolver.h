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

#ifndef F2DOCK_F3DOCK_LOOP_CLOSURE_SOLVER_H
#define F2DOCK_F3DOCK_LOOP_CLOSURE_SOLVER_H

#include <array>
#include <vector>

namespace f3dock {
namespace loop_closure {

using Point3 = std::array<double, 3>;
using Backbone9 = std::array<Point3, 9>; // N1,CA1,C1,N2,CA2,C2,N3,CA3,C3

struct ClosureConstraints {
  Point3 n1;
  Point3 ca1;
  Point3 ca3;
  Point3 c3;
};

struct LoopClosureInput {
  Backbone9 backbone;
  ClosureConstraints constraints;
};

class LoopClosureSolver {
public:
  static constexpr int kMaxSolutions = 16;

  // Returns true if at least one valid closure is found.
  bool solve(const LoopClosureInput &input,
             std::vector<Backbone9> &solutions) const;
};

} // namespace loop_closure
} // namespace f3dock

#endif
