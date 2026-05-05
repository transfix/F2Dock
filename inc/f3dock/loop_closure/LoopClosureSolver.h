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
  bool solve(const LoopClosureInput& input, std::vector<Backbone9>& solutions) const;
};

}  // namespace loop_closure
}  // namespace f3dock

#endif
