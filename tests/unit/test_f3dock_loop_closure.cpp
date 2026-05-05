#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "f3dock/loop_closure/LoopClosureSolver.h"

namespace {

using f3dock::loop_closure::Backbone9;
using f3dock::loop_closure::LoopClosureInput;
using f3dock::loop_closure::LoopClosureSolver;

Backbone9 sample_backbone() {
  // N1, CA1, C1, N2, CA2, C2, N3, CA3, C3
  return Backbone9{{
      {7.773, -9.710, -7.320},
      {6.331, -9.839, -7.259},
      {6.039, -10.924, -6.233},
      {5.295, -10.500, -5.215},
      {4.908, -11.360, -4.113},
      {3.793, -12.330, -4.472},
      {3.822, -13.525, -3.885},
      {2.642, -14.377, -3.863},
      {1.658, -13.856, -2.821},
  }};
}

TEST(F3DockLoopClosure, FindsSolutionsForKnownFragment) {
  LoopClosureInput input;
  input.backbone = sample_backbone();
  input.constraints.n1 = input.backbone[0];
  input.constraints.ca1 = input.backbone[1];
  input.constraints.ca3 = input.backbone[7];
  input.constraints.c3 = input.backbone[8];

  LoopClosureSolver solver;
  std::vector<Backbone9> solutions;
  const bool ok = solver.solve(input, solutions);

  EXPECT_TRUE(ok);
  EXPECT_FALSE(solutions.empty());
  EXPECT_LE(static_cast<int>(solutions.size()), LoopClosureSolver::kMaxSolutions);

  for (const auto& sol : solutions) {
    for (const auto& atom : sol) {
      EXPECT_TRUE(std::isfinite(atom[0]));
      EXPECT_TRUE(std::isfinite(atom[1]));
      EXPECT_TRUE(std::isfinite(atom[2]));
    }
  }
}

TEST(F3DockLoopClosure, RejectsImpossibleClosureConstraints) {
  LoopClosureInput input;
  input.backbone = sample_backbone();
  input.constraints.n1 = input.backbone[0];
  input.constraints.ca1 = input.backbone[1];
  input.constraints.ca3 = {100.0, 100.0, 100.0};
  input.constraints.c3 = {102.0, 98.0, 101.0};

  LoopClosureSolver solver;
  std::vector<Backbone9> solutions;
  const bool ok = solver.solve(input, solutions);

  EXPECT_FALSE(ok);
  EXPECT_TRUE(solutions.empty());
}

}  // namespace
