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

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "f3dock/flex/FlexApply.h"
#include "f3dock/flex/FlexKinematics.h"
#include "f3dock/flex/FlexSampler.h"

namespace {

using f3dock::flex::apply_flex_state;
using f3dock::flex::FlexState;
using f3dock::flex::HingeSpec;
using f3dock::flex::ShearSpec;

constexpr double kEps = 1e-12;

TEST(F3DockFlexApply, EmptyArraysReturnTrue) {
  FlexState state;
  EXPECT_TRUE(apply_flex_state(state, nullptr, nullptr, nullptr, 0));
}

TEST(F3DockFlexApply, NullPointerWithNonZeroNIsRejected) {
  FlexState state;
  double x = 0.0, y = 0.0, z = 0.0;
  EXPECT_FALSE(apply_flex_state(state, nullptr, &y, &z, 1));
  EXPECT_FALSE(apply_flex_state(state, &x, nullptr, &z, 1));
  EXPECT_FALSE(apply_flex_state(state, &x, &y, nullptr, 1));
}

TEST(F3DockFlexApply, IdentityStateIsNoop) {
  // An empty hinge+shear list must not touch the coordinates -- this is
  // what guarantees F3Dock-mode state 0 matches F2Dock rigid-mode output.
  FlexState state;
  double x[3] = {1.0, 2.0, 3.0};
  double y[3] = {-4.0, 5.0, -6.0};
  double z[3] = {7.0, -8.0, 9.0};
  ASSERT_TRUE(apply_flex_state(state, x, y, z, 3));
  EXPECT_DOUBLE_EQ(x[0], 1.0);
  EXPECT_DOUBLE_EQ(x[1], 2.0);
  EXPECT_DOUBLE_EQ(x[2], 3.0);
  EXPECT_DOUBLE_EQ(y[0], -4.0);
  EXPECT_DOUBLE_EQ(y[1], 5.0);
  EXPECT_DOUBLE_EQ(y[2], -6.0);
  EXPECT_DOUBLE_EQ(z[0], 7.0);
  EXPECT_DOUBLE_EQ(z[1], -8.0);
  EXPECT_DOUBLE_EQ(z[2], 9.0);
}

TEST(F3DockFlexApply, HingeNinetyDegreesAroundZ) {
  // Rotate (1, 0, 5) by +90 degrees around the z-axis through origin:
  // expected (0, 1, 5).
  FlexState state;
  HingeSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 1.0};
  h.angle_radians = M_PI / 2.0;
  state.hinges.push_back(h);

  double x[1] = {1.0};
  double y[1] = {0.0};
  double z[1] = {5.0};
  ASSERT_TRUE(apply_flex_state(state, x, y, z, 1));
  EXPECT_NEAR(x[0], 0.0, kEps);
  EXPECT_NEAR(y[0], 1.0, kEps);
  EXPECT_NEAR(z[0], 5.0, kEps);
}

TEST(F3DockFlexApply, HingePreservesPointsOnAxis) {
  // A point lying exactly on the axis is fixed under any rotation.
  FlexState state;
  HingeSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 1.0};
  h.angle_radians = 1.2345;
  state.hinges.push_back(h);

  double x[1] = {0.0};
  double y[1] = {0.0};
  double z[1] = {7.0};
  ASSERT_TRUE(apply_flex_state(state, x, y, z, 1));
  EXPECT_NEAR(x[0], 0.0, kEps);
  EXPECT_NEAR(y[0], 0.0, kEps);
  EXPECT_NEAR(z[0], 7.0, kEps);
}

TEST(F3DockFlexApply, TwoHingesCompose) {
  // Two +90 degree rotations around z compose to +180 degrees:
  // (1, 0, 0) -> (-1, 0, 0).
  FlexState state;
  HingeSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 1.0};
  h.angle_radians = M_PI / 2.0;
  state.hinges.push_back(h);
  state.hinges.push_back(h);

  double x[1] = {1.0};
  double y[1] = {0.0};
  double z[1] = {0.0};
  ASSERT_TRUE(apply_flex_state(state, x, y, z, 1));
  EXPECT_NEAR(x[0], -1.0, kEps);
  EXPECT_NEAR(y[0], 0.0, kEps);
  EXPECT_NEAR(z[0], 0.0, kEps);
}

TEST(F3DockFlexApply, ShearPlanarOffset) {
  // FlexKinematics::applyShear computes p' = p + dot(p, n) * factor * d.
  // With n = (0,0,1), d = (1,0,0), factor = 0.5, the point (0, 0, 2)
  // becomes (0 + 2*0.5*1, 0, 2) = (1, 0, 2).
  FlexState state;
  ShearSpec s;
  s.plane_normal = {0.0, 0.0, 1.0};
  s.shear_direction = {1.0, 0.0, 0.0};
  s.shear_factor = 0.5;
  state.shears.push_back(s);

  double x[1] = {0.0};
  double y[1] = {0.0};
  double z[1] = {2.0};
  ASSERT_TRUE(apply_flex_state(state, x, y, z, 1));
  EXPECT_NEAR(x[0], 1.0, kEps);
  EXPECT_NEAR(y[0], 0.0, kEps);
  EXPECT_NEAR(z[0], 2.0, kEps);
}

TEST(F3DockFlexApply, ShearLeavesInPlanePointAlone) {
  // A point with zero projection on the plane normal is fixed under
  // the planar shear formula above.
  FlexState state;
  ShearSpec s;
  s.plane_normal = {0.0, 0.0, 1.0};
  s.shear_direction = {1.0, 0.0, 0.0};
  s.shear_factor = 0.9;
  state.shears.push_back(s);

  double x[1] = {3.0};
  double y[1] = {-4.0};
  double z[1] = {0.0};
  ASSERT_TRUE(apply_flex_state(state, x, y, z, 1));
  EXPECT_NEAR(x[0], 3.0, kEps);
  EXPECT_NEAR(y[0], -4.0, kEps);
  EXPECT_NEAR(z[0], 0.0, kEps);
}

TEST(F3DockFlexApply, HingeThenShearAppliedInOrder) {
  // Order matters: hinge first, then shear. Rotate (1, 0, 0) by +90
  // around z -> (0, 1, 0), then a shear with n = (0,1,0), d = (0,0,1),
  // factor = 1 maps (0, 1, 0) -> (0, 1, 1).
  FlexState state;

  HingeSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 1.0};
  h.angle_radians = M_PI / 2.0;
  state.hinges.push_back(h);

  ShearSpec s;
  s.plane_normal = {0.0, 1.0, 0.0};
  s.shear_direction = {0.0, 0.0, 1.0};
  s.shear_factor = 1.0;
  state.shears.push_back(s);

  double x[1] = {1.0};
  double y[1] = {0.0};
  double z[1] = {0.0};
  ASSERT_TRUE(apply_flex_state(state, x, y, z, 1));
  EXPECT_NEAR(x[0], 0.0, kEps);
  EXPECT_NEAR(y[0], 1.0, kEps);
  EXPECT_NEAR(z[0], 1.0, kEps);
}

TEST(F3DockFlexApply, DegenerateHingeAxisIsRejected) {
  FlexState state;
  HingeSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 0.0}; // zero length
  h.angle_radians = 1.0;
  state.hinges.push_back(h);

  double x[1] = {1.0};
  double y[1] = {0.0};
  double z[1] = {0.0};
  EXPECT_FALSE(apply_flex_state(state, x, y, z, 1));
}

TEST(F3DockFlexApply, DegenerateShearNormalIsRejected) {
  FlexState state;
  ShearSpec s;
  s.plane_normal = {0.0, 0.0, 0.0}; // zero length
  s.shear_direction = {1.0, 0.0, 0.0};
  s.shear_factor = 0.5;
  state.shears.push_back(s);

  double x[1] = {0.0};
  double y[1] = {0.0};
  double z[1] = {1.0};
  EXPECT_FALSE(apply_flex_state(state, x, y, z, 1));
}

TEST(F3DockFlexApply, MultipleAtomsAreTransformedIndependently) {
  // The hinge must rotate every supplied atom, not just the first.
  FlexState state;
  HingeSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 1.0};
  h.angle_radians = M_PI / 2.0;
  state.hinges.push_back(h);

  double x[3] = {1.0, 0.0, 2.0};
  double y[3] = {0.0, 1.0, 3.0};
  double z[3] = {0.0, 0.0, 0.0};
  ASSERT_TRUE(apply_flex_state(state, x, y, z, 3));
  EXPECT_NEAR(x[0], 0.0, kEps);
  EXPECT_NEAR(y[0], 1.0, kEps);
  EXPECT_NEAR(x[1], -1.0, kEps);
  EXPECT_NEAR(y[1], 0.0, kEps);
  EXPECT_NEAR(x[2], -3.0, kEps);
  EXPECT_NEAR(y[2], 2.0, kEps);
}

} // namespace
