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

#include "f3dock/flex/FlexKinematics.h"

namespace {

using f3dock::flex::FlexKinematics;
using f3dock::flex::HingeSpec;
using f3dock::flex::Point3;
using f3dock::flex::ShearSpec;

TEST(F3DockFlex, AppliesHingeRotationAroundZAxis) {
  const double kHalfPi = std::acos(-1.0) * 0.5;

  const std::vector<Point3> input = {
      {1.0, 0.0, 0.0},
      {0.0, 2.0, 1.0},
  };

  const HingeSpec hinge = {
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0},
      kHalfPi,
  };

  std::vector<Point3> output;
  ASSERT_TRUE(FlexKinematics::applyHinge(input, hinge, &output));
  ASSERT_EQ(output.size(), input.size());

  EXPECT_NEAR(output[0][0], 0.0, 1e-12);
  EXPECT_NEAR(output[0][1], 1.0, 1e-12);
  EXPECT_NEAR(output[0][2], 0.0, 1e-12);

  EXPECT_NEAR(output[1][0], -2.0, 1e-12);
  EXPECT_NEAR(output[1][1], 0.0, 1e-12);
  EXPECT_NEAR(output[1][2], 1.0, 1e-12);
}

TEST(F3DockFlex, RejectsInvalidHingeAxis) {
  const std::vector<Point3> input = {{1.0, 2.0, 3.0}};
  const HingeSpec hinge = {
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
      0.1,
  };

  std::vector<Point3> output = {{9.0, 9.0, 9.0}};
  EXPECT_FALSE(FlexKinematics::applyHinge(input, hinge, &output));
  EXPECT_TRUE(output.empty());
}

TEST(F3DockFlex, AppliesPlanarShearAlongProjectedDirection) {
  const std::vector<Point3> input = {
      {0.0, 0.0, 2.0},
      {1.0, -1.0, -2.0},
      {3.0, 1.0, 0.0},
  };

  const ShearSpec shear = {
      {0.0, 0.0, 1.0},
      {2.0, 0.0, 0.5},
      0.5,
  };

  std::vector<Point3> output;
  ASSERT_TRUE(FlexKinematics::applyShear(input, shear, &output));
  ASSERT_EQ(output.size(), input.size());

  // z=2 shifts by +1 along +x.
  EXPECT_NEAR(output[0][0], 1.0, 1e-12);
  EXPECT_NEAR(output[0][1], 0.0, 1e-12);
  EXPECT_NEAR(output[0][2], 2.0, 1e-12);

  // z=-2 shifts by -1 along +x.
  EXPECT_NEAR(output[1][0], 0.0, 1e-12);
  EXPECT_NEAR(output[1][1], -1.0, 1e-12);
  EXPECT_NEAR(output[1][2], -2.0, 1e-12);

  // z=0 remains unchanged.
  EXPECT_NEAR(output[2][0], 3.0, 1e-12);
  EXPECT_NEAR(output[2][1], 1.0, 1e-12);
  EXPECT_NEAR(output[2][2], 0.0, 1e-12);
}

TEST(F3DockFlex, RejectsDegenerateShearDirection) {
  const std::vector<Point3> input = {{0.0, 0.0, 1.0}};
  const ShearSpec shear = {
      {0.0, 0.0, 1.0},
      {0.0, 0.0, 5.0},
      1.0,
  };

  std::vector<Point3> output = {{7.0, 7.0, 7.0}};
  EXPECT_FALSE(FlexKinematics::applyShear(input, shear, &output));
  EXPECT_TRUE(output.empty());
}

} // namespace
