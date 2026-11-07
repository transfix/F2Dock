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

#include <array>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "f3dock/icp/IcpAligner.h"

namespace {

using f3dock::icp::IcpAligner;
using f3dock::icp::IcpAlignmentResult;
using f3dock::icp::Point3;

std::vector<Point3> makeModelCloud() {
  return {
      {0.0, 0.0, 0.0},  {1.0, 0.0, 0.0},  {0.0, 1.0, 0.0},  {0.0, 0.0, 1.0},
      {1.0, 2.0, 3.0},  {-2.0, 1.0, 0.5}, {0.7, -1.2, 2.3}, {-1.8, -0.4, 1.1},
      {2.5, -1.0, 0.2}, {1.2, 1.7, -0.8},
  };
}

Point3 applyTransform(const std::array<std::array<double, 3>, 3> &r,
                      const Point3 &t, const Point3 &p) {
  return {
      r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2] + t[0],
      r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2] + t[1],
      r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2] + t[2],
  };
}

double rmsd(const std::vector<Point3> &a, const std::vector<Point3> &b) {
  if (a.size() != b.size() || a.empty()) {
    return 0.0;
  }

  double sum_sq = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    const double dx = a[i][0] - b[i][0];
    const double dy = a[i][1] - b[i][1];
    const double dz = a[i][2] - b[i][2];
    sum_sq += dx * dx + dy * dy + dz * dz;
  }

  return std::sqrt(sum_sq / static_cast<double>(a.size()));
}

TEST(F3DockIcp, AlignsRigidlyTransformedPointCloud) {
  const std::vector<Point3> model = makeModelCloud();

  const double theta = 0.2;
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const std::array<std::array<double, 3>, 3> rotation = {{
      {c, -s, 0.0},
      {s, c, 0.0},
      {0.0, 0.0, 1.0},
  }};
  const Point3 translation = {0.25, -0.15, 0.1};

  std::vector<Point3> moving;
  moving.reserve(model.size());
  for (const auto &p : model) {
    moving.push_back(applyTransform(rotation, translation, p));
  }

  IcpAligner aligner;
  std::vector<Point3> aligned;
  IcpAlignmentResult result;
  const bool ok = aligner.alignPointToPoint(model, moving, aligned, &result);

  ASSERT_TRUE(ok);
  ASSERT_EQ(aligned.size(), model.size());

  const double before = rmsd(model, moving);
  const double after = rmsd(model, aligned);
  EXPECT_GT(before, 0.1);
  EXPECT_LT(after, 1e-3);

  EXPECT_NEAR(result.rmsd_before, before, 1e-9);
  EXPECT_NEAR(result.rmsd_after, after, 1e-9);
  EXPECT_LT(result.rmsd_after, result.rmsd_before);
}

TEST(F3DockIcp, RejectsTooFewPoints) {
  const std::vector<Point3> model = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };
  const std::vector<Point3> moving = model;

  IcpAligner aligner;
  std::vector<Point3> aligned;
  const bool ok = aligner.alignPointToPoint(model, moving, aligned);

  EXPECT_FALSE(ok);
  EXPECT_TRUE(aligned.empty());
}

// ── Point-to-plane tests ────────────────────────────────────────────────

// Build a curved-surface point cloud (graph z = 0.1*(x^2 + y^2) on a grid).
// Dense enough that local PCA recovers reasonable surface normals, which is
// what the point-to-plane variant requires.
std::vector<Point3> makeSurfaceCloud() {
  std::vector<Point3> cloud;
  for (int ix = -3; ix <= 3; ++ix) {
    for (int iy = -3; iy <= 3; ++iy) {
      const double x = 0.5 * ix;
      const double y = 0.5 * iy;
      const double z = 0.1 * (x * x + y * y);
      cloud.push_back({x, y, z});
    }
  }
  return cloud;
}

TEST(F3DockIcp, PointToPlaneAlignsRigidlyTransformedSurface) {
  const std::vector<Point3> model = makeSurfaceCloud();

  const double theta = 0.15;
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const std::array<std::array<double, 3>, 3> rotation = {{
      {c, -s, 0.0},
      {s, c, 0.0},
      {0.0, 0.0, 1.0},
  }};
  const Point3 translation = {0.2, -0.1, 0.05};

  std::vector<Point3> moving;
  moving.reserve(model.size());
  for (const auto &p : model) {
    moving.push_back(applyTransform(rotation, translation, p));
  }

  IcpAligner aligner;
  std::vector<Point3> aligned;
  IcpAlignmentResult result;
  const bool ok = aligner.alignPointToPlane(model, moving, aligned, &result);

  ASSERT_TRUE(ok);
  ASSERT_EQ(aligned.size(), model.size());

  EXPECT_GT(result.rmsd_before, 0.05);
  // Point-to-plane converges to a small residual even when normals are
  // estimated from finite samples; allow some slack vs. point-to-point.
  EXPECT_LT(result.rmsd_after, 0.05);
  EXPECT_LT(result.rmsd_after, result.rmsd_before);
}

TEST(F3DockIcp, PointToPlaneRejectsTooFewPoints) {
  const std::vector<Point3> tiny = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };

  IcpAligner aligner;
  std::vector<Point3> aligned;
  EXPECT_FALSE(aligner.alignPointToPlane(tiny, tiny, aligned));
  EXPECT_TRUE(aligned.empty());
}

TEST(F3DockIcp, PointToPlaneRejectsTinyNeighborhood) {
  const std::vector<Point3> model = makeSurfaceCloud();
  const std::vector<Point3> moving = model;

  IcpAligner aligner;
  std::vector<Point3> aligned;
  // num_neighbors < 3 cannot fit a plane.
  EXPECT_FALSE(
      aligner.alignPointToPlane(model, moving, aligned, nullptr, 0.0, 2));
  EXPECT_TRUE(aligned.empty());
}

} // namespace
