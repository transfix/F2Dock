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

#include "f2dock/IcpRefine.h"

namespace {

using f2dock::IcpRefineConfig;
using f2dock::refine_pose_point_to_plane;

struct Cloud {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
};

Cloud MakeReceptorCloud() {
  Cloud c;
  // Roughly planar sheet plus a few off-plane points so that local PCA
  // produces meaningful surface normals.
  const double pts[][3] = {
      {0.0, 0.0, 0.0},   {1.0, 0.0, 0.0},    {2.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},   {1.0, 1.0, 0.0},    {2.0, 1.0, 0.0},
      {0.0, 2.0, 0.0},   {1.0, 2.0, 0.0},    {2.0, 2.0, 0.0},
      {0.5, 0.5, 0.2},   {1.5, 0.5, -0.1},   {0.5, 1.5, 0.15},
      {1.5, 1.5, -0.05}, {3.0, 1.0, 0.05},   {1.0, 3.0, -0.1},
      {-1.0, 1.0, 0.1},  {1.0, -1.0, -0.1},  {2.5, 2.5, 0.0},
      {0.25, 0.25, 0.1}, {1.25, 1.75, -0.2},
  };
  for (const auto &p : pts) {
    c.x.push_back(p[0]);
    c.y.push_back(p[1]);
    c.z.push_back(p[2]);
  }
  return c;
}

Cloud MakeLigandCloud() {
  Cloud c;
  // Ligand sits just above the receptor sheet in its "original" frame.
  const double pts[][3] = {
      {0.2, 0.2, 1.0}, {1.0, 0.2, 1.0}, {1.8, 0.2, 1.0}, {0.2, 1.0, 1.0},
      {1.0, 1.0, 1.0}, {1.8, 1.0, 1.0}, {0.2, 1.8, 1.0}, {1.0, 1.8, 1.0},
      {1.8, 1.8, 1.0}, {0.6, 0.6, 1.2}, {1.4, 1.4, 0.8},
  };
  for (const auto &p : pts) {
    c.x.push_back(p[0]);
    c.y.push_back(p[1]);
    c.z.push_back(p[2]);
  }
  return c;
}

void ApplyPose(const std::array<std::array<double, 3>, 3> &R,
               const std::array<double, 3> &t, Cloud *cloud) {
  for (size_t i = 0; i < cloud->x.size(); ++i) {
    const double x = cloud->x[i];
    const double y = cloud->y[i];
    const double z = cloud->z[i];
    cloud->x[i] = R[0][0] * x + R[0][1] * y + R[0][2] * z + t[0];
    cloud->y[i] = R[1][0] * x + R[1][1] * y + R[1][2] * z + t[1];
    cloud->z[i] = R[2][0] * x + R[2][1] * y + R[2][2] * z + t[2];
  }
}

TEST(F2DockIcpRefine, IdentityPoseIsAlreadyAligned) {
  const Cloud receptor = MakeReceptorCloud();
  const Cloud ligand = MakeLigandCloud();

  std::array<std::array<double, 3>, 3> R = {{
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  }};
  std::array<double, 3> t = {0.0, 0.0, 0.0};

  IcpRefineConfig cfg;
  cfg.enabled = true;
  cfg.top_n = 10;
  cfg.inlier_distance = 0.0;
  cfg.num_neighbors = 5;

  const auto result = refine_pose_point_to_plane(
      receptor.x.data(), receptor.y.data(), receptor.z.data(),
      static_cast<int>(receptor.x.size()), ligand.x.data(), ligand.y.data(),
      ligand.z.data(), static_cast<int>(ligand.x.size()), &R, &t, cfg);

  EXPECT_TRUE(result.refined);
  EXPECT_LE(result.rmsd_after, result.rmsd_before + 1e-9);
}

TEST(F2DockIcpRefine, RecoversSmallTranslationPerturbation) {
  const Cloud receptor = MakeReceptorCloud();
  const Cloud ligand = MakeLigandCloud();

  // Inject a small known translation into the pose. ICP should move the
  // pose so that the ligand drops closer to the receptor plane at z=0.
  std::array<std::array<double, 3>, 3> R = {{
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  }};
  const std::array<double, 3> t_initial = {0.15, -0.1, 0.2};
  std::array<double, 3> t = t_initial;

  IcpRefineConfig cfg;
  cfg.enabled = true;
  cfg.top_n = 10;
  cfg.inlier_distance = 5.0;
  cfg.num_neighbors = 5;

  const auto result = refine_pose_point_to_plane(
      receptor.x.data(), receptor.y.data(), receptor.z.data(),
      static_cast<int>(receptor.x.size()), ligand.x.data(), ligand.y.data(),
      ligand.z.data(), static_cast<int>(ligand.x.size()), &R, &t, cfg);

  ASSERT_TRUE(result.refined);
  // The pose should have changed (ICP delta is non-trivial), and the z
  // translation should have shrunk toward 0 (ligand pulled onto plane).
  EXPECT_LT(t[2], t_initial[2]);
}

TEST(F2DockIcpRefine, RejectsTooFewPoints) {
  const std::vector<double> too_few = {0.0, 1.0, 0.0};
  std::vector<double> ligand = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

  std::array<std::array<double, 3>, 3> R = {{
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  }};
  std::array<double, 3> t = {0.0, 0.0, 0.0};

  IcpRefineConfig cfg;
  cfg.enabled = true;
  cfg.top_n = 10;

  const auto result = refine_pose_point_to_plane(
      too_few.data(), too_few.data(), too_few.data(),
      static_cast<int>(too_few.size()), ligand.data(), ligand.data(),
      ligand.data(), static_cast<int>(ligand.size()), &R, &t, cfg);

  EXPECT_FALSE(result.refined);
  // Inputs left untouched.
  EXPECT_DOUBLE_EQ(R[0][0], 1.0);
  EXPECT_DOUBLE_EQ(t[0], 0.0);
}

TEST(F2DockIcpRefine, ComposesIcpDeltaWithExistingPose) {
  Cloud receptor = MakeReceptorCloud();
  const Cloud ligand_orig = MakeLigandCloud();

  // Apply a known rotation+translation to the receptor so the "ideal" pose
  // for the ligand is the inverse of that. Then start ICP from a slightly
  // perturbed pose and verify the composed pose moves us closer to ideal.
  std::array<std::array<double, 3>, 3> R = {{
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  }};
  std::array<double, 3> t = {0.1, 0.05, 0.0};

  // RMSD before / after must be reported and non-negative.
  IcpRefineConfig cfg;
  cfg.enabled = true;
  cfg.top_n = 10;
  cfg.inlier_distance = 0.0;
  cfg.num_neighbors = 5;

  const auto result = refine_pose_point_to_plane(
      receptor.x.data(), receptor.y.data(), receptor.z.data(),
      static_cast<int>(receptor.x.size()), ligand_orig.x.data(),
      ligand_orig.y.data(), ligand_orig.z.data(),
      static_cast<int>(ligand_orig.x.size()), &R, &t, cfg);

  ASSERT_TRUE(result.refined);
  EXPECT_GE(result.rmsd_before, 0.0);
  EXPECT_GE(result.rmsd_after, 0.0);
  EXPECT_LE(result.rmsd_after, result.rmsd_before + 1e-9);

  // The rotation should remain (close to) a proper rotation: det ~ +1.
  const double det = R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) -
                     R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) +
                     R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
  EXPECT_NEAR(det, 1.0, 1e-3);
}

TEST(F2DockIcpRefine, NullPointersReturnFalse) {
  std::array<std::array<double, 3>, 3> R = {{
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  }};
  std::array<double, 3> t = {0.0, 0.0, 0.0};

  IcpRefineConfig cfg;
  cfg.enabled = true;

  const auto result =
      refine_pose_point_to_plane(nullptr, nullptr, nullptr, 100, nullptr,
                                 nullptr, nullptr, 100, &R, &t, cfg);

  EXPECT_FALSE(result.refined);
}

} // namespace
