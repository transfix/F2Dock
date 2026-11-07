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

// Unit tests for the Phase 4 task-2 receptor-domain atom partition.
// Mirrors the apply()/state pattern used by f3dock::flex.

#include "f3dock/domain/DomainGraph.h"
#include "f3dock/domain/DomainPartition.h"
#include "f3dock/domain/DomainSpec.h"

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using f3dock::domain::DomainGraph;
using f3dock::domain::DomainPartition;
using f3dock::domain::DomainSpec;
using f3dock::domain::JointSpec;
using f3dock::domain::JointType;
using f3dock::domain::ReceptorDomainConfig;
using f3dock::domain::RigidTransform;

namespace {

constexpr double kEps = 1e-9;
constexpr double kPi = 3.14159265358979323846;

// Build a simple two-domain configuration: domain 1 owns atoms
// [0, split), domain 2 owns atoms [split, total). Joint is a hinge
// about the z-axis through the origin, rest angle `rest_angle`.
ReceptorDomainConfig MakeTwoDomainHinge(int split, int total,
                                        double rest_angle) {
  ReceptorDomainConfig cfg;
  DomainSpec a;
  a.id = 1;
  a.first_atom_index = 0;
  a.last_atom_index = split;
  cfg.domains.push_back(a);
  DomainSpec b;
  b.id = 2;
  b.first_atom_index = split;
  b.last_atom_index = total;
  cfg.domains.push_back(b);
  JointSpec j;
  j.parent_id = 1;
  j.child_id = 2;
  j.type = JointType::Hinge;
  j.axis_point = {0.0, 0.0, 0.0};
  j.axis_direction = {0.0, 0.0, 1.0};
  j.initial_angle_radians = rest_angle;
  cfg.joints.push_back(j);
  return cfg;
}

} // namespace

// ─────────────────────────────────────────────────────────────────────
// Build / validation
// ─────────────────────────────────────────────────────────────────────

TEST(DomainPartition, EmptyConfigBuildsDisabled) {
  ReceptorDomainConfig cfg;
  DomainPartition p;
  std::string err;
  ASSERT_TRUE(p.build(cfg, 10, &err));
  EXPECT_FALSE(p.enabled());
  EXPECT_EQ(p.size(), 10u);
  for (std::size_t i = 0; i < p.size(); ++i) {
    EXPECT_EQ(p.domain_for_atom(i), -1);
  }
}

TEST(DomainPartition, BuildsDenseLookup) {
  ReceptorDomainConfig cfg = MakeTwoDomainHinge(/*split=*/3, /*total=*/5, 0.0);
  DomainPartition p;
  std::string err;
  ASSERT_TRUE(p.build(cfg, 5, &err));
  EXPECT_TRUE(p.enabled());
  EXPECT_EQ(p.domain_for_atom(0), 1);
  EXPECT_EQ(p.domain_for_atom(1), 1);
  EXPECT_EQ(p.domain_for_atom(2), 1);
  EXPECT_EQ(p.domain_for_atom(3), 2);
  EXPECT_EQ(p.domain_for_atom(4), 2);
}

TEST(DomainPartition, UnclaimedAtomsStayUnclaimed) {
  ReceptorDomainConfig cfg;
  DomainSpec d;
  d.id = 1;
  d.first_atom_index = 2;
  d.last_atom_index = 4;
  cfg.domains.push_back(d);
  DomainPartition p;
  std::string err;
  ASSERT_TRUE(p.build(cfg, 6, &err));
  EXPECT_EQ(p.domain_for_atom(0), -1);
  EXPECT_EQ(p.domain_for_atom(1), -1);
  EXPECT_EQ(p.domain_for_atom(2), 1);
  EXPECT_EQ(p.domain_for_atom(3), 1);
  EXPECT_EQ(p.domain_for_atom(4), -1);
  EXPECT_EQ(p.domain_for_atom(5), -1);
}

TEST(DomainPartition, RejectsOutOfRange) {
  ReceptorDomainConfig cfg;
  DomainSpec d;
  d.id = 1;
  d.first_atom_index = 0;
  d.last_atom_index = 11;
  cfg.domains.push_back(d);
  DomainPartition p;
  std::string err;
  EXPECT_FALSE(p.build(cfg, 10, &err));
  EXPECT_NE(err.find("exceeds"), std::string::npos);
  EXPECT_FALSE(p.enabled());
}

TEST(DomainPartition, RejectsNegativeFirstAtom) {
  ReceptorDomainConfig cfg;
  DomainSpec d;
  d.id = 7;
  d.first_atom_index = -1;
  d.last_atom_index = 3;
  cfg.domains.push_back(d);
  DomainPartition p;
  std::string err;
  EXPECT_FALSE(p.build(cfg, 10, &err));
  EXPECT_FALSE(p.enabled());
}

TEST(DomainPartition, RejectsInvertedRange) {
  ReceptorDomainConfig cfg;
  DomainSpec d;
  d.id = 9;
  d.first_atom_index = 5;
  d.last_atom_index = 3;
  cfg.domains.push_back(d);
  DomainPartition p;
  std::string err;
  EXPECT_FALSE(p.build(cfg, 10, &err));
}

TEST(DomainPartition, RejectsOverlap) {
  ReceptorDomainConfig cfg;
  DomainSpec a;
  a.id = 1;
  a.first_atom_index = 0;
  a.last_atom_index = 5;
  cfg.domains.push_back(a);
  DomainSpec b;
  b.id = 2;
  b.first_atom_index = 3;
  b.last_atom_index = 7;
  cfg.domains.push_back(b);
  DomainPartition p;
  std::string err;
  EXPECT_FALSE(p.build(cfg, 10, &err));
  EXPECT_NE(err.find("claimed"), std::string::npos);
}

// ─────────────────────────────────────────────────────────────────────
// apply()
// ─────────────────────────────────────────────────────────────────────

TEST(DomainPartition, ApplyIdentityRestPoseIsBitForBit) {
  // Hinge rest angle 0 -> all world transforms identity, output ==
  // input regardless of partition labels.
  ReceptorDomainConfig cfg = MakeTwoDomainHinge(/*split=*/2, /*total=*/4,
                                                /*rest_angle=*/0.0);
  DomainPartition p;
  ASSERT_TRUE(p.build(cfg, 4, nullptr));
  DomainGraph graph;
  ASSERT_TRUE(cfg.build_graph(&graph, nullptr));

  const double x_in[] = {1.0, 2.0, 3.0, 4.0};
  const double y_in[] = {0.5, 1.5, 2.5, 3.5};
  const double z_in[] = {-1.0, -2.0, -3.0, -4.0};
  double x_out[4], y_out[4], z_out[4];
  std::string err;
  ASSERT_TRUE(p.apply(graph, x_in, y_in, z_in, x_out, y_out, z_out, &err));
  for (int i = 0; i < 4; ++i) {
    EXPECT_DOUBLE_EQ(x_out[i], x_in[i]);
    EXPECT_DOUBLE_EQ(y_out[i], y_in[i]);
    EXPECT_DOUBLE_EQ(z_out[i], z_in[i]);
  }
}

TEST(DomainPartition, ApplyHingeOnlyRotatesChildAtoms) {
  // Hinge rotates domain 2 by +pi/2 about z through origin.
  // Domain 1 atoms (the parent / root) should sit unchanged.
  ReceptorDomainConfig cfg = MakeTwoDomainHinge(/*split=*/2, /*total=*/4,
                                                /*rest_angle=*/kPi / 2.0);
  DomainPartition p;
  ASSERT_TRUE(p.build(cfg, 4, nullptr));
  DomainGraph graph;
  ASSERT_TRUE(cfg.build_graph(&graph, nullptr));

  const double x_in[] = {1.0, 2.0, 1.0, 0.0};
  const double y_in[] = {0.0, 0.0, 0.0, 1.0};
  const double z_in[] = {0.0, 0.0, 5.0, 5.0};
  double x_out[4], y_out[4], z_out[4];
  ASSERT_TRUE(p.apply(graph, x_in, y_in, z_in, x_out, y_out, z_out, nullptr));

  // Parent (domain 1) atoms unchanged.
  EXPECT_NEAR(x_out[0], 1.0, kEps);
  EXPECT_NEAR(y_out[0], 0.0, kEps);
  EXPECT_NEAR(z_out[0], 0.0, kEps);
  EXPECT_NEAR(x_out[1], 2.0, kEps);
  EXPECT_NEAR(y_out[1], 0.0, kEps);
  EXPECT_NEAR(z_out[1], 0.0, kEps);
  // Child (domain 2) atoms rotated +90 deg about z: (x,y) -> (-y, x).
  EXPECT_NEAR(x_out[2], 0.0, kEps);
  EXPECT_NEAR(y_out[2], 1.0, kEps);
  EXPECT_NEAR(z_out[2], 5.0, kEps);
  EXPECT_NEAR(x_out[3], -1.0, kEps);
  EXPECT_NEAR(y_out[3], 0.0, kEps);
  EXPECT_NEAR(z_out[3], 5.0, kEps);
}

TEST(DomainPartition, ApplySupportsAliasedOutput) {
  // Pass the same pointers as input and output; per-atom writes are
  // independent so the result must match a non-aliased call.
  ReceptorDomainConfig cfg = MakeTwoDomainHinge(/*split=*/1, /*total=*/2,
                                                /*rest_angle=*/kPi / 2.0);
  DomainPartition p;
  ASSERT_TRUE(p.build(cfg, 2, nullptr));
  DomainGraph graph;
  ASSERT_TRUE(cfg.build_graph(&graph, nullptr));

  double x[] = {3.0, 1.0};
  double y[] = {0.0, 0.0};
  double z[] = {7.0, 7.0};
  ASSERT_TRUE(p.apply(graph, x, y, z, x, y, z, nullptr));
  EXPECT_NEAR(x[0], 3.0, kEps);
  EXPECT_NEAR(y[0], 0.0, kEps);
  EXPECT_NEAR(x[1], 0.0, kEps);
  EXPECT_NEAR(y[1], 1.0, kEps);
}

TEST(DomainPartition, ApplyFailsWhenGraphMissingDomain) {
  ReceptorDomainConfig cfg = MakeTwoDomainHinge(/*split=*/1, /*total=*/2,
                                                /*rest_angle=*/0.0);
  DomainPartition p;
  ASSERT_TRUE(p.build(cfg, 2, nullptr));
  // Build an empty graph that does not contain domain 1 or 2.
  DomainGraph graph;
  double x[2] = {0.0, 0.0};
  double y[2] = {0.0, 0.0};
  double z[2] = {0.0, 0.0};
  std::string err;
  EXPECT_FALSE(p.apply(graph, x, y, z, x, y, z, &err));
  EXPECT_NE(err.find("missing"), std::string::npos);
}

TEST(DomainPartition, ApplyCopiesUnclaimedAtomsUnchanged) {
  // Mostly unclaimed receptor (single mid-range domain). Unclaimed
  // atoms must come out bit-for-bit equal.
  ReceptorDomainConfig cfg;
  DomainSpec d;
  d.id = 1;
  d.first_atom_index = 2;
  d.last_atom_index = 4;
  cfg.domains.push_back(d);
  DomainPartition p;
  ASSERT_TRUE(p.build(cfg, 5, nullptr));
  DomainGraph graph;
  ASSERT_TRUE(cfg.build_graph(&graph, nullptr));

  const double x_in[] = {1.0, 2.0, 3.0, 4.0, 5.0};
  const double y_in[] = {10.0, 20.0, 30.0, 40.0, 50.0};
  const double z_in[] = {-1.0, -2.0, -3.0, -4.0, -5.0};
  double x_out[5], y_out[5], z_out[5];
  ASSERT_TRUE(p.apply(graph, x_in, y_in, z_in, x_out, y_out, z_out, nullptr));
  for (int i : {0, 1, 4}) {
    EXPECT_DOUBLE_EQ(x_out[i], x_in[i]);
    EXPECT_DOUBLE_EQ(y_out[i], y_in[i]);
    EXPECT_DOUBLE_EQ(z_out[i], z_in[i]);
  }
  // Domain 1 is the root with no joints (rest pose identity).
  for (int i : {2, 3}) {
    EXPECT_DOUBLE_EQ(x_out[i], x_in[i]);
    EXPECT_DOUBLE_EQ(y_out[i], y_in[i]);
    EXPECT_DOUBLE_EQ(z_out[i], z_in[i]);
  }
}
