// Unit tests for the Phase 4 task-3 receptor-domain hinge sampler.
// Mirrors the pattern used by FlexSampler / DomainPartition tests.

#include "f3dock/domain/DomainGraph.h"
#include "f3dock/domain/DomainPartition.h"
#include "f3dock/domain/DomainSampler.h"
#include "f3dock/domain/DomainSpec.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using f3dock::domain::build_graph_with_state;
using f3dock::domain::DomainGraph;
using f3dock::domain::DomainHingeSweep;
using f3dock::domain::DomainPartition;
using f3dock::domain::DomainSampler;
using f3dock::domain::DomainSamplingConfig;
using f3dock::domain::DomainSpec;
using f3dock::domain::DomainState;
using f3dock::domain::HingeJointAngle;
using f3dock::domain::JointSpec;
using f3dock::domain::JointType;
using f3dock::domain::ReceptorDomainConfig;

namespace {

constexpr double kEps = 1e-9;
constexpr double kPi = 3.14159265358979323846;

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

TEST(DomainSampler, EmptyConfigEmitsSingleState) {
  DomainSamplingConfig cfg;
  EXPECT_FALSE(cfg.enabled());
  EXPECT_EQ(cfg.state_count(), 1u);
  std::vector<DomainState> states;
  ASSERT_TRUE(DomainSampler::enumerate(cfg, &states));
  ASSERT_EQ(states.size(), 1u);
  EXPECT_EQ(states[0].state_id, 0u);
  EXPECT_TRUE(states[0].hinge_angles.empty());
}

TEST(DomainSampler, ZeroStepHingeIsDisabled) {
  DomainSamplingConfig cfg;
  DomainHingeSweep h;
  h.parent_id = 1;
  h.child_id = 2;
  h.min_angle_radians = 0.0;
  h.step_radians = 0.5;
  h.num_steps = 0;
  cfg.hinges.push_back(h);
  EXPECT_FALSE(cfg.enabled());
  EXPECT_EQ(cfg.state_count(), 1u);
  std::vector<DomainState> states;
  ASSERT_TRUE(DomainSampler::enumerate(cfg, &states));
  ASSERT_EQ(states.size(), 1u);
  EXPECT_TRUE(states[0].hinge_angles.empty());
}

TEST(DomainSampler, SingleHingeSweep) {
  DomainSamplingConfig cfg;
  DomainHingeSweep h;
  h.parent_id = 1;
  h.child_id = 2;
  h.min_angle_radians = -1.0;
  h.step_radians = 0.5;
  h.num_steps = 5;
  cfg.hinges.push_back(h);
  EXPECT_TRUE(cfg.enabled());
  EXPECT_EQ(cfg.state_count(), 5u);
  std::vector<DomainState> states;
  ASSERT_TRUE(DomainSampler::enumerate(cfg, &states));
  ASSERT_EQ(states.size(), 5u);
  for (std::size_t i = 0; i < states.size(); ++i) {
    EXPECT_EQ(states[i].state_id, i);
    ASSERT_EQ(states[i].hinge_angles.size(), 1u);
    EXPECT_EQ(states[i].hinge_angles[0].parent_id, 1);
    EXPECT_EQ(states[i].hinge_angles[0].child_id, 2);
    EXPECT_NEAR(states[i].hinge_angles[0].angle_radians,
                -1.0 + 0.5 * static_cast<double>(i), kEps);
  }
}

TEST(DomainSampler, CrossProductTwoHinges) {
  DomainSamplingConfig cfg;
  DomainHingeSweep a;
  a.parent_id = 1;
  a.child_id = 2;
  a.min_angle_radians = 0.0;
  a.step_radians = 1.0;
  a.num_steps = 3;
  cfg.hinges.push_back(a);
  DomainHingeSweep b;
  b.parent_id = 1;
  b.child_id = 3;
  b.min_angle_radians = -0.5;
  b.step_radians = 0.25;
  b.num_steps = 4;
  cfg.hinges.push_back(b);
  EXPECT_EQ(cfg.state_count(), 12u);
  std::vector<DomainState> states;
  ASSERT_TRUE(DomainSampler::enumerate(cfg, &states));
  ASSERT_EQ(states.size(), 12u);

  // Last hinge varies fastest. For state_id k, k = i*4 + j with i in
  // [0,3), j in [0,4). hinge[0].angle = i; hinge[1].angle = -0.5 + 0.25*j.
  for (std::size_t k = 0; k < states.size(); ++k) {
    const std::size_t i = k / 4u;
    const std::size_t j = k % 4u;
    ASSERT_EQ(states[k].hinge_angles.size(), 2u);
    EXPECT_NEAR(states[k].hinge_angles[0].angle_radians, static_cast<double>(i),
                kEps);
    EXPECT_NEAR(states[k].hinge_angles[1].angle_radians,
                -0.5 + 0.25 * static_cast<double>(j), kEps);
  }
}

TEST(DomainSampler, ParserRoundTrip) {
  DomainSamplingConfig cfg;
  bool ok = false;
  EXPECT_TRUE(DomainSampler::parse_param("receptorJointHingeSweep",
                                         "1 2 -1.5 0.5 7", &cfg, &ok));
  EXPECT_TRUE(ok);
  ASSERT_EQ(cfg.hinges.size(), 1u);
  EXPECT_EQ(cfg.hinges[0].parent_id, 1);
  EXPECT_EQ(cfg.hinges[0].child_id, 2);
  EXPECT_NEAR(cfg.hinges[0].min_angle_radians, -1.5, kEps);
  EXPECT_NEAR(cfg.hinges[0].step_radians, 0.5, kEps);
  EXPECT_EQ(cfg.hinges[0].num_steps, 7u);
}

TEST(DomainSampler, ParserUnknownKey) {
  DomainSamplingConfig cfg;
  bool ok = true;
  EXPECT_FALSE(DomainSampler::parse_param("somethingElse", "x", &cfg, &ok));
  EXPECT_TRUE(cfg.hinges.empty());
}

TEST(DomainSampler, ParserRejectsBadInput) {
  DomainSamplingConfig cfg;
  bool ok = true;
  EXPECT_TRUE(DomainSampler::parse_param("receptorJointHingeSweep", "1 2 0.0",
                                         &cfg, &ok));
  EXPECT_FALSE(ok);
  ok = true;
  EXPECT_TRUE(DomainSampler::parse_param("receptorJointHingeSweep",
                                         "1 2 0.0 0.1 -3", &cfg, &ok));
  EXPECT_FALSE(ok);
}

TEST(DomainSampler, BuildGraphWithStateOverridesHingeAngle) {
  // Rest pose: rotate child by 0; state pose: rotate by pi/2 about z.
  ReceptorDomainConfig rc = MakeTwoDomainHinge(/*split=*/2, /*total=*/4,
                                               /*rest_angle=*/0.0);
  DomainState state;
  state.state_id = 7;
  HingeJointAngle ov;
  ov.parent_id = 1;
  ov.child_id = 2;
  ov.angle_radians = kPi / 2.0;
  state.hinge_angles.push_back(ov);

  DomainGraph graph;
  std::string err;
  ASSERT_TRUE(build_graph_with_state(rc, state, &graph, &err)) << err;

  f3dock::domain::RigidTransform world;
  ASSERT_TRUE(graph.worldTransform(2, &world));

  // For a rotation of pi/2 about z through the origin, the point
  // (1,0,0) in child-local coords should map to (0,1,0) in world.
  const double rx = world.rotation[0][0] * 1.0 + world.rotation[0][1] * 0.0 +
                    world.rotation[0][2] * 0.0 + world.translation[0];
  const double ry = world.rotation[1][0] * 1.0 + world.rotation[1][1] * 0.0 +
                    world.rotation[1][2] * 0.0 + world.translation[1];
  const double rz = world.rotation[2][0] * 1.0 + world.rotation[2][1] * 0.0 +
                    world.rotation[2][2] * 0.0 + world.translation[2];
  EXPECT_NEAR(rx, 0.0, 1e-9);
  EXPECT_NEAR(ry, 1.0, 1e-9);
  EXPECT_NEAR(rz, 0.0, 1e-9);
}

TEST(DomainSampler, BuildGraphRejectsUnknownJoint) {
  ReceptorDomainConfig rc = MakeTwoDomainHinge(2, 4, 0.0);
  DomainState state;
  HingeJointAngle ov;
  ov.parent_id = 99;
  ov.child_id = 100;
  ov.angle_radians = 0.0;
  state.hinge_angles.push_back(ov);
  DomainGraph graph;
  std::string err;
  EXPECT_FALSE(build_graph_with_state(rc, state, &graph, &err));
  EXPECT_NE(err.find("unknown joint"), std::string::npos);
}

TEST(DomainSampler, BuildGraphRejectsFixedJointOverride) {
  ReceptorDomainConfig rc;
  DomainSpec a;
  a.id = 1;
  a.first_atom_index = 0;
  a.last_atom_index = 2;
  rc.domains.push_back(a);
  DomainSpec b;
  b.id = 2;
  b.first_atom_index = 2;
  b.last_atom_index = 4;
  rc.domains.push_back(b);
  JointSpec j;
  j.parent_id = 1;
  j.child_id = 2;
  j.type = JointType::Fixed;
  rc.joints.push_back(j);

  DomainState state;
  HingeJointAngle ov;
  ov.parent_id = 1;
  ov.child_id = 2;
  ov.angle_radians = 0.5;
  state.hinge_angles.push_back(ov);
  DomainGraph graph;
  std::string err;
  EXPECT_FALSE(build_graph_with_state(rc, state, &graph, &err));
  EXPECT_NE(err.find("not a hinge"), std::string::npos);
}

TEST(DomainSampler, EndToEndPartitionWithSampledState) {
  // 4-atom receptor: domain 1 = atoms [0,2), domain 2 = atoms [2,4).
  // Domain 2 is hinged about z through the origin. Rest angle 0
  // leaves atom 2 == (1,0,0) and atom 3 == (0,1,0) unchanged; the
  // state's pi/2 override rotates them to (0,1,0) and (-1,0,0).
  ReceptorDomainConfig rc = MakeTwoDomainHinge(2, 4, 0.0);
  DomainPartition part;
  std::string err;
  ASSERT_TRUE(part.build(rc, /*total_atoms=*/4, &err)) << err;

  DomainState state;
  HingeJointAngle ov;
  ov.parent_id = 1;
  ov.child_id = 2;
  ov.angle_radians = kPi / 2.0;
  state.hinge_angles.push_back(ov);

  DomainGraph graph;
  ASSERT_TRUE(build_graph_with_state(rc, state, &graph, &err)) << err;

  const double xin[4] = {0.5, -0.5, 1.0, 0.0};
  const double yin[4] = {0.0, 0.0, 0.0, 1.0};
  const double zin[4] = {0.0, 0.0, 0.0, 0.0};
  double xout[4] = {0}, yout[4] = {0}, zout[4] = {0};

  ASSERT_TRUE(part.apply(graph, xin, yin, zin, xout, yout, zout, &err)) << err;

  // Domain 1 atoms are unchanged.
  EXPECT_NEAR(xout[0], 0.5, kEps);
  EXPECT_NEAR(yout[0], 0.0, kEps);
  EXPECT_NEAR(xout[1], -0.5, kEps);
  EXPECT_NEAR(yout[1], 0.0, kEps);

  // Domain 2 atoms are rotated +90 deg about z through the origin.
  EXPECT_NEAR(xout[2], 0.0, 1e-9);
  EXPECT_NEAR(yout[2], 1.0, 1e-9);
  EXPECT_NEAR(zout[2], 0.0, kEps);
  EXPECT_NEAR(xout[3], -1.0, 1e-9);
  EXPECT_NEAR(yout[3], 0.0, 1e-9);
  EXPECT_NEAR(zout[3], 0.0, kEps);
}
