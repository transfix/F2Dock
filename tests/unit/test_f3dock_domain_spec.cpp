// Unit tests for the Phase 4 task-1 receptor domain-spec parser and
// validator. Parser is a static class mirroring f3dock::flex::FlexSampler.
//
// Tests focus on the *configuration* layer: PARAMS_IN-style key/value
// parsing, semantic validation (duplicate ids, unknown joints, multiple
// roots, cycles), and `build_graph()` producing a DomainGraph whose
// `worldTransform` reflects the joints' rest-pose transforms.

#include "f3dock/domain/DomainGraph.h"
#include "f3dock/domain/DomainSpec.h"

#include <cmath>
#include <string>

#include <gtest/gtest.h>

using f3dock::domain::DomainGraph;
using f3dock::domain::DomainSpec;
using f3dock::domain::DomainSpecParser;
using f3dock::domain::JointSpec;
using f3dock::domain::JointType;
using f3dock::domain::Point3;
using f3dock::domain::ReceptorDomainConfig;
using f3dock::domain::RigidTransform;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEps = 1e-9;

void ExpectPointNear(const Point3 &a, const Point3 &b, double eps = kEps) {
  EXPECT_NEAR(a[0], b[0], eps);
  EXPECT_NEAR(a[1], b[1], eps);
  EXPECT_NEAR(a[2], b[2], eps);
}

} // namespace

// ─────────────────────────────────────────────────────────────────────
// Defaults / enablement
// ─────────────────────────────────────────────────────────────────────

TEST(DomainSpec, DefaultConfigIsDisabled) {
  ReceptorDomainConfig cfg;
  EXPECT_FALSE(cfg.enabled());
  std::string err;
  EXPECT_EQ(cfg.find_root(&err), -1);
  EXPECT_FALSE(err.empty());
}

// ─────────────────────────────────────────────────────────────────────
// Parser: receptorDomain
// ─────────────────────────────────────────────────────────────────────

TEST(DomainSpecParser, ParsesReceptorDomain) {
  ReceptorDomainConfig cfg;
  bool ok = false;
  EXPECT_TRUE(
      DomainSpecParser::parse_param("receptorDomain", "7 0 1500", &cfg, &ok));
  EXPECT_TRUE(ok);
  ASSERT_EQ(cfg.domains.size(), 1u);
  EXPECT_EQ(cfg.domains[0].id, 7);
  EXPECT_EQ(cfg.domains[0].first_atom_index, 0);
  EXPECT_EQ(cfg.domains[0].last_atom_index, 1500);
  EXPECT_TRUE(cfg.enabled());
}

TEST(DomainSpecParser, ReceptorDomainIsCaseInsensitive) {
  ReceptorDomainConfig cfg;
  bool ok = false;
  EXPECT_TRUE(
      DomainSpecParser::parse_param("ReCePtOrDoMaIn", "1 0 10", &cfg, &ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(cfg.domains.size(), 1u);
}

TEST(DomainSpecParser, RejectsWrongTokenCount) {
  ReceptorDomainConfig cfg;
  bool ok = true;
  EXPECT_TRUE(
      DomainSpecParser::parse_param("receptorDomain", "1 2", &cfg, &ok));
  EXPECT_FALSE(ok);
  EXPECT_TRUE(cfg.domains.empty());
}

TEST(DomainSpecParser, RejectsInvertedRange) {
  ReceptorDomainConfig cfg;
  bool ok = true;
  EXPECT_TRUE(
      DomainSpecParser::parse_param("receptorDomain", "1 100 50", &cfg, &ok));
  EXPECT_FALSE(ok);
}

// ─────────────────────────────────────────────────────────────────────
// Parser: receptorJointFixed / receptorJointHinge
// ─────────────────────────────────────────────────────────────────────

TEST(DomainSpecParser, ParsesReceptorJointFixed) {
  ReceptorDomainConfig cfg;
  bool ok = false;
  EXPECT_TRUE(
      DomainSpecParser::parse_param("receptorJointFixed", "1 2", &cfg, &ok));
  EXPECT_TRUE(ok);
  ASSERT_EQ(cfg.joints.size(), 1u);
  EXPECT_EQ(cfg.joints[0].parent_id, 1);
  EXPECT_EQ(cfg.joints[0].child_id, 2);
  EXPECT_EQ(cfg.joints[0].type, JointType::Fixed);
}

TEST(DomainSpecParser, ParsesReceptorJointHinge) {
  ReceptorDomainConfig cfg;
  bool ok = false;
  // parent=1 child=2 axis_point=(1,0,0) axis_dir=(0,0,1) initAngle=0
  EXPECT_TRUE(DomainSpecParser::parse_param("receptorJointHinge",
                                            "1 2 1 0 0 0 0 1 0", &cfg, &ok));
  EXPECT_TRUE(ok);
  ASSERT_EQ(cfg.joints.size(), 1u);
  const JointSpec &j = cfg.joints[0];
  EXPECT_EQ(j.parent_id, 1);
  EXPECT_EQ(j.child_id, 2);
  EXPECT_EQ(j.type, JointType::Hinge);
  ExpectPointNear(j.axis_point, {1.0, 0.0, 0.0});
  ExpectPointNear(j.axis_direction, {0.0, 0.0, 1.0});
  EXPECT_DOUBLE_EQ(j.initial_angle_radians, 0.0);
}

TEST(DomainSpecParser, IgnoresUnrelatedKeys) {
  ReceptorDomainConfig cfg;
  bool ok = true;
  EXPECT_FALSE(
      DomainSpecParser::parse_param("flexHingeAxis", "1 2 3", &cfg, &ok));
  EXPECT_TRUE(ok); // untouched
  EXPECT_TRUE(cfg.domains.empty());
  EXPECT_TRUE(cfg.joints.empty());
}

TEST(DomainSpecParser, RejectsNullInputs) {
  ReceptorDomainConfig cfg;
  bool ok = true;
  EXPECT_FALSE(DomainSpecParser::parse_param(nullptr, "x", &cfg, &ok));
  EXPECT_FALSE(
      DomainSpecParser::parse_param("receptorDomain", "x", nullptr, &ok));
  EXPECT_FALSE(
      DomainSpecParser::parse_param("receptorDomain", "x", &cfg, nullptr));
}

// ─────────────────────────────────────────────────────────────────────
// Validation rules
// ─────────────────────────────────────────────────────────────────────

TEST(ReceptorDomainConfig, FindsSingleRoot) {
  ReceptorDomainConfig cfg;
  cfg.domains.push_back(DomainSpec{1, 0, 10});
  cfg.domains.push_back(DomainSpec{2, 10, 20});
  JointSpec j;
  j.parent_id = 1;
  j.child_id = 2;
  j.type = JointType::Fixed;
  cfg.joints.push_back(j);

  std::string err;
  EXPECT_EQ(cfg.find_root(&err), 1);
  EXPECT_TRUE(err.empty());
}

TEST(ReceptorDomainConfig, DuplicateDomainIdRejected) {
  ReceptorDomainConfig cfg;
  cfg.domains.push_back(DomainSpec{1, 0, 10});
  cfg.domains.push_back(DomainSpec{1, 10, 20});

  std::string err;
  EXPECT_EQ(cfg.find_root(&err), -1);
  EXPECT_NE(err.find("duplicate"), std::string::npos);
}

TEST(ReceptorDomainConfig, JointWithUnknownDomainRejected) {
  ReceptorDomainConfig cfg;
  cfg.domains.push_back(DomainSpec{1, 0, 10});
  JointSpec j;
  j.parent_id = 1;
  j.child_id = 42; // unknown
  cfg.joints.push_back(j);

  std::string err;
  EXPECT_EQ(cfg.find_root(&err), -1);
  EXPECT_NE(err.find("unknown"), std::string::npos);
}

TEST(ReceptorDomainConfig, MultiParentRejected) {
  ReceptorDomainConfig cfg;
  cfg.domains.push_back(DomainSpec{1, 0, 10});
  cfg.domains.push_back(DomainSpec{2, 10, 20});
  cfg.domains.push_back(DomainSpec{3, 20, 30});
  JointSpec j1{1, 3, JointType::Fixed, {0, 0, 0}, {0, 0, 1}, 0.0};
  JointSpec j2{2, 3, JointType::Fixed, {0, 0, 0}, {0, 0, 1}, 0.0};
  cfg.joints.push_back(j1);
  cfg.joints.push_back(j2);

  std::string err;
  EXPECT_EQ(cfg.find_root(&err), -1);
  EXPECT_NE(err.find("more than one parent"), std::string::npos);
}

TEST(ReceptorDomainConfig, MultiRootRejected) {
  ReceptorDomainConfig cfg;
  cfg.domains.push_back(DomainSpec{1, 0, 10});
  cfg.domains.push_back(DomainSpec{2, 10, 20});
  // No joints -> both domains are roots.
  std::string err;
  EXPECT_EQ(cfg.find_root(&err), -1);
  EXPECT_NE(err.find("more than one root"), std::string::npos);
}

// ─────────────────────────────────────────────────────────────────────
// build_graph() produces a DomainGraph whose rest-pose transforms match
// ─────────────────────────────────────────────────────────────────────

TEST(ReceptorDomainConfig, BuildGraphFixedJointIsIdentity) {
  ReceptorDomainConfig cfg;
  cfg.domains.push_back(DomainSpec{10, 0, 5});
  cfg.domains.push_back(DomainSpec{20, 5, 10});
  JointSpec j;
  j.parent_id = 10;
  j.child_id = 20;
  j.type = JointType::Fixed;
  cfg.joints.push_back(j);

  DomainGraph g;
  std::string err;
  ASSERT_TRUE(cfg.build_graph(&g, &err)) << err;

  RigidTransform t;
  ASSERT_TRUE(g.worldTransform(20, &t));
  // Fixed joint with parent at world identity -> child world is identity.
  ExpectPointNear(t.translation, {0.0, 0.0, 0.0});
  EXPECT_NEAR(t.rotation[0][0], 1.0, kEps);
  EXPECT_NEAR(t.rotation[1][1], 1.0, kEps);
  EXPECT_NEAR(t.rotation[2][2], 1.0, kEps);
}

TEST(ReceptorDomainConfig, BuildGraphHingeRotatesAboutAxisLine) {
  // Hinge axis: line through (1,0,0) with direction +z, initial angle = pi/2.
  // Child-local point (1,0,0) lies on the axis line -> should map to (1,0,0)
  // in the parent frame regardless of angle.
  // Child-local point (2,0,0) is one unit off the axis along +x; after
  // pi/2 about +z it should be at (1,1,0) in the parent frame.
  ReceptorDomainConfig cfg;
  cfg.domains.push_back(DomainSpec{1, 0, 5});
  cfg.domains.push_back(DomainSpec{2, 5, 10});

  JointSpec j;
  j.parent_id = 1;
  j.child_id = 2;
  j.type = JointType::Hinge;
  j.axis_point = {1.0, 0.0, 0.0};
  j.axis_direction = {0.0, 0.0, 1.0};
  j.initial_angle_radians = kPi / 2.0;
  cfg.joints.push_back(j);

  DomainGraph g;
  std::string err;
  ASSERT_TRUE(cfg.build_graph(&g, &err)) << err;

  RigidTransform t;
  ASSERT_TRUE(g.worldTransform(2, &t));

  // (1,0,0) on the axis -> fixed.
  Point3 on_axis = t.apply({1.0, 0.0, 0.0});
  ExpectPointNear(on_axis, {1.0, 0.0, 0.0}, 1e-9);

  // (2,0,0) -> (1,1,0).
  Point3 off_axis = t.apply({2.0, 0.0, 0.0});
  ExpectPointNear(off_axis, {1.0, 1.0, 0.0}, 1e-9);
}

TEST(ReceptorDomainConfig, BuildGraphRejectsZeroLengthHingeAxis) {
  ReceptorDomainConfig cfg;
  cfg.domains.push_back(DomainSpec{1, 0, 5});
  cfg.domains.push_back(DomainSpec{2, 5, 10});
  JointSpec j;
  j.parent_id = 1;
  j.child_id = 2;
  j.type = JointType::Hinge;
  j.axis_direction = {0.0, 0.0, 0.0};
  cfg.joints.push_back(j);

  DomainGraph g;
  std::string err;
  EXPECT_FALSE(cfg.build_graph(&g, &err));
  EXPECT_NE(err.find("axis direction"), std::string::npos);
}
