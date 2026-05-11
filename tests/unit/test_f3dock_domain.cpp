#include <cmath>

#include <gtest/gtest.h>

#include "f3dock/domain/DomainGraph.h"

namespace {

using f3dock::domain::DomainGraph;
using f3dock::domain::Point3;
using f3dock::domain::RigidTransform;

RigidTransform makeTranslation(double tx, double ty, double tz) {
  RigidTransform t = RigidTransform::identity();
  t.translation = {tx, ty, tz};
  return t;
}

RigidTransform makeRotationZ(double radians) {
  const double c = std::cos(radians);
  const double s = std::sin(radians);
  RigidTransform t = RigidTransform::identity();
  t.rotation = {{{c, -s, 0.0}, {s, c, 0.0}, {0.0, 0.0, 1.0}}};
  return t;
}

TEST(F3DockDomain, ComputesWorldTransformAcrossTwoEdges) {
  DomainGraph graph;
  ASSERT_TRUE(graph.addDomain(0));
  ASSERT_TRUE(graph.addDomain(1));
  ASSERT_TRUE(graph.addDomain(2));

  // Domain 2 -> Domain 1: rotate around Z by +90 degrees.
  ASSERT_TRUE(graph.setParent(2, 1, makeRotationZ(std::acos(-1.0) * 0.5)));
  // Domain 1 -> Domain 0: translate +10 in x.
  ASSERT_TRUE(graph.setParent(1, 0, makeTranslation(10.0, 0.0, 0.0)));

  RigidTransform world;
  ASSERT_TRUE(graph.worldTransform(2, &world));

  const Point3 p_local = {1.0, 0.0, 0.0};
  const Point3 p_world = world.apply(p_local);

  // Rotated point is (0,1,0), then translated by (+10,0,0).
  EXPECT_NEAR(p_world[0], 10.0, 1e-12);
  EXPECT_NEAR(p_world[1], 1.0, 1e-12);
  EXPECT_NEAR(p_world[2], 0.0, 1e-12);
}

TEST(F3DockDomain, RejectsCycles) {
  DomainGraph graph;
  ASSERT_TRUE(graph.addDomain(10));
  ASSERT_TRUE(graph.addDomain(11));
  ASSERT_TRUE(graph.addDomain(12));

  ASSERT_TRUE(graph.setParent(11, 10, makeTranslation(1.0, 0.0, 0.0)));
  ASSERT_TRUE(graph.setParent(12, 11, makeTranslation(1.0, 0.0, 0.0)));

  // Would create cycle: 10 -> 12 -> 11 -> 10.
  EXPECT_FALSE(graph.setParent(10, 12, makeTranslation(0.0, 0.0, 0.0)));
}

TEST(F3DockDomain, RejectsMissingNodesAndNullOutput) {
  DomainGraph graph;
  ASSERT_TRUE(graph.addDomain(1));
  ASSERT_TRUE(graph.addDomain(2));

  EXPECT_FALSE(graph.setParent(2, 99, makeTranslation(0.0, 0.0, 0.0)));

  RigidTransform world;
  EXPECT_FALSE(graph.worldTransform(99, &world));
  EXPECT_FALSE(graph.worldTransform(1, nullptr));
}

} // namespace
