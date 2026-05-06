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
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
      {1.0, 2.0, 3.0},
      {-2.0, 1.0, 0.5},
      {0.7, -1.2, 2.3},
      {-1.8, -0.4, 1.1},
      {2.5, -1.0, 0.2},
      {1.2, 1.7, -0.8},
  };
}

Point3 applyTransform(const std::array<std::array<double, 3>, 3> &r,
                      const Point3 &t,
                      const Point3 &p) {
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
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
  };
  const std::vector<Point3> moving = model;

  IcpAligner aligner;
  std::vector<Point3> aligned;
  const bool ok = aligner.alignPointToPoint(model, moving, aligned);

  EXPECT_FALSE(ok);
  EXPECT_TRUE(aligned.empty());
}

} // namespace
