#pragma once

#include <array>
#include <vector>

namespace f3dock::icp {

using Point3 = std::array<double, 3>;

struct IcpAlignmentResult {
  std::array<std::array<double, 3>, 3> rotation{};
  Point3 translation{};
  double rmsd_before = 0.0;
  double rmsd_after = 0.0;
};

class IcpAligner {
public:
  static constexpr int kMinPoints = 5;

  bool alignPointToPoint(const std::vector<Point3> &model,
                         const std::vector<Point3> &moving,
                         std::vector<Point3> &aligned,
                         IcpAlignmentResult *result = nullptr,
                         double inlier_distance = 0.0) const;
};

} // namespace f3dock::icp
