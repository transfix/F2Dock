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

  // Point-to-plane ICP. The model surface normals are estimated from the
  // model point cloud via local PCA over the `num_neighbors` nearest
  // neighbors of each model point; `flatness` is reserved for future
  // anisotropy weighting and accepted for libicp parity.
  // Requires at least kMinPoints in each cloud and num_neighbors >= 3.
  bool alignPointToPlane(const std::vector<Point3> &model,
                         const std::vector<Point3> &moving,
                         std::vector<Point3> &aligned,
                         IcpAlignmentResult *result = nullptr,
                         double inlier_distance = 0.0, int num_neighbors = 10,
                         double flatness = 5.0) const;
};

} // namespace f3dock::icp
