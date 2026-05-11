#include "f3dock/icp/IcpAligner.h"

#include <cmath>
#include <vector>

#include "libicp/icpPointToPlane.h"
#include "libicp/icpPointToPoint.h"
#include "libicp/matrix.h"

namespace f3dock::icp {
namespace {

double computeRmsd(const std::vector<Point3> &a, const std::vector<Point3> &b) {
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

std::vector<double> flatten(const std::vector<Point3> &points) {
  std::vector<double> data(points.size() * 3);
  for (size_t i = 0; i < points.size(); ++i) {
    data[i * 3 + 0] = points[i][0];
    data[i * 3 + 1] = points[i][1];
    data[i * 3 + 2] = points[i][2];
  }
  return data;
}

void applyTransform(const libicp::Matrix &rotation,
                    const libicp::Matrix &translation,
                    const std::vector<Point3> &moving,
                    std::vector<Point3> &aligned) {
  aligned.resize(moving.size());
  for (size_t i = 0; i < moving.size(); ++i) {
    const double x = moving[i][0];
    const double y = moving[i][1];
    const double z = moving[i][2];

    aligned[i][0] = rotation.val[0][0] * x + rotation.val[0][1] * y +
                    rotation.val[0][2] * z + translation.val[0][0];
    aligned[i][1] = rotation.val[1][0] * x + rotation.val[1][1] * y +
                    rotation.val[1][2] * z + translation.val[1][0];
    aligned[i][2] = rotation.val[2][0] * x + rotation.val[2][1] * y +
                    rotation.val[2][2] * z + translation.val[2][0];
  }
}

void fillResult(const libicp::Matrix &rotation,
                const libicp::Matrix &translation,
                const std::vector<Point3> &model,
                const std::vector<Point3> &moving,
                const std::vector<Point3> &aligned,
                IcpAlignmentResult *result) {
  if (result == nullptr) {
    return;
  }
  result->rmsd_before = computeRmsd(model, moving);
  result->rmsd_after = computeRmsd(model, aligned);
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      result->rotation[row][col] = rotation.val[row][col];
    }
    result->translation[row] = translation.val[row][0];
  }
}

} // namespace

bool IcpAligner::alignPointToPoint(const std::vector<Point3> &model,
                                   const std::vector<Point3> &moving,
                                   std::vector<Point3> &aligned,
                                   IcpAlignmentResult *result,
                                   double inlier_distance) const {
  aligned.clear();

  if (model.size() < static_cast<size_t>(kMinPoints) ||
      moving.size() < static_cast<size_t>(kMinPoints)) {
    return false;
  }

  std::vector<double> model_data = flatten(model);
  std::vector<double> moving_data = flatten(moving);

  libicp::Matrix rotation = libicp::Matrix::eye(3);
  libicp::Matrix translation(3, 1);
  libicp::IcpPointToPoint icp(model_data.data(),
                              static_cast<int32_t>(model.size()), 3);
  icp.fit(moving_data.data(), static_cast<int32_t>(moving.size()), rotation,
          translation, inlier_distance);

  applyTransform(rotation, translation, moving, aligned);
  fillResult(rotation, translation, model, moving, aligned, result);
  return true;
}

bool IcpAligner::alignPointToPlane(const std::vector<Point3> &model,
                                   const std::vector<Point3> &moving,
                                   std::vector<Point3> &aligned,
                                   IcpAlignmentResult *result,
                                   double inlier_distance, int num_neighbors,
                                   double flatness) const {
  aligned.clear();

  if (model.size() < static_cast<size_t>(kMinPoints) ||
      moving.size() < static_cast<size_t>(kMinPoints) || num_neighbors < 3) {
    return false;
  }
  // Cap neighbor count to the available model size.
  if (static_cast<size_t>(num_neighbors) > model.size()) {
    num_neighbors = static_cast<int>(model.size());
  }

  std::vector<double> model_data = flatten(model);
  std::vector<double> moving_data = flatten(moving);

  libicp::Matrix rotation = libicp::Matrix::eye(3);
  libicp::Matrix translation(3, 1);
  libicp::IcpPointToPlane icp(model_data.data(),
                              static_cast<int32_t>(model.size()), 3,
                              num_neighbors, flatness);
  icp.fit(moving_data.data(), static_cast<int32_t>(moving.size()), rotation,
          translation, inlier_distance);

  applyTransform(rotation, translation, moving, aligned);
  fillResult(rotation, translation, model, moving, aligned, result);
  return true;
}

} // namespace f3dock::icp
