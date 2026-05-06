#include "f3dock/flex/FlexKinematics.h"

#include <cmath>

namespace {

using f3dock::flex::Point3;

double dot(const Point3 &a, const Point3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Point3 sub(const Point3 &a, const Point3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Point3 add(const Point3 &a, const Point3 &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Point3 mul(const Point3 &a, double s) {
  return {a[0] * s, a[1] * s, a[2] * s};
}

Point3 cross(const Point3 &a, const Point3 &b) {
  return {
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
  };
}

bool normalize(const Point3 &in, Point3 *out) {
  const double n = std::sqrt(dot(in, in));
  if (n <= 1e-12) {
    return false;
  }
  *out = {in[0] / n, in[1] / n, in[2] / n};
  return true;
}

} // namespace

namespace f3dock {
namespace flex {

bool FlexKinematics::applyHinge(const std::vector<Point3> &input,
                                const HingeSpec &hinge,
                                std::vector<Point3> *output) {
  if (output == nullptr) {
    return false;
  }

  Point3 axis = {0.0, 0.0, 0.0};
  if (!normalize(hinge.axis_direction, &axis)) {
    output->clear();
    return false;
  }

  const double c = std::cos(hinge.angle_radians);
  const double s = std::sin(hinge.angle_radians);

  output->clear();
  output->reserve(input.size());

  for (const Point3 &p : input) {
    const Point3 v = sub(p, hinge.axis_point);
    const Point3 term1 = mul(v, c);
    const Point3 term2 = mul(cross(axis, v), s);
    const Point3 term3 = mul(axis, dot(axis, v) * (1.0 - c));
    output->push_back(add(hinge.axis_point, add(add(term1, term2), term3)));
  }

  return true;
}

bool FlexKinematics::applyShear(const std::vector<Point3> &input,
                                const ShearSpec &shear,
                                std::vector<Point3> *output) {
  if (output == nullptr) {
    return false;
  }

  Point3 normal = {0.0, 0.0, 0.0};
  if (!normalize(shear.plane_normal, &normal)) {
    output->clear();
    return false;
  }

  // Keep shear direction tangential to the shear plane.
  const Point3 projected = sub(shear.shear_direction,
                               mul(normal, dot(shear.shear_direction, normal)));
  Point3 direction = {0.0, 0.0, 0.0};
  if (!normalize(projected, &direction)) {
    output->clear();
    return false;
  }

  output->clear();
  output->reserve(input.size());

  for (const Point3 &p : input) {
    const double signed_distance = dot(p, normal);
    const Point3 displacement =
        mul(direction, shear.shear_factor * signed_distance);
    output->push_back(add(p, displacement));
  }

  return true;
}

} // namespace flex
} // namespace f3dock
