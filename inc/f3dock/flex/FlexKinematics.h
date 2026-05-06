#pragma once

#include <array>
#include <vector>

namespace f3dock {
namespace flex {

using Point3 = std::array<double, 3>;

struct HingeSpec {
  Point3 axis_point;
  Point3 axis_direction;
  double angle_radians;
};

struct ShearSpec {
  Point3 plane_normal;
  Point3 shear_direction;
  double shear_factor;
};

class FlexKinematics {
 public:
  // Applies a rigid rotation around an arbitrary axis using Rodrigues formula.
  static bool applyHinge(const std::vector<Point3> &input,
                         const HingeSpec &hinge,
                         std::vector<Point3> *output);

  // Applies planar shear: p' = p + dot(p, n) * factor * d.
  static bool applyShear(const std::vector<Point3> &input,
                         const ShearSpec &shear,
                         std::vector<Point3> *output);
};

} // namespace flex
} // namespace f3dock
