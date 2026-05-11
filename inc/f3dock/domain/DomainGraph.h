#pragma once

#include <array>
#include <map>
#include <vector>

namespace f3dock {
namespace domain {

using Point3 = std::array<double, 3>;
using Matrix3 = std::array<std::array<double, 3>, 3>;

struct RigidTransform {
  Matrix3 rotation;
  Point3 translation;

  static RigidTransform identity();

  Point3 apply(const Point3 &point) const;
  RigidTransform compose(const RigidTransform &other) const;
};

class DomainGraph {
public:
  bool addDomain(int id);

  // Sets parent relationship and child-local transform (parent <- child).
  bool setParent(int child_id, int parent_id,
                 const RigidTransform &child_to_parent);

  bool worldTransform(int id, RigidTransform *out) const;

private:
  struct Node {
    int id = -1;
    int parent_id = -1;
    bool has_parent = false;
    RigidTransform child_to_parent = RigidTransform::identity();
  };

  bool hasCycleFrom(int id) const;

  std::map<int, Node> nodes_;
};

} // namespace domain
} // namespace f3dock
