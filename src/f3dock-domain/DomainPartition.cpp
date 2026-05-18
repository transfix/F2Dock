#include "f3dock/domain/DomainPartition.h"

#include <algorithm>
#include <sstream>

namespace f3dock {
namespace domain {

bool DomainPartition::build(const ReceptorDomainConfig &config,
                            std::size_t total_atoms, std::string *error) {
  enabled_ = false;
  atom_to_domain_.assign(total_atoms, -1);

  if (!config.enabled()) {
    return true;
  }

  for (const auto &d : config.domains) {
    if (d.first_atom_index < 0 || d.last_atom_index < d.first_atom_index) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "receptor domain " << d.id << " has invalid range ["
           << d.first_atom_index << ", " << d.last_atom_index << ")";
        *error = os.str();
      }
      atom_to_domain_.clear();
      return false;
    }
    if (static_cast<std::size_t>(d.last_atom_index) > total_atoms) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "receptor domain " << d.id << " range [" << d.first_atom_index
           << ", " << d.last_atom_index << ") exceeds receptor atom count "
           << total_atoms;
        *error = os.str();
      }
      atom_to_domain_.clear();
      return false;
    }
    for (int i = d.first_atom_index; i < d.last_atom_index; ++i) {
      const std::size_t idx = static_cast<std::size_t>(i);
      if (atom_to_domain_[idx] != -1) {
        if (error != nullptr) {
          std::ostringstream os;
          os << "receptor atom " << i << " is claimed by both domain "
             << atom_to_domain_[idx] << " and domain " << d.id;
          *error = os.str();
        }
        atom_to_domain_.clear();
        return false;
      }
      atom_to_domain_[idx] = d.id;
    }
  }

  enabled_ = true;
  return true;
}

int DomainPartition::domain_for_atom(std::size_t atom_index) const {
  if (atom_index >= atom_to_domain_.size()) {
    return -1;
  }
  return atom_to_domain_[atom_index];
}

bool DomainPartition::apply(const DomainGraph &graph, const double *x_in,
                            const double *y_in, const double *z_in,
                            double *x_out, double *y_out, double *z_out,
                            std::string *error) const {
  if (x_in == nullptr || y_in == nullptr || z_in == nullptr ||
      x_out == nullptr || y_out == nullptr || z_out == nullptr) {
    if (error != nullptr) {
      *error = "DomainPartition::apply: null pointer";
    }
    return false;
  }

  // Cache one RigidTransform per distinct domain id (-1 == identity)
  // so we do not re-query the graph per atom.
  std::vector<int> seen_ids;
  seen_ids.reserve(atom_to_domain_.size());
  for (int id : atom_to_domain_) {
    if (id >= 0) {
      seen_ids.push_back(id);
    }
  }
  std::sort(seen_ids.begin(), seen_ids.end());
  seen_ids.erase(std::unique(seen_ids.begin(), seen_ids.end()), seen_ids.end());

  std::vector<RigidTransform> world(seen_ids.empty() ? 0 : seen_ids.size());
  for (std::size_t k = 0; k < seen_ids.size(); ++k) {
    if (!graph.worldTransform(seen_ids[k], &world[k])) {
      if (error != nullptr) {
        std::ostringstream os;
        os << "domain " << seen_ids[k]
           << " is referenced by the partition but missing from the graph";
        *error = os.str();
      }
      return false;
    }
  }

  for (std::size_t i = 0; i < atom_to_domain_.size(); ++i) {
    const int id = atom_to_domain_[i];
    if (id < 0) {
      x_out[i] = x_in[i];
      y_out[i] = y_in[i];
      z_out[i] = z_in[i];
      continue;
    }
    auto it = std::lower_bound(seen_ids.begin(), seen_ids.end(), id);
    const RigidTransform &t = world[it - seen_ids.begin()];
    const Point3 p = {x_in[i], y_in[i], z_in[i]};
    const Point3 q = t.apply(p);
    x_out[i] = q[0];
    y_out[i] = q[1];
    z_out[i] = q[2];
  }
  return true;
}

} // namespace domain
} // namespace f3dock
