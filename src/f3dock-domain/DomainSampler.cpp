#include "f3dock/domain/DomainSampler.h"

#include "f3dock/util/StringUtil.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>

namespace f3dock {
namespace domain {

namespace {

// Parse `n` space-separated doubles from `value`. Mirrors the helper
// in DomainSpec.cpp so the two files stay independent.
bool parse_doubles(const char *value, double *out, std::size_t n) {
  if (value == nullptr) {
    return false;
  }
  const char *p = value;
  for (std::size_t i = 0; i < n; ++i) {
    while (*p == ' ' || *p == '\t') {
      ++p;
    }
    if (*p == '\0') {
      return false;
    }
    char *end = nullptr;
    out[i] = std::strtod(p, &end);
    if (end == p) {
      return false;
    }
    p = end;
  }
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') {
    ++p;
  }
  return *p == '\0';
}

} // namespace

bool DomainSamplingConfig::enabled() const {
  for (const auto &h : hinges) {
    if (h.num_steps > 0) {
      return true;
    }
  }
  return false;
}

std::size_t DomainSamplingConfig::state_count() const {
  std::size_t n = 1;
  bool any = false;
  for (const auto &h : hinges) {
    if (h.num_steps > 0) {
      n *= h.num_steps;
      any = true;
    }
  }
  return any ? n : 1;
}

bool DomainSampler::enumerate(const DomainSamplingConfig &config,
                              std::vector<DomainState> *out) {
  if (out == nullptr) {
    return false;
  }
  out->clear();

  // Collect enabled hinges and their step counts.
  std::vector<const DomainHingeSweep *> active;
  active.reserve(config.hinges.size());
  for (const auto &h : config.hinges) {
    if (h.num_steps > 0) {
      active.push_back(&h);
    }
  }

  if (active.empty()) {
    DomainState s;
    s.state_id = 0;
    out->push_back(std::move(s));
    return true;
  }

  std::vector<std::size_t> idx(active.size(), 0);
  std::size_t sid = 0;
  while (true) {
    DomainState s;
    s.state_id = sid++;
    s.hinge_angles.reserve(active.size());
    for (std::size_t i = 0; i < active.size(); ++i) {
      HingeJointAngle a;
      a.parent_id = active[i]->parent_id;
      a.child_id = active[i]->child_id;
      a.angle_radians = active[i]->min_angle_radians +
                        static_cast<double>(idx[i]) * active[i]->step_radians;
      s.hinge_angles.push_back(a);
    }
    out->push_back(std::move(s));

    // Advance lexicographically: last index varies fastest.
    std::size_t k = active.size();
    while (k > 0) {
      --k;
      if (++idx[k] < active[k]->num_steps) {
        break;
      }
      idx[k] = 0;
      if (k == 0) {
        return true;
      }
    }
  }
}

bool DomainSampler::parse_param(const char *key, const char *value,
                                DomainSamplingConfig *config, bool *ok) {
  if (key == nullptr || config == nullptr || ok == nullptr) {
    return false;
  }

  const std::string_view key_sv(key);

  if (f3dock::util::iequals(key_sv, "receptorJointHingeSweep")) {
    double v[5] = {0};
    if (!parse_doubles(value, v, 5)) {
      std::printf("Error: receptorJointHingeSweep expects 5 values: "
                  "parentId childId minAngleRad stepRad numSteps\n");
      *ok = false;
      return true;
    }
    DomainHingeSweep s;
    s.parent_id = static_cast<int>(v[0]);
    s.child_id = static_cast<int>(v[1]);
    s.min_angle_radians = v[2];
    s.step_radians = v[3];
    const double ns = v[4];
    if (ns < 0.0) {
      std::printf("Error: receptorJointHingeSweep numSteps must be "
                  ">= 0 (got %g)\n",
                  ns);
      *ok = false;
      return true;
    }
    s.num_steps = static_cast<std::size_t>(ns);
    config->hinges.push_back(s);
    *ok = true;
    return true;
  }

  return false;
}

} // namespace domain
} // namespace f3dock
