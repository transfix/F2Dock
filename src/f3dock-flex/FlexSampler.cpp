/*
  Copyright 2026 The University of Texas at Austin

        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of F2Dock.

  F2Dock is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  F2Dock is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/

#include "f3dock/flex/FlexSampler.h"

#include "f3dock/util/StringUtil.h"

#include <cmath>
#include <cstdio>
#include <cstring>

namespace f3dock {
namespace flex {

bool FlexSamplingConfig::enabled() const {
  for (const auto &h : hinges) {
    if (h.num_steps > 0) {
      return true;
    }
  }
  for (const auto &s : shears) {
    if (s.num_steps > 0) {
      return true;
    }
  }
  return false;
}

std::size_t FlexSamplingConfig::state_count() const {
  std::size_t total = 1;
  for (const auto &h : hinges) {
    if (h.num_steps > 0) {
      total *= h.num_steps;
    }
  }
  for (const auto &s : shears) {
    if (s.num_steps > 0) {
      total *= s.num_steps;
    }
  }
  return total;
}

namespace {

bool axis_is_nondegenerate(const Point3 &v) {
  const double n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  return n2 > 1e-24;
}

} // namespace

bool FlexSampler::enumerate(const FlexSamplingConfig &config,
                            std::vector<FlexState> *out) {
  if (out == nullptr) {
    return false;
  }

  // Validate enabled knobs before allocating output.
  for (const auto &h : config.hinges) {
    if (h.num_steps > 0 && !axis_is_nondegenerate(h.axis_direction)) {
      out->clear();
      return false;
    }
  }
  for (const auto &s : config.shears) {
    if (s.num_steps > 0 && !axis_is_nondegenerate(s.plane_normal)) {
      out->clear();
      return false;
    }
  }

  // Build the list of enabled "knobs" with their step counts so we can
  // iterate a mixed-radix counter across them.
  struct Knob {
    bool is_hinge;
    std::size_t index; // index into config.hinges or config.shears
    std::size_t radix; // num_steps
  };
  std::vector<Knob> knobs;
  knobs.reserve(config.hinges.size() + config.shears.size());
  for (std::size_t i = 0; i < config.hinges.size(); ++i) {
    if (config.hinges[i].num_steps > 0) {
      knobs.push_back({true, i, config.hinges[i].num_steps});
    }
  }
  for (std::size_t i = 0; i < config.shears.size(); ++i) {
    if (config.shears[i].num_steps > 0) {
      knobs.push_back({false, i, config.shears[i].num_steps});
    }
  }

  const std::size_t total = config.state_count();
  out->clear();
  out->reserve(total);

  // Mixed-radix counter: digit[i] in [0, knobs[i].radix). Increment from
  // the last digit so emitted order is lexicographic over (hinge0,
  // hinge1, ..., shear0, shear1, ...).
  std::vector<std::size_t> digits(knobs.size(), 0);

  for (std::size_t s = 0; s < total; ++s) {
    FlexState state;
    state.state_id = s;
    for (std::size_t k = 0; k < knobs.size(); ++k) {
      const Knob &kn = knobs[k];
      const std::size_t step = digits[k];
      if (kn.is_hinge) {
        const HingeAxisSpec &h = config.hinges[kn.index];
        HingeSpec hs;
        hs.axis_point = h.axis_point;
        hs.axis_direction = h.axis_direction;
        hs.angle_radians =
            h.min_angle_radians + static_cast<double>(step) * h.step_radians;
        state.hinges.push_back(hs);
      } else {
        const ShearPlaneSpec &sp = config.shears[kn.index];
        ShearSpec ss;
        ss.plane_normal = sp.plane_normal;
        ss.shear_direction = sp.shear_direction;
        ss.shear_factor =
            sp.min_shear_factor + static_cast<double>(step) * sp.step;
        state.shears.push_back(ss);
      }
    }
    out->push_back(std::move(state));

    // Increment mixed-radix counter (from least-significant = last).
    for (std::size_t k = knobs.size(); k-- > 0;) {
      if (++digits[k] < knobs[k].radix) {
        break;
      }
      digits[k] = 0;
    }
  }

  return true;
}

namespace {

// Parse exactly N space-separated doubles from `value`. Returns true
// on success. Extra whitespace is tolerated; trailing garbage after
// N numbers is an error.
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
    const double v = std::strtod(p, &end);
    if (end == p) {
      return false;
    }
    out[i] = v;
    p = end;
  }
  // Trailing whitespace only.
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') {
    ++p;
  }
  return *p == '\0';
}

} // namespace

bool FlexSampler::parse_param(const char *key, const char *value,
                              FlexSamplingConfig *config, bool *ok) {
  if (key == nullptr || config == nullptr || ok == nullptr) {
    return false;
  }

  if (f3dock::util::iequals(key, "flexHingeAxis")) {
    double v[9] = {0};
    if (!parse_doubles(value, v, 9)) {
      std::printf("Error: flexHingeAxis expects 9 numbers: "
                  "ax ay az dx dy dz minAngleRad stepRad numSteps\n");
      *ok = false;
      return true;
    }
    if (v[8] < 0.0) {
      std::printf("Error: flexHingeAxis numSteps must be >= 0\n");
      *ok = false;
      return true;
    }
    HingeAxisSpec h;
    h.axis_point = {v[0], v[1], v[2]};
    h.axis_direction = {v[3], v[4], v[5]};
    h.min_angle_radians = v[6];
    h.step_radians = v[7];
    h.num_steps = static_cast<std::size_t>(v[8]);
    config->hinges.push_back(h);
    *ok = true;
    return true;
  }

  if (f3dock::util::iequals(key, "flexShearPlane")) {
    double v[9] = {0};
    if (!parse_doubles(value, v, 9)) {
      std::printf("Error: flexShearPlane expects 9 numbers: "
                  "nx ny nz dx dy dz minFactor step numSteps\n");
      *ok = false;
      return true;
    }
    if (v[8] < 0.0) {
      std::printf("Error: flexShearPlane numSteps must be >= 0\n");
      *ok = false;
      return true;
    }
    ShearPlaneSpec s;
    s.plane_normal = {v[0], v[1], v[2]};
    s.shear_direction = {v[3], v[4], v[5]};
    s.min_shear_factor = v[6];
    s.step = v[7];
    s.num_steps = static_cast<std::size_t>(v[8]);
    config->shears.push_back(s);
    *ok = true;
    return true;
  }

  return false;
}

} // namespace flex
} // namespace f3dock
