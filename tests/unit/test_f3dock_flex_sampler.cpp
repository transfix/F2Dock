#include <cmath>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "f3dock/flex/FlexSampler.h"

namespace {

using f3dock::flex::FlexSampler;
using f3dock::flex::FlexSamplingConfig;
using f3dock::flex::FlexState;
using f3dock::flex::HingeAxisSpec;
using f3dock::flex::ShearPlaneSpec;

TEST(F3DockFlexSampler, EmptyConfigIsDisabledAndYieldsOneIdentityState) {
  FlexSamplingConfig cfg;
  EXPECT_FALSE(cfg.enabled());
  EXPECT_EQ(cfg.state_count(), 1u);

  std::vector<FlexState> states;
  ASSERT_TRUE(FlexSampler::enumerate(cfg, &states));
  ASSERT_EQ(states.size(), 1u);
  EXPECT_EQ(states[0].state_id, 0u);
  EXPECT_TRUE(states[0].hinges.empty());
  EXPECT_TRUE(states[0].shears.empty());
}

TEST(F3DockFlexSampler, DisabledHinge_ZeroSteps_StaysIdentity) {
  FlexSamplingConfig cfg;
  HingeAxisSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 1.0};
  h.min_angle_radians = 1.0;
  h.step_radians = 0.5;
  h.num_steps = 0;
  cfg.hinges.push_back(h);

  EXPECT_FALSE(cfg.enabled());
  EXPECT_EQ(cfg.state_count(), 1u);

  std::vector<FlexState> states;
  ASSERT_TRUE(FlexSampler::enumerate(cfg, &states));
  ASSERT_EQ(states.size(), 1u);
  EXPECT_TRUE(states[0].hinges.empty());
}

TEST(F3DockFlexSampler, SingleHingeAngleSweep) {
  FlexSamplingConfig cfg;
  HingeAxisSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 1.0};
  h.min_angle_radians = 0.0;
  h.step_radians = 0.25;
  h.num_steps = 4; // 0.00, 0.25, 0.50, 0.75
  cfg.hinges.push_back(h);

  EXPECT_TRUE(cfg.enabled());
  EXPECT_EQ(cfg.state_count(), 4u);

  std::vector<FlexState> states;
  ASSERT_TRUE(FlexSampler::enumerate(cfg, &states));
  ASSERT_EQ(states.size(), 4u);

  for (std::size_t k = 0; k < 4u; ++k) {
    EXPECT_EQ(states[k].state_id, k);
    ASSERT_EQ(states[k].hinges.size(), 1u);
    EXPECT_TRUE(states[k].shears.empty());
    EXPECT_DOUBLE_EQ(states[k].hinges[0].angle_radians,
                     static_cast<double>(k) * 0.25);
  }

  // State 0 must be the identity (angle 0).
  EXPECT_DOUBLE_EQ(states[0].hinges[0].angle_radians, 0.0);
}

TEST(F3DockFlexSampler, CrossProductOfHingeAndShear) {
  FlexSamplingConfig cfg;
  HingeAxisSpec h;
  h.axis_point = {0.0, 0.0, 0.0};
  h.axis_direction = {0.0, 0.0, 1.0};
  h.min_angle_radians = 0.0;
  h.step_radians = 1.0;
  h.num_steps = 3; // 0, 1, 2
  cfg.hinges.push_back(h);

  ShearPlaneSpec s;
  s.plane_normal = {0.0, 0.0, 1.0};
  s.shear_direction = {1.0, 0.0, 0.0};
  s.min_shear_factor = 0.0;
  s.step = 0.5;
  s.num_steps = 2; // 0.0, 0.5
  cfg.shears.push_back(s);

  EXPECT_EQ(cfg.state_count(), 6u);

  std::vector<FlexState> states;
  ASSERT_TRUE(FlexSampler::enumerate(cfg, &states));
  ASSERT_EQ(states.size(), 6u);

  // Lexicographic: (h=0, s=0), (h=0, s=1), (h=1, s=0), ...
  const double kExpected[6][2] = {
      {0.0, 0.0}, {0.0, 0.5}, {1.0, 0.0}, {1.0, 0.5}, {2.0, 0.0}, {2.0, 0.5},
  };
  for (std::size_t i = 0; i < 6u; ++i) {
    ASSERT_EQ(states[i].hinges.size(), 1u);
    ASSERT_EQ(states[i].shears.size(), 1u);
    EXPECT_DOUBLE_EQ(states[i].hinges[0].angle_radians, kExpected[i][0]);
    EXPECT_DOUBLE_EQ(states[i].shears[0].shear_factor, kExpected[i][1]);
    EXPECT_EQ(states[i].state_id, i);
  }
}

TEST(F3DockFlexSampler, DegenerateAxisDirectionIsRejected) {
  FlexSamplingConfig cfg;
  HingeAxisSpec h;
  h.axis_direction = {0.0, 0.0, 0.0};
  h.num_steps = 1;
  cfg.hinges.push_back(h);

  std::vector<FlexState> states;
  EXPECT_FALSE(FlexSampler::enumerate(cfg, &states));
  EXPECT_TRUE(states.empty());
}

TEST(F3DockFlexSampler, EnumerateRejectsNullOutput) {
  FlexSamplingConfig cfg;
  EXPECT_FALSE(FlexSampler::enumerate(cfg, nullptr));
}

TEST(F3DockFlexSampler, ParsesHingeAxisParameter) {
  FlexSamplingConfig cfg;
  bool ok = false;
  const bool consumed = FlexSampler::parse_param(
      "flexHingeAxis", "1 2 3  0 0 1  0.0 0.25 4", &cfg, &ok);
  EXPECT_TRUE(consumed);
  EXPECT_TRUE(ok);
  ASSERT_EQ(cfg.hinges.size(), 1u);
  EXPECT_DOUBLE_EQ(cfg.hinges[0].axis_point[0], 1.0);
  EXPECT_DOUBLE_EQ(cfg.hinges[0].axis_point[1], 2.0);
  EXPECT_DOUBLE_EQ(cfg.hinges[0].axis_point[2], 3.0);
  EXPECT_DOUBLE_EQ(cfg.hinges[0].axis_direction[2], 1.0);
  EXPECT_DOUBLE_EQ(cfg.hinges[0].step_radians, 0.25);
  EXPECT_EQ(cfg.hinges[0].num_steps, 4u);
}

TEST(F3DockFlexSampler, ParsesShearPlaneParameter) {
  FlexSamplingConfig cfg;
  bool ok = false;
  const bool consumed = FlexSampler::parse_param(
      "flexShearPlane", "0 0 1  1 0 0  0.0 0.1 3", &cfg, &ok);
  EXPECT_TRUE(consumed);
  EXPECT_TRUE(ok);
  ASSERT_EQ(cfg.shears.size(), 1u);
  EXPECT_DOUBLE_EQ(cfg.shears[0].plane_normal[2], 1.0);
  EXPECT_DOUBLE_EQ(cfg.shears[0].shear_direction[0], 1.0);
  EXPECT_DOUBLE_EQ(cfg.shears[0].step, 0.1);
  EXPECT_EQ(cfg.shears[0].num_steps, 3u);
}

TEST(F3DockFlexSampler, ParseParamReturnsFalseForUnknownKey) {
  FlexSamplingConfig cfg;
  bool ok = true;
  const bool consumed = FlexSampler::parse_param("bandwidth", "1.0", &cfg, &ok);
  EXPECT_FALSE(consumed);
  // `ok` is intentionally untouched for unknown keys.
  EXPECT_TRUE(ok);
}

TEST(F3DockFlexSampler, ParseParamRejectsMalformedValue) {
  FlexSamplingConfig cfg;
  bool ok = true;
  const bool consumed =
      FlexSampler::parse_param("flexHingeAxis", "1 2 3 not enough", &cfg, &ok);
  EXPECT_TRUE(consumed);
  EXPECT_FALSE(ok);
  EXPECT_TRUE(cfg.hinges.empty());
}

TEST(F3DockFlexSampler, ParseParamRejectsNegativeNumSteps) {
  FlexSamplingConfig cfg;
  bool ok = true;
  const bool consumed = FlexSampler::parse_param(
      "flexHingeAxis", "0 0 0 0 0 1 0 0.1 -2", &cfg, &ok);
  EXPECT_TRUE(consumed);
  EXPECT_FALSE(ok);
}

} // namespace
