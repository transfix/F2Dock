// SPDX-License-Identifier: LGPL-2.1-or-later
//
// Unit tests for f3dock::DockMode parsing and stringification.

#include <gtest/gtest.h>

#include "f3dock/DockMode.h"

using f3dock::DockMode;
using f3dock::parse_dock_mode;
using f3dock::to_string;

TEST(DockMode, ParsesF2Synonyms) {
  DockMode m = DockMode::kF3Dock; // start non-default to detect writes
  ASSERT_TRUE(parse_dock_mode("f2", &m));
  EXPECT_EQ(m, DockMode::kF2Dock);

  ASSERT_TRUE(parse_dock_mode("F2", &m));
  EXPECT_EQ(m, DockMode::kF2Dock);

  ASSERT_TRUE(parse_dock_mode("f2dock", &m));
  EXPECT_EQ(m, DockMode::kF2Dock);

  ASSERT_TRUE(parse_dock_mode("rigid", &m));
  EXPECT_EQ(m, DockMode::kF2Dock);
}

TEST(DockMode, ParsesF3Synonyms) {
  DockMode m = DockMode::kF2Dock;
  ASSERT_TRUE(parse_dock_mode("f3", &m));
  EXPECT_EQ(m, DockMode::kF3Dock);

  ASSERT_TRUE(parse_dock_mode("F3Dock", &m));
  EXPECT_EQ(m, DockMode::kF3Dock);

  ASSERT_TRUE(parse_dock_mode("flex", &m));
  EXPECT_EQ(m, DockMode::kF3Dock);

  ASSERT_TRUE(parse_dock_mode("FLEXIBLE", &m));
  EXPECT_EQ(m, DockMode::kF3Dock);
}

TEST(DockMode, RejectsUnknownAndNullInputs) {
  DockMode m = DockMode::kF2Dock;
  EXPECT_FALSE(parse_dock_mode("", &m));
  EXPECT_FALSE(parse_dock_mode("rigid-ish", &m));
  EXPECT_FALSE(parse_dock_mode("f4", &m));
  EXPECT_FALSE(parse_dock_mode(nullptr, &m));
  EXPECT_FALSE(parse_dock_mode("f2", nullptr));
}

TEST(DockMode, ToStringIsHumanReadable) {
  EXPECT_STREQ(to_string(DockMode::kF2Dock), "F2Dock (rigid)");
  EXPECT_STREQ(to_string(DockMode::kF3Dock), "F3Dock (flexible)");
}

TEST(DockMode, IntCastIsStable) {
  // PARAMS_IN stores the mode as int (C ABI). Ensure the values are
  // the documented constants so existing input files / serialized
  // params keep working as we evolve the enum.
  EXPECT_EQ(static_cast<int>(DockMode::kF2Dock), 0);
  EXPECT_EQ(static_cast<int>(DockMode::kF3Dock), 1);
}
