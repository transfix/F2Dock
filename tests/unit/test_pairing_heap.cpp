/*
 * Unit tests for PairingHeap.
 *
 * PairingHeap is used internally by TopValues to maintain the top-K docking
 * solutions efficiently.  It supports Insert, Find_Min, Delete_Min, and
 * Decrease_Key operations.
 */

#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <climits>
#include "fast-PQ/PairingHeap.h"

// ── Basic operations ───────────────────────────────────────────────────────────

TEST(PairingHeap, EmptyOnConstruction) {
    PairingHeap ph;
    EXPECT_TRUE(ph.isEmpty());
}

TEST(PairingHeap, InsertSingle) {
    PairingHeap ph;
    ph.Insert(42, 3.14);
    EXPECT_FALSE(ph.isEmpty());
}

TEST(PairingHeap, FindMinSingle) {
    PairingHeap ph;
    ph.Insert(1, 10.0);
    int mx;
    VAL_TYPE mk;
    ph.Find_Min(mx, mk);
    EXPECT_EQ(mx, 1);
    EXPECT_DOUBLE_EQ(mk, 10.0);
}

TEST(PairingHeap, DeleteMinSingle) {
    PairingHeap ph;
    ph.Insert(1, 5.0);
    int mx;
    VAL_TYPE mk;
    ph.Delete_Min(mx, mk);
    EXPECT_EQ(mx, 1);
    EXPECT_DOUBLE_EQ(mk, 5.0);
    EXPECT_TRUE(ph.isEmpty());
}

// ── Ordering ───────────────────────────────────────────────────────────────────

TEST(PairingHeap, ExtractsInOrder) {
    PairingHeap ph;
    // Insert in random order
    ph.Insert(3, 30.0);
    ph.Insert(1, 10.0);
    ph.Insert(5, 50.0);
    ph.Insert(2, 20.0);
    ph.Insert(4, 40.0);

    // Extract min repeatedly — should come out in ascending key order
    std::vector<double> extracted;
    while (!ph.isEmpty()) {
        int mx;
        VAL_TYPE mk;
        ph.Delete_Min(mx, mk);
        extracted.push_back(mk);
    }

    ASSERT_EQ(extracted.size(), 5u);
    EXPECT_DOUBLE_EQ(extracted[0], 10.0);
    EXPECT_DOUBLE_EQ(extracted[1], 20.0);
    EXPECT_DOUBLE_EQ(extracted[2], 30.0);
    EXPECT_DOUBLE_EQ(extracted[3], 40.0);
    EXPECT_DOUBLE_EQ(extracted[4], 50.0);
}

TEST(PairingHeap, DuplicateKeys) {
    PairingHeap ph;
    ph.Insert(1, 5.0);
    ph.Insert(2, 5.0);
    ph.Insert(3, 5.0);

    int mx;
    VAL_TYPE mk;
    ph.Delete_Min(mx, mk);
    EXPECT_DOUBLE_EQ(mk, 5.0);

    ph.Delete_Min(mx, mk);
    EXPECT_DOUBLE_EQ(mk, 5.0);

    ph.Delete_Min(mx, mk);
    EXPECT_DOUBLE_EQ(mk, 5.0);

    EXPECT_TRUE(ph.isEmpty());
}

// ── Decrease key ───────────────────────────────────────────────────────────────

TEST(PairingHeap, DecreaseKey) {
    PairingHeap ph;
    int node1 = ph.Insert(1, 100.0);
    int node2 = ph.Insert(2, 200.0);

    // Decrease node2's key to be smaller than node1
    ph.Decrease_Key(node2, 50.0);

    int mx;
    VAL_TYPE mk;
    ph.Find_Min(mx, mk);
    EXPECT_EQ(mx, 2);
    EXPECT_DOUBLE_EQ(mk, 50.0);
}

// ── Multi-pass and auxiliary tree modes ─────────────────────────────────────────

TEST(PairingHeap, MultiPassMode) {
    PairingHeap ph(true);
    ph.Insert(1, 30.0);
    ph.Insert(2, 10.0);
    ph.Insert(3, 20.0);

    int mx;
    VAL_TYPE mk;
    ph.Delete_Min(mx, mk);
    EXPECT_DOUBLE_EQ(mk, 10.0);
}

TEST(PairingHeap, AuxTreeMode) {
    PairingHeap ph(false, true);
    ph.Insert(1, 30.0);
    ph.Insert(2, 10.0);
    ph.Insert(3, 20.0);

    int mx;
    VAL_TYPE mk;
    ph.Delete_Min(mx, mk);
    EXPECT_DOUBLE_EQ(mk, 10.0);
}

// ── Stress test ────────────────────────────────────────────────────────────────

TEST(PairingHeap, LargeInsertExtract) {
    PairingHeap ph;
    const int N = 1000;
    for (int i = N; i > 0; --i) {
        ph.Insert(i, static_cast<double>(i));
    }

    // Should extract 1, 2, 3, ..., N
    double prev = -1.0;
    for (int i = 0; i < N; ++i) {
        int mx;
        VAL_TYPE mk;
        ph.Delete_Min(mx, mk);
        EXPECT_GT(mk, prev);
        prev = mk;
    }
    EXPECT_TRUE(ph.isEmpty());
}
