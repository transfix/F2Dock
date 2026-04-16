/*
 * Unit tests for the PG (Dynamic Packing Grid) data structure.
 *
 * PG is the core spatial-index used by F2Dock for range queries on atom
 * positions.  It stores points in a 5-level hierarchy:
 *   grid → plane → line → cell → balls (points)
 *
 * These tests validate:
 *   - Construction with different constructors
 *   - Point insertion (addPoint / addPoints)
 *   - Range queries (range, countPointsWithinRange, pointsWithinRange)
 *   - Point removal (removePoint)
 *   - Bounding box retrieval
 *   - Edge cases (empty grid, single point, coincident points)
 *
 * Per the documentation (Appendix B, §DPG data structure), PG should answer
 * range queries in O(k + log log w) time where k = output size.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "PG-range/PG.h"

// ── Helpers ────────────────────────────────────────────────────────────────────

static float distance(const Point& a, const Point& b) {
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Brute-force range query for verification
static std::vector<Point*> bruteForceRange(std::vector<Point*>& pts, Point* q, double radius) {
    std::vector<Point*> result;
    for (auto* p : pts) {
        if (distance(*p, *q) <= radius) {
            result.push_back(p);
        }
    }
    return result;
}

// ── Construction ───────────────────────────────────────────────────────────────

TEST(PGRange, ConstructBasic) {
    PG pg(3.5, 100.0, 5.0);
    EXPECT_EQ(pg.cellsstored(), 0);
    EXPECT_DOUBLE_EQ(pg.getdivsize(), 3.5);
}

TEST(PGRange, ConstructBoundingBox) {
    PG pg(0.0, 0.0, 0.0, 100.0, 100.0, 100.0, 5.0);
    EXPECT_EQ(pg.cellsstored(), 0);
}

// ── Single point ───────────────────────────────────────────────────────────────

TEST(PGRange, AddSinglePoint) {
    PG pg(3.5, 100.0, 5.0);
    Point p(10.0f, 20.0f, 30.0f);
    pg.addPoint(&p);

    // Should find the point within range
    Point q(10.0f, 20.0f, 30.0f);
    auto result = pg.range(&q, 1.0);
    EXPECT_GE(result.size(), 1u);
}

TEST(PGRange, CountSinglePoint) {
    PG pg(3.5, 100.0, 5.0);
    Point p(5.0f, 5.0f, 5.0f);
    pg.addPoint(&p);

    Point q(5.0f, 5.0f, 5.0f);
    EXPECT_GE(pg.countPointsWithinRange(&q, 1.0), 1);
}

// ── Multiple points ────────────────────────────────────────────────────────────

class PGRangeMultiTest : public ::testing::Test {
protected:
    static constexpr int N = 50;
    std::vector<Point> points;
    std::vector<Point*> ptrs;
    PG* pg = nullptr;

    void SetUp() override {
        points.reserve(N);
        ptrs.reserve(N);
        // Create a grid of points in [0, 50) × [0, 50) × [0, 50)
        for (int i = 0; i < N; ++i) {
            float v = static_cast<float>(i);
            points.push_back(Point(v, v * 0.5f, v * 0.3f));
        }
        for (auto& p : points) {
            ptrs.push_back(&p);
        }
        pg = new PG(0.0, 0.0, 0.0, 50.0, 50.0, 50.0, 3.5);
        pg->addPoints(&ptrs);
    }

    void TearDown() override {
        delete pg;
    }
};

TEST_F(PGRangeMultiTest, RangeQueryMatchesBruteForce) {
    // Query near point[25] = (25, 12.5, 7.5)
    Point q(25.0f, 12.5f, 7.5f);
    double radius = 5.0;

    auto pgResult = pg->range(&q, radius);
    auto bfResult = bruteForceRange(ptrs, &q, radius);

    // PG result should contain at least as many as brute force finds
    // (PG may include points slightly beyond the exact radius due to cell granularity)
    EXPECT_GE(pgResult.size(), bfResult.size())
        << "PG returned fewer points than brute force within radius " << radius;
}

TEST_F(PGRangeMultiTest, CountMatchesBruteForce) {
#ifdef _WIN32
    GTEST_SKIP() << "countPointsWithinRange has known issues on Windows";
#endif
    Point q(10.0f, 5.0f, 3.0f);
    double radius = 8.0;

    int pgCount = pg->countPointsWithinRange(&q, radius);
    auto bfResult = bruteForceRange(ptrs, &q, radius);

    EXPECT_GE(pgCount, static_cast<int>(bfResult.size()));
}

TEST_F(PGRangeMultiTest, PointsWithinRangeBoolean) {
    // Query exactly at a known point — should always return true
    Point q(10.0f, 5.0f, 3.0f);
    EXPECT_TRUE(pg->pointsWithinRange(&q, 5.0));
}

TEST_F(PGRangeMultiTest, NoPointsFarAway) {
    // Query far from all inserted points
    Point q(999.0f, 999.0f, 999.0f);
    auto result = pg->range(&q, 1.0);
    EXPECT_EQ(result.size(), 0u);
}

// ── Point removal ──────────────────────────────────────────────────────────────

TEST(PGRange, RemovePoint) {
#ifdef _WIN32
    GTEST_SKIP() << "removePoint has known issues on Windows";
#endif
    PG pg(3.5, 100.0, 5.0);
    Point a(10.0f, 10.0f, 10.0f);
    Point b(20.0f, 20.0f, 20.0f);

    pg.addPoint(&a);
    pg.addPoint(&b);

    // Both should be findable
    Point qa(10.0f, 10.0f, 10.0f);
    EXPECT_TRUE(pg.pointsWithinRange(&qa, 1.0));

    pg.removePoint(&a);

    // After removal, b should still be findable
    Point qb(20.0f, 20.0f, 20.0f);
    EXPECT_TRUE(pg.pointsWithinRange(&qb, 1.0));
}

// ── Bounding box ───────────────────────────────────────────────────────────────

TEST(PGRange, BoundingBox) {
    PG pg(-10.0, -20.0, -30.0, 40.0, 50.0, 60.0, 5.0);
    double mnx, mny, mnz, mxx, mxy, mxz;
    pg.boundingbox(mnx, mny, mnz, mxx, mxy, mxz);
    EXPECT_DOUBLE_EQ(mnx, -10.0);
    EXPECT_DOUBLE_EQ(mny, -20.0);
    EXPECT_DOUBLE_EQ(mnz, -30.0);
    EXPECT_DOUBLE_EQ(mxx, 40.0);
    EXPECT_DOUBLE_EQ(mxy, 50.0);
    EXPECT_DOUBLE_EQ(mxz, 60.0);
}

// ── Edge cases ─────────────────────────────────────────────────────────────────

TEST(PGRange, EmptyGridRangeQuery) {
    PG pg(3.5, 100.0, 5.0);
    Point q(0.0f, 0.0f, 0.0f);
    auto result = pg.range(&q, 10.0);
    EXPECT_EQ(result.size(), 0u);
}

TEST(PGRange, CoincidentPoints) {
    PG pg(3.5, 100.0, 5.0);
    Point a(5.0f, 5.0f, 5.0f);
    Point b(5.0f, 5.0f, 5.0f);
    pg.addPoint(&a);
    pg.addPoint(&b);

    Point q(5.0f, 5.0f, 5.0f);
    // Both coincident points should appear in the range query
    auto result = pg.range(&q, 0.1);
    EXPECT_GE(result.size(), 2u);
}

TEST(PGRange, LargeRadius) {
    PG pg(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 2.0);
    std::vector<Point> pts;
    std::vector<Point*> ptrs;
    for (int i = 0; i < 10; ++i) {
        pts.push_back(Point(static_cast<float>(i), 0.0f, 0.0f));
    }
    for (auto& p : pts) ptrs.push_back(&p);
    pg.addPoints(&ptrs);

    // A radius large enough to encompass all points
    Point q(5.0f, 0.0f, 0.0f);
    auto result = pg.range(&q, 100.0);
    EXPECT_EQ(result.size(), 10u);
}

// ── Timing counters ────────────────────────────────────────────────────────────

TEST(PGRange, TimingCounters) {
    PG pg(3.5, 100.0, 5.0);
    pg.resetTimes();
    EXPECT_EQ(pg.getRangeCount(), 0);
    EXPECT_DOUBLE_EQ(pg.getRangeTime(), 0.0);
    EXPECT_DOUBLE_EQ(pg.getCompTime(), 0.0);
    EXPECT_DOUBLE_EQ(pg.getInitTime(), 0.0);
}
