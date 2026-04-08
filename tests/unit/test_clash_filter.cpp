/*
 * Unit tests for clashFilter.
 *
 * The clash filter detects steric clashes between static and moving molecule
 * atoms using an octree-based approximation scheme (per Appendix B §Clash Filter).
 *
 * Two atoms clash when:
 *   distance(center_i, center_j) < eqmDistFrac * r_eqm_XY
 *
 * Default parameters from the documentation:
 *   - eqmDistFrac = 0.5
 *   - minRadius = 2.0
 *   - maxLeafSize = 10
 *   - epsilon = 0.5 (approximation factor)
 *
 * clashTolerance varies by complexType:
 *   - A (antibody): 2
 *   - E (enzyme):   9
 *   - G (generic):  10
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "fast-clash/clashFilter.h"

// Helper: create an atom array in the format expected by clashFilter:
//   5 doubles per atom: x, y, z, charge, radius
static std::vector<double> makeAtomArray(
    const std::vector<std::array<double, 5>>& atoms) {
    std::vector<double> flat;
    flat.reserve(atoms.size() * 5);
    for (auto& a : atoms) {
        flat.push_back(a[0]); // x
        flat.push_back(a[1]); // y
        flat.push_back(a[2]); // z
        flat.push_back(a[3]); // charge
        flat.push_back(a[4]); // radius
    }
    return flat;
}

// ── Construction ───────────────────────────────────────────────────────────────

TEST(ClashFilter, ConstructWithAtoms) {
    // Two well-separated atoms — no clash expected
    auto staticAtoms = makeAtomArray({
        {0.0, 0.0, 0.0, 0.0, 1.5}
    });
    auto movingAtoms = makeAtomArray({
        {100.0, 100.0, 100.0, 0.0, 1.5}
    });

    clashFilter cf(1, staticAtoms.data(), 1, movingAtoms.data(), false);
    // Just verify it doesn't crash
    SUCCEED();
}

// ── No clashes when atoms are far apart ─────────────────────────────────────

TEST(ClashFilter, NoClashesWhenFarApart) {
    auto staticAtoms = makeAtomArray({
        {0.0, 0.0, 0.0, 0.0, 1.8},
        {3.0, 0.0, 0.0, 0.0, 1.8},
        {6.0, 0.0, 0.0, 0.0, 1.8}
    });
    auto movingAtoms = makeAtomArray({
        {0.0, 50.0, 0.0, 0.0, 1.8},
        {3.0, 50.0, 0.0, 0.0, 1.8}
    });

    clashFilter cf(3, staticAtoms.data(), 2, movingAtoms.data(), false);
    cf.setEpsilon(0.5);

    int nClashes = 0, nSevere = 0;
    double interaction = 0.0;

    // Identity transform — atoms are 50 Å apart
    Matrix identity;
    identity.reset();
    bool ok = cf.computeInteractions(identity, &nClashes, &nSevere, &interaction);
    EXPECT_TRUE(ok);
    EXPECT_EQ(nClashes, 0);
    EXPECT_EQ(nSevere, 0);
}

// ── Clashes detected for overlapping atoms ──────────────────────────────────

TEST(ClashFilter, DetectsOverlapNaively) {
    // Place atoms on top of each other
    auto staticAtoms = makeAtomArray({
        {0.0, 0.0, 0.0, 0.0, 1.8}
    });
    auto movingAtoms = makeAtomArray({
        {0.0, 0.0, 0.0, 0.0, 1.8}
    });

    clashFilter cf(1, staticAtoms.data(), 1, movingAtoms.data(), false);

    int nClashes = 0, nSevere = 0;
    double interaction = 0.0;

    Matrix identity;
    identity.reset();
    bool ok = cf.computeInteractionsNaively(identity, &nClashes, &nSevere, &interaction);
    EXPECT_TRUE(ok);
    // Perfectly overlapping atoms should produce at least 1 clash
    EXPECT_GE(nClashes, 1);
}

// ── Approximation vs naive comparison ───────────────────────────────────────

TEST(ClashFilter, ApproximateMatchesNaive) {
    // Create a small system with known overlaps
    auto staticAtoms = makeAtomArray({
        {0.0, 0.0, 0.0, 0.0, 1.8},
        {5.0, 0.0, 0.0, 0.0, 1.8},
        {10.0, 0.0, 0.0, 0.0, 1.8}
    });
    auto movingAtoms = makeAtomArray({
        {0.5, 0.0, 0.0, 0.0, 1.8},  // close to static[0]
        {50.0, 0.0, 0.0, 0.0, 1.8}  // far away
    });

    clashFilter cf(3, staticAtoms.data(), 2, movingAtoms.data(), false);

    Matrix identity;
    identity.reset();

    int nClashNaive = 0, nSevereNaive = 0;
    double interNaive = 0.0;
    cf.computeInteractionsNaively(identity, &nClashNaive, &nSevereNaive, &interNaive);

    int nClashApprox = 0, nSevereApprox = 0;
    double interApprox = 0.0;
    cf.computeInteractions(identity, &nClashApprox, &nSevereApprox, &interApprox);

    // The approximate algorithm may over-count but should not miss clashes
    EXPECT_GE(nClashApprox, nClashNaive)
        << "Approximate method found fewer clashes than naive";
}

// ── Transformation matrix affects clash detection ───────────────────────────

TEST(ClashFilter, TransformMovesAtomsApart) {
    // Atoms overlap at origin
    auto staticAtoms = makeAtomArray({
        {0.0, 0.0, 0.0, 0.0, 1.8}
    });
    auto movingAtoms = makeAtomArray({
        {0.0, 0.0, 0.0, 0.0, 1.8}
    });

    clashFilter cf(1, staticAtoms.data(), 1, movingAtoms.data(), false);

    // Translate the moving atom far away
    Matrix trans = Matrix::translation(100.0f, 0.0f, 0.0f);

    int nClashes = 0, nSevere = 0;
    double interaction = 0.0;
    cf.computeInteractionsNaively(trans, &nClashes, &nSevere, &interaction);

    EXPECT_EQ(nClashes, 0);
}

// ── Parameter setters ──────────────────────────────────────────────────────────

TEST(ClashFilter, SetEpsilon) {
    auto staticAtoms = makeAtomArray({{0, 0, 0, 0, 1.8}});
    auto movingAtoms = makeAtomArray({{0, 0, 0, 0, 1.8}});
    clashFilter cf(1, staticAtoms.data(), 1, movingAtoms.data(), false);

    EXPECT_TRUE(cf.setEpsilon(0.3));
    EXPECT_TRUE(cf.setMinRadius(1.5));
    EXPECT_TRUE(cf.setMaxLeafSize(20));
}

// ── Documentation: clashWeight per complex type ─────────────────────────────
// These are integration-level assertions about expected parameter values.
// The actual PARAMS_IN defaults are verified in test_file_formats.cpp.

TEST(ClashFilter, DocumentedDefaults) {
    // Per documentation Table params.pdf:
    // eqmDistFrac = 0.5 (all types)
    // clashTolerance: A=2, E=9, G=10
    // clashWeight: A=-30, E/G=-0.5

    // We can't directly read PARAMS_IN here, but we verify the
    // clash filter accepts the documented epsilon (eqmDistFrac)
    auto staticAtoms = makeAtomArray({{0, 0, 0, 0, 1.8}});
    auto movingAtoms = makeAtomArray({{0, 0, 0, 0, 1.8}});
    clashFilter cf(1, staticAtoms.data(), 1, movingAtoms.data(), false);

    // epsilon = 0.5 per documentation
    EXPECT_TRUE(cf.setEpsilon(0.5));
}
