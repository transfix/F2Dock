/*
 * Unit tests for CCVOpenGLMath::Matrix, CCVOpenGLMath::Vector, CCVOpenGLMath::Tuple
 *
 * Validates the linear-algebra primitives used throughout F2Dock for
 * rotation matrices, transformation matrices, and 3-D vector operations.
 */

#include <gtest/gtest.h>
#include <cmath>
#include "math/Matrix.h"
#include "math/Vector.h"
#include "math/Tuple.h"

using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;
using CCVOpenGLMath::Tuple;

static constexpr float kEps = 1e-5f;

// ── Tuple ──────────────────────────────────────────────────────────────────────

TEST(Tuple, DefaultConstructor) {
    Tuple t;
    // Default-constructed Tuple should be zero (implementation dependent, but
    // we verify the subscript operator works)
    (void)t[0];
    (void)t[3];
}

TEST(Tuple, SetAndAccess) {
    Tuple t;
    t.set(1.0f, 2.0f, 3.0f, 4.0f);
    EXPECT_FLOAT_EQ(t[0], 1.0f);
    EXPECT_FLOAT_EQ(t[1], 2.0f);
    EXPECT_FLOAT_EQ(t[2], 3.0f);
    EXPECT_FLOAT_EQ(t[3], 4.0f);
}

TEST(Tuple, CopyConstructorAndAssignment) {
    Tuple a;
    a.set(5.0f, 6.0f, 7.0f, 8.0f);

    Tuple b(a);
    EXPECT_FLOAT_EQ(b[0], 5.0f);
    EXPECT_FLOAT_EQ(b[3], 8.0f);

    Tuple c;
    c = a;
    EXPECT_FLOAT_EQ(c[2], 7.0f);
}

// ── Vector ─────────────────────────────────────────────────────────────────────

TEST(Vector, Construction) {
    Vector v(1.0f, 0.0f, 0.0f, 0.0f);
    EXPECT_FLOAT_EQ(v[0], 1.0f);
    EXPECT_FLOAT_EQ(v[1], 0.0f);
}

TEST(Vector, Norm) {
    Vector v(3.0f, 4.0f, 0.0f, 0.0f);
    EXPECT_NEAR(v.norm(), 5.0f, kEps);
}

TEST(Vector, Normalize) {
    Vector v(0.0f, 0.0f, 5.0f, 0.0f);
    v.normalize();
    EXPECT_NEAR(v.norm(), 1.0f, kEps);
    EXPECT_NEAR(v[2], 1.0f, kEps);
}

TEST(Vector, DotProduct) {
    Vector a(1.0f, 0.0f, 0.0f, 0.0f);
    Vector b(0.0f, 1.0f, 0.0f, 0.0f);
    EXPECT_NEAR(a.dot(b), 0.0f, kEps);

    Vector c(1.0f, 2.0f, 3.0f, 0.0f);
    Vector d(4.0f, 5.0f, 6.0f, 0.0f);
    // dot = 4 + 10 + 18 = 32
    EXPECT_NEAR(c.dot(d), 32.0f, kEps);
}

TEST(Vector, CrossProduct) {
    // x × y = z
    Vector x(1.0f, 0.0f, 0.0f, 0.0f);
    Vector y(0.0f, 1.0f, 0.0f, 0.0f);
    Vector z = x.cross(y);
    EXPECT_NEAR(z[0], 0.0f, kEps);
    EXPECT_NEAR(z[1], 0.0f, kEps);
    EXPECT_NEAR(z[2], 1.0f, kEps);
}

TEST(Vector, ArithmeticOperators) {
    Vector a(1.0f, 2.0f, 3.0f, 0.0f);
    Vector b(4.0f, 5.0f, 6.0f, 0.0f);

    Vector sum = a + b;
    EXPECT_FLOAT_EQ(sum[0], 5.0f);
    EXPECT_FLOAT_EQ(sum[1], 7.0f);
    EXPECT_FLOAT_EQ(sum[2], 9.0f);

    Vector diff = b - a;
    EXPECT_FLOAT_EQ(diff[0], 3.0f);

    Vector scaled = a * 2.0f;
    EXPECT_FLOAT_EQ(scaled[0], 2.0f);
    EXPECT_FLOAT_EQ(scaled[2], 6.0f);
}

TEST(Vector, Negation) {
    Vector v(1.0f, -2.0f, 3.0f, 0.0f);
    Vector neg = -v;
    EXPECT_FLOAT_EQ(neg[0], -1.0f);
    EXPECT_FLOAT_EQ(neg[1], 2.0f);
    EXPECT_FLOAT_EQ(neg[2], -3.0f);
}

TEST(Vector, BadVector) {
    Vector bad = Vector::badVector();
    EXPECT_TRUE(bad.isBad());
}

// ── Matrix ─────────────────────────────────────────────────────────────────────

TEST(Matrix, IdentityDefault) {
    Matrix m;
    m.reset();
    // reset() sets identity
    EXPECT_FLOAT_EQ(m.get(0, 0), 1.0f);
    EXPECT_FLOAT_EQ(m.get(1, 1), 1.0f);
    EXPECT_FLOAT_EQ(m.get(2, 2), 1.0f);
    EXPECT_FLOAT_EQ(m.get(3, 3), 1.0f);
    EXPECT_FLOAT_EQ(m.get(0, 1), 0.0f);
    EXPECT_FLOAT_EQ(m.get(1, 0), 0.0f);
}

TEST(Matrix, IdentityDeterminant) {
    Matrix m;
    m.reset();
    EXPECT_NEAR(m.determinant(), 1.0f, kEps);
}

TEST(Matrix, TranslationMatrix) {
    Matrix t = Matrix::translation(3.0f, 4.0f, 5.0f);
    // Translation matrix should be:
    // 1 0 0 3
    // 0 1 0 4
    // 0 0 1 5
    // 0 0 0 1
    EXPECT_FLOAT_EQ(t.get(0, 3), 3.0f);
    EXPECT_FLOAT_EQ(t.get(1, 3), 4.0f);
    EXPECT_FLOAT_EQ(t.get(2, 3), 5.0f);
    EXPECT_FLOAT_EQ(t.get(0, 0), 1.0f);
    EXPECT_FLOAT_EQ(t.get(3, 3), 1.0f);
}

TEST(Matrix, TranslationAppliedToVector) {
    Matrix t = Matrix::translation(10.0f, 20.0f, 30.0f);
    Vector v(1.0f, 2.0f, 3.0f, 1.0f); // w=1 for point
    Vector result = t * v;
    EXPECT_NEAR(result[0], 11.0f, kEps);
    EXPECT_NEAR(result[1], 22.0f, kEps);
    EXPECT_NEAR(result[2], 33.0f, kEps);
}

TEST(Matrix, RotationX90) {
    // F2Dock convention: rotationX uses transposed rotation matrix.
    // rotationX(π/2) * (0,1,0) = (0, cos(π/2), -sin(π/2)) = (0, 0, -1)
    float angle = static_cast<float>(M_PI / 2.0);
    Matrix rx = Matrix::rotationX(angle);

    Vector v(0.0f, 1.0f, 0.0f, 0.0f);
    Vector result = rx * v;
    EXPECT_NEAR(result[0], 0.0f, kEps);
    EXPECT_NEAR(result[1], 0.0f, kEps);
    EXPECT_NEAR(result[2], -1.0f, kEps);
}

TEST(Matrix, RotationY90) {
    // F2Dock convention: rotationY uses transposed rotation matrix.
    // rotationY(π/2) * (0,0,1) = (-sin(π/2), 0, cos(π/2)) = (-1, 0, 0)
    float angle = static_cast<float>(M_PI / 2.0);
    Matrix ry = Matrix::rotationY(angle);

    Vector v(0.0f, 0.0f, 1.0f, 0.0f);
    Vector result = ry * v;
    EXPECT_NEAR(result[0], -1.0f, kEps);
    EXPECT_NEAR(result[1], 0.0f, kEps);
    EXPECT_NEAR(result[2], 0.0f, kEps);
}

TEST(Matrix, RotationZ90) {
    // F2Dock convention: rotationZ uses transposed rotation matrix.
    // rotationZ(π/2) * (1,0,0) = (cos(π/2), -sin(π/2), 0) = (0, -1, 0)
    float angle = static_cast<float>(M_PI / 2.0);
    Matrix rz = Matrix::rotationZ(angle);

    Vector v(1.0f, 0.0f, 0.0f, 0.0f);
    Vector result = rz * v;
    EXPECT_NEAR(result[0], 0.0f, kEps);
    EXPECT_NEAR(result[1], -1.0f, kEps);
    EXPECT_NEAR(result[2], 0.0f, kEps);
}

TEST(Matrix, RotationPreservesLength) {
    // A rotation should not change the magnitude of a vector.
    float angle = 1.234f;
    Matrix r = Matrix::rotationX(angle) * Matrix::rotationY(0.567f);

    Vector v(1.0f, 2.0f, 3.0f, 0.0f);
    Vector rv = r * v;

    float origLen = v.norm();
    float rotLen  = rv.norm();
    EXPECT_NEAR(rotLen, origLen, kEps);
}

TEST(Matrix, InverseOfIdentity) {
    Matrix m;
    m.reset();
    Matrix inv = m.inverse();
    // Should still be identity
    EXPECT_NEAR(inv.get(0, 0), 1.0f, kEps);
    EXPECT_NEAR(inv.get(1, 2), 0.0f, kEps);
}

TEST(Matrix, InverseOfTranslation) {
    Matrix t = Matrix::translation(5.0f, -3.0f, 7.0f);
    Matrix inv = t.inverse();
    // T^-1 * T should be identity
    Matrix prod = inv * t;
    EXPECT_NEAR(prod.get(0, 0), 1.0f, kEps);
    EXPECT_NEAR(prod.get(0, 3), 0.0f, kEps);
    EXPECT_NEAR(prod.get(1, 3), 0.0f, kEps);
    EXPECT_NEAR(prod.get(2, 3), 0.0f, kEps);
}

TEST(Matrix, InverseOfRotation) {
    float angle = 0.7f;
    Matrix r = Matrix::rotationZ(angle);
    Matrix inv = r.inverse();
    Matrix prod = inv * r;

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            float expected = (i == j) ? 1.0f : 0.0f;
            EXPECT_NEAR(prod.get(i, j), expected, kEps)
                << "at (" << i << "," << j << ")";
        }
}

TEST(Matrix, Transpose) {
    Matrix m(
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    );
    Matrix mt = m.transpose();
    EXPECT_FLOAT_EQ(mt.get(0, 1), m.get(1, 0));
    EXPECT_FLOAT_EQ(mt.get(2, 3), m.get(3, 2));
}

TEST(Matrix, ScaleMatrix) {
    Matrix s = Matrix::scale(2.0f, 3.0f, 4.0f);
    Vector v(1.0f, 1.0f, 1.0f, 1.0f);
    Vector result = s * v;
    EXPECT_NEAR(result[0], 2.0f, kEps);
    EXPECT_NEAR(result[1], 3.0f, kEps);
    EXPECT_NEAR(result[2], 4.0f, kEps);
}

TEST(Matrix, Determinant3x3Subblock) {
    // Known determinant: rotation matrices have det = 1
    float angle = 1.5f;
    Matrix r = Matrix::rotationX(angle);
    EXPECT_NEAR(r.determinant(), 1.0f, kEps);
}

TEST(Matrix, MultiplicationAssociativity) {
    Matrix a = Matrix::rotationX(0.5f);
    Matrix b = Matrix::rotationY(0.7f);
    Matrix c = Matrix::translation(1.0f, 2.0f, 3.0f);

    Matrix ab_c = (a * b) * c;
    Matrix a_bc = a * (b * c);

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            EXPECT_NEAR(ab_c.get(i, j), a_bc.get(i, j), kEps)
                << "at (" << i << "," << j << ")";
}

TEST(Matrix, CopyAndAssignment) {
    Matrix orig = Matrix::translation(1.0f, 2.0f, 3.0f);
    Matrix copy(orig);
    EXPECT_FLOAT_EQ(copy.get(0, 3), 1.0f);

    Matrix assigned;
    assigned = orig;
    EXPECT_FLOAT_EQ(assigned.get(1, 3), 2.0f);
}

// ── Docking-specific: transformation matrix applied to atoms ───────────────────
// Per the documentation, F2Dock output contains a 3×4 transformation matrix
// (mat1..mat12) that transforms the moving molecule. This test verifies
// that a composed rotation+translation can round-trip through inverse.

TEST(Matrix, DockingTransformRoundTrip) {
    // Simulate a docking result: rotation then translation
    Matrix rot = Matrix::rotationX(0.3f) * Matrix::rotationY(0.6f) * Matrix::rotationZ(0.9f);
    Matrix trans = Matrix::translation(12.5f, -3.7f, 8.2f);
    Matrix dockTransform = trans * rot;

    // Apply to an atom position
    Vector atom(10.0f, 20.0f, 30.0f, 1.0f);
    Vector moved = dockTransform * atom;

    // Inverse should recover original
    Matrix inv = dockTransform.inverse();
    Vector recovered = inv * moved;

    EXPECT_NEAR(recovered[0], 10.0f, 1e-3f);
    EXPECT_NEAR(recovered[1], 20.0f, 1e-3f);
    EXPECT_NEAR(recovered[2], 30.0f, 1e-3f);
}
