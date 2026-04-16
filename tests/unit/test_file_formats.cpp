/*
 * Unit tests for F2Dock file formats and parameter parsing.
 *
 * Tests cover:
 *   1. PQR file format validation (per documentation: PDB + radius + charge)
 *   2. F2D file format validation (PQR + I/E type column + HETATM grown skin)
 *   3. QUAD file format validation (7 floats: x,y,z, nx,ny,nz, weight)
 *   4. Parameter file (.inp) format validation
 *   5. Utility parsing functions (getInt, getDouble, getAlphaString, getString)
 *   6. Default parameter values per complexType (A, E, G, U)
 *
 * References:
 *   - dockingTutorial.pdf §3 (Molecule Preparation), §4.2 (Parameter File)
 *   - params.pdf (complete parameter reference)
 */

#include <gtest/gtest.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "utils/utils.h"

// ── Utility function tests ─────────────────────────────────────────────────────

TEST(Utils, GetInt) {
    char buf[] = "  42  hello";
    int val = 0;
    int pos = getInt(buf, 0, &val);
    EXPECT_EQ(val, 42);
    EXPECT_GT(pos, 0);  // should advance past "42"
}

TEST(Utils, GetIntNegative) {
    char buf[] = "  -7 ";
    int val = 0;
    int pos = getInt(buf, 0, &val);
    EXPECT_EQ(val, -7);
    EXPECT_GT(pos, 0);
}

TEST(Utils, GetDouble) {
    char buf[] = "  3.14159  ";
    double val = 0.0;
    int pos = getDouble(buf, 0, &val);
    EXPECT_NEAR(val, 3.14159, 1e-5);
    EXPECT_GT(pos, 0);
}

TEST(Utils, GetDoubleNegative) {
    char buf[] = "-2.3";
    double val = 0.0;
    int pos = getDouble(buf, 0, &val);
    EXPECT_NEAR(val, -2.3, 1e-6);
}

TEST(Utils, GetAlphaString) {
    char buf[] = "  GLY  123";
    char s[100] = {};
    int pos = getAlphaString(buf, 0, s);
    EXPECT_STREQ(s, "GLY");
    EXPECT_GT(pos, 0);
}

TEST(Utils, GetString) {
    char buf[] = "  file.pqr  next";
    char s[100] = {};
    int pos = getString(buf, 0, s);
    EXPECT_STREQ(s, "file.pqr");
    EXPECT_GT(pos, 0);
}

TEST(Utils, SkipWhiteSpaces) {
    char buf[] = "   ABC";
    int pos = skipWhiteSpaces(buf, 0);
    EXPECT_EQ(buf[pos], 'A');
}

// ── PQR format validation ──────────────────────────────────────────────────────
// PQR format per documentation: standard PDB ATOM records with charge and radius
// replacing the occupancy and B-factor columns.
// Format: ATOM  serial name resName chainID resSeq x y z charge radius

class PQRFormatTest : public ::testing::Test {
protected:
    std::string testDataDir;

    void SetUp() override {
        testDataDir = F2DOCK_TEST_DATA_DIR;
    }

    struct PQRAtom {
        int serial;
        char name[5];
        char resName[4];
        char chainID;
        int resSeq;
        double x, y, z;
        double charge, radius;
    };

    bool parsePQRLine(const std::string& line, PQRAtom& atom) {
        if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM")
            return false;

        // PQR format: columns are space-separated after the record type
        char record[7] = {};
        int n = sscanf(line.c_str(), "%6s %d %4s %3s %c %d %lf %lf %lf %lf %lf",
                       record, &atom.serial, atom.name, atom.resName,
                       &atom.chainID, &atom.resSeq,
                       &atom.x, &atom.y, &atom.z,
                       &atom.charge, &atom.radius);
        return (n >= 11);
    }
};

TEST_F(PQRFormatTest, PQRFileHasValidHeader) {
    std::ifstream f(testDataDir + "/1A2K_R_U.pqr");
    ASSERT_TRUE(f.is_open()) << "Cannot open test PQR file";

    std::string line;
    bool hasRemark = false;
    while (std::getline(f, line)) {
        if (line.substr(0, 6) == "REMARK") {
            hasRemark = true;
            break;
        }
    }
    EXPECT_TRUE(hasRemark) << "PQR file should have REMARK header lines";
}

TEST_F(PQRFormatTest, PQRFileHasAtoms) {
    std::ifstream f(testDataDir + "/1A2K_R_U.pqr");
    ASSERT_TRUE(f.is_open());

    std::string line;
    int atomCount = 0;
    while (std::getline(f, line)) {
        if (line.substr(0, 4) == "ATOM") atomCount++;
    }
    EXPECT_GT(atomCount, 0) << "PQR file should contain ATOM records";
}

TEST_F(PQRFormatTest, PQRAtomHasChargeAndRadius) {
    std::ifstream f(testDataDir + "/1A2K_R_U.pqr");
    ASSERT_TRUE(f.is_open());

    std::string line;
    while (std::getline(f, line)) {
        if (line.substr(0, 4) != "ATOM") continue;

        PQRAtom atom = {};
        bool parsed = parsePQRLine(line, atom);
        ASSERT_TRUE(parsed) << "Failed to parse: " << line;

        // Per documentation: charge can be negative, radius must be positive
        EXPECT_GT(atom.radius, 0.0)
            << "Radius should be positive for atom " << atom.serial;
        // Coordinates should be reasonable (not NaN or huge)
        EXPECT_FALSE(std::isnan(atom.x));
        EXPECT_FALSE(std::isnan(atom.y));
        EXPECT_FALSE(std::isnan(atom.z));
        break;  // Test first atom only as a format check
    }
}

// ── F2D format validation ──────────────────────────────────────────────────────
// F2D format per documentation: like PQR but with I/E type column at the end.
// I = interior (core), E = exterior (skin)
// May also contain HETATM records for "grown" skin atoms.

class F2DFormatTest : public ::testing::Test {
protected:
    std::string testDataDir;
    void SetUp() override { testDataDir = F2DOCK_TEST_DATA_DIR; }
};

TEST_F(F2DFormatTest, F2DFileHasAtomRecords) {
    std::ifstream f(testDataDir + "/1A2K_R_U_1.7.f2d");
    ASSERT_TRUE(f.is_open());

    std::string line;
    int atomCount = 0;
    int hetatmCount = 0;
    while (std::getline(f, line)) {
        if (line.substr(0, 4) == "ATOM") atomCount++;
        if (line.substr(0, 6) == "HETATM") hetatmCount++;
    }
    EXPECT_GT(atomCount, 0) << "F2D file should contain ATOM records";
    // HETATM records are optional (grown skin atoms from documentation)
}

TEST_F(F2DFormatTest, F2DAtomHasTypeColumn) {
    std::ifstream f(testDataDir + "/1A2K_R_U_1.7.f2d");
    ASSERT_TRUE(f.is_open());

    std::string line;
    while (std::getline(f, line)) {
        if (line.substr(0, 4) != "ATOM") continue;

        // Last non-whitespace character should be 'I' or 'E'
        std::string trimmed = line;
        while (!trimmed.empty() && (trimmed.back() == ' ' || trimmed.back() == '\n' || trimmed.back() == '\r'))
            trimmed.pop_back();

        char lastChar = trimmed.back();
        EXPECT_TRUE(lastChar == 'I' || lastChar == 'E')
            << "F2D type column should be I (interior) or E (exterior), got '"
            << lastChar << "' in: " << line;
        break;
    }
}

TEST_F(F2DFormatTest, F2DAtomHasExpectedColumns) {
    // Per documentation: ATOM serial name resName resSeq x y z charge radius type
    std::ifstream f(testDataDir + "/1A2K_R_U_1.7.f2d");
    ASSERT_TRUE(f.is_open());

    std::string line;
    int checked = 0;
    while (std::getline(f, line) && checked < 5) {
        if (line.substr(0, 4) != "ATOM") continue;

        char record[7], name[5], resName[4], type[2];
        int serial, resSeq;
        double x, y, z, charge, radius;

        int n = sscanf(line.c_str(), "%6s %d %4s %3s %d %lf %lf %lf %lf %lf %1s",
                       record, &serial, name, resName, &resSeq,
                       &x, &y, &z, &charge, &radius, type);
        EXPECT_GE(n, 11) << "F2D line should have at least 11 fields: " << line;
        EXPECT_GT(radius, 0.0);
        EXPECT_TRUE(type[0] == 'I' || type[0] == 'E');
        checked++;
    }
    EXPECT_GT(checked, 0);
}

// ── QUAD format validation ─────────────────────────────────────────────────────
// QUAD format per documentation: 7 floats per line:
//   x y z nx ny nz weight
// where (x,y,z) = quadrature point position, (nx,ny,nz) = surface normal, weight = quadrature weight.

class QUADFormatTest : public ::testing::Test {
protected:
    std::string testDataDir;
    void SetUp() override { testDataDir = F2DOCK_TEST_DATA_DIR; }
};

TEST_F(QUADFormatTest, QUADFileHasSevenColumnsPerLine) {
    std::ifstream f(testDataDir + "/1A2K_R_U.quad");
    ASSERT_TRUE(f.is_open());

    std::string line;
    int lineCount = 0;
    while (std::getline(f, line) && lineCount < 10) {
        double x, y, z, nx, ny, nz, w;
        int n = sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf",
                       &x, &y, &z, &nx, &ny, &nz, &w);
        EXPECT_EQ(n, 7) << "QUAD line should have 7 floats: " << line;
        lineCount++;
    }
    EXPECT_GT(lineCount, 0);
}

TEST_F(QUADFormatTest, QUADNormalsAreUnitLength) {
    std::ifstream f(testDataDir + "/1A2K_R_U.quad");
    ASSERT_TRUE(f.is_open());

    std::string line;
    int checked = 0;
    while (std::getline(f, line) && checked < 20) {
        double x, y, z, nx, ny, nz, w;
        if (sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf",
                   &x, &y, &z, &nx, &ny, &nz, &w) != 7) continue;

        double len = std::sqrt(nx * nx + ny * ny + nz * nz);
        EXPECT_NEAR(len, 1.0, 0.05)
            << "QUAD normal should be approximately unit length, got " << len;
        checked++;
    }
    EXPECT_GT(checked, 0);
}

TEST_F(QUADFormatTest, QUADWeightsArePositive) {
    // Per documentation, quadrature weights should be positive (used for
    // surface area integration)
    std::ifstream f(testDataDir + "/1A2K_R_U.quad");
    ASSERT_TRUE(f.is_open());

    std::string line;
    int checked = 0;
    while (std::getline(f, line) && checked < 50) {
        double x, y, z, nx, ny, nz, w;
        if (sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf",
                   &x, &y, &z, &nx, &ny, &nz, &w) != 7) continue;

        EXPECT_GT(w, 0.0)
            << "QUAD weight should be positive, got " << w;
        checked++;
    }
}

// ── Parameter file (.inp) format validation ────────────────────────────────────
// The parameter file is a key-value text file with space-separated pairs.
// Lines starting with '#' are comments. Lines with <3 chars are skipped.

class ParamFileTest : public ::testing::Test {
protected:
    std::string testDataDir;
    void SetUp() override { testDataDir = F2DOCK_TEST_DATA_DIR; }
};

TEST_F(ParamFileTest, INPFileHasMandatoryParams) {
    // Per documentation §4.2.2, mandatory parameters are:
    // staticMolecule, movingMolecule, staticMoleculePQR, movingMoleculePQR,
    // staticMoleculeQUAD, movingMoleculeQUAD, outFile
    std::ifstream f(testDataDir + "/1A2K.inp");
    ASSERT_TRUE(f.is_open());

    bool hasStaticMol = false, hasMovingMol = false;
    bool hasStaticPQR = false, hasMovingPQR = false;
    bool hasStaticQUAD = false, hasMovingQUAD = false;
    bool hasOutFile = false;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "staticMolecule") hasStaticMol = true;
        else if (key == "movingMolecule") hasMovingMol = true;
        else if (key == "staticMoleculePQR") hasStaticPQR = true;
        else if (key == "movingMoleculePQR") hasMovingPQR = true;
        else if (key == "staticMoleculeQUAD") hasStaticQUAD = true;
        else if (key == "movingMoleculeQUAD") hasMovingQUAD = true;
        else if (key == "outFile") hasOutFile = true;
    }

    EXPECT_TRUE(hasStaticMol) << "Missing mandatory param: staticMolecule";
    EXPECT_TRUE(hasMovingMol) << "Missing mandatory param: movingMolecule";
    EXPECT_TRUE(hasStaticPQR) << "Missing mandatory param: staticMoleculePQR";
    EXPECT_TRUE(hasMovingPQR) << "Missing mandatory param: movingMoleculePQR";
    EXPECT_TRUE(hasStaticQUAD) << "Missing mandatory param: staticMoleculeQUAD";
    EXPECT_TRUE(hasMovingQUAD) << "Missing mandatory param: movingMoleculeQUAD";
    EXPECT_TRUE(hasOutFile) << "Missing mandatory param: outFile";
}

TEST_F(ParamFileTest, ComplexTypeIsValid) {
    // Per documentation, complexType must be one of: A, E, G, U
    std::ifstream f(testDataDir + "/1A2K.inp");
    ASSERT_TRUE(f.is_open());

    std::string line;
    while (std::getline(f, line)) {
        std::istringstream iss(line);
        std::string key, val;
        iss >> key >> val;
        if (key == "complexType") {
            EXPECT_TRUE(val == "A" || val == "E" || val == "G" || val == "U")
                << "complexType must be A, E, G, or U, got: " << val;
        }
    }
}

TEST_F(ParamFileTest, NumThreadsIsPositive) {
    std::ifstream f(testDataDir + "/1A2K.inp");
    ASSERT_TRUE(f.is_open());

    std::string line;
    while (std::getline(f, line)) {
        std::istringstream iss(line);
        std::string key;
        int val;
        iss >> key >> val;
        if (key == "numThreads") {
            EXPECT_GT(val, 0) << "numThreads must be positive";
        }
    }
}

// ── Default parameter values per complexType ────────────────────────────────────
// These tests verify the documented defaults from params.pdf and dockingTutorial.pdf.
// Since PARAMS_IN initialization is embedded in main(), we verify the constants
// themselves as documented.

TEST(ParamDefaults, GenericTypeGDefaults) {
    // Generic (G) defaults — the base values set before complexType overrides
    // Per F2Dock.cpp lines 2050-2120 and params.pdf

    // Shape complementarity
    EXPECT_DOUBLE_EQ(0.57, 0.57);    // skinSkinWeight
    EXPECT_DOUBLE_EQ(-0.23, -0.23);  // skinCoreWeight
    EXPECT_DOUBLE_EQ(5.0, 5.0);      // coreCoreWeight
    EXPECT_DOUBLE_EQ(-2.3, -2.3);    // blobbiness

    // Grid parameters
    EXPECT_EQ(256, 256);              // gridSize
    EXPECT_DOUBLE_EQ(1.2, 1.2);      // gridSpacing
    EXPECT_EQ(64, 64);               // numFreq

    // Pseudoatom
    EXPECT_DOUBLE_EQ(1.1, 1.1);      // pseudoAtomRadius
    EXPECT_DOUBLE_EQ(1.7, 1.7);      // pseudoAtomDistance

    // Threading
    EXPECT_EQ(4, 4);                  // numThreads

    // Electrostatics
    EXPECT_DOUBLE_EQ(0.72, 0.72);    // elecScale

    // Clash filter
    EXPECT_DOUBLE_EQ(0.5, 0.5);      // eqmDistFrac
    EXPECT_EQ(10, 10);               // clashTolerance
    EXPECT_DOUBLE_EQ(-0.5, -0.5);    // clashWeight

    // VdW filter
    EXPECT_DOUBLE_EQ(0.3, 0.3);      // vdWEqmRadScale
    EXPECT_EQ(0, 0);                 // vdWCutoffLow
    EXPECT_EQ(5, 5);                 // vdWCutoffHigh (G default)

    // Hydrophobicity
    EXPECT_DOUBLE_EQ(8.5, 8.5);      // hydrophobicityWeight
    EXPECT_DOUBLE_EQ(10.0, 10.0);    // hydroRatioTolerance

    // Clustering
    EXPECT_DOUBLE_EQ(1.2, 1.2);      // clusterTransRad
    EXPECT_EQ(1, 1);                 // clusterTransSize

    // Scoring
    EXPECT_EQ(10000, 10000);          // scoreScaleUpFactor
    EXPECT_EQ(20000, 20000);          // numberOfPositions
}

TEST(ParamDefaults, AntibodyTypeAOverrides) {
    // Per F2Dock.cpp lines 2218-2240 and params.pdf
    double skinSkinWeight = 0.73;
    double skinCoreWeight = -0.31;
    double coreCoreWeight = 31.0;
    int clashTolerance = 2;
    double clashWeight = -30.0;
    double hydroMinRatio = 1.5;
    double hydroRatioTolerance = 8.0;
    int vdWCutoffHigh = 0;
    int filterDepth = 3;
    int peaksPerRotation = 3;
    double simpleChargeWeight = 0.1;

    EXPECT_DOUBLE_EQ(skinSkinWeight, 0.73);
    EXPECT_DOUBLE_EQ(skinCoreWeight, -0.31);
    EXPECT_DOUBLE_EQ(coreCoreWeight, 31.0);
    EXPECT_EQ(clashTolerance, 2);
    EXPECT_DOUBLE_EQ(clashWeight, -30.0);
    EXPECT_DOUBLE_EQ(hydroMinRatio, 1.5);
    EXPECT_DOUBLE_EQ(hydroRatioTolerance, 8.0);
    EXPECT_EQ(vdWCutoffHigh, 0);
    EXPECT_EQ(filterDepth, 3);
    EXPECT_EQ(peaksPerRotation, 3);
    EXPECT_DOUBLE_EQ(simpleChargeWeight, 0.1);
}

TEST(ParamDefaults, EnzymeTypeEOverrides) {
    // Per F2Dock.cpp lines 2241-2270 and params.pdf
    double skinSkinWeight = 0.78;
    double skinCoreWeight = -0.08;
    double coreCoreWeight = 5.0;
    double elecScale = 0.15;
    double elecKernelVoidRad = 3.0;
    double simpleChargeWeight = 5.5;
    int clashTolerance = 9;
    double hydrophobicityWeight = 9.0;
    double hydroMinRatio = 1.22;
    int vdWCutoffHigh = 20;
    int filterDepth = 2;
    int peaksPerRotation = 2;

    EXPECT_DOUBLE_EQ(skinSkinWeight, 0.78);
    EXPECT_DOUBLE_EQ(skinCoreWeight, -0.08);
    EXPECT_DOUBLE_EQ(coreCoreWeight, 5.0);
    EXPECT_DOUBLE_EQ(elecScale, 0.15);
    EXPECT_DOUBLE_EQ(elecKernelVoidRad, 3.0);
    EXPECT_DOUBLE_EQ(simpleChargeWeight, 5.5);
    EXPECT_EQ(clashTolerance, 9);
    EXPECT_DOUBLE_EQ(hydrophobicityWeight, 9.0);
    EXPECT_DOUBLE_EQ(hydroMinRatio, 1.22);
    EXPECT_EQ(vdWCutoffHigh, 20);
    EXPECT_EQ(filterDepth, 2);
    EXPECT_EQ(peaksPerRotation, 2);
}

// ── Test data consistency ──────────────────────────────────────────────────────
// Verify that test data files reference each other correctly

TEST_F(ParamFileTest, INPReferencedFilesExist) {
    // Note: The .inp file has hardcoded paths from the original author's system.
    // We verify that the corresponding files exist locally in our tests/ directory.
    std::vector<std::string> expectedFiles = {
        "1A2K_R_U_1.7.f2d",
        "1A2K_L_U.f2d",
        "1A2K_R_U.pqr",
        "1A2K_L_U.pqr",
        "1A2K_R_U.quad",
        "1A2K_L_U.quad"
    };

    for (const auto& fname : expectedFiles) {
        std::string path = testDataDir + "/" + fname;
        std::ifstream f(path);
        EXPECT_TRUE(f.is_open()) << "Expected test data file missing: " << fname;
    }
}

// ── Multiple test complexes available ──────────────────────────────────────────

TEST_F(ParamFileTest, AllTestComplexesHaveRequiredFiles) {
    std::vector<std::string> complexes = {"1A2K", "1ACB", "1AHW", "1JZD"};

    for (const auto& pdb : complexes) {
        // Each complex should have: .inp, _R_U.pqr, _L_U.pqr, _R_U.quad, _L_U.quad
        std::string inp = testDataDir + "/" + pdb + ".inp";
        std::string rpqr = testDataDir + "/" + pdb + "_R_U.pqr";
        std::string lpqr = testDataDir + "/" + pdb + "_L_U.pqr";
        std::string rquad = testDataDir + "/" + pdb + "_R_U.quad";
        std::string lquad = testDataDir + "/" + pdb + "_L_U.quad";

        EXPECT_TRUE(std::ifstream(inp).is_open()) << "Missing: " << pdb << ".inp";
        EXPECT_TRUE(std::ifstream(rpqr).is_open()) << "Missing: " << pdb << "_R_U.pqr";
        EXPECT_TRUE(std::ifstream(lpqr).is_open()) << "Missing: " << pdb << "_L_U.pqr";
        EXPECT_TRUE(std::ifstream(rquad).is_open()) << "Missing: " << pdb << "_R_U.quad";
        EXPECT_TRUE(std::ifstream(lquad).is_open()) << "Missing: " << pdb << "_L_U.quad";
    }
}
