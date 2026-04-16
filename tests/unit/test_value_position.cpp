/*
 * Unit tests for ValuePosition3D.
 *
 * ValuePosition3D stores a single docking solution with all its scoring
 * components.  Per the documentation output format (§4.2.4), each solution has:
 *   rank, score, shape, ssr, ccr, scr, elec, hbond, hydrophobicity,
 *   smplcomp, vdw, clashes, pgsol, pgsolh, deldispe,
 *   mat1-mat12 (transformation matrix), conf, rmsd
 */

#include <gtest/gtest.h>
#include "f2dock/ValuePosition3D.h"

TEST(ValuePosition3D, DefaultConstruction) {
    ValuePosition3D vp;
    // After construction, verify reset state.
    // Per ValuePosition3D.cpp: reset() uses -100000000 as a sentinel for
    // "no valid solution" on primary score fields.
    vp.reset();
    EXPECT_DOUBLE_EQ(vp.m_Value, -100000000.0);
    EXPECT_DOUBLE_EQ(vp.m_RealValue, -100000000.0);
    EXPECT_DOUBLE_EQ(vp.m_ImaginaryValue, -100000000.0);
    EXPECT_DOUBLE_EQ(vp.m_elecValue, 0.0);
    EXPECT_DOUBLE_EQ(vp.m_hbondValue, 0.0);
    EXPECT_DOUBLE_EQ(vp.m_hydrophobicityValue, 0.0);
    EXPECT_DOUBLE_EQ(vp.m_vdWPotential, 0.0);
    EXPECT_DOUBLE_EQ(vp.m_simpComp, 0.0);
    EXPECT_DOUBLE_EQ(vp.m_pGsol, 0.0);
    EXPECT_DOUBLE_EQ(vp.m_pGsolH, 0.0);
    EXPECT_DOUBLE_EQ(vp.m_delDispE, 0.0);
    EXPECT_EQ(vp.m_nClashes, 0);
    EXPECT_EQ(vp.m_RotationIndex, 0);
    EXPECT_EQ(vp.m_FineRotationIndex, 0);
    EXPECT_EQ(vp.m_ConformationIndex, -1);
    EXPECT_DOUBLE_EQ(vp.m_Translation[0], 0.0);
    EXPECT_DOUBLE_EQ(vp.m_Translation[1], 0.0);
    EXPECT_DOUBLE_EQ(vp.m_Translation[2], 0.0);
    // origScore and rerankerScore reset to m_Value sentinel
    EXPECT_DOUBLE_EQ(vp.m_origScore, -100000000.0);
    EXPECT_DOUBLE_EQ(vp.m_rerankerScore, -100000000.0);
    EXPECT_EQ(vp.m_origRank, -1);
    EXPECT_EQ(vp.m_rerankerRank, -1);
}

TEST(ValuePosition3D, SetAndCopy) {
    ValuePosition3D a;
    a.m_Value = 100.0;
    a.m_elecValue = -5.5;
    a.m_nClashes = 3;
    a.m_Translation[0] = 1.0;
    a.m_Translation[1] = 2.0;
    a.m_Translation[2] = 3.0;
    a.m_RotationIndex = 42;
    a.m_hydrophobicityValue = 7.3;
    a.m_pGsol = -12.4;

    ValuePosition3D b;
    b = a;

    EXPECT_DOUBLE_EQ(b.m_Value, 100.0);
    EXPECT_DOUBLE_EQ(b.m_elecValue, -5.5);
    EXPECT_EQ(b.m_nClashes, 3);
    EXPECT_DOUBLE_EQ(b.m_Translation[0], 1.0);
    EXPECT_DOUBLE_EQ(b.m_Translation[1], 2.0);
    EXPECT_DOUBLE_EQ(b.m_Translation[2], 3.0);
    EXPECT_EQ(b.m_RotationIndex, 42);
    EXPECT_DOUBLE_EQ(b.m_hydrophobicityValue, 7.3);
    EXPECT_DOUBLE_EQ(b.m_pGsol, -12.4);
}

TEST(ValuePosition3D, SetMethod) {
    ValuePosition3D a;
    a.m_Value = 55.5;
    a.m_vdWPotential = -8.1;
    a.m_SkinSkinRealValue = 3.3;
    a.m_CoreCoreRealValue = 1.1;
    a.m_SkinCoreRealValue = -2.2;

    ValuePosition3D b;
    b.set(a);

    EXPECT_DOUBLE_EQ(b.m_Value, 55.5);
    EXPECT_DOUBLE_EQ(b.m_vdWPotential, -8.1);
    EXPECT_DOUBLE_EQ(b.m_SkinSkinRealValue, 3.3);
    EXPECT_DOUBLE_EQ(b.m_CoreCoreRealValue, 1.1);
    EXPECT_DOUBLE_EQ(b.m_SkinCoreRealValue, -2.2);
}

TEST(ValuePosition3D, ResetClearsAll) {
    ValuePosition3D vp;
    vp.m_Value = 999.0;
    vp.m_nClashes = 100;
    vp.m_Translation[0] = 50.0;
    vp.m_pGsolH = -20.0;
    vp.m_delDispE = -15.0;

    vp.reset();

    // Primary score fields reset to -100000000 sentinel
    EXPECT_DOUBLE_EQ(vp.m_Value, -100000000.0);
    // Filter/energy fields reset to 0
    EXPECT_EQ(vp.m_nClashes, 0);
    EXPECT_DOUBLE_EQ(vp.m_Translation[0], 0.0);
    EXPECT_DOUBLE_EQ(vp.m_pGsolH, 0.0);
    EXPECT_DOUBLE_EQ(vp.m_delDispE, 0.0);
}

// ── Output format field coverage ───────────────────────────────────────────────
// Per documentation §4.2.4, the output format has these columns:
// rank score shape ssr ccr scr elec hbond hydrophobicity smplcomp vdw
// clashes pgsol pgsolh deldispe mat1..mat12 conf rmsd
//
// ValuePosition3D should store all of these.

TEST(ValuePosition3D, AllOutputFieldsAccessible) {
    ValuePosition3D vp;

    // Score components
    vp.m_Value = 1000.0;              // score
    vp.m_RealValue = 800.0;           // shape
    vp.m_SkinSkinRealValue = 500.0;   // ssr
    vp.m_CoreCoreRealValue = 200.0;   // ccr
    vp.m_SkinCoreRealValue = 100.0;   // scr
    vp.m_elecValue = -10.0;           // elec
    vp.m_hbondValue = -5.0;           // hbond
    vp.m_hydrophobicityValue = 3.0;   // hydrophobicity
    vp.m_simpComp = 2.0;              // smplcomp
    vp.m_vdWPotential = -20.0;        // vdw
    vp.m_nClashes = 5;                // clashes
    vp.m_pGsol = -15.0;              // pgsol
    vp.m_pGsolH = -12.0;             // pgsolh
    vp.m_delDispE = -8.0;            // deldispe

    // Transformation (mat1..mat12 = 3×4 rotation+translation matrix)
    vp.m_Translation[0] = 1.5;
    vp.m_Translation[1] = 2.5;
    vp.m_Translation[2] = 3.5;
    vp.m_RotationIndex = 10;
    vp.m_FineRotationIndex = 20;
    vp.m_ConformationIndex = 0;       // conf

    // Verify all are stored correctly
    EXPECT_DOUBLE_EQ(vp.m_Value, 1000.0);
    EXPECT_DOUBLE_EQ(vp.m_SkinSkinRealValue, 500.0);
    EXPECT_DOUBLE_EQ(vp.m_elecValue, -10.0);
    EXPECT_EQ(vp.m_nClashes, 5);
    EXPECT_DOUBLE_EQ(vp.m_pGsol, -15.0);
    EXPECT_DOUBLE_EQ(vp.m_delDispE, -8.0);
}

// ── Reranking fields ───────────────────────────────────────────────────────────
// Per documentation §5, GB-rerank produces a reranked score

TEST(ValuePosition3D, RerankingFields) {
    ValuePosition3D vp;
    vp.m_origScore = 500.0;
    vp.m_origRank = 1;
    vp.m_rerankerScore = 750.0;
    vp.m_rerankerRank = 3;
    vp.m_clusterPenalty = -10.0;

    EXPECT_DOUBLE_EQ(vp.m_origScore, 500.0);
    EXPECT_EQ(vp.m_origRank, 1);
    EXPECT_DOUBLE_EQ(vp.m_rerankerScore, 750.0);
    EXPECT_EQ(vp.m_rerankerRank, 3);
    EXPECT_DOUBLE_EQ(vp.m_clusterPenalty, -10.0);
}
