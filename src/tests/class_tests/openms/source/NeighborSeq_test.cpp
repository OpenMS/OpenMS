// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/ID/NeighborSeq.h>
#include <vector>


using namespace OpenMS;
using namespace std;

START_TEST(NeighborSeq, "$Id$")

// NS()=delete;


// Test section for the generateSpectrum function
// The spectra were generated via TOPPView and contained b-and y-ion
START_SECTION(MSSpectrum generateSpectrum(const String& peptide_sequence))
{
  NeighborSeq ns({AASequence::fromString("TEST")});
  MSSpectrum spec_1 = ns.generateSpectrum(AASequence::fromString("PEPT"));
  MSSpectrum spec_2 = ns.generateSpectrum(AASequence::fromString("AR"));
  MSSpectrum spec_3 = ns.generateSpectrum(AASequence::fromString("VGLPINQR"));

  // Test "PEPT"
  TEST_REAL_SIMILAR(spec_1[0].getPosition()[0], 120.0655);
  TEST_REAL_SIMILAR(spec_1[1].getPosition()[0], 217.1182);
  TEST_REAL_SIMILAR(spec_1[2].getPosition()[0], 227.1026);
  TEST_REAL_SIMILAR(spec_1[3].getPosition()[0], 324.1553);
  TEST_REAL_SIMILAR(spec_1[4].getPosition()[0], 346.1608);
  

  // Test "AR"
  TEST_REAL_SIMILAR(spec_2[0].getPosition()[0], 175.1189);

  // Test "VGLPINQR"
  TEST_REAL_SIMILAR(spec_3[0].getPosition()[0], 157.0971);
  TEST_REAL_SIMILAR(spec_3[1].getPosition()[0], 175.1189);
  TEST_REAL_SIMILAR(spec_3[2].getPosition()[0], 270.1812);
  TEST_REAL_SIMILAR(spec_3[3].getPosition()[0], 303.1775);
  TEST_REAL_SIMILAR(spec_3[4].getPosition()[0], 367.2339);
  TEST_REAL_SIMILAR(spec_3[5].getPosition()[0], 417.2204);
  TEST_REAL_SIMILAR(spec_3[6].getPosition()[0], 480.3180);
  TEST_REAL_SIMILAR(spec_3[7].getPosition()[0], 530.3045);
  TEST_REAL_SIMILAR(spec_3[8].getPosition()[0], 594.3609);
  TEST_REAL_SIMILAR(spec_3[9].getPosition()[0], 627.3578);
  TEST_REAL_SIMILAR(spec_3[10].getPosition()[0], 722.4195);
  TEST_REAL_SIMILAR(spec_3[11].getPosition()[0], 740.4413);
  TEST_REAL_SIMILAR(spec_3[12].getPosition()[0], 797.4628);

}
END_SECTION

// Test section for the compareSpectra function
START_SECTION(
  static bool isNeighborSpectrum(const MSSpectrum& spec1, const MSSpectrum& spec2, const double min_shared_ion_fraction, const double mz_bin_size))
{
  MSSpectrum spec1({Peak1D(100.00, 1.0),
                    Peak1D(200.00, 1.0),
                    Peak1D(300.00, 1.0)});

  MSSpectrum spec2({Peak1D(100.05, 1.0),
                    Peak1D(200.05, 1.0),
                    Peak1D(300.05, 1.0)});

  MSSpectrum spec3({Peak1D(101.00, 1.0),
                    Peak1D(201.00, 1.0),
                    Peak1D(301.00, 1.0)});

  MSSpectrum spec4({Peak1D(100.05, 1.0),
                    Peak1D(201.00, 1.0),
                    Peak1D(300.05, 1.0),
                    Peak1D(301.00, 1.0)});

  double min_shared_ion_fraction = 0.5; 
  
  NeighborSeq ns({AASequence::fromString("TEST")});

  // bin interval is from [a,b[
  TEST_TRUE (ns.isNeighborSpectrum(spec1, spec2, min_shared_ion_fraction, 1.0))
  TEST_FALSE(ns.isNeighborSpectrum(spec1, spec3, min_shared_ion_fraction, 1.0))
  TEST_TRUE (ns.isNeighborSpectrum(spec1, spec4, min_shared_ion_fraction, 1.0))
  TEST_FALSE(ns.isNeighborSpectrum(spec2, spec3, min_shared_ion_fraction, 1.0))
  TEST_TRUE (ns.isNeighborSpectrum(spec2, spec4, min_shared_ion_fraction, 1.0))
  TEST_TRUE (ns.isNeighborSpectrum(spec3, spec4, min_shared_ion_fraction, 1.0))

  TEST_FALSE(ns.isNeighborSpectrum(spec1, spec2, min_shared_ion_fraction, 0.05))
  TEST_FALSE(ns.isNeighborSpectrum(spec1, spec3, min_shared_ion_fraction, 0.05))
  TEST_FALSE(ns.isNeighborSpectrum(spec1, spec4, min_shared_ion_fraction, 0.05))
  TEST_FALSE(ns.isNeighborSpectrum(spec2, spec3, min_shared_ion_fraction, 0.05))
  TEST_TRUE(ns.isNeighborSpectrum(spec2, spec4, min_shared_ion_fraction, 0.05))
  TEST_TRUE(ns.isNeighborSpectrum(spec3, spec4, min_shared_ion_fraction, 0.05))
}
END_SECTION


// Test section for the findCandidatePositions function
START_SECTION(static int computeSharedIonCount(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& mz_bin_size))
{
  MSSpectrum spec1({Peak1D(100.00, 1.0),
                    Peak1D(200.00, 1.0),
                    Peak1D(300.00, 1.0)});

  MSSpectrum spec2({Peak1D(100.05, 1.0),
                    Peak1D(200.05, 1.0),
                    Peak1D(300.05, 1.0)});

  MSSpectrum spec3({Peak1D(101.00, 1.0),
                    Peak1D(201.00, 1.0),
                    Peak1D(301.00, 1.0)});

  MSSpectrum spec4({Peak1D(100.05, 1.0),
                    Peak1D(201.00, 1.0),
                    Peak1D(300.05, 1.0),
                    Peak1D(301.00, 1.0)});

  NeighborSeq ns({AASequence::fromString("TEST")});

  // bin interval is from [a,b[
  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec2, 2.0), 3)
  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec3, 2.0), 3)
  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec4, 2.0), 3)
  TEST_EQUAL(ns.computeSharedIonCount(spec2, spec3, 2.0), 3)
  TEST_EQUAL(ns.computeSharedIonCount(spec2, spec4, 2.0), 3)
  TEST_EQUAL(ns.computeSharedIonCount(spec3, spec4, 2.0), 3)

  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec2, 1.0), 3)
  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec3, 1.0), 0)
  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec4, 1.0), 2)
  TEST_EQUAL(ns.computeSharedIonCount(spec2, spec3, 1.0), 0)
  TEST_EQUAL(ns.computeSharedIonCount(spec2, spec4, 1.0), 2)
  TEST_EQUAL(ns.computeSharedIonCount(spec3, spec4, 1.0), 2)

  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec2, 0.1), 3)
  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec3, 0.1), 0)
  TEST_EQUAL(ns.computeSharedIonCount(spec1, spec4, 0.1), 2)
  TEST_EQUAL(ns.computeSharedIonCount(spec2, spec3, 0.1), 0)
  TEST_EQUAL(ns.computeSharedIonCount(spec2, spec4, 0.1), 2)
  TEST_EQUAL(ns.computeSharedIonCount(spec3, spec4, 0.1), 2)
}
END_SECTION

START_SECTION(bool isNeighborPeptide(const AASequence& neighbor_candidate,
                             const double mass_tolerance_pc,
                             const bool mass_tolerance_pc_ppm,
                             const double min_shared_ion_fraction,
                             const double mz_bin_size))
{
  const AASequence AA_VELQSK = AASequence::fromString("VELQSK");
  const AASequence AA_VESQLK = AASequence::fromString("VESQLK");
  std::vector<AASequence> seqs = {AASequence::fromString("VELQSK"), AASequence::fromString("SVQELK"), AASequence::fromString("TVDQLK"), AASequence::fromString("VGEFK")};
  NeighborSeq ns(std::move(seqs));
  // VELQSK has neighbor VESQLK
  // SVQELK has neighbor VESQLK
  // TVDQLK has neighbor VESQLK
  // VGEFK  has neighbor GLDFK
  // GIDFK  has neighbor GLDFK

  const double pc_tolerance = 0.01;
  const double mz_bin_size = 0.05;
  TEST_TRUE(std::abs(AA_VELQSK.getMonoWeight() - AA_VESQLK.getMonoWeight()) < pc_tolerance)
  TEST_TRUE(ns.computeSharedIonCount(ns.generateSpectrum(AA_VELQSK), ns.generateSpectrum(AA_VESQLK), mz_bin_size) == 5)

  TEST_TRUE(ns.isNeighborPeptide(AASequence::fromString("VESQLK"), pc_tolerance, false, 0.25, mz_bin_size))
  TEST_FALSE(ns.isNeighborPeptide(AASequence::fromString("VESQLK"), pc_tolerance, false, 5*2.0/( (6-1)*2*2) + 0.1, mz_bin_size))
  TEST_TRUE(ns.isNeighborPeptide(AASequence::fromString("GLDFK"), pc_tolerance, false, 0.25, mz_bin_size))

  auto stats = ns.getNeighborStats();
  TEST_EQUAL(stats.unfindable_peptides, 0)
  TEST_EQUAL(stats.findable_no_neighbors, 0)
  TEST_EQUAL(stats.findable_one_neighbor, 4)
  TEST_EQUAL(stats.findable_multiple_neighbors, 0)

  // test one neighbor again, which causes 3 candidates to have multiple neighbors
  TEST_TRUE(ns.isNeighborPeptide(AASequence::fromString("VESQLK"), 0.01, false, 0.25, 0.05))
  auto stats2 = ns.getNeighborStats();
  TEST_EQUAL(stats2.unfindable_peptides, 0)
  TEST_EQUAL(stats2.findable_no_neighbors, 0)
  TEST_EQUAL(stats2.findable_one_neighbor, 1)
  TEST_EQUAL(stats2.findable_multiple_neighbors, 3)
}
END_SECTION

START_SECTION(NeighborStats getNeighborStats() const)
{
  NOT_TESTABLE // tested above
}
END_SECTION

END_TEST
