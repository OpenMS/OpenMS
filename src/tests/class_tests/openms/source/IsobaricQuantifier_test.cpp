// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTThirtyTwoPlexQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricQuantifier, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricQuantifier* ptr = nullptr;
IsobaricQuantifier* null_ptr = nullptr;

ItraqFourPlexQuantitationMethod quant_meth;

START_SECTION((IsobaricQuantifier(const IsobaricQuantitationMethod *const quant_method)))
{
	ptr = new IsobaricQuantifier(&quant_meth);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IsobaricQuantifier())
{
	delete ptr;
}
END_SECTION

START_SECTION((IsobaricQuantifier(const IsobaricQuantifier &other)))
{
  IsobaricQuantifier quantifier(&quant_meth);
  IsobaricQuantifier * quantifier2 = new IsobaricQuantifier(quantifier);

	TEST_NOT_EQUAL(quantifier2, null_ptr)
  delete quantifier2;

  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((IsobaricQuantifier& operator=(const IsobaricQuantifier &rhs)))
{
  IsobaricQuantifier quantifier(&quant_meth);
  IsobaricQuantifier quantifier2(&quant_meth);

  quantifier2 = quantifier;

  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void quantify(const ConsensusMap &consensus_map_in, ConsensusMap &consensus_map_out)))
{
	ConsensusXMLFile cm_file;
	ConsensusMap cm_in, cm_out;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("IsobaricQuantifier.consensusXML"),cm_in);

	IsobaricQuantifier iq(&quant_meth);
	Param p;
	p.setValue("normalization", "true");
  p.setValue("isotope_correction", "true");
	iq.setParameters(p);
	iq.quantify(cm_in,cm_out);

	String cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_out);

	WHITELIST("<?xml-stylesheet");
	// WHITELIST("<?xml-stylesheet,consensusElement id=");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("IsobaricQuantifier_out.consensusXML"));
}
END_SECTION

// TMT 32-plex: verify that the identity matrix is always used (user-supplied correction_matrix is ignored)
// and that quantify() preserves intensities when isotope correction is disabled.
START_SECTION(([EXTRA] TMT 32plex identity matrix enforced in quantify))
{
  TMTThirtyTwoPlexQuantitationMethod tmt32;

  // Verify getIsotopeCorrectionMatrix() returns identity even after setting non-identity values
  {
    Param p = tmt32.getParameters();
    // Try to set a non-identity correction matrix (with some non-zero off-diagonal values)
    std::vector<std::string> fake_correction(32, "1.0/0.5/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA");
    p.setValue("correction_matrix", fake_correction);
    tmt32.setParameters(p);

    Matrix<double> m = tmt32.getIsotopeCorrectionMatrix();
    // Must still be identity
    for (Size r = 0; r < m.rows(); ++r)
    {
      for (Size c = 0; c < m.cols(); ++c)
      {
        if (r == c)
        {
          TEST_REAL_SIMILAR(m(r, c), 1.0)
        }
        else
        {
          TEST_REAL_SIMILAR(m(r, c), 0.0)
        }
      }
    }
  }

  // Build synthetic ConsensusMap with known intensities and run quantify()
  {
    const Size n_channels = 32;
    TMTThirtyTwoPlexQuantitationMethod tmt32_q;
    auto channels = tmt32_q.getChannelInformation();

    ConsensusMap cm_in;
    // Set up column headers (one per channel)
    for (Size i = 0; i < n_channels; ++i)
    {
      ConsensusMap::ColumnHeader header;
      header.setMetaValue("channel_name", channels[i].name);
      cm_in.getColumnHeaders()[i] = header;
    }

    // Create one consensus feature with known intensities per channel
    ConsensusFeature cf;
    cf.setRT(100.0);
    cf.setMZ(500.0);
    double total_intensity = 0.0;
    for (Size i = 0; i < n_channels; ++i)
    {
      FeatureHandle fh;
      fh.setMapIndex(i);
      fh.setIntensity(1000.0 * (i + 1));
      fh.setRT(100.0);
      fh.setMZ(channels[i].center);
      fh.setUniqueId(i + 1);
      cf.insert(fh);
      total_intensity += 1000.0 * (i + 1);
    }
    cf.setIntensity(total_intensity);
    cm_in.push_back(cf);

    // Run quantify with isotope_correction disabled (identity matrix would throw in corrector)
    IsobaricQuantifier iq(&tmt32_q);
    Param p;
    p.setValue("isotope_correction", "false");
    p.setValue("normalization", "false");
    iq.setParameters(p);

    ConsensusMap cm_out;
    iq.quantify(cm_in, cm_out);

    TEST_EQUAL(cm_out.size(), 1)
    ABORT_IF(cm_out.size() != 1)

    // Intensities must be preserved (no correction applied)
    auto cf_it = cm_out[0].begin();
    for (Size i = 0; i < n_channels; ++i)
    {
      TEST_REAL_SIMILAR(cf_it->getIntensity(), 1000.0 * (i + 1))
      ++cf_it;
    }
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
