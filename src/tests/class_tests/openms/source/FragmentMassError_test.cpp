// --------------------------------------------------------------------------
//           OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/QC/FragmentMassError.h>

#include <vector>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

//////////////////////////

using namespace OpenMS;

START_TEST(FragmentMassError, "$Id$")

  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////

  //--------------------------------------------------------------------
  // Input Data
  //--------------------------------------------------------------------

  //-------------------------------------------------
  // create MSExperiment
  //-------------------------------------------------

  // Precursor
  Precursor precursor_cid;
  std::set<Precursor::ActivationMethod> am_cid;
  am_cid.insert(Precursor::ActivationMethod::CID);
  precursor_cid.setActivationMethods(am_cid);

  Precursor precursor_ecd;
  std::set<Precursor::ActivationMethod> am_ecd;
  am_ecd.insert(Precursor::ActivationMethod::ECD);
  precursor_ecd.setActivationMethods(am_ecd);

  Precursor precursor_sori;
  std::set<Precursor::ActivationMethod> am_sori;
  am_sori.insert(Precursor::ActivationMethod::SORI);
  precursor_sori.setActivationMethods(am_sori);


  // create spectrum for HIMALAYA
  // generate b- and y-ion spectrum of peptide seq with charge 1
  // set precursor, RT and MSLevel
  PeakSpectrum ms_spec_2_himalaya;
  TheoreticalSpectrumGenerator theo_gen_hi;
  theo_gen_hi.getSpectrum(ms_spec_2_himalaya, AASequence::fromString("HIMALAYA"), 1, 1);
  ms_spec_2_himalaya.setPrecursors({precursor_cid});
  ms_spec_2_himalaya.setRT(3.7);
  ms_spec_2_himalaya.setMSLevel(2);

  // create spectrum for ALABAMA
  // generate c- and z-ion spectrum of peptide seq with charge 2
  // set precursor, RT and MSLevel
  PeakSpectrum ms_spec_2_alabama;
  TheoreticalSpectrumGenerator theo_gen_al;
  Param theo_gen_settings_al = theo_gen_al.getParameters();
  theo_gen_settings_al.setValue("add_c_ions", "true");
  theo_gen_settings_al.setValue("add_z_ions", "true");
  theo_gen_settings_al.setValue("add_b_ions", "false");
  theo_gen_settings_al.setValue("add_y_ions", "false");
  theo_gen_al.setParameters(theo_gen_settings_al);
  theo_gen_al.getSpectrum(ms_spec_2_alabama, AASequence::fromString("ALABAMA"), 2, 2);
  ms_spec_2_alabama.setPrecursors({precursor_ecd});
  ms_spec_2_alabama.setRT(2);
  ms_spec_2_alabama.setMSLevel(2);


  //shift every peak by 0.001 mz
  for(Peak1D& peak : ms_spec_2_himalaya)
  {
    double mz = peak.getMZ() + 0.001;
    peak.setMZ(mz);
  };

  for(Peak1D& peak : ms_spec_2_alabama)
  {
    double mz = peak.getMZ() + 0.001;
    peak.setMZ(mz);
  };

  // MSSpectrum with fragmentation method SORI (not supported)
  MSSpectrum ms_spec_2_sori;
  ms_spec_2_sori.setPrecursors({precursor_sori});
  ms_spec_2_sori.setRT(7);
  ms_spec_2_sori.setMSLevel(2);

  //MSSpectrum with no given fragmentation method given
  MSSpectrum ms_spec_2_no_precursor;
  ms_spec_2_no_precursor.setRT(8);
  ms_spec_2_no_precursor.setMSLevel(2);

  // MS1Spectrum
  MSSpectrum ms_spec_1;
  ms_spec_1.setRT(5);
  ms_spec_1.setMSLevel(1);
  ms_spec_1.setPrecursors({precursor_cid});

  // MSSpectrum without peaks
  MSSpectrum ms_spec_2_excp;
  ms_spec_2_excp.setRT(4);
  ms_spec_2_excp.setMSLevel(2);
  ms_spec_2_excp.setPrecursors({precursor_cid});

  // MS_Empty
  MSSpectrum ms_spec_empty;

  // MSExperiment with unsorted spectrum
  MSExperiment exp_unsort;
  exp_unsort.setSpectra({ ms_spec_1, ms_spec_2_himalaya, ms_spec_2_alabama, ms_spec_2_excp, ms_spec_empty, ms_spec_2_sori, ms_spec_2_no_precursor});

  // sorted MSSpectrum
  MSExperiment exp(exp_unsort);
  exp.sortSpectra();

  //-------------------------------------------------
  // create FeatureMap
  //-------------------------------------------------

  // PepHit Himalaya
  PeptideHit pep_hit_hi;
  pep_hit_hi.setSequence(AASequence::fromString("HIMALAYA"));
  pep_hit_hi.setCharge(1);

  // PepHit Alabama
  PeptideHit pep_hit_al;
  pep_hit_al.setSequence(AASequence::fromString("ALABAMA"));
  pep_hit_al.setCharge(2);

  // PeptideHit Peptide
  PeptideHit pep_hit_pe;
  pep_hit_pe.setSequence(AASequence::fromString("PEPTIDE"));
  pep_hit_pe.setCharge(3);

  // PepID empty
  PeptideIdentification pep_id_empty;
  pep_id_empty.setRT(6);

  // PepID himalaya
  PeptideIdentification pep_id_hi;
  pep_id_hi.setRT(3.72);
  pep_id_hi.setMZ(888);
  pep_id_hi.setHits({pep_hit_hi});

  // PepID alabama
  PeptideIdentification pep_id_al;
  pep_id_al.setRT(2);
  pep_id_al.setMZ(264);
  pep_id_al.setHits({pep_hit_al});

  // PepID with RT out of tolerance
  PeptideIdentification pep_id_tol_out;
  pep_id_tol_out.setRT(2.1);
  pep_id_tol_out.setMZ(266);
  pep_id_tol_out.setHits({pep_hit_pe});

  // PepID matches with ms1 spectrum
  PeptideIdentification pep_id_ms1;
  pep_id_ms1.setRT(5);
  pep_id_ms1.setMZ(266);
  pep_id_ms1.setHits({pep_hit_pe});

  // PepID peak_RT does not exist in msExp
  PeptideIdentification pep_id_notExist;
  pep_id_notExist.setRT(10);
  pep_id_notExist.setMZ(266);
  pep_id_notExist.setHits({pep_hit_pe});

  // PepID matches with spectrum without peaks
  PeptideIdentification pep_id_excp;
  pep_id_excp.setRT(4);
  pep_id_excp.setMZ(266);
  pep_id_excp.setHits({pep_hit_pe});

  // PepID matches spectrum which have no Precursor
  PeptideIdentification pep_id_no_precursor;
  pep_id_no_precursor.setRT(8);
  pep_id_no_precursor.setMZ(266);
  pep_id_no_precursor.setHits({pep_hit_pe});

  // PepID matches spectrum with fragmentation method sori
  PeptideIdentification pep_id_sori;
  pep_id_sori.setRT(7);
  pep_id_sori.setMZ(266);
  pep_id_sori.setHits({pep_hit_pe});

  // Feature valid data
  Feature feat_valid;
  feat_valid.setPeptideIdentifications({pep_id_hi, pep_id_al});

  // Feature empty
  Feature feat_empty;

  // FeatureMap valid data
  FeatureMap fmap;
  fmap.setUnassignedPeptideIdentifications({pep_id_empty});
  fmap.push_back(feat_empty);
  fmap.push_back(feat_valid);

  // FeatureMap tol_out
  FeatureMap fmap_tol_out;
  fmap_tol_out.setUnassignedPeptideIdentifications({pep_id_tol_out});

  // FeatureMap ms1
  FeatureMap fmap_ms1;
  fmap_ms1.setUnassignedPeptideIdentifications({pep_id_ms1});

  // FeatureMap not exist
  FeatureMap fmap_notExist;
  fmap_notExist.setUnassignedPeptideIdentifications({pep_id_notExist});

  // FeatureMap for no peak
  FeatureMap fmap_excp;
  fmap_excp.setUnassignedPeptideIdentifications({pep_id_excp});

  // FeatureMap no precursor
  FeatureMap fmap_no_precursor;
  fmap_no_precursor.setUnassignedPeptideIdentifications({pep_id_no_precursor});

  // FeatureMap with fragmentation method sori
  FeatureMap fmap_sori;
  fmap_sori.setUnassignedPeptideIdentifications({pep_id_sori});


  //--------------------------------------------------------------------
  // Tests
  //--------------------------------------------------------------------

  //////////////////////////////////////////////////////////////////
  // start Section
  /////////////////////////////////////////////////////////////////

  FragmentMassError* ptr = nullptr;
  FragmentMassError* nulpt = nullptr;
  START_SECTION(FragmentMassError())
  {
    ptr = new FragmentMassError();
    TEST_NOT_EQUAL(ptr, nulpt)
  }
  END_SECTION

  START_SECTION(~FragmentMassError())
  {
    delete ptr;
  }
  END_SECTION


  FragmentMassError frag_ma_err;
  FragmentMassError frag_ma_err_flag;
  FragmentMassError frag_ma_err_tol_out;
  FragmentMassError frag_ma_err_ms1;
  FragmentMassError frag_ma_err_notExits;
  FragmentMassError frag_ma_err_excp;
  FragmentMassError frag_ma_err_no_precusor;
  FragmentMassError frag_ma_err_sori;

  //tests compute function
  START_SECTION(void compute(FeatureMap& fmap, const MSExperiment& exp, const double tolerance = 20, const String tolerance_unit = "ppm"))
  {
    //test with valid input
    frag_ma_err.compute(fmap, exp);
    std::vector<FragmentMassError::FMEStatistics> result;
    result = frag_ma_err.getResults();

    TEST_REAL_SIMILAR(result[0].average_ppm, 4.698439) // mz: 0.001
    TEST_REAL_SIMILAR(result[0].variance_ppm, 10.93094) //mz: 5.915844


    //test with valid input and flags
    frag_ma_err_flag.compute(fmap, exp, 1, "mz");
    std::vector<FragmentMassError::FMEStatistics> result_flag;
    result_flag = frag_ma_err_flag.getResults();

    TEST_REAL_SIMILAR(result_flag[0].average_ppm, 5.938193)
    TEST_REAL_SIMILAR(result_flag[0].variance_ppm, 36.45247
    )

    // test if MSExperiment is not sorted
    TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, frag_ma_err_ms1.compute(fmap, exp_unsort),"MSExperiment is not sorted by ascending RT")

    // test with matching ms1 spectrum
    TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, frag_ma_err_ms1.compute(fmap_ms1, exp), "The matching retention time of the mzML is not a MS2 Spectrum.")

    // test if RT from FeatureMap does not match to any RT in MSExp
    TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, frag_ma_err_notExits.compute(fmap_notExist, exp);, "The retention time of the mzML and featureXML file does not match.")

    // test with RT out of tolerance
    TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, frag_ma_err_tol_out.compute(fmap_tol_out, exp), "PeptideID with RT 2.1 s does not have a matching MS2 Spectrum. Closest RT was 3.7, which seems to far off.")

    // test with no fragmentation method given
    TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, frag_ma_err_no_precusor.compute(fmap_no_precursor, exp), "No fragmentation method given.")

    // test with fragmentation method sori, which is not supported
    TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, frag_ma_err_sori.compute(fmap_sori, exp), "Fragmentation method is not supported.")

    // test if spectrum has no peaks
    frag_ma_err_excp.compute(fmap_excp, exp);
    std::vector<FragmentMassError::FMEStatistics> result_excp;
    result_excp = frag_ma_err_excp.getResults();

    TEST_REAL_SIMILAR(result_excp[0].average_ppm, 0)
    TEST_REAL_SIMILAR(result_excp[0].variance_ppm, 0)


  }
  END_SECTION


  START_SECTION(QCBase::Status requires() const override)
  {
    QCBase::Status stat = QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
    TEST_EQUAL(stat, frag_ma_err.requires())
  }
  END_SECTION


  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
END_TEST


