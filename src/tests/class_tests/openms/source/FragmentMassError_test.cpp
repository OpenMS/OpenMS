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

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/QC/FragmentMassError.h>
#include <vector>
//////////////////////////

using namespace OpenMS;

// Functions to create input data

// create a MSSpectrum with Precursor, MSLevel and RT
const MSSpectrum createMSSpectrum(UInt ms_level, double rt, const String& id, Precursor::ActivationMethod precursor_method = Precursor::ActivationMethod::CID)
{
  Precursor precursor;
  std::set<Precursor::ActivationMethod> am;
  am.insert(precursor_method);
  precursor.setActivationMethods(am);

  MSSpectrum ms_spec;
  ms_spec.setRT(rt);
  ms_spec.setMSLevel(ms_level);
  ms_spec.setPrecursors({precursor});
  ms_spec.setNativeID(id);

  return ms_spec;
}

// create a PeptideIdentifiaction with a PeptideHit (sequence, charge), rt and mz
// default values for sequence PEPTIDE
const PeptideIdentification createPeptideIdentification(const String& id, const String& sequence = String("PEPTIDE"), Int charge = 3, double mz = 266)
{
  PeptideHit peptide_hit;
  peptide_hit.setSequence(AASequence::fromString(sequence));
  peptide_hit.setCharge(charge);

  PeptideIdentification peptide_id;
  peptide_id.setMetaValue("spectrum_reference",id);
  peptide_id.setMZ(mz);
  peptide_id.setHits({peptide_hit});

  return peptide_id;
}

START_TEST(FragmentMassError, "$Id$")

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

  //tests compute function
  START_SECTION(void compute(FeatureMap& fmap, const MSExperiment& exp, const std::map<String,UInt64>& map_to_spectrum, const ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, const double tolerance = 20))
  {
    //--------------------------------------------------------------------
    // create valid input data
    //--------------------------------------------------------------------
    // FeatureMap
    FeatureMap fmap;

    // empty PeptideIdentification
    PeptideIdentification pep_id_empty;
    pep_id_empty.setRT(6);

    // empty Feature
    Feature feat_empty;

    // put valid data in fmap
    fmap.setUnassignedPeptideIdentifications({createPeptideIdentification("XTandem::1", "HIMALAYA", 1, 888), createPeptideIdentification("XTandem::2", "ALABAMA", 2, 264), pep_id_empty});
    fmap.push_back(feat_empty);
    // set ProteinIdentifications
    ProteinIdentification protId;
    ProteinIdentification::SearchParameters param;
    param.fragment_mass_tolerance_ppm = false;
    param.fragment_mass_tolerance = 0.3;
    protId.setSearchParameters(param);
    fmap.setProteinIdentifications( {protId} );

    // MSExperiment
    MSExperiment exp;

    // create b- and y-ion spectrum of peptide sequence HIMALAYA with charge 1
    // shift every peak by 5 ppm
    PeakSpectrum ms_spec_2_himalaya = createMSSpectrum(2, 3.7, "XTandem::1");
    TheoreticalSpectrumGenerator theo_gen_hi;
    theo_gen_hi.getSpectrum(ms_spec_2_himalaya, AASequence::fromString("HIMALAYA"), 1, 1);
    for(Peak1D& peak : ms_spec_2_himalaya) peak.setMZ(Math::ppmToMass(peak.getMZ(), 5.0) + peak.getMZ());

    // create c- and z-ion spectrum of peptide sequence ALABAMA with charge 2
    // shift every peak by 5 ppm
    PeakSpectrum ms_spec_2_alabama = createMSSpectrum(2, 2, "XTandem::2", Precursor::ActivationMethod::ECD);
    TheoreticalSpectrumGenerator theo_gen_al;
    Param theo_gen_settings_al = theo_gen_al.getParameters();
    theo_gen_settings_al.setValue("add_c_ions", "true");
    theo_gen_settings_al.setValue("add_z_ions", "true");
    theo_gen_settings_al.setValue("add_b_ions", "false");
    theo_gen_settings_al.setValue("add_y_ions", "false");
    theo_gen_al.setParameters(theo_gen_settings_al);
    theo_gen_al.getSpectrum(ms_spec_2_alabama, AASequence::fromString("ALABAMA"), 2, 2);
    for(Peak1D& peak : ms_spec_2_alabama) peak.setMZ(Math::ppmToMass(peak.getMZ(), 5.0) + peak.getMZ());
    
    // empty MSSpectrum
    MSSpectrum ms_spec_empty;

    // put valid data in exp
    exp.setSpectra({ms_spec_empty, ms_spec_2_alabama, ms_spec_2_himalaya});
    
    // map the MSExperiment
    QCBase::SpectraMap spectra_map(exp);

    //--------------------------------------------------------------------
    // test with valid input - default parameter
    //--------------------------------------------------------------------
    frag_ma_err.compute(fmap, exp, spectra_map);
    std::vector<FragmentMassError::FMEStatistics> result = frag_ma_err.getResults();

    TEST_REAL_SIMILAR(result[0].average_ppm, 5.0)
    TEST_REAL_SIMILAR(result[0].variance_ppm, 0.0)  // offset is constant, i.e. no variance

    //--------------------------------------------------------------------
    // test with valid input - ToleranceUnit PPM
    //--------------------------------------------------------------------

    FragmentMassError frag_ma_err_ppm;
    frag_ma_err_ppm.compute(fmap, exp, spectra_map, FragmentMassError::ToleranceUnit::PPM, 6);
    std::vector<FragmentMassError::FMEStatistics> result_ppm = frag_ma_err_ppm.getResults();

    TEST_REAL_SIMILAR(result_ppm[0].average_ppm, 5.0)
    TEST_REAL_SIMILAR(result_ppm[0].variance_ppm, 0.0) // offset is constant, i.e. no variance

    //--------------------------------------------------------------------
    // test with valid input and flags
    //--------------------------------------------------------------------
    FragmentMassError frag_ma_err_flag_da;
    frag_ma_err_flag_da.compute(fmap, exp, spectra_map, FragmentMassError::ToleranceUnit::DA, 1);
    std::vector<FragmentMassError::FMEStatistics> result_flag_da = frag_ma_err_flag_da.getResults();

    TEST_REAL_SIMILAR(result_flag_da[0].average_ppm, 5.0)
    TEST_REAL_SIMILAR(result_flag_da[0].variance_ppm, 0.0)  // offset is constant, i.e. no variance

    //--------------------------------------------------------------------
    // test with missing toleranceUnit and toleranceValue in featureMap
    //--------------------------------------------------------------------

    // featureMap with missing ProteinIdentifications
    FeatureMap fmap_auto;

    TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, frag_ma_err.compute(fmap_auto, exp, spectra_map, FragmentMassError::ToleranceUnit::AUTO), "There is no information about fragment mass tolerance given in the FeatureXML. Please choose a fragment_mass_unit")

    //--------------------------------------------------------------------
    // test with matching ms1 spectrum
    //--------------------------------------------------------------------

    // fmap with PeptideIdentification with RT matching to a MS1 Spectrum
    fmap.setUnassignedPeptideIdentifications({createPeptideIdentification("XTandem::3")});

    // set MS1 Spectrum to exp
    exp.setSpectra({createMSSpectrum(1, 5, "XTandem::3")});
    spectra_map.calculateMap(exp);

    TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, frag_ma_err.compute(fmap, exp, spectra_map), "The matching spectrum of the mzML is not an MS2 Spectrum.")

    //--------------------------------------------------------------------
    // test with no given fragmentation method
    //--------------------------------------------------------------------

    // put PeptideIdentification with RT matching to MSSpectrum without given fragmentation method to fmap
    fmap.setUnassignedPeptideIdentifications({createPeptideIdentification("XTandem::4")});

    // create MSExperiment with no given fragmentation method
    MSSpectrum ms_spec_2_no_precursor;
    ms_spec_2_no_precursor.setRT(8);
    ms_spec_2_no_precursor.setMSLevel(2);
    ms_spec_2_no_precursor.setNativeID("XTandem::4");
    exp.setSpectra({ms_spec_2_no_precursor});
    spectra_map.calculateMap(exp);

    TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, frag_ma_err.compute(fmap, exp, spectra_map), "No fragmentation method given.")

    //--------------------------------------------------------------------
    // test with fragmentation method SORI, which is not supported
    //--------------------------------------------------------------------

    // put PeptideIdentification with RT matching to MSSpectrum with fragmentation method SORI to fmap
    FeatureMap fmap_sori;
    fmap_sori.setProteinIdentifications( {protId} );
    fmap_sori.setUnassignedPeptideIdentifications({createPeptideIdentification("XTandem::5")});

    // MSExperiment with fragmentation method SORI (not supported)
    exp.setSpectra({createMSSpectrum(2, 7, "XTandem::5", Precursor::ActivationMethod::SORI)});
    spectra_map.calculateMap(exp);

    TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, frag_ma_err.compute(fmap_sori, exp, spectra_map), "Fragmentation method is not supported.")

    //--------------------------------------------------------------------
    // test if spectrum has no peaks
    //--------------------------------------------------------------------

    // put PeptideIdentification with RT matching to MSSpectrum with no peaks to fmap
    fmap.setUnassignedPeptideIdentifications({createPeptideIdentification("XTandem::6")});

    // MSExperiment without peaks
    exp.setSpectra({createMSSpectrum(2, 4, "XTandem::6")});
    spectra_map.calculateMap(exp);

    FragmentMassError frag_ma_err_excp;
    frag_ma_err_excp.compute(fmap, exp, spectra_map);
    std::vector<FragmentMassError::FMEStatistics> result_excp;
    result_excp = frag_ma_err_excp.getResults();

    TEST_REAL_SIMILAR(result_excp[0].average_ppm, 0)
    TEST_REAL_SIMILAR(result_excp[0].variance_ppm, 0)
  }
  END_SECTION


  START_SECTION(const String& getName() const override)
  {
    TEST_EQUAL(frag_ma_err.getName(), "FragmentMassError")
  }
  END_SECTION


  START_SECTION(QCBase::Status requires() const override)
  {
    QCBase::Status stat = QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
    TEST_EQUAL(stat, frag_ma_err.requires())
  }
  END_SECTION

END_TEST


