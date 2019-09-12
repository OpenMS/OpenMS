// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
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
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#include <OpenMS/QC/FragmentMassError.h>

#include <assert.h>
#include <string>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

namespace OpenMS
{
  const std::string FragmentMassError::names_of_toleranceUnit[] = {"ppm", "da", "auto"};

  template <typename MIV>
  void twoSpecErrors (MIV& mi, std::vector<double>& ppms, std::vector<double>& dalton, double& accumulator_ppm, UInt32& counter_ppm)
  {
    while (mi != mi.end())
    {
      // difference between peaks
      auto dalt_diff = mi->getMZ() - mi.ref().getMZ();
      auto ppm_diff = Math::getPPM(mi->getMZ(), mi.ref().getMZ());

      ppms.push_back(ppm_diff);
      dalton.push_back(dalt_diff);

      // for statistics
      accumulator_ppm += ppm_diff;
      ++counter_ppm;
      ++mi;
    }
  }


  PeakSpectrum getTheoSpec_(const Precursor::ActivationMethod& fm, const AASequence& seq, const int charge)
  {
    // initialize a TheoreticalSpectrumGenerator
    TheoreticalSpectrumGenerator theo_gen;

    // get current parameters (default)
    // default with b and y ions
    Param theo_gen_settings = theo_gen.getParameters();

    if (fm == Precursor::ActivationMethod::CID || fm == Precursor::ActivationMethod::HCID)
    {
      theo_gen_settings.setValue("add_b_ions", "true");
      theo_gen_settings.setValue("add_y_ions", "true");
    }
    else if (fm == Precursor::ActivationMethod::ECD || fm == Precursor::ActivationMethod::ETD)
    {
      theo_gen_settings.setValue("add_c_ions", "true");
      theo_gen_settings.setValue("add_z_ions", "true");
      theo_gen_settings.setValue("add_b_ions", "false");
      theo_gen_settings.setValue("add_y_ions", "false");
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Fragmentation method is not supported.");
    }

    // set changed parameters
    theo_gen.setParameters(theo_gen_settings);

    // generate b/y or c/z-ion spectrum of peptide seq
    PeakSpectrum theo_spectrum;
    theo_gen.getSpectrum(theo_spectrum, seq, 1, charge <= 2 ? 1 : 2);

    return theo_spectrum;
  }


  void FragmentMassError::compute(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit, double tolerance)
  {
    FMEStatistics result;

    // accumulates ppm errors over all first PeptideHits
    double accumulator_ppm{};

    // counts number of ppm errors
    UInt32 counter_ppm{};

    //---------------------------------------------------------------------
    // Prepare MSExperiment
    //---------------------------------------------------------------------

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 6, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    //-------------------------------------------------------------------
    // find tolerance unit and value
    //------------------------------------------------------------------
    if (tolerance_unit == ToleranceUnit::AUTO)
    {
      if (fmap.getProteinIdentifications().empty() )
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There is no information about fragment mass tolerance given in the FeatureXML. Please choose a fragment_mass_unit");
      }
      tolerance_unit = fmap.getProteinIdentifications()[0].getSearchParameters().fragment_mass_tolerance_ppm ? ToleranceUnit::PPM : ToleranceUnit::DA;
      tolerance = fmap.getProteinIdentifications()[0].getSearchParameters().fragment_mass_tolerance;
    }

    // computes the FragmentMassError
    auto lamCompPPM = [&exp, &map_to_spectrum, tolerance, tolerance_unit, &accumulator_ppm, &counter_ppm, &window_mower_filter](PeptideIdentification& pep_id)
    {
      if (pep_id.getHits().empty())
      {
        OPENMS_LOG_WARN << "PeptideHits of PeptideIdentification with RT: " << pep_id.getRT() << " and MZ: " << pep_id.getMZ() << " is empty.";
        return;
      }

      //---------------------------------------------------------------------
      // FIND DATA FOR THEORETICAL SPECTRUM
      //---------------------------------------------------------------------

      // sequence
      const AASequence& seq = pep_id.getHits()[0].getSequence();

      // charge: re-calulated from masses since much more robust this way (PepID annotation of pep_id.getHits()[0].getCharge() could be wrong)
      Int charge = static_cast<Int>(round(seq.getMonoWeight() / pep_id.getMZ()));

      //-----------------------------------------------------------------------
      // GET EXPERIMENTAL SPECTRUM MATCHING TO PEPTIDEIDENTIFICTION
      //-----------------------------------------------------------------------

      if (!pep_id.metaValueExists("spectrum_reference"))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectrum reference annotated at peptide identifiction!");
      }
      const MSSpectrum& exp_spectrum = exp[map_to_spectrum.at(pep_id.getMetaValue("spectrum_reference").toString())];

      if (exp_spectrum.getMSLevel() != 2)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The matching spectrum of the mzML is not an MS2 Spectrum.");
      }
      if (exp_spectrum.getPrecursors().empty() || exp_spectrum.getPrecursors()[0].getActivationMethods().empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No fragmentation method given.");
      }

      //---------------------------------------------------------------------
      // CREATE THEORETICAL SPECTRUM
      //---------------------------------------------------------------------
      PeakSpectrum theo_spectrum = getTheoSpec_((*exp_spectrum.getPrecursors()[0].getActivationMethods().begin()), seq, charge);

      //-----------------------------------------------------------------------
      // COMPARE THEORETICAL AND EXPERIMENTAL SPECTRUM
      //-----------------------------------------------------------------------
      if (exp_spectrum.empty() || theo_spectrum.empty())
      {
        OPENMS_LOG_WARN << "The spectrum with RT: " + String(exp_spectrum.getRT()) + " is empty." << "\n";
        return;
      }

      auto exp_spectrum_filtered(exp_spectrum);
      window_mower_filter.filterPeakSpectrum(exp_spectrum_filtered);

      // stores ppms for one spectrum
      DoubleList ppms{};
      DoubleList dalton{};

      // iterator, finds nearest peak of a target container to a given peak in a reference container
      if (tolerance_unit == ToleranceUnit::DA)
      {
        using MIV = MatchedIterator<MSSpectrum, DaTrait, true>;
        MIV mi(theo_spectrum, exp_spectrum_filtered, tolerance);
        twoSpecErrors(mi, ppms, dalton, accumulator_ppm, counter_ppm);
      }
      else
      {
        using MIV = MatchedIterator<MSSpectrum, PpmTrait, true>;
        MIV mi(theo_spectrum, exp_spectrum_filtered, tolerance);
        twoSpecErrors(mi, ppms, dalton, accumulator_ppm, counter_ppm);
      }

      //-----------------------------------------------------------------------
      // WRITE PPM ERROR IN PEPTIDEHIT
      //-----------------------------------------------------------------------
      pep_id.getHits()[0].setMetaValue("fragment_mass_error_ppm", ppms);
      pep_id.getHits()[0].setMetaValue("fragment_mass_error_da", dalton);
    };

    auto lamVar = [&result](const PeptideIdentification& pep_id)
    {
      if (pep_id.getHits().empty())
      {
        OPENMS_LOG_WARN << "There is a Peptideidentification(RT: " << pep_id.getRT() << ", MZ: " << pep_id.getMZ() <<  ") without PeptideHits. " << "\n";
        return;
      }
      for (const auto& ppm : (pep_id.getHits()[0].getMetaValue("fragment_mass_error_ppm")).toDoubleList())
      {
        result.variance_ppm += pow((ppm - result.average_ppm), 2);
      }
    };

    // computation of ppms
    QCBase::iterateFeatureMap(fmap, lamCompPPM);
    // if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (counter_ppm == 0)
    {
      results_.push_back(result);
      return;
    }

    // computes average
    result.average_ppm = accumulator_ppm / counter_ppm;

    // computes variance
    QCBase::iterateFeatureMap(fmap, lamVar);

    result.variance_ppm = result.variance_ppm / counter_ppm;

    results_.push_back(result);

  }

  const String& FragmentMassError::getName() const
  {
    return name_;
  }

  const std::vector<FragmentMassError::FMEStatistics>& FragmentMassError::getResults() const
  {
    return results_;
  }


  QCBase::Status FragmentMassError::requires() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS



