// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/QC/FragmentMassError.h>
#include <cassert>
#include <string>

namespace OpenMS
{
  // Using matched iterator for aligned spectra calculate mz errors
  template<typename MIV>
  void twoSpecErrors(MIV& mi, std::vector<double>& ppms, std::vector<double>& dalton, double& accumulator_ppm, UInt32& counter_ppm)
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

  void FragmentMassError::calculateFME_(PeptideIdentification& pep_id, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, bool& print_warning, double tolerance,
                                        FragmentMassError::ToleranceUnit tolerance_unit, double& accumulator_ppm, UInt32& counter_ppm, WindowMower& window_mower_filter)
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

    // charge: re-calculated from masses since much more robust this way (PepID annotation of pep_id.getHits()[0].getCharge() could be wrong)
    Int charge = static_cast<Int>(round(seq.getMonoWeight() / pep_id.getMZ()));

    //-----------------------------------------------------------------------
    // GET EXPERIMENTAL SPECTRUM MATCHING TO PEPTIDEIDENTIFICATION
    //-----------------------------------------------------------------------

    if (!pep_id.metaValueExists("spectrum_reference"))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectrum reference annotated at peptide identifiction!");
    }
    const MSSpectrum& exp_spectrum = exp[map_to_spectrum.at(pep_id.getSpectrumReference())];

    if (exp_spectrum.getMSLevel() != 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Spectrum with wrong MS level provided. MS2 expected.");
    }
    Precursor::ActivationMethod act_method;
    if (exp_spectrum.getPrecursors().empty())
    {
      if (print_warning)
      {
        OPENMS_LOG_WARN << "No MS2 activation method provided. Using CID as fallback to compute fragment mass errors." << std::endl;
      }
      print_warning = false; // only print it once
      act_method = Precursor::ActivationMethod::CID;
    }
    else
    {
      if (exp_spectrum.getPrecursors()[0].getActivationMethods().empty())
      {
        if (print_warning)
        {
          OPENMS_LOG_WARN << "No MS2 activation method provided. Using CID as fallback to compute fragment mass errors." << std::endl;
        }
        print_warning = false; // only print it once
        act_method = Precursor::ActivationMethod::CID;
      }
      act_method = *exp_spectrum.getPrecursors()[0].getActivationMethods().begin();
    }

    //---------------------------------------------------------------------
    // CREATE THEORETICAL SPECTRUM
    //---------------------------------------------------------------------
    PeakSpectrum theo_spectrum = TheoreticalSpectrumGenerator::generateSpectrum(act_method, seq, charge);

    //-----------------------------------------------------------------------
    // COMPARE THEORETICAL AND EXPERIMENTAL SPECTRUM
    //-----------------------------------------------------------------------
    if (exp_spectrum.empty() || theo_spectrum.empty())
    {
      OPENMS_LOG_WARN << "The spectrum with RT: " + String(exp_spectrum.getRT()) + " is empty."
                      << "\n";
      return;
    }

    auto exp_spectrum_filtered(exp_spectrum);
    window_mower_filter.filterPeakSpectrum(exp_spectrum_filtered);

    // stores ppms for one spectrum
    DoubleList ppms {};
    DoubleList dalton {};

    // iterator, finds nearest peak of a target container to a given peak in a reference container
    if (tolerance_unit == FragmentMassError::ToleranceUnit::DA)
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
    pep_id.getHits()[0].setMetaValue(Constants::UserParam::FRAGMENT_ERROR_PPM_USERPARAM, ppms);
    pep_id.getHits()[0].setMetaValue(Constants::UserParam::FRAGMENT_ERROR_DA_USERPARAM, dalton);
    if (ppms.size() > 1)
    {
      pep_id.getHits()[0].setMetaValue(Constants::UserParam::FRAGMENT_ERROR_PPM_USERPARAM + "_variance", Math::variance(ppms.begin(), ppms.end()));
    }
    if (dalton.size() > 1)
    {
      pep_id.getHits()[0].setMetaValue(Constants::UserParam::FRAGMENT_ERROR_DA_USERPARAM + "_variance", Math::variance(dalton.begin(), dalton.end()));
    }
  }

  void FragmentMassError::calculateVariance_(FragmentMassError::Statistics& result, const PeptideIdentification& pep_id, const UInt num_ppm)
  {
    if (pep_id.getHits().empty())
    {
      OPENMS_LOG_WARN << "There is a Peptideidentification(RT: " << pep_id.getRT() << ", MZ: " << pep_id.getMZ() << ") without PeptideHits. "
                      << "\n";
      return;
    }
    for (const auto& ppm : (pep_id.getHits()[0].getMetaValue("fragment_mass_error_ppm")).toDoubleList())
    {
      double tmp = ppm - result.average_ppm;
      result.variance_ppm += (tmp * tmp / num_ppm);
    }
  }

  void FragmentMassError::compute(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit, double tolerance)
  {
    Statistics result;

    bool has_pepIDs = QCBase::hasPepID(fmap);
    // if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (!has_pepIDs)
    {
      results_.push_back(result);
      return;
    }
    // accumulates ppm errors over all first PeptideHits
    double accumulator_ppm {};

    // counts number of ppm errors
    UInt32 counter_ppm {};

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
      if (fmap.getProteinIdentifications().empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No information about fragment mass tolerance given in the FeatureMap. Please choose a fragment_mass_unit and tolerance manually.");
      }
      tolerance_unit = fmap.getProteinIdentifications()[0].getSearchParameters().fragment_mass_tolerance_ppm ? ToleranceUnit::PPM : ToleranceUnit::DA;
      tolerance = fmap.getProteinIdentifications()[0].getSearchParameters().fragment_mass_tolerance;
      if (tolerance <= 0.0)
      { // some engines, e.g. MSGF+ have no fragment tolerance parameter. It will be 0.0.
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No information about fragment mass tolerance given in the FeatureMap. Please choose a fragment_mass_unit and tolerance manually.");
      }
    }

    bool print_warning {false};

    // computes the FragmentMassError
    std::function<void(PeptideIdentification&)> fCompPPM = [&exp, &map_to_spectrum, &print_warning, tolerance, tolerance_unit, &accumulator_ppm, &counter_ppm,
                                                            &window_mower_filter](PeptideIdentification& pep_id) {
      calculateFME_(pep_id, exp, map_to_spectrum, print_warning, tolerance, tolerance_unit, accumulator_ppm, counter_ppm, window_mower_filter);
    };

    auto fVar = [&result, &counter_ppm](const PeptideIdentification& pep_id) { calculateVariance_(result, pep_id, counter_ppm); };

    // computation of ppms
    fmap.applyFunctionOnPeptideIDs(fCompPPM);
    // if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (counter_ppm == 0)
    {
      results_.push_back(result);
      return;
    }

    // computes average
    result.average_ppm = accumulator_ppm / counter_ppm;

    // computes variance
    fmap.applyFunctionOnPeptideIDs(fVar);

    results_.push_back(result);
  }

  void FragmentMassError::compute(std::vector<PeptideIdentification>& pep_ids, const ProteinIdentification::SearchParameters& search_params, const MSExperiment& exp,
                                  const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit, double tolerance)
  {
    Statistics result;

    if (pep_ids.empty())
    {
      results_.push_back(result);
      return;
    }
    // accumulates ppm errors over all first PeptideHits
    double accumulator_ppm {};

    // counts number of ppm errors
    UInt32 counter_ppm {};

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
      tolerance_unit = search_params.fragment_mass_tolerance_ppm ? ToleranceUnit::PPM : ToleranceUnit::DA;
      tolerance = search_params.fragment_mass_tolerance;
      if (tolerance <= 0.0)
      { // some engines, e.g. MSGF+ have no fragment tolerance parameter. It will be 0.0.
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No information about fragment mass tolerance given. Please choose a fragment_mass_unit and tolerance manually.");
      }
    }

    bool print_warning {false};

    // computation of ppms
    // computes the FragmentMassError
    for (auto& pep_id : pep_ids)
    {
      calculateFME_(pep_id, exp, map_to_spectrum, print_warning, tolerance, tolerance_unit, accumulator_ppm, counter_ppm, window_mower_filter);

      // if there are no matching peaks, the counter is zero and it is not possible to find ppms
      if (counter_ppm == 0)
      {
        results_.push_back(result);
        return;
      }
      // computes average
      result.average_ppm = accumulator_ppm / counter_ppm;

      calculateVariance_(result, pep_id, counter_ppm);
    }

    results_.push_back(result);
  }

  const String& FragmentMassError::getName() const
  {
    static const String& name = "FragmentMassError";
    return name;
  }

  const std::vector<FragmentMassError::Statistics>& FragmentMassError::getResults() const
  {
    return results_;
  }


  QCBase::Status FragmentMassError::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS
