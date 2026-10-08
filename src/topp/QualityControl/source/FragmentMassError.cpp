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
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include "FragmentMassError.h"
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

  void FragmentMassError::calculateFME_(QCBase::AnnotatedIdentification& id, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, bool& print_warning, double tolerance,
                                        FragmentMassError::ToleranceUnit tolerance_unit, double& accumulator_ppm, UInt32& counter_ppm, WindowMower& window_mower_filter)
  {
    if (id.top == nullptr)
    {
      OPENMS_LOG_WARN << "PeptideHits of PeptideIdentification with RT: " << id.rt << " and MZ: " << id.mz << " is empty.";
      return;
    }

    //---------------------------------------------------------------------
    // FIND DATA FOR THEORETICAL SPECTRUM
    //---------------------------------------------------------------------

    // sequence
    const AASequence& seq = id.top->getSequence();

    // charge: re-calculated from masses since much more robust this way (PepID annotation of the top hit's charge could be wrong)
    Int charge = static_cast<Int>(round(seq.getMonoWeight() / id.mz));

    //-----------------------------------------------------------------------
    // GET EXPERIMENTAL SPECTRUM MATCHING TO PEPTIDEIDENTIFICATION
    //-----------------------------------------------------------------------

    if (id.spectrum_reference.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectrum reference annotated at peptide identifiction!");
    }
    const MSSpectrum& exp_spectrum = exp[map_to_spectrum.at(id.spectrum_reference)];

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
      OPENMS_LOG_WARN << "The spectrum with RT: " + StringUtils::toStr(exp_spectrum.getRT()) + " is empty."
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
    id.top->setMetaValue(Constants::UserParam::FRAGMENT_ERROR_PPM_USERPARAM, ppms);
    id.top->setMetaValue(Constants::UserParam::FRAGMENT_ERROR_DA_USERPARAM, dalton);
    if (ppms.size() > 1)
    {
      id.top->setMetaValue(Constants::UserParam::FRAGMENT_ERROR_PPM_USERPARAM + "_variance", Math::variance(ppms.begin(), ppms.end()));
    }
    if (dalton.size() > 1)
    {
      id.top->setMetaValue(Constants::UserParam::FRAGMENT_ERROR_DA_USERPARAM + "_variance", Math::variance(dalton.begin(), dalton.end()));
    }
  }

  void FragmentMassError::calculateVariance_(FragmentMassError::Statistics& result, const QCBase::AnnotatedIdentification& id, const UInt num_ppm)
  {
    if (id.top == nullptr)
    {
      OPENMS_LOG_WARN << "There is a Peptideidentification(RT: " << id.rt << ", MZ: " << id.mz << ") without PeptideHits. "
                      << "\n";
      return;
    }
    for (const auto& ppm : (id.top->getMetaValue("fragment_mass_error_ppm")).toDoubleList())
    {
      double tmp = ppm - result.average_ppm;
      result.variance_ppm += (tmp * tmp / num_ppm);
    }
  }

  void FragmentMassError::compute(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit, double tolerance)
  {
    IdentificationDataConverter::editAsIdentificationData(fmap, [&](FeatureMap& map) { computeNative_(map, exp, map_to_spectrum, tolerance_unit, tolerance); });
  }

  void FragmentMassError::computeNative_(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit, double tolerance)
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
      const auto* search = QCBase::searchParameters(fmap.getIdentificationData());
      if (search == nullptr)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No information about fragment mass tolerance given in the FeatureMap. Please choose a fragment_mass_unit and tolerance manually.");
      }
      tolerance_unit = search->fragment_mass_tolerance_ppm ? ToleranceUnit::PPM : ToleranceUnit::DA;
      tolerance = search->fragment_mass_tolerance;
      if (tolerance <= 0.0)
      { // some engines, e.g. MSGF+ have no fragment tolerance parameter. It will be 0.0.
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No information about fragment mass tolerance given in the FeatureMap. Please choose a fragment_mass_unit and tolerance manually.");
      }
    }

    bool print_warning {false};

    // computation of ppms, of the identifications of the features and the unassigned ones
    QCBase::annotateIdentifications(fmap, [&](Feature*, std::vector<QCBase::AnnotatedIdentification>& identifications) {
      for (auto& id : identifications)
      {
        calculateFME_(id, exp, map_to_spectrum, print_warning, tolerance, tolerance_unit, accumulator_ppm, counter_ppm, window_mower_filter);
      }
    });
    // if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (counter_ppm == 0)
    {
      results_.push_back(result);
      return;
    }

    // computes average
    result.average_ppm = accumulator_ppm / counter_ppm;

    // computes variance
    QCBase::visitIdentifications(fmap, [&](const Feature*, const std::vector<QCBase::AnnotatedIdentification>& identifications) {
      for (const auto& id : identifications)
      {
        calculateVariance_(result, id, counter_ppm);
      }
    });

    results_.push_back(result);
  }

  void FragmentMassError::compute(PeptideIdentificationList& pep_ids, const ProteinIdentification::SearchParameters& search_params, const MSExperiment& exp,
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
    // Both overloads use the same two-pass computation so their average/variance agree:
    // variance must be accumulated against the FINAL average over all ppm errors, not a
    // moving/partial average inside the loop.
    // first pass: accumulate all ppm errors
    for (auto& pep_id : pep_ids)
    {
      auto id = QCBase::annotated(pep_id);
      calculateFME_(id, exp, map_to_spectrum, print_warning, tolerance, tolerance_unit, accumulator_ppm, counter_ppm, window_mower_filter);
    }

    // if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (counter_ppm == 0)
    {
      results_.push_back(result);
      return;
    }

    // computes average
    result.average_ppm = accumulator_ppm / counter_ppm;

    // computes variance (second pass: against the final average)
    for (auto& pep_id : pep_ids)
    {
      calculateVariance_(result, QCBase::annotated(pep_id), counter_ppm);
    }

    results_.push_back(result);
  }

  const std::string& FragmentMassError::getName() const
  {
    static const std::string& name = "FragmentMassError";
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
