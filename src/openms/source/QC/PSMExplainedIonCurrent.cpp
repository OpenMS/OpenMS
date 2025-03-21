// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Waschischeck$
// $Authors: Tom Waschischeck$
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/QC/PSMExplainedIonCurrent.h>
#include <cfloat>
#include <numeric>

namespace OpenMS
{
  template<typename MIV>
  double sumOfMatchedIntensities(MIV& mi)
  {
    double sum = 0;
    while (mi != mi.end())
    {
      sum += mi->getIntensity();
      ++mi;
    }
    return sum;
  }

  double PSMExplainedIonCurrent::annotatePSMExplainedIonCurrent_(PeptideIdentification& pep_id, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, WindowMower& filter,
                                                                 PSMExplainedIonCurrent::ToleranceUnit tolerance_unit, double tolerance)
  {
    if (pep_id.getHits().empty())
    {
      OPENMS_LOG_DEBUG << "PeptideHits of PeptideIdentification with RT: " << pep_id.getRT() << " and MZ: " << pep_id.getMZ() << " is empty.";
      return DBL_MAX;
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
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The matching spectrum of the mzML is not an MS2 Spectrum.");
    }
    Precursor::ActivationMethod act_method;
    if (exp_spectrum.getPrecursors().empty() || exp_spectrum.getPrecursors()[0].getActivationMethods().empty())
    {
      OPENMS_LOG_DEBUG << "No MS2 activation method provided. Using CID as fallback to compute fragment mass errors." << std::endl;
      act_method = Precursor::ActivationMethod::CID;
    }
    else
    {
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
      return DBL_MAX;
    }

    // filter the spectrum
    PeakSpectrum filtered_spec(exp_spectrum);
    filter.filterPeakSpectrum(filtered_spec);

    double sum_of_intensities = 0;
    for (const auto& peak : filtered_spec)
    {
      sum_of_intensities += peak.getIntensity();
    }

    if (sum_of_intensities <= 0)
    {
      OPENMS_LOG_WARN << "The spectrum with RT: " + String(exp_spectrum.getRT()) + " has only peaks with intensity 0."
                      << "\n";
      return DBL_MAX;
    }

    double correctness = 0;
    // iterator, finds nearest peak of a target container to a given peak in a reference container
    if (tolerance_unit == PSMExplainedIonCurrent::ToleranceUnit::DA)
    {
      using MIV = MatchedIterator<MSSpectrum, DaTrait, true>;
      MIV mi(theo_spectrum, filtered_spec, tolerance);
      correctness = sumOfMatchedIntensities(mi) / sum_of_intensities;
    }
    else
    {
      using MIV = MatchedIterator<MSSpectrum, PpmTrait, true>;
      MIV mi(theo_spectrum, filtered_spec, tolerance);
      correctness = sumOfMatchedIntensities(mi) / sum_of_intensities;
    }

    pep_id.getHits()[0].setMetaValue(Constants::UserParam::PSM_EXPLAINED_ION_CURRENT_USERPARAM, correctness);

    return correctness;
  }

  void PSMExplainedIonCurrent::compute(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit, double tolerance)
  {
    Statistics result;

    bool has_pepIDs = QCBase::hasPepID(fmap);
    // if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (!has_pepIDs)
    {
      results_.push_back(result);
      return;
    }

    //---------------------------------------------------------------------
    // Prepare Spectrum Filter
    //---------------------------------------------------------------------

    WindowMower wm_filter;
    Param filter_param = wm_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 6, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    wm_filter.setParameters(filter_param);

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

    std::vector<double> correctnesses;

    std::function<void(PeptideIdentification&)> fCorrectness = [&exp, &map_to_spectrum, &correctnesses, &wm_filter, tolerance, tolerance_unit](PeptideIdentification& pep_id) {
      double correctness = annotatePSMExplainedIonCurrent_(pep_id, exp, map_to_spectrum, wm_filter, tolerance_unit, tolerance);
      if (correctness != DBL_MAX)
      {
        correctnesses.push_back(correctness);
      }
    };

    fmap.applyFunctionOnPeptideIDs(fCorrectness);

    if (correctnesses.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Couldn't calculate PSM correctness for any spectra! Check log for more information.");
    }

    result.average_correctness = Math::mean(correctnesses.begin(), correctnesses.end());
    result.variance_correctness = Math::variance(correctnesses.begin(), correctnesses.end(), result.average_correctness);

    results_.push_back(result);
  }

  void PSMExplainedIonCurrent::compute(std::vector<PeptideIdentification>& pep_ids, const ProteinIdentification::SearchParameters& search_params, const MSExperiment& exp,
                                       const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit, double tolerance)
  {
    Statistics result;

    if (pep_ids.empty())
    {
      results_.push_back(result);
      return;
    }

    //---------------------------------------------------------------------
    // Prepare Spectrum Filter
    //---------------------------------------------------------------------

    WindowMower wm_filter;
    Param filter_param = wm_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 6, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    wm_filter.setParameters(filter_param);

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

    std::vector<double> correctnesses;

    for (auto& pep_id : pep_ids)
    {
      double correctness = annotatePSMExplainedIonCurrent_(pep_id, exp, map_to_spectrum, wm_filter, tolerance_unit, tolerance);
      if (correctness != DBL_MAX)
      {
        correctnesses.push_back(correctness);
      }
    }

    if (correctnesses.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Couldn't calculate PSM correctness for any spectra! Check log for more information.");
    }

    result.average_correctness = Math::mean(correctnesses.begin(), correctnesses.end());
    result.variance_correctness = Math::variance(correctnesses.begin(), correctnesses.end(), result.average_correctness);

    results_.push_back(result);
  }

  const String& PSMExplainedIonCurrent::getName() const
  {
    static const String& name = "PSMExplainedIonCurrent";
    return name;
  }

  const std::vector<PSMExplainedIonCurrent::Statistics>& PSMExplainedIonCurrent::getResults() const
  {
    return results_;
  }


  QCBase::Status PSMExplainedIonCurrent::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }

}; // namespace OpenMS
