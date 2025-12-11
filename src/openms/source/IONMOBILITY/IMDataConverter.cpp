// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/IMDataConverter.h>

#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/KERNEL/MSExperiment.h>


#include <cstddef>
#include <cmath>
#include <map>
#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>

namespace OpenMS
{
  std::vector<std::pair<double, MSExperiment>> IMDataConverter::splitByFAIMSCV(PeakMap&& exp)
  {
    std::vector<std::pair<double, MSExperiment>> result;

    std::set<double> CVs = FAIMSHelper::getCompensationVoltages(exp);

    if (CVs.empty())
    {
      OPENMS_LOG_INFO << "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN." << std::endl;
      // Represent "no FAIMS" as a single-element group that carries the original experiment.
      result.emplace_back(std::numeric_limits<double>::quiet_NaN(), std::move(exp));
      return result;
    }

    // use a map internally to accumulate per-CV experiments in ascending CV order
    std::map<double, MSExperiment> split_peakmap;
    for (double cv : CVs)
    {
      MSExperiment pm;
      pm.getExperimentalSettings() = exp.getExperimentalSettings();
      split_peakmap.emplace(cv, std::move(pm));
    }

    // fill up the PeakMaps by moving spectra from the input PeakMap
    // Keep MS2 spectra which might not carry a FAIMS CV (typically they do) by assigning them to the last seen FAIMS CV (nearest previous in run order)
    double last_faims_cv = std::numeric_limits<double>::quiet_NaN();
    for (MSSpectrum& spec : exp)
    {
      if (spec.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
      {
        last_faims_cv = spec.getDriftTime();
        if (const auto map_it = split_peakmap.find(last_faims_cv); map_it == split_peakmap.end())
        {
          OPENMS_LOG_WARN << "Encountered spectrum with unexpected FAIMS CV (not in detected set): " << last_faims_cv << std::endl;
          continue;
        }
        else
        {
          map_it->second.addSpectrum(std::move(spec));
        }
        continue;
      }

      // No FAIMS CV on spectrum: if it's MS2+ and we have a prior FAIMS CV context, assign it to that bin
      if (!std::isnan(last_faims_cv) && spec.getMSLevel() > 1)
      {
        if (const auto map_it = split_peakmap.find(last_faims_cv); map_it != split_peakmap.end())
        {
          map_it->second.addSpectrum(std::move(spec));
          continue;
        }
      }

      // Otherwise skip with a warning
      OPENMS_LOG_WARN << "Skipping spectrum without FAIMS CV (no prior FAIMS CV context or unexpected layout)." << std::endl;
    }
    
    exp.clear(true);

    // move from map (sorted by CV) into vector of (CV, MSExperiment) pairs
    result.reserve(split_peakmap.size());
    for (auto& kv : split_peakmap)
    {
      result.emplace_back(kv.first, std::move(kv.second));
    }
    return result;
  }

  MSExperiment IMDataConverter::reshapeIMFrameToMany(MSSpectrum im_frame)
  {
    MSExperiment out;

    if (im_frame.empty())
    {// nothing to split (we do not even check for IM data, for robustness)
      return out;
    }

    // check if data is sorted by IM... if not, sort
    if (! im_frame.isSortedByIM())
    { // sorts the spectrum (and its binary data arrays) according to IM
      im_frame.sortByIonMobility();
    }

    // can throw if IM float data array is missing
    const auto [im_data_index, im_unit] = im_frame.getIMData();
    // Capture IM array by Ref, because .getIMData() is expensive to call for every peak!
    const auto& im_data = im_frame.getFloatDataArrays()[im_data_index];

    // copy meta data (RT, name, ...) without the raw data and without the IM array
    MSSpectrum prototype = im_frame;
    prototype.clear(false);

    // adds a new spectrum with drift time to `out`
    auto addSpectrum = [&out, im_unit = im_unit, &prototype](double drift_time_avg) {
      // keeps RT identical for all scans, since they are from the same IM-frame
      // keeps MSlevel
      out.addSpectrum(prototype);
      auto& spec = out.getSpectra().back();
      // copy drift-time unit from parent scan
      spec.setDriftTime(drift_time_avg);
      spec.setDriftTimeUnit(im_unit);
      return &spec;
    };

    MSSpectrum* last_spec{};
    // Separate spec for each IM value:
    OPENMS_PRECONDITION(std::is_sorted(im_data.begin(), im_data.end()), "we sorted it... what happened???");
    using IMV_t = MSSpectrum::FloatDataArray::value_type;
    IMV_t im_last = std::numeric_limits<IMV_t>::max();
    for (Size i = 0; i < im_data.size(); ++i)// is sorted now!
    {
      const IMV_t im = im_data[i];
      if (im != im_last)
      {
        im_last = im;
        last_spec = addSpectrum(im);
      }
      last_spec->push_back(im_frame[i]);// copy the m/z of the peak
    }
    out.sortSpectra(true);
    out.updateRanges();
    return out;
  }

  std::tuple<std::vector<MSExperiment>, Math::BinContainer> IMDataConverter::splitExperimentByIonMobility(MSExperiment&& in,
                                                                                                          UInt number_of_bins,
                                                                                                          double bin_extension_abs,
                                                                                                          double mz_binning_width,
                                                                                                          MZ_UNITS mz_binning_width_unit)
  {
    if (number_of_bins == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot split into 0 bins.", String(number_of_bins));
    }
    if (bin_extension_abs < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Overlap must not be negative.", String(bin_extension_abs));
    }
    std::vector<MSExperiment> results(number_of_bins);
    in.updateRanges();
    // find the IM range
    const auto range_IM = RangeMobility(in.spectrumRanges());
    if (range_IM.getSpan() / number_of_bins < bin_extension_abs * 2)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Bin size (") + String(range_IM.getSpan() / number_of_bins) + ") is smaller than the overlap.", String(bin_extension_abs*2));
    }

    // compute the bins
    const auto bins = Math::createBins(range_IM.getMin(), range_IM.getMax(), number_of_bins, bin_extension_abs);

    // results for each IM-frame: all spectra per bin, to get merged
    MSExperiment binned_spectra;

    SpectraMerger merger;
    auto p = merger.getParameters();
    const auto ms_levels = in.getMSLevels();
    p.setValue("block_method:ms_levels", IntList(ms_levels.begin(), ms_levels.end())); // merge all MS levels
    p.setValue("mz_binning_width", mz_binning_width);
    p.setValue("mz_binning_width_unit", String(MZ_UNIT_NAMES[(int)mz_binning_width_unit]));
    p.setValue("block_method:rt_block_size", INT_MAX);
    p.setValue("block_method:rt_max_length", 10e10);


    for (auto& frame : in)
    {
      // For data without ion mobility, simply append the result (only
      // collapse for scans that actually have a float data array).
      if (! frame.containsIMData())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Spectrum does not contain 'wide' IM data.", frame.getNativeID());
      }
      
      MSExperiment frame_melt = IMDataConverter::reshapeIMFrameToMany(std::move(frame));
      for (size_t i = 0; i < bins.size(); ++i)
      {
        binned_spectra.clear(false);
        // check if spectrum goes into this bin
        for (auto&& spec : frame_melt)
        {
          if (bins[i].contains(spec.getDriftTime()))
          { // spectrum goes into this bin
            binned_spectra.addSpectrum(std::move(spec));
          }
        }
        // collapse spectra in this bin
        if (!binned_spectra.empty())
        {
          merger.setParameters(p);
          merger.mergeSpectraBlockWise(binned_spectra);
          assert(binned_spectra.size() == 1);
          results[i].addSpectrum(std::move(binned_spectra.getSpectra().back()));
          results[i].getSpectra().back().setDriftTime(bins[i].center());
        }
      }
    }
    for (auto& result : results)
    {
      result.ExperimentalSettings::operator=(in);
      result.updateRanges();
    }
    return {std::move(results), std::move(bins)};
  }

  void annotateAsIM(OpenMS::DataArrays::FloatDataArray& fda, const DriftTimeUnit unit)
  {
    const auto& cv = ControlledVocabulary::getPSIMSCV();
    const ControlledVocabulary::CVTerm* term;
    switch (unit)
    {
      case DriftTimeUnit::MILLISECOND:
        term = &cv.getTerm("MS:1002816");
        break;
      case DriftTimeUnit::VSSC:
         term = &cv.getTerm("MS:1003008");
        break;
      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unit cannot be converted into CV term.", toString(unit));
    }
    fda.setName(term->name);
  }

  
  /// private: Process a stack of drift time spectra
  void processDriftTimeStack(std::vector<const MSSpectrum*>& stack, MSExperiment& result)
  {
    if (stack.empty()) return;

    // copy meta data without the raw data and without the IM array
    MSSpectrum new_spec = *stack[0];
    new_spec.clear(false);

    // create new FDA
    OpenMS::DataArrays::FloatDataArray& fda = new_spec.getFloatDataArrays().emplace_back();
    IMDataConverter::setIMUnit(fda, new_spec.getDriftTimeUnit());
    for (const auto& s : stack)
    {
      new_spec.insert(new_spec.end(), s->begin(), s->end()); // append data
      fda.insert(fda.end(), s->size(), s->getDriftTime());   // create IM array
    }
    new_spec.setDriftTime(IMTypes::DRIFTTIME_NOT_SET);// drift time is now encoded in the FloatDataArray
    new_spec.setDriftTimeUnit(DriftTimeUnit::NONE);   // drift time is now encoded in the FloatDataArray
    
    // append to PeakMap
    result.getSpectra().push_back(std::move(new_spec));

    stack.clear();
  }



  MSExperiment IMDataConverter::reshapeIMFrameToSingle(const MSExperiment& exp)
  {
    MSExperiment result;

    if (exp.empty())
    {
      return result;
    }      

    std::vector<const MSSpectrum*> stack;
    double curr_rt = std::numeric_limits<double>::max();
    for (const auto& spec : exp)
    {
      // copy non-IM or already framed spectra
      // throws Exception if spec has mixed IM format
      if (IMTypes::determineIMFormat(spec) != IMFormat::MULTIPLE_SPECTRA)
      {
        processDriftTimeStack(stack, result); // clear current stack
        result.getSpectra().push_back(spec);
        continue;
      }
      
      // new stack starts here. Process all previous spectra
      if (spec.getRT() != curr_rt)
      {
        processDriftTimeStack(stack, result);
        curr_rt = spec.getRT();
      }
      stack.push_back(&spec);
    }
     
     processDriftTimeStack(stack, result);
     
     return result;
  }

  void IMDataConverter::setIMUnit(DataArrays::FloatDataArray& fda, const DriftTimeUnit unit)
  {
    const auto& cv = ControlledVocabulary::getPSIMSCV();
    switch (unit)
    {
      case DriftTimeUnit::MILLISECOND: 
        fda.setName(cv.getTerm("MS:1002816").name); // MS:1002816 ! mean ion mobility array
        return;
      case DriftTimeUnit::VSSC:
        fda.setName(cv.getTerm("MS:1003008").name); // MS:1003008 ! raw inverse reduced ion mobility array
        return;
      default:
        // invalid enum ...
        // There is no CV term which can be used to describe the FDA
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unit is not a valid IM unit for float data arrays", toString(unit));
    }
  }

  bool IMDataConverter::getIMUnit(const DataArrays::FloatDataArray& fda, DriftTimeUnit& unit)
  {
    const auto& cv = ControlledVocabulary::getPSIMSCV();
    if (fda.getName().hasPrefix(Constants::UserParam::ION_MOBILITY))
    { // fallback for non-standard IM arrays (as created by Mobi-DIK, or "Ion Mobility Centroid" from PeakPickerIM)
      if (fda.getName().hasSubstring("MS:1002815"))
      {
        unit = DriftTimeUnit::VSSC;
      }
      else
      {
        unit = DriftTimeUnit::MILLISECOND;
      }
      return true;
    }
    try
    {
      const auto& cv_term = cv.getTermByName(fda.getName()); // may throw if term is unknown

      if (cv.isChildOf(cv_term.id, "MS:1002893")) // is child of generic 'ion mobility array'?
      {
        if (cv_term.units.find("MS:1002814") != cv_term.units.end())
        { // MS:1002814 ! volt-second per square centimeter
          unit = DriftTimeUnit::VSSC;
        }
        else if (cv_term.units.find("UO:0000028") != cv_term.units.end())
        { // UO:0000028 ! millisecond
          unit = DriftTimeUnit::MILLISECOND;
        }
        else
        { // fallback
          OPENMS_LOG_WARN << "Warning: FloatDataArray for IonMobility data '" << cv_term.id << " " << cv_term.name << "' does not contain proper units!" << std::endl;
          unit = DriftTimeUnit::NONE;
        }
        return true;
      }
    }
    catch (...)
    {
    }
    return false;
  }

  double IMDataConverter::convertVSSCToCCS(double IM, double mz, int charge)
  {
    // Constants for the Mason-Schamp equation (confirmed values with alpha and MaxQuant)
    constexpr double bruker_CCS_coef = 1059.62245; // constant coefficient for Bruker in the Mason-Schamp equation
    constexpr double IM_N2_gas_mass = 28.0; // N2 gas mass (like in alpha code)

    const double mass = mz * charge;
    const double reduced_mass = mass * IM_N2_gas_mass / (mass + IM_N2_gas_mass);
    return IM * charge * bruker_CCS_coef / std::sqrt(reduced_mass); // Mason-Schamp equation
  }

  void IMDataConverter::convertVSSCToCCS(MSExperiment& spectra)
  {
    OPENMS_LOG_INFO << "Converting 1/k0 to CCS values." << std::endl;

    for (auto& s : spectra)
    {
      if (s.getPrecursors().empty())
      {
        continue; // skip spectra without precursor info
      }
      const double mz = s.getPrecursors()[0].getMZ();
      const int charge = s.getPrecursors()[0].getCharge();
      if (charge == 0)
      {
        continue; // skip spectra with unknown charge
      }

      // Convert spectrum-level drift time if valid
      const double IM = s.getDriftTime();
      if (IM > 0 && !std::isnan(IM))
      {
        s.setDriftTime(convertVSSCToCCS(IM, mz, charge));
      }

      // Also convert IM values in float data arrays (per-peak ion mobility)
      for (auto& fda : s.getFloatDataArrays())
      {
        DriftTimeUnit unit;
        if (getIMUnit(fda, unit) && unit == DriftTimeUnit::VSSC)
        {
          for (auto& im_value : fda)
          {
            if (im_value > 0 && !std::isnan(im_value))
            {
              im_value = static_cast<float>(convertVSSCToCCS(im_value, mz, charge));
            }
          }
        }
      }
    }
  }

}  //end namespace OpenMS
