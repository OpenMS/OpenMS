// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

using namespace std;

namespace OpenMS
{
  void SpectrumMetaDataLookup::getSpectrumMetaData(Size index,
                                                   SpectrumMetaData& meta) const
  {
    if (index >= n_spectra_)
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     index, n_spectra_);
    }
    meta = metadata_[index];
  }


  void SpectrumMetaDataLookup::getSpectrumMetaData(
    const MSSpectrum& spectrum, SpectrumMetaData& meta,
    const boost::regex& scan_regexp, const map<Size, double>& precursor_rts)
  {
    meta.native_id = spectrum.getNativeID();
    meta.rt = spectrum.getRT();
    meta.ms_level = spectrum.getMSLevel();
    if (!scan_regexp.empty())
    {
      meta.scan_number = extractScanNumber(meta.native_id, scan_regexp, true);
      if (meta.scan_number < 0)
      {
        OPENMS_LOG_ERROR << "Error: Could not extract scan number from spectrum native ID '" + meta.native_id + "' using regular expression '" + scan_regexp.str() + "'." << endl;
      }
    }
    if (!spectrum.getPrecursors().empty())
    {
      meta.precursor_mz = spectrum.getPrecursors()[0].getMZ();
      meta.precursor_charge = spectrum.getPrecursors()[0].getCharge();
      if (!precursor_rts.empty())
      {
        // precursor RT is RT of previous spectrum with lower MS level:
        map<Size, double>::const_iterator pos = 
          precursor_rts.find(meta.ms_level - 1);
        if (pos != precursor_rts.end()) // found
        {
          meta.precursor_rt = pos->second;
        }
        else
        {
          OPENMS_LOG_ERROR << "Error: Could not set precursor RT for spectrum with native ID '" + meta.native_id + "' - precursor spectrum not found." << endl;
        }
      }
    } 
  }


  void SpectrumMetaDataLookup::getSpectrumMetaData(const String& spectrum_ref,
                                                   SpectrumMetaData& meta,
                                                   MetaDataFlags flags) const
  {
    for (std::vector<boost::regex>::const_iterator it = 
           reference_formats.begin(); it != reference_formats.end(); ++it)
    {
      boost::smatch match;
      bool found = boost::regex_search(spectrum_ref, match, *it);
      if (found)
      {
        // first try to extract the requested meta data from the reference:
        if (((flags & MDF_RT) == MDF_RT) && match["RT"].matched)
        {
          String value = match["RT"].str();
          if (!value.empty())
          {
            meta.rt = value.toDouble();
            flags &= ~MDF_RT; // unset flag
          }
        }
        if (((flags & MDF_PRECURSORRT) == MDF_PRECURSORRT) && 
            match["PRECRT"].matched)
        {
          String value = match["PRECRT"].str();
          if (!value.empty())
          {
            meta.precursor_rt = value.toDouble();
            flags &= ~MDF_PRECURSORRT; // unset flag
          }
        }
        if (((flags & MDF_PRECURSORMZ) == MDF_PRECURSORMZ) && 
            match["MZ"].matched)
        {
          String value = match["MZ"].str();
          if (!value.empty())
          {
            meta.precursor_mz = value.toDouble();
            flags &= ~MDF_PRECURSORMZ; // unset flag
          }
        }
        if (((flags & MDF_PRECURSORCHARGE) == MDF_PRECURSORCHARGE) &&
            match["CHARGE"].matched)
        {
          String value = match["CHARGE"].str();
          if (!value.empty())
          {
            meta.precursor_charge = value.toDouble();
            flags &= ~MDF_PRECURSORCHARGE; // unset flag
          }
        }
        if (((flags & MDF_MSLEVEL) == MDF_MSLEVEL) && match["LEVEL"].matched)
        {
          String value = match["LEVEL"].str();
          if (!value.empty())
          {
            meta.ms_level = value.toInt();
            flags &= ~MDF_MSLEVEL; // unset flag
          }
        }
        if (((flags & MDF_SCANNUMBER) == MDF_SCANNUMBER) && 
            match["SCAN"].matched)
        {
          String value = match["SCAN"].str();
          if (!value.empty())
          {
            meta.scan_number = value.toInt();
            flags &= ~MDF_SCANNUMBER; // unset flag
          }
        }
        if (((flags & MDF_NATIVEID) == MDF_NATIVEID) && match["ID"].matched)
        {
          meta.native_id = match["ID"].str();
          if (!meta.native_id.empty())
          {
            flags &= ~MDF_NATIVEID; // unset flag
          }
        }
        if (flags) // not all requested values have been found -> look them up
        {
          Size index = findByRegExpMatch_(spectrum_ref, it->str(), match);
          meta = metadata_[index];
        }
        return; // use the first reference format that matches
      }
    }
  }


bool SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(vector<PeptideIdentification>& peptides,
                                                        const MSExperiment& exp)
{
    // Check if the experiment has spectra
    if (exp.getSpectra().empty())
    {
        OPENMS_LOG_INFO << "No spectra found in the experiment. Skipping RT annotation." << endl;
        return false;
    }

    SpectrumLookup lookup;
    lookup.readSpectra(exp.getSpectra());
    bool success = true;

    // Iterate over peptide IDs and annotate missing RT values
    for (auto& pep : peptides)
    {
        if (std::isnan(pep.getRT())) // Only annotate peptides with missing RT
        {
            String native_id = pep.getSpectrumReference();
            try
            {
                // Look up spectrum index by native ID and assign RT
                Size index = lookup.findByNativeID(native_id);
                pep.setRT(exp.getSpectra()[index].getRT());
            }
            catch (Exception::ElementNotFound&)
            {
                // Log error if the spectrum reference is not found
                OPENMS_LOG_ERROR << "Error: Failed to look up retention time for peptide identification with spectrum reference '"
                                 << native_id << "' - no spectrum with corresponding native ID found." << endl;
                success = false;
            }
        }
    }

    return success;
}

  bool SpectrumMetaDataLookup::addMissingIMToPeptideIDs(vector<PeptideIdentification>& peptides,
                                const MSExperiment& exp)
  { 
    // Check if the experiment has spectra
    if (exp.getSpectra().empty())
    {
        OPENMS_LOG_INFO << "No spectra found in the experiment. Skipping IM annotation." << endl;
        return false;
    }
    SpectrumLookup lookup;
    bool all_ids_have_im = true;
    lookup.readSpectra(exp.getSpectra());
    // Iterate over peptide_ids and annotate IM values stored in MSExperiment
    for (auto& pep : peptides)
    {
      String native_id = pep.getSpectrumReference();
      Size index = lookup.findByNativeID(native_id);
      const MSSpectrum& spec =  exp.getSpectra()[index];
      // Check if desired IM format is present
	  if (IMTypes::determineIMFormat(spec) == IMFormat::MULTIPLE_SPECTRA)
	  {
        pep.setMetaValue(Constants::UserParam::IM, spec.getDriftTime());
	  }
	  else
	  {
		all_ids_have_im = false;
	  }

    }
    return all_ids_have_im;
  }

  bool SpectrumMetaDataLookup::addMissingSpectrumReferences(vector<PeptideIdentification>& peptides, const String& filename,
    bool stop_on_error, 
    bool override_spectra_data, 
    bool override_spectra_references,
    vector<ProteinIdentification> proteins)
  {
    bool success = true;
    PeakMap exp;
    SpectrumMetaDataLookup lookup;
    if (lookup.empty())
    {
      FileHandler fh;
      auto opts = fh.getOptions();
      opts.setFillData(false);
      opts.setSkipXMLChecks(true);
      fh.setOptions(opts);
      fh.loadExperiment(filename, exp, {FileTypes::MZXML, FileTypes::MZML, FileTypes::MZDATA, FileTypes::MGF}, OpenMS::ProgressLogger::NONE, true, true);
      lookup.readSpectra(exp.getSpectra());
      lookup.setSpectraDataRef(filename);
    }
    if (override_spectra_data)
    {
      vector<String> spectra_data(1);
      spectra_data[0] = "file://" + lookup.spectra_data_ref;
      for (auto& prot : proteins)
      {
        prot.setMetaValue("spectra_data", spectra_data);
      }
    }
    for (auto& pep : peptides)
    {
      // spectrum reference already set? skip if we don't want to overwrite
      if (!override_spectra_references && pep.metaValueExists("spectrum_reference"))
      {
        continue;
      }

      try
      {
        Size index = lookup.findByRT(pep.getRT());
        SpectrumMetaDataLookup::SpectrumMetaData meta;
        lookup.getSpectrumMetaData(index, meta);
        pep.setSpectrumReference( meta.native_id);
      }
      catch (Exception::ElementNotFound&)
      {
        OPENMS_LOG_ERROR << "Error: Failed to look up spectrum native ID for peptide identification with retention time '" + String(pep.getRT()) + "'." << endl;
        success = false;
        if (stop_on_error)
        {
          break;
        }
      }
    }

    return success;
  }


} // namespace OpenMS

