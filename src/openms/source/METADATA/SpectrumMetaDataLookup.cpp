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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

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
        LOG_ERROR << "Error: Could not extract scan number from spectrum native ID '" + meta.native_id + "' using regular expression '" + scan_regexp.str() + "'." << endl;
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
          LOG_ERROR << "Error: Could not set precursor RT for spectrum with native ID '" + meta.native_id + "' - precursor spectrum not found." << endl;
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


  bool SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(vector<PeptideIdentification>& peptides, const String& filename,
    bool stop_on_error)
  {
    PeakMap exp;
    SpectrumLookup lookup;
    bool success = true;
    for (vector<PeptideIdentification>::iterator it = peptides.begin();
            it != peptides.end(); ++it)
    {
      if (boost::math::isnan(it->getRT()))
      {
        if (lookup.empty())
        {
          FileHandler().loadExperiment(filename, exp);
          lookup.readSpectra(exp.getSpectra());
        }
        String spectrum_id = it->getMetaValue("spectrum_reference");
        try
        {
          Size index = lookup.findByNativeID(spectrum_id);
          it->setRT(exp[index].getRT());
        }
        catch (Exception::ElementNotFound&)
        {
          LOG_ERROR << "Error: Failed to look up retention time for peptide identification with spectrum reference '" + spectrum_id + "' - no spectrum with corresponding native ID found." << endl;
          success = false;
          if (stop_on_error) break;
        }
      }
    }
    return success;
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
      FileHandler().loadExperiment(filename, exp);
      lookup.readSpectra(exp.getSpectra());
      lookup.setSpectraDataRef(filename);
    }
    if (override_spectra_data)
    {
      vector<String> spectra_data(1);
      spectra_data[0] = "file://" + lookup.spectra_data_ref;
      for (vector<ProteinIdentification>::iterator it =
            proteins.begin(); it !=
            proteins.end(); ++it)
      {
        it->setMetaValue("spectra_data", spectra_data);
      }
    }
    for (auto it = peptides.begin(); it != peptides.end(); ++it)
    {
      // spectrum reference already set? skip if we don't want to overwrite
      if (!override_spectra_references 
	&& it->metaValueExists("spectrum_reference"))
      {
        continue;
      }

      try
      {
        Size index = lookup.findByRT(it->getRT());
        SpectrumMetaDataLookup::SpectrumMetaData meta;
        lookup.getSpectrumMetaData(index, meta);
        it->setMetaValue("spectrum_reference", meta.native_id);
      }
      catch (Exception::ElementNotFound&)
      {
        LOG_ERROR << "Error: Failed to look up spectrum native ID for peptide identification with retention time '" + String(it->getRT()) + "'." << endl;
        success = false;
        if (stop_on_error) break;
      }
    }

    return success;
  }


} // namespace OpenMS

