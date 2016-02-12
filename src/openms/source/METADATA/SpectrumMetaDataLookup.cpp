// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <algorithm>
#include <cmath>
#include <list>
#include <numeric>
#include <boost/regex.hpp>

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
    const MSSpectrum<>& spectrum, SpectrumMetaData& meta, 
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

  bool SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(
    vector<PeptideIdentification>& peptides, const String& filename,
    bool stop_on_error, bool reset_basename)
  {
    PeakMap exp;
    SpectrumLookup lookup;
    bool success = true;
    String bn = basename(filename.c_str());
    for (vector<PeptideIdentification>::iterator it = peptides.begin();
         it != peptides.end(); ++it)
    {
      if (reset_basename)
      {
        it->setBaseName(bn);
      }
      if (boost::math::isnan(it->getRT()))
      {
        if (lookup.empty()) // load raw data only if we have to
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
          LOG_ERROR << "Error: Failed to look up retention time for peptide ID with spectrum reference '" + spectrum_id + "' - no spectrum with corresponding native ID found." << endl;
          success = false;
          if (stop_on_error) break;
        }
      }
    }
    return success;
  }

  bool SpectrumMetaDataLookup::addMissingSpectrumReferencestoPeptideIDs(
    vector<PeptideIdentification>& peptides, const String& filename,
    bool stop_on_error, bool add_ionmatches, bool reset_basename)
  {
    MSExperiment<> exp;
    SpectrumLookup lookup;
    bool success = true;
    String bn = basename(filename.c_str());
    for (vector<PeptideIdentification>::iterator it = peptides.begin();
         it != peptides.end(); ++it)
    {
      if (reset_basename)
      {
        it->setBaseName(bn);
      }
      if (!it->metaValueExists("spectrum_reference"))
      {
        if (lookup.empty()) // load raw data only if we have to
        {
          FileHandler().loadExperiment(filename, exp);
          lookup.readSpectra(exp.getSpectra());
        }
        try
        {
          Size index = lookup.findByRT(it->getRT());
          it->setMetaValue("spectrum_reference", exp[index].getNativeID());

          if (add_ionmatches)
          {
            addIonMatches_(*it, exp[index]);
          }
        }
        catch (Exception::ElementNotFound&)
        {
          LOG_ERROR << "Error: Failed to look up spectrum reference by retention time for peptide ID with RT'" + String(it->getRT()) + "' - no spectrum with corresponding RT found." << endl;
          success = false;
          if (stop_on_error) break;
        }
      }
      else
      {
        if (add_ionmatches)
        {
          if (lookup.empty()) // load raw data only if we have to
          {
            FileHandler().loadExperiment(filename, exp);
            lookup.readSpectra(exp.getSpectra());
          }
          String spectrum_id = it->getMetaValue("spectrum_reference");
          try
          {
            Size index = lookup.findByNativeID(spectrum_id);
            addIonMatches_(*it, exp[index]);
          }
          catch (Exception::ElementNotFound&)
          {
            LOG_ERROR << "Error: Failed to look up spectrum reference '" + spectrum_id + "' - for adding IonMatches; no spectrum with corresponding native ID found." << endl;
            success = false;
            if (stop_on_error) break;
          }
        }
      }
    }
    return success;
  }

  void SpectrumMetaDataLookup::addIonMatches_(PeptideIdentification& pi, const MSSpectrum<Peak1D>& spec)
  {
    if (!spec.empty())
    {
        TheoreticalSpectrumGenerator tg = TheoreticalSpectrumGenerator();
        Param tgp(tg.getDefaults());
        tgp.setValue("add_metainfo", "true");
        tgp.setValue("add_losses", "true");
        tgp.setValue("add_precursor_peaks", "true");
        tgp.setValue("add_abundant_immonium_ions", "true");
        tgp.setValue("add_first_prefix_ion", "true");
        tgp.setValue("add_y_ions", "true");
        tgp.setValue("add_b_ions", "true");
        tgp.setValue("add_a_ions", "true");
        tgp.setValue("add_c_ions", "true");
        tgp.setValue("add_x_ions", "true");
        tgp.setValue("add_z_ions", "true");

        tg.setParameters(tgp);
        SpectrumAlignment sa;
        Param sap = sa.getDefaults();
        sap.setValue("tolerance", 0.5, "...");
        sa.setParameters(sap);
        for (vector<PeptideHit>::iterator ph = pi.getHits().begin(); ph != pi.getHits().end(); ++ph)
        {
          RichPeakSpectrum rich_spec;
          vector<pair<Size, Size> > al;
          tg.getSpectrum(rich_spec, ph->getSequence(), 2); //will get y5++ or b2+
          // convert rich spectrum to simple spectrum
          MSSpectrum<Peak1D> gen_spec, new_spec(spec);
          for (RichPeakSpectrum::Iterator it = rich_spec.begin(); it != rich_spec.end(); ++it)
          {
            gen_spec.push_back(static_cast<Peak1D>(*it));
          }
          if (!new_spec.isSorted())
              new_spec.sortByPosition();
          sa.getSpectrumAlignment(al, gen_spec, new_spec); //peaks from theor. may be matched to none or one in spec!
          StringList ions;
          double match_intensity = 0;
          //TODO ionseries max length
          StringList allowed_types = ListUtils::create<String>("y,b,a,c,x,z");
          map<String, vector<bool> > ion_series;
          for (StringList::iterator st = allowed_types.begin(); st != allowed_types.end(); ++st)
          {
            ion_series.insert(make_pair(*st, vector<bool>(ph->getSequence().size()-1, false)));
          }
          map<double,double> fragmenterrors; // mapped from intensity to error for topN statistics
          double nint = 0;
          double cint = 0;
          for (vector<pair<Size, Size > >::const_iterator it = al.begin(); it != al.end(); ++it)
          {
            match_intensity += new_spec[it->second].getIntensity();
            String ion_name = rich_spec[it->first].getMetaValue("IonName");
            const boost::regex nt_regex("[a,b,c][[:digit:]]+[+]+");
            const boost::regex ct_regex("[x,y,z][[:digit:]]+[+]+");
            if (boost::regex_match(ion_name, nt_regex))
            {
              nint += new_spec[it->second].getIntensity();
            }
            else if (boost::regex_match(ion_name, ct_regex))
            {
              cint += new_spec[it->second].getIntensity();
            }

            ions.push_back(ion_name);
            fragmenterrors.insert(make_pair(new_spec[it->second].getIntensity(),std::fabs(new_spec[it->second].getMZ() - rich_spec[it->first].getMZ())));
            String ion_type = ion_name.prefix(1);
            if (ListUtils::contains(allowed_types, ion_type))
            {
              try
              {
                int i = ion_name.substr(1).remove('+').toInt() - 1;
                ion_series[ion_type].at(i) = true;
              }
              catch (std::out_of_range)
              {
                LOG_WARN << ion_type << ion_name.substr(1).remove('+').toInt() -1
                         << " ions not foreseen for" << ph->getSequence().toString() << endl;
                continue;
              }
            }
          }
          ph->setMetaValue("NTermIonCurrentRatio", nint/match_intensity);
          ph->setMetaValue("CTermIonCurrentRatio", cint/match_intensity);

          if(fragmenterrors.empty())
          {
              ph->setMetaValue("median_fragment_error", 0);
              ph->setMetaValue("IQR_fragment_error", 0);
              ph->setMetaValue("top7_meanfragmenterror", 0);
              ph->setMetaValue("top7_MSEfragmenterror", 0);
              ph->setMetaValue("top7_stddevfragmenterror", 0);
          }
          else
          {
            vector<double> fe;
            fe.reserve(fragmenterrors.size());
            for (map<double,double>::const_iterator fit = fragmenterrors.begin(); fit != fragmenterrors.end(); ++fit)
            {
              fe.push_back(fit->second);
            }
            std::size_t mid = fe.size()/2;
            std::size_t lq = fe.size()/4;
            std::size_t uq = lq + mid;
            std::nth_element(fe.begin(), fe.begin()+mid, fe.end());
            if(fe.size() % 2 != 0)
            {
              ph->setMetaValue("median_fragment_error", fe[mid]);
            }
            else
            {
              double right2mid = fe[mid];
              std::nth_element(fe.begin(), fe.begin()+mid-1, fe.end());
              ph->setMetaValue("median_fragment_error", (right2mid+fe[mid-1])/2.0);
            }
            std::nth_element(fe.begin(),          fe.begin() + lq, fe.end());
            std::nth_element(fe.begin() + lq + 1, fe.begin() + mid, fe.end());
            std::nth_element(fe.begin() + mid + 1, fe.begin() + uq, fe.end());
            ph->setMetaValue("IQR_fragment_error", fe[uq]-fe[lq]);

            vector<double> topn_fe;
            topn_fe.reserve(7);
            //apparently i cannot do this: map<double,double>::reverse_iterator rtop7 = fragmenterrors.rbegin() + 8;
            for (map<double,double>::reverse_iterator fit = fragmenterrors.rbegin(); fit != fragmenterrors.rend() && topn_fe.size() < 7; ++fit)
            {
              topn_fe.push_back(fit->second);
            }
            double sum = std::accumulate(topn_fe.begin(), topn_fe.end(), 0.0);
            double mean = sum / topn_fe.size();

            std::vector<double> diff(topn_fe.size());
            std::transform(topn_fe.begin(), topn_fe.end(), diff.begin(),
                           std::bind2nd(std::minus<double>(), mean));
            double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            double m_sq_sum = (sq_sum / topn_fe.size());
            double stdev = std::sqrt(sq_sum / topn_fe.size());

            ph->setMetaValue("top7_meanfragmenterror", mean);
            ph->setMetaValue("top7_MSEfragmenterror", m_sq_sum);
            ph->setMetaValue("top7_stddevfragmenterror", stdev);
          }

          ph->setMetaValue("matched_ions", ListUtils::concatenate(ions, ","));
          ph->setMetaValue("matched_intensity", match_intensity);
          ph->setMetaValue("matched_ion_number", ions.size());

          String max_series;
          int series_stretch = 0;
          for (map<String, vector<bool> >::iterator tt = ion_series.begin(); tt != ion_series.end(); ++tt)
          {
            int stretch = 0;
            for (vector<bool>::iterator it = tt->second.begin(); it != tt->second.end(); ++it)
            {
              if (*it)
              {
                ++stretch;
              }
              else
              {
                stretch = 0;
              }
            }
            if (stretch > series_stretch)
            {
              series_stretch = stretch;
              max_series = tt->first;
            }
          }
          ph->setMetaValue("max_series_type", max_series);
          ph->setMetaValue("max_series_size", series_stretch);

          float sum_intensity = 0;
          int peak_number = 0;
          for (MSSpectrum<Peak1D>::const_iterator pt = spec.begin(); pt != spec.end(); ++pt)
          {
            sum_intensity += pt->getIntensity();
            //TODO parent peak intensity complement pairs number
            ++peak_number;
          }
          ph->setMetaValue("peak_number", peak_number);
          ph->setMetaValue("sum_intensity", sum_intensity);

          float sn_by_matched_intensity = (match_intensity/ions.size())/
                  ((sum_intensity-match_intensity)/(peak_number-ions.size()));
          ph->setMetaValue("sn_by_matched_intensity", sn_by_matched_intensity);

          bool precursor = false;
          Precursor p = spec.getPrecursors().front();
          new_spec.sortByPosition();
          //TODO precursor_H2O_loss and precursor_NH3_loss
          if (new_spec.findNearest(p.getMZ(),sap.getValue("tolerance"),
                               sap.getValue("tolerance")) > -1)
          {
            precursor = true;
          }
          ph->setMetaValue("precursor_in_ms2", precursor);

          float median = 0;
          new_spec.sortByIntensity();
          if(new_spec.size() % 2 == 0)
            median = (new_spec[new_spec.size()/2 - 1].getIntensity() + new_spec[new_spec.size()/2].getIntensity()) / 2;
          else
            median = new_spec[new_spec.size()/2].getIntensity();
          float sign_int= 0;
          float nois_int = 0;
          size_t sign_cnt= 0;
          size_t nois_cnt = 0;
          for (MSSpectrum<Peak1D>::const_iterator pt = spec.begin(); pt != spec.end(); ++pt)
          {
            if (pt->getIntensity() <= median)
            {
              ++nois_cnt;
              nois_int += pt->getIntensity();
            }
            else
            {
              ++sign_cnt;
              sign_int += pt->getIntensity();
            }
          }
          float sn_by_median_intensity = (sign_int/sign_cnt)/(nois_int/nois_cnt);
          ph->setMetaValue("sn_by_median_intensity", sn_by_median_intensity);

        }


    }
  }

} // namespace OpenMS
