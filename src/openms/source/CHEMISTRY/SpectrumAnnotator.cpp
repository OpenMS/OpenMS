// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrumHelper.h>

#include <algorithm>
#include <numeric>
#include <boost/regex.hpp>

using namespace std;

namespace OpenMS
{
  const boost::regex SpectrumAnnotator::nt_regex_("[a,b,c][[:digit:]]+[+]+");
  const boost::regex SpectrumAnnotator::ct_regex_("[x,y,z][[:digit:]]+[+]+");

  SpectrumAnnotator::SpectrumAnnotator() :
    DefaultParamHandler("SpectrumAnnotator")
  {
    defaults_.setValue("basic_statistics", "true", "If set, meta values for peak_number, sum_intensity, matched_ion_number, matched_intensity are added");
    defaults_.setValidStrings("basic_statistics", ListUtils::create<String>("true,false"));
    /** adds meta values
     peak_number
     sum_intensity
     matched_ion_number
     matched_intensity
    */
    defaults_.setValue("list_of_ions_matched", "true", "If set, meta values for matched_ions are added");
    defaults_.setValidStrings("list_of_ions_matched", ListUtils::create<String>("true,false"));
    /** adds meta values
     matched_ions
    */
    defaults_.setValue("max_series", "true", "If set, meta values for max_series_type, max_series_size are added");
    defaults_.setValidStrings("max_series", ListUtils::create<String>("true,false"));
    /** adds meta values
     max_series_type
     max_series_size
    */
    defaults_.setValue("S/N_statistics", "true", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("S/N_statistics", ListUtils::create<String>("true,false"));
    /** adds meta values
     sn_by_matched_intensity
     sn_by_median_intensity
    */
    defaults_.setValue("precursor_statistics", "true", "If set, meta values for precursor_in_ms2 are added");
    defaults_.setValidStrings("precursor_statistics", ListUtils::create<String>("true,false"));
    /** adds meta values
     precursor_in_ms2
    */
    defaults_.setValue("topNmatch_fragmenterrors", 7, "If set n > 0,  meta values for topN_meanfragmenterror, topN_MSEfragmenterror, topN_stddevfragmenterror are added");
    defaults_.setValidStrings("topNmatch_fragmenterrors", ListUtils::create<String>("true,false"));
    /** adds meta values
     topN_meanfragmenterror
     topN_MSEfragmenterror
     topN_stddevfragmenterror
    */
    // TODO right now topN matched, but additional information would be gained if the topN would be assembled from NOT ONLY the matched (high intense peak not identified?)

    defaults_.setValue("fragmenterror_statistics", "true", "If set, meta values for median_fragment_error, IQR_fragment_error are added");
    defaults_.setValidStrings("fragmenterror_statistics", ListUtils::create<String>("true,false"));
    /** adds meta values
     median_fragment_error
     IQR_fragment_error
    */
    defaults_.setValue("terminal_series_match_ratio", "true", "If set, meta values for NTermIonCurrentRatio, CTermIonCurrentRatio are added");
    defaults_.setValidStrings("terminal_series_match_ratio", ListUtils::create<String>("true,false"));
    /** adds meta values
     NTermIonCurrentRatio
     CTermIonCurrentRatio
    */

    defaultsToParam_();
  }

  SpectrumAnnotator::SpectrumAnnotator(const SpectrumAnnotator & rhs) :
    DefaultParamHandler(rhs)
  {
  }

  SpectrumAnnotator & SpectrumAnnotator::operator=(const SpectrumAnnotator & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
    }
    return *this;
  }

  SpectrumAnnotator::~SpectrumAnnotator()
  {
  }

  RichPeakSpectrum SpectrumAnnotator::annotateMatches(const PeptideHit& ph, const MSSpectrum<Peak1D>& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const
  {
    MSSpectrum<RichPeak1D> rich_theoretical_spec;
    MSSpectrum<RichPeak1D> rich_input_spectrum = MSSpectrumHelper::clone(spec);

    vector<pair<Size, Size> > al;
    tg.getSpectrum(rich_theoretical_spec, ph.getSequence(), ph.getCharge()-1); //will get y5++, b2+, ... for precursor+++
    if (!rich_theoretical_spec.isSorted())
    {
      rich_theoretical_spec.sortByPosition();
    }
    sa.getSpectrumAlignment(al, rich_theoretical_spec, rich_input_spectrum); //peaks from theor. may be matched to none or one in spec!
    for (vector<pair<Size, Size > >::const_iterator it = al.begin(); it != al.end(); ++it)
    {
      rich_input_spectrum[it->second].setMetaValue("IonMatchError", std::fabs(rich_input_spectrum[it->second].getMZ() - rich_theoretical_spec[it->first].getMZ()));
    }
    Param sap = sa.getParameters();
    rich_input_spectrum.setMetaValue("fragment_match_tolerance", sap.getValue("tolerance"));
    return rich_input_spectrum;
  }

  void SpectrumAnnotator::addIonMatchStatistics(PeptideIdentification& pi, const MSSpectrum<Peak1D>& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const
  {
    if (!spec.empty())
    {
      for (vector<PeptideHit>::iterator ph = pi.getHits().begin(); ph != pi.getHits().end(); ++ph)
      {
        RichPeakSpectrum rich_spec = annotateMatches(*ph, spec, tg, sa);
        rich_spec.sortByIntensity();

        StringList ions;
        double sum_intensity = 0;
        double match_intensity = 0;
        vector<double> fragmenterrors, intensities, mzs; // sorted by ascending intensity via rich_spec.sortByIntensity for topN statistics
        fragmenterrors.reserve(rich_spec.size());
        intensities.reserve(rich_spec.size());
        mzs.reserve(rich_spec.size());

        double nint = 0;
        double cint = 0;

        StringList allowed_types = ListUtils::create<String>("y,b,a,c,x,z");
        map<String, vector<bool> > ion_series;
        for (StringList::iterator st = allowed_types.begin(); st != allowed_types.end(); ++st)
        {
          ion_series.insert(make_pair(*st, vector<bool>(ph->getSequence().size()-1, false)));
        }

        for (RichPeakSpectrum::Iterator it = rich_spec.begin(); it != rich_spec.end(); ++it)
        {
          sum_intensity += it->getIntensity();
          if (it->metaValueExists("IonMatchError"))
          {
            fragmenterrors.push_back(it->getMetaValue("IonMatchError"));
            intensities.push_back(it->getIntensity());
            match_intensity += it->getIntensity();
            mzs.push_back(it->getMZ());

            String ion_name = it->getMetaValue("IonName");
            ions.push_back(ion_name);
            if (terminal_series_match_ratio_)
            {
              if (boost::regex_match(ion_name, nt_regex_))
              {
                nint += it->getIntensity();
              }
              else if (boost::regex_match(ion_name, ct_regex_))
              {
                cint += it->getIntensity();
              }
            }
            if (max_series_)
            {
              String ion_type = ion_name.prefix(1);
              if (ListUtils::contains(allowed_types, ion_type)) //case sensitive?
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
          }
        }
        if (basic_statistics_)
        {
          ph->setMetaValue("matched_ions", ListUtils::concatenate(ions, ","));
          ph->setMetaValue("matched_intensity", match_intensity);
          ph->setMetaValue("matched_ion_number", ions.size());
          ph->setMetaValue("peak_number", spec.size());
          ph->setMetaValue("sum_intensity", sum_intensity);
        }
        if (terminal_series_match_ratio_)
        {
          ph->setMetaValue("NTermIonCurrentRatio", nint/match_intensity);
          ph->setMetaValue("CTermIonCurrentRatio", cint/match_intensity);
        }
        if (topNmatch_fragmenterrors_)
        {
          if (fragmenterrors.empty())
          {
            ph->setMetaValue("median_fragment_error", 0);
            ph->setMetaValue("IQR_fragment_error", 0);
            ph->setMetaValue("topN_meanfragmenterror", 0);
            ph->setMetaValue("topN_MSEfragmenterror", 0);
            ph->setMetaValue("topN_stddevfragmenterror", 0);
          }
          else
          {
            vector<double> fe(fragmenterrors);
            fe.reserve(fragmenterrors.size());
            std::size_t mid = fe.size()/2;
            std::size_t lq = fe.size()/4;
            std::size_t uq = lq + mid;
            std::nth_element(fe.begin(), fe.begin()+mid, fe.end());
            if (fe.size() % 2 != 0)
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
            std::reverse_copy(fragmenterrors.begin(), fragmenterrors.end(), topn_fe.begin());
            topn_fe.resize(topNmatch_fragmenterrors_);
            double sum = std::accumulate(topn_fe.begin(), topn_fe.end(), 0.0);
            double mean = sum / topn_fe.size();

            std::vector<double> diff(topn_fe.size());
            std::transform(topn_fe.begin(), topn_fe.end(), diff.begin(),
                           std::bind2nd(std::minus<double>(), mean));
            double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            double m_sq_sum = (sq_sum / topn_fe.size());
            double stdev = std::sqrt(sq_sum / topn_fe.size());

            ph->setMetaValue("topN_meanfragmenterror", mean);
            ph->setMetaValue("topN_MSEfragmenterror", m_sq_sum);
            ph->setMetaValue("topN_stddevfragmenterror", stdev);
          }
        }
        if (max_series_)
        {
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
        }
        //TODO parent peak intensity complement pairs number
        if (SN_statistics_)
        {
          float sn_by_matched_intensity = (match_intensity/ions.size())/
                  ((sum_intensity-match_intensity)/(spec.size()-ions.size()));
          ph->setMetaValue("sn_by_matched_intensity", sn_by_matched_intensity);

          float median = 0;
          //rich_spec is already in sorted order of intensity
          if (rich_spec.size() % 2 == 0)
            median = (rich_spec[rich_spec.size()/2 - 1].getIntensity() + rich_spec[rich_spec.size()/2].getIntensity()) / 2;
          else
            median = rich_spec[rich_spec.size()/2].getIntensity();
          float sign_int= 0;
          float nois_int = 0;
          size_t sign_count= 0;
          size_t nois_count = 0;
          for (MSSpectrum<Peak1D>::const_iterator pt = spec.begin(); pt != spec.end(); ++pt)
          {
            if (pt->getIntensity() <= median)
            {
              ++nois_count;
              nois_int += pt->getIntensity();
            }
            else
            {
              ++sign_count;
              sign_int += pt->getIntensity();
            }
          }
          float sn_by_median_intensity = (sign_int/sign_count)/(nois_int/nois_count);

          ph->setMetaValue("sn_by_median_intensity", sn_by_median_intensity);
        }
        if (precursor_statistics_)
        {
          bool precursor = false;
          for (std::vector<Precursor>::const_iterator pit = spec.getPrecursors().begin(); pit != spec.getPrecursors().end(); ++pit)
          {
            rich_spec.sortByPosition();
            //TODO what about precursor_H2O_loss and precursor_NH3_loss
            if (rich_spec.findNearest(pit->getMZ(),sa.getParameters().getValue("tolerance"),
                                 sa.getParameters().getValue("tolerance")) > -1)
            {
              precursor = true;
            }
          }
          ph->setMetaValue("precursor_in_ms2", precursor);
        }
        //TODO add "FragmentArray"s
      }
    }
  }

  void SpectrumAnnotator::updateMembers_()
  {
    basic_statistics_ = param_.getValue("basic_statistics").toBool();
    list_of_ions_matched_ = param_.getValue("list_of_ions_matched").toBool();
    max_series_ = param_.getValue("max_series").toBool();
    SN_statistics_ = param_.getValue("S/N_statistics").toBool();
    precursor_statistics_ = param_.getValue("precursor_statistics").toBool();
    topNmatch_fragmenterrors_ = (uint)param_.getValue("topNmatch_fragmenterrors");
    fragmenterror_statistics_ = param_.getValue("fragmenterror_statistics").toBool();
    terminal_series_match_ratio_ = param_.getValue("terminal_series_match_ratio").toBool();
  }

}
