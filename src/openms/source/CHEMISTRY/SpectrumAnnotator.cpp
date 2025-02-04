// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/MATH/MathFunctions.h>

using namespace std;

namespace OpenMS
{
  const boost::regex SpectrumAnnotator::nt_regex_("[a,b,c][[:digit:]]+[+]+");
  const boost::regex SpectrumAnnotator::ct_regex_("[x,y,z][[:digit:]]+[+]+");
  const boost::regex SpectrumAnnotator::noloss_regex_("[a,b,c,x,y,z][[:digit:]]+[+]+");
  const boost::regex SpectrumAnnotator::seriesposition_regex_("[a,b,c,x,y,z]([[:digit:]]+)[+,-]+[[:word:]]*[+]*");

  SpectrumAnnotator::SpectrumAnnotator() :
    DefaultParamHandler("SpectrumAnnotator")
  {
    defaults_.setValue("basic_statistics", "true", "If set, meta values for peak_number, sum_intensity, matched_ion_number, matched_intensity are added");
    defaults_.setValidStrings("basic_statistics", {"true","false"});
    /** adds meta values
     peak_number
     sum_intensity
     matched_ion_number
     matched_intensity
    */
    defaults_.setValue("list_of_ions_matched", "true", "If set, meta values for matched_ions are added");
    defaults_.setValidStrings("list_of_ions_matched", {"true","false"});
    /** adds meta values
     matched_ions
    */
    defaults_.setValue("max_series", "true", "If set, meta values for max_series_type, max_series_size are added");
    defaults_.setValidStrings("max_series", {"true","false"});
    /** adds meta values
     max_series_type
     max_series_size
    */
    defaults_.setValue("S/N_statistics", "true", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("S/N_statistics", {"true","false"});
    /** adds meta values
     sn_by_matched_intensity
     sn_by_median_intensity
    */
    defaults_.setValue("precursor_statistics", "true", "If set, meta values for precursor_in_ms2 are added");
    defaults_.setValidStrings("precursor_statistics", {"true","false"});
    /** adds meta values
     precursor_in_ms2
    */
    defaults_.setValue("topNmatch_fragmenterrors", unsigned(7), "If set n > 0,  meta values for topN_meanfragmenterror, topN_MSEfragmenterror, topN_stddevfragmenterror are added");
    /** adds meta values
     topN_meanfragmenterror
     topN_MSEfragmenterror
     topN_stddevfragmenterror
    */
    // TODO right now topN matched, but additional information would be gained if the topN would be assembled from NOT ONLY the matched (high intense peak not identified?)

    defaults_.setValue("fragmenterror_statistics", "true", "If set, meta values for median_fragment_error, IQR_fragment_error are added");
    defaults_.setValidStrings("fragmenterror_statistics", {"true","false"});
    /** adds meta values
     median_fragment_error
     IQR_fragment_error
    */
    defaults_.setValue("terminal_series_match_ratio", "true", "If set, meta values for NTermIonCurrentRatio, CTermIonCurrentRatio are added");
    defaults_.setValidStrings("terminal_series_match_ratio", {"true","false"});
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

  SpectrumAnnotator::~SpectrumAnnotator() = default;

  void SpectrumAnnotator::annotateMatches(PeakSpectrum& spec, const PeptideHit& ph, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const
  {
    PeakSpectrum theoretical_spec;
    vector<pair<Size, Size> > al;
    Int zmin = 1;
    Int zmax = 2;
    tg.getSpectrum(theoretical_spec, ph.getSequence(), zmin, min(ph.getCharge(), zmax));
    OPENMS_PRECONDITION(theoretical_spec.isSorted(), "TheoreticalSpectrumGenerator::getSpectrum did not yield a sorted spectrum!")

    if (!spec.isSorted())
    {
      spec.sortByPosition();
    }
    sa.getSpectrumAlignment(al, theoretical_spec, spec);  // peaks from theor. may be matched to none or one in spec!

    PeakSpectrum::StringDataArray theo_annot = theoretical_spec.getStringDataArrays().front();
    PeakSpectrum::StringDataArray type_annotations = PeakSpectrum::StringDataArray();
    PeakSpectrum::FloatDataArray error_annotations = PeakSpectrum::FloatDataArray();
    type_annotations.setName(Constants::UserParam::IonNames);
    error_annotations.setName("IonMatchError");
    type_annotations.resize(spec.size());
    error_annotations.resize(spec.size());
    for (vector<pair<Size, Size > >::const_iterator it = al.begin(); it != al.end(); ++it)
    {
        error_annotations[it->second] = std::fabs(spec[it->second].getMZ() - theoretical_spec[it->first].getMZ());
        type_annotations[it->second] = theo_annot[it->first];
    }
    const Param& sap = sa.getParameters();
    spec.setMetaValue("fragment_mass_tolerance", sap.getValue("tolerance"));
    spec.setMetaValue("fragment_mass_tolerance_ppm", sap.getValue("is_relative_tolerance").toBool());
    spec.setStringDataArrays(PeakSpectrum::StringDataArrays(1, type_annotations));
    spec.setFloatDataArrays(PeakSpectrum::FloatDataArrays(1, error_annotations));
  }

  void SpectrumAnnotator::addIonMatchStatistics(PeptideIdentification& pi, MSSpectrum& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const
  {
    if (!spec.empty())
    {
      for (vector<PeptideHit>::iterator ph = pi.getHits().begin(); ph != pi.getHits().end(); ++ph)
      {
        annotateMatches(spec, *ph, tg, sa);
        spec.sortByIntensity();

        StringList ions;
        double sum_intensity = 0;
        double match_intensity = 0;
        vector<double> fragmenterrors, intensities, mzs;  // sorted by ascending intensity via spec.sortByIntensity for topN statistics
        fragmenterrors.reserve(spec.size());
        intensities.reserve(spec.size());
        mzs.reserve(spec.size());

        double nint = 0;
        double cint = 0;

        StringList allowed_types = ListUtils::create<String>("y,b,a,c,x,z");
        map<String, vector<bool> > ion_series;
        for (StringList::iterator st = allowed_types.begin(); st != allowed_types.end(); ++st)
        {
          ion_series.insert(make_pair(*st, vector<bool>(ph->getSequence().size()-1, false)));
        }

        PeakSpectrum::StringDataArray type_annotations = PeakSpectrum::StringDataArray();
        PeakSpectrum::FloatDataArray error_annotations = PeakSpectrum::FloatDataArray();
        for (PeakSpectrum::StringDataArrays::iterator it = spec.getStringDataArrays().begin(); it != spec.getStringDataArrays().end(); ++it)
        {
          if (it->getName() == Constants::UserParam::IonNames)
            type_annotations = *it;
        }
        for (PeakSpectrum::FloatDataArrays::iterator it = spec.getFloatDataArrays().begin(); it != spec.getFloatDataArrays().end(); ++it)
        {
          if (it->getName() == "IonMatchError")
            error_annotations = *it;
        }

        for (size_t i = 0; i < spec.size(); ++i)
        {
          sum_intensity += spec[i].getIntensity();
          if (!type_annotations.at(i).empty())  // implies error_annotations is set, too.
          {
            fragmenterrors.push_back(error_annotations.at(i));
            intensities.push_back(spec[i].getIntensity());
            match_intensity += spec[i].getIntensity();
            mzs.push_back(spec[i].getMZ());
            const String& ion_name = type_annotations.at(i);
            {
              ions.push_back(ion_name);
              if (terminal_series_match_ratio_)
              {
                if (boost::regex_match(ion_name, nt_regex_))
                {
                  nint += spec[i].getIntensity();
                }
                else if (boost::regex_match(ion_name, ct_regex_))
                {
                  cint += spec[i].getIntensity();
                }
              }
              if (max_series_)  // without loss max series is sometimes pretty crummy
              {
                const String& ion_type = ion_name.prefix(1);
                boost::cmatch what;
                if (boost::regex_match(ion_name.c_str(), what, seriesposition_regex_) &&
                        ListUtils::contains(allowed_types, ion_type))
                {
                  // what[0] contains the whole string
                  // what[1] contains the response code
                  try
                  {
                    int i = std::atoi(what[1].first);
                    ion_series[ion_type].at(i-1) = true;
                  }
                  catch (std::out_of_range&)
                  {
                    OPENMS_LOG_WARN << "Note: Ions of " << ion_type << ion_name.substr(1).remove('+').toInt()
                             << " will be ignored for max_series " << ph->getSequence().toString() << endl;
                    continue;
                  }
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
            std::size_t mid = fe.size() / 2;
            std::size_t lq = fe.size() / 4;
            std::size_t uq = lq + mid;
            std::nth_element(fe.begin(), fe.begin()+mid, fe.end());
            if (fe.size() % 2 != 0)
            {
              ph->setMetaValue("median_fragment_error", fe[mid]);
            }
            else
            {
              double right2mid = fe[mid];
              std::nth_element(fe.begin(), fe.begin() + mid-1, fe.end());
              ph->setMetaValue("median_fragment_error", (right2mid + fe[mid-1]) / 2.0);
            }
            std::nth_element(fe.begin(),          fe.begin() + lq, fe.end());
            std::nth_element(fe.begin() + lq + 1, fe.begin() + mid, fe.end());
            std::nth_element(fe.begin() + mid + 1, fe.begin() + uq, fe.end());
            ph->setMetaValue("IQR_fragment_error", fe[uq]-fe[lq]);

            vector<double> topn_fe;
            topn_fe.resize(fragmenterrors.size());
            std::reverse_copy(fragmenterrors.begin(), fragmenterrors.end(), topn_fe.begin());  // fragmenterrors is sortByIntensity before, get TopN from the back of the vector
            topn_fe.resize(topNmatch_fragmenterrors_);

            double mean = Math::mean(topn_fe.begin(), topn_fe.end());
            double stdev = Math::sd(topn_fe.begin(), topn_fe.end(), mean);

            double sq_sum = 0;
            for (std::vector<double>::iterator it = topn_fe.begin(); it != topn_fe.end(); ++it)
            {
              sq_sum += *it * *it;
            }
            double m_sq_sum = (sq_sum / topn_fe.size());

            ph->setMetaValue("topN_meanfragmenterror", mean);
            ph->setMetaValue("topN_MSEfragmenterror", m_sq_sum);
            ph->setMetaValue("topN_stddevfragmenterror", stdev);
          }
        }
        if (max_series_)
        {
          String max_series;
          int max_stretch = 0;
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
              if (stretch > max_stretch)
              {
                max_stretch = stretch;
                max_series = tt->first;
              }
            }
          }
          ph->setMetaValue("max_series_type", max_series);
          ph->setMetaValue("max_series_size", max_stretch);
        }
        //TODO parent peak intensity complement pairs number
        if (SN_statistics_)
        {
          float sn_by_matched_intensity = (match_intensity / ions.size()) /
                  ((sum_intensity-match_intensity) / (spec.size()-ions.size()));
          if (spec.size() - ions.size() == 0)
          {
            sn_by_matched_intensity = 0;
          }
          ph->setMetaValue("sn_by_matched_intensity", sn_by_matched_intensity);

          float median = 0;
          // spec is already in sorted order of intensity
          if (spec.size() % 2 == 0)
            median = (spec[spec.size() / 2 - 1].getIntensity() + spec[spec.size() / 2].getIntensity()) / 2;
          else
            median = spec[spec.size() / 2].getIntensity();
          float sign_int= 0;
          float nois_int = 0;
          size_t sign_count= 0;
          size_t nois_count = 0;
          for (MSSpectrum::const_iterator pt = spec.begin(); pt != spec.end(); ++pt)
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
          float sn_by_median_intensity = (sign_int / sign_count) / (nois_int / nois_count);
          if (nois_count == 0 || sign_count == 0)
          {
            sn_by_median_intensity = 0;
          }
          ph->setMetaValue("sn_by_median_intensity", sn_by_median_intensity);
        }
        //TODO charge related features might be worth looking at in the future
        if (precursor_statistics_)
        {
          bool precursor = false;
          for (std::vector<Precursor>::const_iterator pit = spec.getPrecursors().begin(); pit != spec.getPrecursors().end(); ++pit)
          {
            spec.sortByPosition();
            //TODO what about precursor_H2O_loss and precursor_NH3_loss
            if (spec.findNearest(pit->getMZ(),sa.getParameters().getValue("tolerance"),
                                 sa.getParameters().getValue("tolerance")) > -1)
            {
              precursor = true;
            }
          }
          ph->setMetaValue("precursor_in_ms2", precursor);
        }
        //TODO add "FragmentArray"s

        const Param& sap = sa.getParameters();
        pi.setMetaValue("fragment_match_tolerance", (double)sap.getValue("tolerance"));
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
    topNmatch_fragmenterrors_ = (unsigned)param_.getValue("topNmatch_fragmenterrors");
    fragmenterror_statistics_ = param_.getValue("fragmenterror_statistics").toBool();
    terminal_series_match_ratio_ = param_.getValue("terminal_series_match_ratio").toBool();
  }

}
