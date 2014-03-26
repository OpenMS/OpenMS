// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>

#include <cmath>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

namespace OpenMS
{

  typedef std::iterator_traits< std::vector<double>::const_iterator >::difference_type DiffType;

  Size findNearest(const std::vector<double>& tmp, double mz)
  {
    // no peak => no search
    if (tmp.empty()) throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "There must be at least one peak to determine the nearest peak!");

    // search for position for inserting
    std::vector<double>::const_iterator it = std::lower_bound(tmp.begin(), tmp.end(), mz);
    // border cases
    if (it == tmp.begin()) return 0;

    if (it == tmp.end()) return tmp.size() - 1;

    // the peak before or the current peak are closest
    std::vector<double>::const_iterator it2 = it;
    --it2;
    if (std::fabs(*it - mz) < std::fabs(*it2 - mz))
    {
      return Size(it - tmp.begin());
    }
    else
    {
      return Size(it2 - tmp.begin());
    }
  }

  bool hasActivationMethod(OpenMS::OMSInterfaces::SpectrumPtr s, std::string method)
  {
    for (std::vector<OpenMS::OMSInterfaces::Precursor>::const_iterator it = s->getPrecursors().begin(); it != s->getPrecursors().end(); ++it)
    {
      for (std::set<std::string>::const_iterator it_a = it->activation_methods.begin(); it_a != it->activation_methods.end(); ++it_a)
      {
        if (method == *it_a)
        {
          // found matching activation method
          return true;
        }
      }
    }
    return false;
  }

  IsobaricChannelExtractor::IsobaricChannelExtractor(const IsobaricQuantitationMethod* const quant_method) :
    DefaultParamHandler("IsobaricChannelExtractor"),
    quant_method_(quant_method),
    selected_activation_(""),
    reporter_mass_shift_(0.1),
    min_precursor_intensity_(1.0),
    keep_unannotated_precursor_(true),
    min_reporter_intensity_(0.0),
    remove_low_intensity_quantifications_(false),
    min_precursor_purity_(0.0),
    max_precursor_isotope_deviation_(0.02)
  {
    setDefaultParams_();
  }

  IsobaricChannelExtractor::IsobaricChannelExtractor(const IsobaricChannelExtractor& other) :
    DefaultParamHandler(other),
    quant_method_(other.quant_method_),
    selected_activation_(other.selected_activation_),
    reporter_mass_shift_(other.reporter_mass_shift_),
    min_precursor_intensity_(other.min_precursor_intensity_),
    keep_unannotated_precursor_(other.keep_unannotated_precursor_),
    min_reporter_intensity_(other.min_reporter_intensity_),
    remove_low_intensity_quantifications_(other.remove_low_intensity_quantifications_),
    min_precursor_purity_(other.min_precursor_purity_),
    max_precursor_isotope_deviation_(other.max_precursor_isotope_deviation_)
  {
  }

  IsobaricChannelExtractor& IsobaricChannelExtractor::operator=(const IsobaricChannelExtractor& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    quant_method_ = rhs.quant_method_;
    selected_activation_ = rhs.selected_activation_;
    reporter_mass_shift_ = rhs.reporter_mass_shift_;
    min_precursor_intensity_ = rhs.min_precursor_intensity_;
    keep_unannotated_precursor_ = rhs.keep_unannotated_precursor_;
    min_reporter_intensity_ = rhs.min_reporter_intensity_;
    remove_low_intensity_quantifications_ = rhs.remove_low_intensity_quantifications_;
    min_precursor_purity_ = rhs.min_precursor_purity_;
    max_precursor_isotope_deviation_ = rhs.max_precursor_isotope_deviation_;

    return *this;
  }

  void IsobaricChannelExtractor::setDefaultParams_()
  {
    defaults_.setValue("select_activation", Precursor::NamesOfActivationMethod[Precursor::HCID], "Operate only on MSn scans where any of its precursors features a certain activation method (e.g., usually HCD for iTRAQ). Set to empty string if you want to disable filtering.");
    StringList activation_list;
    activation_list.insert(activation_list.begin(), Precursor::NamesOfActivationMethod, Precursor::NamesOfActivationMethod + Precursor::SIZE_OF_ACTIVATIONMETHOD - 1);
    activation_list.push_back(""); // allow disabling this
    
    defaults_.setValidStrings("select_activation", activation_list);

    defaults_.setValue("reporter_mass_shift", 0.1, "Allowed shift (left to right) in Da from the expected position.");
    defaults_.setMinFloat("reporter_mass_shift", 0.00000001);
    defaults_.setMaxFloat("reporter_mass_shift", 0.5);

    defaults_.setValue("min_precursor_intensity", 1.0, "Minimum intensity of the precursor to be extracted. MS/MS scans having a precursor with a lower intensity will not be considered for quantitation.");
    defaults_.setMinFloat("min_precursor_intensity", 0.0);

    defaults_.setValue("keep_unannotated_precursor", "true", "Flag if precursor with missing intensity value or missing precursor spectrum should be included or not.");
    defaults_.setValidStrings("keep_unannotated_precursor", ListUtils::create<String>("true,false"));

    defaults_.setValue("min_reporter_intensity", 0.0, "Minimum intenesity of the individual reporter ions to be used extracted.");
    defaults_.setMinFloat("min_reporter_intensity", 0.0);

    defaults_.setValue("discard_low_intensity_quantifications", "false", "Remove all reporter intensities if a single reporter is below the threshold given in min_reporter_intensity.");
    defaults_.setValidStrings("discard_low_intensity_quantifications", ListUtils::create<String>("true,false"));

    defaults_.setValue("min_precursor_purity", 0.0, "Minimum fraction of the total intensity in the isolation window of the precursor spectrum attributable to the selected precursor.");
    defaults_.setMinFloat("min_precursor_purity", 0.0);
    defaults_.setMaxFloat("min_precursor_purity", 1.0);

    defaults_.setValue("precursor_isotope_deviation", 0.02, "Maximum allowed deviation between theoretical and observed isotopic peaks of the precursor peak in the isolation window to be counted as part of the precursor.");
    defaults_.setMinFloat("precursor_isotope_deviation", 0.0);
    defaults_.addTag("precursor_isotope_deviation", "advanced");

    defaultsToParam_();
  }

  void IsobaricChannelExtractor::updateMembers_()
  {
    selected_activation_ = getParameters().getValue("select_activation");
    reporter_mass_shift_ = getParameters().getValue("reporter_mass_shift");
    min_precursor_intensity_ = getParameters().getValue("min_precursor_intensity");
    keep_unannotated_precursor_ = getParameters().getValue("keep_unannotated_precursor") == "true";
    min_reporter_intensity_ = getParameters().getValue("min_reporter_intensity");
    remove_low_intensity_quantifications_ = getParameters().getValue("discard_low_intensity_quantifications") == "true";
    min_precursor_purity_ = getParameters().getValue("min_precursor_purity");
    max_precursor_isotope_deviation_ = getParameters().getValue("precursor_isotope_deviation");
  }

  bool IsobaricChannelExtractor::isValidPrecursor_(const OpenMS::OMSInterfaces::Precursor& precursor) const
  {
    return (precursor.intensity == 0.0 && keep_unannotated_precursor_) || !(precursor.intensity < min_precursor_intensity_);
  }

  bool IsobaricChannelExtractor::hasLowIntensityReporter_(const ConsensusFeature& cf) const
  {
    for (ConsensusFeature::const_iterator cf_it = cf.begin();
         cf_it != cf.end();
         ++cf_it)
    {
      if (cf_it->getIntensity() == 0.0)
      {
        return true;
      }
    }

    return false;
  }

  DoubleReal IsobaricChannelExtractor::sumPotentialIsotopePeaks_(const OpenMS::OMSInterfaces::SpectrumPtr precursor,
                                                                 const double lower_mz_bound,
                                                                 const double upper_mz_bound,
                                                                 double theoretical_mz,
                                                                 const double isotope_offset) const
  {
    DoubleReal intensity_contribution = 0.0;

    // move theoretical_mz to first potential isotopic peak
    theoretical_mz += isotope_offset;

    // check if we are still in the isolation window
    while (theoretical_mz > lower_mz_bound && theoretical_mz < upper_mz_bound)
    {
      Size potential_peak = findNearest(precursor->getPos(), theoretical_mz);

      // is isotopic ?
      if (fabs(theoretical_mz - precursor->getPos()[potential_peak]) < max_precursor_isotope_deviation_)
      {
        intensity_contribution += precursor->getIntensity()[potential_peak];
      }
      else
      {
        // we abort in case of missing peaks
        break;
      }

      // update mz with the defined offset
      theoretical_mz += isotope_offset;
    }

    return intensity_contribution;
  }


  DoubleReal IsobaricChannelExtractor::computePrecursorPurity_(const OpenMS::OMSInterfaces::SpectrumPtr ms2_spec, const OpenMS::OMSInterfaces::SpectrumPtr precursor) const
  {
    // we cannot analyze precursors without a charge
    if (ms2_spec->getPrecursors()[0].charge == 0)
      return 1.0;

    // compute boundaries
    std::vector<double>::const_iterator isolation_lower_mz_ = std::lower_bound(precursor->getPos().begin(), precursor->getPos().end(), 
      ms2_spec->getPrecursors()[0].mz - ms2_spec->getPrecursors()[0].lower_offset);
    std::vector<double>::const_iterator isolation_upper_mz_ = std::upper_bound(precursor->getPos().begin(), precursor->getPos().end(), 
      ms2_spec->getPrecursors()[0].mz + ms2_spec->getPrecursors()[0].lower_offset);
    std::vector<double>::const_iterator int_it_ = precursor->getIntensity().begin();
    DiffType iterator_pos = std::distance((std::vector<double>::const_iterator)precursor->getPos().begin(), isolation_lower_mz_);
    std::advance(int_it_, iterator_pos);

    // get total intensity
    double total_intensity = 0;

    for (std::vector<double>::const_iterator mz_it_ = isolation_lower_mz_; mz_it_ != isolation_upper_mz_; ++mz_it_, ++int_it_)
    {
      total_intensity += *int_it_;
    }

    // now get the intensity of the precursor .. we assume everything in the distance of 1/c to belong to the precursor
    // for c == charge of precursor
    
    // precursor mz
    Size precursor_peak_idx = findNearest(precursor->getPos(), ms2_spec->getPrecursors()[0].mz);
    double precursor_intensity = precursor->getIntensity()[precursor_peak_idx];

    // compute the
    double charge_dist = Constants::NEUTRON_MASS_U / (double) ms2_spec->getPrecursors()[0].charge;

    // TODO : not tested in class test! 0-> int always stays equal !
    // search left of precursor for isotopic peaks
    precursor_intensity += sumPotentialIsotopePeaks_(precursor, *isolation_lower_mz_, *isolation_upper_mz_, precursor->getPos()[precursor_peak_idx], -1 * charge_dist);
    // search right of precursor for isotopic peaks
    precursor_intensity += sumPotentialIsotopePeaks_(precursor, *isolation_lower_mz_, *isolation_upper_mz_, precursor->getPos()[precursor_peak_idx], charge_dist);

    return precursor_intensity / total_intensity;
  }

  void IsobaricChannelExtractor::extractChannels(const OMSInterfaces::MSRunIF& ms_data, ConsensusMap& consensus_map)
  {
    if (ms_data.getSpectraNr() == 0)
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.\n";
      throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Experiment has no scans!");
    }

    // clear the output map
    consensus_map.clear(false);
    consensus_map.setExperimentType("labeled_MS2");

    // create predicate for spectrum checking
    LOG_INFO << "Selecting scans with activation mode: " << (selected_activation_ == "" ? "any" : selected_activation_) << "\n";

    // now we have picked data
    // --> assign peaks to channels
    UInt64 element_index(0);

    // remember the current precusor spectrum
    int prec_spec = -1;

    for (Size i = 0; i < ms_data.getSpectraNr(); i++)
    {
      const OpenMS::OMSInterfaces::SpectrumPtr it = ms_data.get_Spectrum(i);

      // remember the last MS1 spectra as we assume it to be the precursor spectrum
      if (ms_data.get_Spectrum(i)->getMSLevel() ==  1)
      {
        // remember potential precursor and continue
        prec_spec = i;
        continue;
      }

      if (selected_activation_ == "" || hasActivationMethod(it, selected_activation_))
      {
        // check if precursor is available
        if (it->getPrecursors().empty())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("No precursor information given for scan native ID ") + it->getNativeID() + " with RT " + String(it->getRT()));
        }

        // check precursor constraints
        if (!isValidPrecursor_(it->getPrecursors()[0]))
        {
          LOG_DEBUG << "Skip spectrum " << it->getNativeID() << ": Precursor doesn't fulfill all constraints." << std::endl;
          continue;
        }

        // check precursor purity if we have a valid precursor ..
        DoubleReal precursor_purity = -1.0;
        if (prec_spec != -1)
        {
          precursor_purity = computePrecursorPurity_(it, ms_data.get_Spectrum(prec_spec));
          // check if purity is high enough
          if (precursor_purity < min_precursor_purity_)
          {
            LOG_DEBUG << "Skip spectrum " << it->getNativeID() << ": Precursor purity is below the threshold. [purity = " << precursor_purity << "]" << std::endl;
            continue;
          }
        }
        else
        {
          LOG_INFO << "No precursor available for spectrum: " << it->getNativeID() << std::endl;
        }

        // store RT&MZ of parent ion as centroid of ConsensusFeature
        ConsensusFeature cf;
        cf.setUniqueId();
        cf.setRT(it->getRT());
        cf.setMZ(it->getPrecursors()[0].mz);

        Peak2D channel_value;
        channel_value.setRT(it->getRT());
        // for each each channel
        UInt64 map_index = 0;
        Peak2D::IntensityType overall_intensity = 0;
        for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator cl_it = quant_method_->getChannelInformation().begin();
             cl_it != quant_method_->getChannelInformation().end();
             ++cl_it)
        {
          // set mz-position of channel
          channel_value.setMZ(cl_it->center);
          // reset intensity
          channel_value.setIntensity(0);

          // compute the iterators
          std::vector<double>::const_iterator mz_end_ = std::upper_bound(it->getPos().begin(), it->getPos().end(), cl_it->center + reporter_mass_shift_);
          std::vector<double>::const_iterator mz_it = std::lower_bound(it->getPos().begin(), it->getPos().end(), cl_it->center - reporter_mass_shift_);
          std::vector<double>::const_iterator int_it = it->getIntensity().begin();
          DiffType iterator_pos = std::distance((std::vector<double>::const_iterator)it->getPos().begin(), mz_it);
          std::advance(int_it, iterator_pos);

          // add up all signals
          for (; mz_it != mz_end_; ++mz_it, ++int_it)
          {
            channel_value.setIntensity(channel_value.getIntensity() + *int_it);
          }

          // discard contribution of this channel as it is below the required intensity threshold
          if (channel_value.getIntensity() < min_reporter_intensity_)
          {
            channel_value.setIntensity(0);
          }

          overall_intensity += channel_value.getIntensity();
          // add channel to ConsensusFeature
          cf.insert(map_index++, channel_value, element_index);
        } // ! channel_iterator

        // check if we keep this feature or if it contains low-intensity quantifications
        if (remove_low_intensity_quantifications_ && hasLowIntensityReporter_(cf))
        {
          continue;
        }

        // check featureHandles are not empty
        if (overall_intensity == 0)
        {
          cf.setMetaValue("all_empty", String("true"));
        }
        // add purity information if we could compute it
        if (precursor_purity != -1.0)
        {
          cf.setMetaValue("precursor_purity", precursor_purity);
        }

        // embed the id of the scan from which the quantitative information was extracted
        cf.setMetaValue("scan_id", it->getNativeID());
        // .. as well as additional meta information
        cf.setMetaValue("precursor_intensity", it->getPrecursors()[0].intensity);
        cf.setMetaValue("precursor_charge", it->getPrecursors()[0].charge);
        
        cf.setIntensity(overall_intensity);
        consensus_map.push_back(cf);

        // the tandem-scan in the order they appear in the experiment
        ++element_index;
      }
    } // ! Experiment iterator

    /// add meta information to the map
    registerChannelsInOutputMap_(consensus_map);
  }

  void IsobaricChannelExtractor::registerChannelsInOutputMap_(ConsensusMap& consensus_map)
  {
    // register the individual channels in the output consensus map
    Int index = 0;
    for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator cl_it = quant_method_->getChannelInformation().begin();
         cl_it != quant_method_->getChannelInformation().end();
         ++cl_it)
    {
      ConsensusMap::FileDescription channel_as_map;
      // label is the channel + description provided in the Params
      channel_as_map.label = quant_method_->getName() + "_" + cl_it->name;

      // TODO(aiche): number of features need to be set later
      channel_as_map.size = consensus_map.size();

      // add some more MetaInfo
      channel_as_map.setMetaValue("channel_name", cl_it->name);
      channel_as_map.setMetaValue("channel_id", cl_it->id);
      channel_as_map.setMetaValue("channel_description", cl_it->description);
      channel_as_map.setMetaValue("channel_center", cl_it->center);
      consensus_map.getFileDescriptions()[index++] = channel_as_map;
    }
  }

} // namespace
