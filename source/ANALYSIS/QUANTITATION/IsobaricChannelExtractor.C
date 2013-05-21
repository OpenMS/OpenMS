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

#include <OpenMS/KERNEL/RangeUtils.h>

namespace OpenMS
{

  IsobaricChannelExtractor::IsobaricChannelExtractor(const IsobaricQuantitationMethod* const quant_method) :
    DefaultParamHandler("IsobaricChannelExtractor"),
    quant_method_(quant_method),
    selected_activation_(""),
    reporter_mass_shift_(0.1),
    min_precursor_intensity_(1.0),
    keep_unannotated_precursor_(true),
    min_reporter_intensity_(0.0),
    remove_low_intensity_quantifications_(false)
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
    remove_low_intensity_quantifications_(other.remove_low_intensity_quantifications_)
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

    return *this;
  }

  void IsobaricChannelExtractor::setDefaultParams_()
  {
    defaults_.setValue("select_activation", Precursor::NamesOfActivationMethod[Precursor::HCID], "Operate only on MSn scans where any of its precursors features a certain activation method (e.g., usually HCD for iTRAQ). Set to empty string if you want to disable filtering.");
    StringList activation_list(std::vector<std::string>(Precursor::NamesOfActivationMethod, &Precursor::NamesOfActivationMethod[Precursor::SIZE_OF_ACTIVATIONMETHOD - 1]));
    activation_list.push_back(""); // allow disabling this
    defaults_.setValidStrings("select_activation", activation_list);

    defaults_.setValue("reporter_mass_shift", 0.1, "Allowed shift (left to right) in Da from the expected position.");
    defaults_.setMinFloat("reporter_mass_shift", 0.00000001);
    defaults_.setMaxFloat("reporter_mass_shift", 0.5);

    defaults_.setValue("min_precursor_intensity", 1.0, "Minimum intensity of the precursor to be extracted. MS/MS scans having a precursor with a lower intensity will not be considered for quantitation.");
    defaults_.setMinFloat("min_precursor_intensity", 0.0);

    defaults_.setValue("keep_unannotated_precursor", "true", "Flag if precursor with missing intensity value or missing precursor spectrum should be included or not.");
    defaults_.setValidStrings("keep_unannotated_precursor", StringList::create("true,false"));

    defaults_.setValue("min_reporter_intensity", 0.0, "Minimum intenesity of the individual reporter ions to be used extracted.");
    defaults_.setMinFloat("min_reporter_intensity", 0.0);

    defaults_.setValue("discard_low_intensity_quantifications", "false", "Remove all reporter intensities if a single reporter is below the threshold given in min_reporter_intensity.");
    defaults_.setValidStrings("discard_low_intensity_quantifications", StringList::create("true,false"));

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
  }

  bool IsobaricChannelExtractor::isValidPrecursor_(const Precursor& precursor) const
  {
    return (precursor.getIntensity() == 0.0 && keep_unannotated_precursor_) || !(precursor.getIntensity() < min_precursor_intensity_);
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

  void IsobaricChannelExtractor::extractChannels(const MSExperiment<Peak1D>& ms_exp_data, ConsensusMap& consensus_map)
  {
    if (ms_exp_data.empty())
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
    HasActivationMethod<MSExperiment<Peak1D>::SpectrumType> activation_predicate(StringList::create(selected_activation_));

    // now we have picked data
    // --> assign peaks to channels
    UInt64 element_index(0);

    for (MSExperiment<>::ConstIterator it = ms_exp_data.begin(); it != ms_exp_data.end(); ++it)
    {
      if (selected_activation_ == "" || activation_predicate(*it))
      {
        // check if precursor is available
        if (it->getPrecursors().empty())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("No precursor information given for scan native ID ") + String(it->getNativeID()) + " with RT " + String(it->getRT()));
        }

        // check precursor constraints
        if (!isValidPrecursor_(it->getPrecursors()[0]))
        {
          LOG_DEBUG << "Skip spectrum " << String(it->getNativeID()) << ": Precursor doesn't fulfill all constraints." << std::endl;
          continue;
        }

        // store RT&MZ of parent ion as centroid of ConsensusFeature
        ConsensusFeature cf;
        cf.setUniqueId();
        cf.setRT(it->getRT());
        cf.setMZ(it->getPrecursors()[0].getMZ());

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

          // add up all signals
          for (MSExperiment<>::SpectrumType::ConstIterator mz_it = it->MZBegin(cl_it->center - reporter_mass_shift_);
               mz_it != it->MZEnd(cl_it->center + reporter_mass_shift_);
               ++mz_it)
          {
            channel_value.setIntensity(channel_value.getIntensity() + mz_it->getIntensity());
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
