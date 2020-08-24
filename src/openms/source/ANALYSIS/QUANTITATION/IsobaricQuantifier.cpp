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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricNormalizer.h>

namespace OpenMS
{
  IsobaricQuantifier::IsobaricQuantifier(const IsobaricQuantitationMethod* const quant_method) :
    DefaultParamHandler("IsobaricQuantifier"),
    quant_method_(quant_method)
  {
    setDefaultParams_();
  }

  IsobaricQuantifier::IsobaricQuantifier(const IsobaricQuantifier& other) :
    DefaultParamHandler(other),
    quant_method_(other.quant_method_)
  {
  }

  IsobaricQuantifier& IsobaricQuantifier::operator=(const IsobaricQuantifier& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    quant_method_ = rhs.quant_method_;

    return *this;
  }

  void IsobaricQuantifier::setDefaultParams_()
  {
    defaults_.setValue("isotope_correction", "true", "Enable isotope correction (highly recommended). "
                                                     "Note that you need to provide a correct isotope correction matrix "
                                                     "otherwise the tool will fail or produce invalid results.");
    defaults_.setValidStrings("isotope_correction", ListUtils::create<String>("true,false"));

    defaults_.setValue("normalization", "false", "Enable normalization of channel intensities with respect to the reference channel. "
                                                 "The normalization is done by using the Median of Ratios (every channel / Reference). "
                                                 "Also the ratio of medians (from any channel and reference) is provided as control measure!");
    defaults_.setValidStrings("normalization", ListUtils::create<String>("true,false"));

    defaultsToParam_();
  }

  void IsobaricQuantifier::updateMembers_()
  {
    isotope_correction_enabled_ = getParameters().getValue("isotope_correction") == "true";
    normalization_enabled_ = getParameters().getValue("normalization") == "true";
  }

  void IsobaricQuantifier::quantify(const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out)
  {
    // precheck incoming map
    if (consensus_map_in.empty())
    {
      OPENMS_LOG_WARN << "Warning: Empty iTRAQ/TMT container. No quantitative information available!" << std::endl;
      return;
    }

    // create output map based on input, we will cleanup the channels while iterating over it
    consensus_map_out = consensus_map_in;

    // init stats
    stats_.reset();
    stats_.channel_count = quant_method_->getNumberOfChannels();

    // apply isotope correction if requested by user
    if (isotope_correction_enabled_)
    {
      stats_ = IsobaricIsotopeCorrector::correctIsotopicImpurities(consensus_map_in, consensus_map_out, quant_method_);
    }
    else
    {
      OPENMS_LOG_WARN << "Warning: Due to deactivated isotope-correction labeling statistics will be based on raw intensities, which might give too optimistic results." << std::endl;
    }

    // compute statistics and embed into output map
    computeLabelingStatistics_(consensus_map_out);

    // apply normalization if requested
    if (normalization_enabled_)
    {
      IsobaricNormalizer normalizer(quant_method_);
      normalizer.normalize(consensus_map_out);
    }
  }

  void IsobaricQuantifier::computeLabelingStatistics_(ConsensusMap& consensus_map_out)
  {
    // number of total quantified spectra
    stats_.number_ms2_total = consensus_map_out.size();

    // Labeling efficiency statistics
    for (size_t i = 0; i < consensus_map_out.size(); ++i)
    {
      // is whole scan empty?!
      if (consensus_map_out[i].getIntensity() == 0) ++stats_.number_ms2_empty;

      // look at single reporters
      for (ConsensusFeature::HandleSetType::const_iterator it_elements = consensus_map_out[i].begin();
           it_elements != consensus_map_out[i].end();
           ++it_elements)
      {
        if (it_elements->getIntensity() == 0)
        {
          String ch_index = consensus_map_out.getColumnHeaders()[it_elements->getMapIndex()].getMetaValue("channel_name");
          ++stats_.empty_channels[ch_index];
        }
      }
    }
    OPENMS_LOG_INFO << "IsobaricQuantifier: skipped " << stats_.number_ms2_empty << " of " << consensus_map_out.size() << " selected scans due to lack of reporter information:\n";
    consensus_map_out.setMetaValue("isoquant:scans_noquant", stats_.number_ms2_empty);
    consensus_map_out.setMetaValue("isoquant:scans_total", consensus_map_out.size());

    OPENMS_LOG_INFO << "IsobaricQuantifier: channels with signal\n";
    for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator cl_it = quant_method_->getChannelInformation().begin();
      cl_it != quant_method_->getChannelInformation().end();
      ++cl_it) // use the same iteration method for printing stats as in IsobaricChannelExtractor which have the same order, so user can make 1:1 comparison
    {
      std::map<String, Size>::const_iterator it_m = stats_.empty_channels.find(cl_it->name);
      if (it_m == stats_.empty_channels.end()) 
      { // should not happen
        OPENMS_LOG_WARN << "Warning: no stats for channel '" << cl_it->name << "'" << std::endl;
        continue;
      }
      OPENMS_LOG_INFO << "  ch " << String(cl_it->name).fillRight(' ', 4) << ": " << (consensus_map_out.size() - it_m->second) << " / " << consensus_map_out.size() << " (" << ((consensus_map_out.size() - it_m->second) * 100 / consensus_map_out.size()) << "%)\n";
      consensus_map_out.setMetaValue(String("isoquant:quantifyable_ch") + it_m->first, (consensus_map_out.size() - it_m->second));
    }
  }

} // namespace
