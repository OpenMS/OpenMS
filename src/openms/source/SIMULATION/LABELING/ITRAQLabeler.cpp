// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/LABELING/ITRAQLabeler.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <Eigen/Dense>

#include <boost/random/uniform_real.hpp>

using std::vector;
using std::pair;
using std::set;

namespace OpenMS
{
  ITRAQLabeler::ITRAQLabeler() :
    BaseLabeler(),
    itraq_type_(),
    channel_map_(),
    isotope_corrections_()
  {
    setName("ITRAQLabeler");
    channel_description_ = "iTRAQ labeling on MS2 level with up to 4 (4plex) or 8 (8plex) channels.";


    // this needs to come first!
    isotope_corrections_.resize(2);
    isotope_corrections_[0].setMatrix<4, 4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
    isotope_corrections_[1].setMatrix<8, 4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);

    // iTRAQ
    defaults_.setValue("iTRAQ", "4plex", "4plex or 8plex iTRAQ?");
    defaults_.setValidStrings("iTRAQ", ListUtils::create<String>("4plex,8plex"));

    defaults_.setValue("reporter_mass_shift", 0.1, "Allowed shift (uniformly distributed - left to right) in Da from the expected position (of e.g. 114.1, 115.1)");
    defaults_.setMinFloat("reporter_mass_shift", 0);
    defaults_.setMaxFloat("reporter_mass_shift", 0.5);

    defaults_.setValue("channel_active_4plex", ListUtils::create<String>("114:myReference"), "Four-plex only: Each channel that was used in the experiment and its description (114-117) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\".");
    defaults_.setValue("channel_active_8plex", ListUtils::create<String>("113:myReference"), "Eight-plex only: Each channel that was used in the experiment and its description (113-121) in format <channel>:<name>, e.g. \"113:myref\",\"115:liver\",\"118:lung\".");

    StringList isotopes = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::FOURPLEX, isotope_corrections_);
    defaults_.setValue("isotope_correction_values_4plex", isotopes, "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' , '116:0.1/0.3/3/0.2' ", ListUtils::create<String>("advanced"));
    isotopes = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::EIGHTPLEX, isotope_corrections_);
    defaults_.setValue("isotope_correction_values_8plex", isotopes, "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '113:0/0.3/4/0' , '116:0.1/0.3/3/0.2' ", ListUtils::create<String>("advanced"));

    defaults_.setValue("Y_contamination", 0.3, "Efficiency of labeling tyrosine ('Y') residues. 0=off, 1=full labeling");
    defaults_.setMinFloat("Y_contamination", 0.0);
    defaults_.setMaxFloat("Y_contamination", 1.0);

    defaultsToParam_();
  }

  ITRAQLabeler::~ITRAQLabeler()
  {
  }

  void ITRAQLabeler::updateMembers_()
  {
    StringList channels_active;

    if (param_.getValue("iTRAQ") == "4plex")
    {
      itraq_type_ = ItraqConstants::FOURPLEX;
      channels_active = param_.getValue("channel_active_4plex");
    }
    else if (param_.getValue("iTRAQ") == "8plex")
    {
      itraq_type_ = ItraqConstants::EIGHTPLEX;
      channels_active = param_.getValue("channel_active_8plex");
    }

    ItraqConstants::initChannelMap(itraq_type_, channel_map_);
    ItraqConstants::updateChannelMap(channels_active, channel_map_);


    // update isotope_corrections_ Matrix with custom values
    StringList channels;
    if (itraq_type_ == ItraqConstants::FOURPLEX)
    {
      channels = param_.getValue("isotope_correction_values_4plex");
    }
    else
    {
      channels = param_.getValue("isotope_correction_values_8plex");
    }
    if (channels.size() > 0)
    {
      ItraqConstants::updateIsotopeMatrixFromStringList(itraq_type_, channels, isotope_corrections_);
    }

    y_labeling_efficiency_ = param_.getValue("Y_contamination");

  }

  void ITRAQLabeler::preCheck(Param& param) const
  {
    // check for valid MS/MS method
    if (!ListUtils::contains(ListUtils::create<String>("disabled,precursor"), param.getValue("RawTandemSignal:status")))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "iTRAQ Labeling does not work with the chosen MS/MS type");
    }
  }

  void ITRAQLabeler::setUpHook(SimTypes::FeatureMapSimVector& features)
  {
    // no action here .. just check for correct # of channels
    Size active_channel_count = 0;
    for (ChannelMapType::ConstIterator it = channel_map_.begin(); it != channel_map_.end(); ++it)
    {
      if (it->second.active)
        ++active_channel_count;
    }
    if (features.size() != active_channel_count)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("iTRAQ Labeling received wrong number of channels: ") + String(active_channel_count) + " defined, but " + String(features.size()) + " given as FASTA files.");
    }
  }

  /// Labeling between digestion and rt simulation
  /// Join all peptides with the same sequence into one feature
  /// channels are retained via metavalues
  /// if a peptide is not present in all channels, then there will be missing meta values! (so don't rely on them being present)
  void ITRAQLabeler::postDigestHook(SimTypes::FeatureMapSimVector& channels)
  {
    // merge channels into a single feature map
    SimTypes::FeatureMapSim final_feature_map = mergeProteinIdentificationsMaps_(channels);

    std::map<String, Size> peptide_to_feature;

    for (Size i = 0; i < channels.size(); ++i)
    {
      for (SimTypes::FeatureMapSim::iterator it_f_o = channels[i].begin();
           it_f_o != channels[i].end();
           ++it_f_o)
      {
        // derive iTRAQ labeled features from original sequence (might be more than one due to partial labeling)
        SimTypes::FeatureMapSim labeled_features;
        labelPeptide_(*it_f_o, labeled_features);
        for (SimTypes::FeatureMapSim::iterator it_f = labeled_features.begin();
             it_f != labeled_features.end();
             ++it_f)
        {
          const String& seq = it_f->getPeptideIdentifications()[0].getHits()[0].getSequence().toString();
          Size f_index;
          //check if we already have a feature for this peptide
          if (peptide_to_feature.count(seq) > 0)
          {
            f_index = peptide_to_feature[seq];
          }
          else // create new feature
          {
            final_feature_map.push_back(*it_f);
            // update map:
            f_index = final_feature_map.size() - 1;
            peptide_to_feature[seq] = f_index;
          }
          // add intensity as metavalue
          final_feature_map[f_index].setMetaValue(getChannelIntensityName(i), it_f->getIntensity());
          // increase overall intensity
          final_feature_map[f_index].setIntensity(final_feature_map[f_index].getIntensity() + it_f->getIntensity());
          mergeProteinAccessions_(final_feature_map[f_index], *it_f);
        }
      }
    }

    channels.clear();
    channels.push_back(final_feature_map);
  }

  /// Labeling between RT and Detectability
  void ITRAQLabeler::postRTHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */)
  {
  }

  /// Labeling between Detectability and Ionization
  void ITRAQLabeler::postDetectabilityHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */)
  {
  }

  /// Labeling between Ionization and RawMS
  void ITRAQLabeler::postIonizationHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */)
  {
  }

  /// Labeling after RawMS
  void ITRAQLabeler::postRawMSHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */)
  {
  }

  void ITRAQLabeler::postRawTandemMSHook(SimTypes::FeatureMapSimVector& fm, SimTypes::MSSimExperiment& exp)
  {
    //std::cout << "Matrix used: \n" << ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_) << "\n\n";

    double rep_shift = param_.getValue("reporter_mass_shift");


    OPENMS_PRECONDITION(fm.size() == 1, "More than one feature map given in ITRAQLabeler::postRawTandemMSHook()!")
    EigenMatrixXdPtr channel_frequency = convertOpenMSMatrix2EigenMatrixXd(ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_));
    Eigen::MatrixXd itraq_intensity_sum(ItraqConstants::CHANNEL_COUNT[itraq_type_], 1);

    std::vector<Matrix<Int> > channel_names(2);
    channel_names[0].setMatrix<4, 1>(ItraqConstants::CHANNELS_FOURPLEX);
    channel_names[1].setMatrix<8, 1>(ItraqConstants::CHANNELS_EIGHTPLEX);

    boost::uniform_real<double> udist(0.0, 1.0);

    // add signal...
    for (SimTypes::MSSimExperiment::iterator it = exp.begin(); it != exp.end(); ++it)
    {
      if (it->getMSLevel() != 2)
        continue;

      // reset sum matrix to 0
      itraq_intensity_sum.setZero();

      // add up signal of all features
      OPENMS_PRECONDITION(it->metaValueExists("parent_feature_ids"), "Meta value 'parent_feature_ids' missing in ITRAQLabeler::postRawTandemMSHook()!")
      IntList parent_fs = it->getMetaValue("parent_feature_ids");
      for (Size i_f = 0; i_f < parent_fs.size(); ++i_f)
      {
        // get RT scaled iTRAQ intensities
        EigenMatrixXdPtr row = getItraqIntensity_(fm[0][parent_fs[i_f]], it->getRT());

        // apply isotope matrix to active channels
        // row * channel_frequencyOld = observed iTRAQ intensities
        Eigen::MatrixXd itraq_intensity_observed = (*channel_frequency) * (*row);
        // add result to sum
        itraq_intensity_sum += itraq_intensity_observed;
      }

      // add signal to MS2 spectrum
      for (Int i_channel = 0; i_channel < ItraqConstants::CHANNEL_COUNT[itraq_type_]; ++i_channel)
      {
        SimTypes::MSSimExperiment::SpectrumType::PeakType p;
        // random shift of +-rep_shift around exact position
        double rnd_shift = udist(rng_->getTechnicalRng()) * 2 * rep_shift - rep_shift;
        p.setMZ(channel_names[itraq_type_].getValue(i_channel, 0) + 0.1 + rnd_shift);
        p.setIntensity(itraq_intensity_sum(i_channel, 0));
        it->push_back(p);
      }
    }

  }

  // CUSTOM FUNCTIONS for iTRAQ:: //

  void ITRAQLabeler::addModificationToPeptideHit_(Feature& feature, const String& modification, const Size& pos) const
  {
    vector<PeptideHit> pep_hits(feature.getPeptideIdentifications()[0].getHits());
    AASequence modified_sequence(pep_hits[0].getSequence());
    modified_sequence.setModification(pos, modification);
    pep_hits[0].setSequence(modified_sequence);
    feature.getPeptideIdentifications()[0].setHits(pep_hits);
  }

  void ITRAQLabeler::labelPeptide_(const Feature& feature, SimTypes::FeatureMapSim& result) const
  {
    // modify with iTRAQ modification (needed for mass calculation and MS/MS signal)
    //site="Y" - low abundance
    //site="N-term"
    //site="K" - lysine
    String modification = (itraq_type_ == ItraqConstants::FOURPLEX ? "iTRAQ4plex" : "iTRAQ8plex");
    vector<PeptideHit> pep_hits(feature.getPeptideIdentifications()[0].getHits());
    AASequence seq(pep_hits[0].getSequence());
    // N-term
    seq.setNTerminalModification(modification);
    // all "K":
    for (Size i = 0; i < seq.size(); ++i)
    {
      if (seq[i] == 'K' && !seq[i].isModified())
        seq.setModification(i, modification);
    }
    result.resize(1);
    result[0] = feature;
    pep_hits[0].setSequence(seq);
    result[0].getPeptideIdentifications()[0].setHits(pep_hits);
    // some "Y":
    // for each "Y" create two new features, depending on labeling efficiency on "Y":
    if (y_labeling_efficiency_ == 0)
      return;

    for (Size i = 0; i < seq.size(); ++i)
    {
      if (seq[i] == 'Y' && !seq[i].isModified())
      {
        if (y_labeling_efficiency_ == 1)
        {
          addModificationToPeptideHit_(result.back(), modification, i);
        }
        else // double number of features:
        {
          Size f_count = result.size();
          for (Size f = 0; f < f_count; ++f)
          {
            // copy feature
            result.push_back(result[f]);
            // modify the copy
            addModificationToPeptideHit_(result.back(), modification, i);
            // adjust intensities:
            result.back().setIntensity(result.back().getIntensity() * y_labeling_efficiency_);
            result[f].setIntensity(result[f].getIntensity() * (1 - y_labeling_efficiency_));
          }
        }
      }
    }


  }

  double ITRAQLabeler::getRTProfileIntensity_(const Feature& f, const double MS2_RT_time) const
  {
    // compute intensity correction factor for feature from RT profile
    const DoubleList& elution_bounds = f.getMetaValue("elution_profile_bounds");
    const DoubleList& elution_ints   = f.getMetaValue("elution_profile_intensities");

    // check that RT is within the elution bound:
    OPENMS_POSTCONDITION(f.getConvexHull().getBoundingBox().encloses(MS2_RT_time, f.getMZ()), "The MS2 spectrum has wrong parent features! The feature does not touch the spectrum's RT!")

    if (MS2_RT_time < elution_bounds[1] || elution_bounds[3] < MS2_RT_time)
    {
      OPENMS_LOG_WARN << "Warn: requesting MS2 RT for " << MS2_RT_time << ", but bounds are only from [" << elution_bounds[1] << "," << elution_bounds[3] << "]\n";
      return 0;
    }

    // do linear interpolation
    double width = elution_bounds[3] - elution_bounds[1];
    double offset = MS2_RT_time - elution_bounds[1];
    Int index = floor(offset / (width / (elution_ints.size() - 1)) + 0.5);

    OPENMS_POSTCONDITION(index < (Int)elution_ints.size(), "Wrong index computation! (Too large)")

    return elution_ints[index];
  }

  EigenMatrixXdPtr ITRAQLabeler::getItraqIntensity_(const Feature& f, const double MS2_RT_time) const
  {

    double factor = getRTProfileIntensity_(f, MS2_RT_time);

    //std::cerr << "\n\nfactor is: " << factor << "\n";
    // fill map with values present (all missing ones remain 0)
    MutableEigenMatrixXdPtr m(new Eigen::MatrixXd(ItraqConstants::CHANNEL_COUNT[itraq_type_], 1));
    m->setZero();
    Size ch(0);
    Size ch_internal(0);
    for (ChannelMapType::ConstIterator it = channel_map_.begin(); it != channel_map_.end(); ++it)
    {
      SimTypes::SimIntensityType intensity(0);
      if (it->second.active && f.metaValueExists(getChannelIntensityName(ch_internal))) // peptide is present in this channel
      {
        intensity = (double) f.getMetaValue(getChannelIntensityName(ch_internal));
        ++ch_internal;
      }
      (* m)(ch, 0) = intensity * factor;
      ++ch;
    }

    return m;
  }

} // namespace OpenMS
