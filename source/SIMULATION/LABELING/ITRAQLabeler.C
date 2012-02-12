// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/LABELING/ITRAQLabeler.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <gsl/gsl_blas.h>

using std::vector;
using std::pair;
using std::set;

namespace OpenMS
{
  ITRAQLabeler::ITRAQLabeler()
    : BaseLabeler(),
		  itraq_type_(),
		  channel_map_(),
		  isotope_corrections_()

  {
    setName("ITRAQLabeler");

		// this needs to come first!
		isotope_corrections_.resize(2);
		isotope_corrections_[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
		isotope_corrections_[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);

		// iTRAQ
		defaults_.setValue("iTRAQ", "4plex", "4plex or 8plex iTRAQ?");
		defaults_.setValidStrings("iTRAQ", StringList::create("4plex,8plex"));

		defaults_.setValue("reporter_mass_shift", 0.1, "Allowed shift (uniformly distributed - left to right) in Da from the expected position (of e.g. 114.1, 115.1)"); 
		defaults_.setMinFloat ("reporter_mass_shift", 0);
		defaults_.setMaxFloat ("reporter_mass_shift", 0.5);

		defaults_.setValue("channel_active_4plex", StringList::create("114:myReference"), "Four-plex only: Each channel that was used in the experiment and its description (114-117) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\"."); 
		defaults_.setValue("channel_active_8plex", StringList::create("113:myReference"), "Eight-plex only: Each channel that was used in the experiment and its description (113-121) in format <channel>:<name>, e.g. \"113:myref\",\"115:liver\",\"118:lung\"."); 

		StringList isotopes = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::FOURPLEX, isotope_corrections_);
		defaults_.setValue("isotope_correction_values_4plex", isotopes, "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' , '116:0.1/0.3/3/0.2' ", StringList::create("advanced"));
    isotopes = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::EIGHTPLEX, isotope_corrections_);
    defaults_.setValue("isotope_correction_values_8plex", isotopes, "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '113:0/0.3/4/0' , '116:0.1/0.3/3/0.2' ", StringList::create("advanced"));

    defaults_.setValue("Y_contamination", 0.3, "Efficiency of labeling tyrosine ('Y') residues. 0=off, 1=full labeling"); 
		defaults_.setMinFloat ("Y_contamination", 0.0);
		defaults_.setMaxFloat ("Y_contamination", 1.0);
		
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
		if (channels.size()>0)
		{
			ItraqConstants::updateIsotopeMatrixFromStringList(itraq_type_, channels, isotope_corrections_);
		}

    y_labeling_efficiency_ = param_.getValue("Y_contamination");

	}
  void ITRAQLabeler::preCheck(Param & param) const
  {
    // check for valid MS/MS method
    if(!StringList::create("disabled,precursor").contains( param.getValue("RawTandemSignal:status")))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "iTRAQ Labeling does not work with the chosen MS/MS type");
    }
  }

  void ITRAQLabeler::setUpHook(FeatureMapSimVector & features)
  {
    // no action here .. just check for correct # of channels
    Size active_channel_count=0;
    for (ChannelMapType::ConstIterator it=channel_map_.begin();it!=channel_map_.end();++it) 
    {
      if (it->second.active) ++active_channel_count;
    }
    if(features.size() != active_channel_count)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("iTRAQ Labeling received wrong number of channels: ") + String(active_channel_count) + " defined, but " + String(features.size()) + " given as FASTA files.");
    }
  }

  /// Labeling between digestion and rt simulation
  /// Join all peptides with the same sequence into one feature
  /// channels are retained via metavalues
  /// if a peptide is not present in all channels, then there will be missing meta values! (so don't rely on them being present)
  void ITRAQLabeler::postDigestHook(FeatureMapSimVector & channels)
  {
    // merge channels into a single feature map
    FeatureMapSim final_feature_map = mergeProteinIdentificationsMaps_(channels);

    std::map<String, Size> peptide_to_feature;

    for (Size i=0; i<channels.size(); ++i)
    {
      for(FeatureMapSim::iterator it_f_o = channels[i].begin() ;
          it_f_o != channels[i].end() ;
          ++it_f_o)
      {
        // derive iTRAQ labeled features from original sequence (might be more than one due to partial labeling)
        FeatureMapSim labeled_features;
        labelPeptide_(*it_f_o, labeled_features);
        for(FeatureMapSim::iterator it_f = labeled_features.begin() ;
            it_f != labeled_features.end() ;
            ++it_f)
        {
          const String & seq = it_f->getPeptideIdentifications()[0].getHits()[0].getSequence().toString ();
          Size f_index;
          //check if we already have a feature for this peptide
          if (peptide_to_feature.count(seq)>0)
          {
            f_index = peptide_to_feature[seq];
          }
          else
          { // create new feature
            final_feature_map.push_back(*it_f);
            // update map:
            f_index=final_feature_map.size()-1;
            peptide_to_feature[seq]=f_index;
          }
          // add intensity as metavalue
          final_feature_map[f_index].setMetaValue(getChannelIntensityName(i), it_f->getIntensity());
          // increase overall intensity
          final_feature_map[f_index].setIntensity( final_feature_map[f_index].getIntensity() + it_f->getIntensity());
          mergeProteinAccessions_(final_feature_map[f_index], *it_f);
        }
      }
    }

    channels.clear();
    channels.push_back(final_feature_map);
  }

  /// Labeling between RT and Detectability
  void ITRAQLabeler::postRTHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Detectability and Ionization
  void ITRAQLabeler::postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Ionization and RawMS
  void ITRAQLabeler::postIonizationHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling after RawMS
  void ITRAQLabeler::postRawMSHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  void ITRAQLabeler::postRawTandemMSHook(FeatureMapSimVector & fm, MSSimExperiment & exp)
  {
    //std::cout << "Matrix used: \n" << ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_) << "\n\n";
		
    DoubleReal rep_shift = param_.getValue("reporter_mass_shift");

    OPENMS_PRECONDITION(fm.size()==1, "More than one feature map given in ITRAQLabeler::postRawTandemMSHook()!")
		gsl_matrix* channel_frequency = ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_).toGslMatrix();
		gsl_matrix* itraq_intensity_observed = Matrix<SimIntensityType>(ItraqConstants::CHANNEL_COUNT[itraq_type_],1).toGslMatrix();
		gsl_matrix* itraq_intensity_sum = Matrix<SimIntensityType>(ItraqConstants::CHANNEL_COUNT[itraq_type_],1).toGslMatrix();
		
		std::vector< Matrix<Int> > channel_names(2);
		channel_names[0].setMatrix<4,1>(ItraqConstants::CHANNELS_FOURPLEX);
		channel_names[1].setMatrix<8,1>(ItraqConstants::CHANNELS_EIGHTPLEX);

		// add signal...
		for (MSSimExperiment::iterator it=exp.begin(); it!=exp.end(); ++it)
		{
      if (it->getMSLevel()!=2) continue;

			// reset sum matrix to 0
			gsl_matrix_scale (itraq_intensity_sum, 0);
			
			// add up signal of all features
      OPENMS_PRECONDITION(it->metaValueExists("parent_feature_ids"),"Meta value 'parent_feature_ids' missing in ITRAQLabeler::postRawTandemMSHook()!")
			IntList parent_fs = (IntList) it->getMetaValue("parent_feature_ids");
			for (Size i_f=0; i_f < parent_fs.size(); ++i_f)
			{
        // get RT scaled iTRAQ intensities
        gsl_matrix* row = getItraqIntensity_(fm[0][parent_fs[i_f]], it->getRT()).toGslMatrix();
				// apply isotope matrix to active channels
				// row * channel_frequency = observed iTRAQ intensities
				gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
												1.0, channel_frequency, row,
												0.0, itraq_intensity_observed);
        // add result to sum
        gsl_matrix_add (itraq_intensity_sum, itraq_intensity_observed);
				gsl_matrix_free (row);
			}
			
			// add signal to MS2 spectrum
			for (Int i_channel=0; i_channel< ItraqConstants::CHANNEL_COUNT[itraq_type_]; ++i_channel)
			{
				MSSimExperiment::SpectrumType::PeakType p;
				// random shift of +-rep_shift around exact position
        DoubleReal rnd_shift = gsl_rng_uniform(rng_->technical_rng) * 2 * rep_shift - rep_shift;
				p.setMZ(channel_names[itraq_type_].getValue(i_channel,0) + 0.1 + rnd_shift);
				p.setIntensity(gsl_matrix_get(itraq_intensity_sum, i_channel, 0));
				//std::cout << "inserted iTRAQ peak: " << p << "\n";
				it->push_back(p);
			}
		}
		
		gsl_matrix_free (channel_frequency);
		gsl_matrix_free (itraq_intensity_observed);
		gsl_matrix_free (itraq_intensity_sum);
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

  void ITRAQLabeler::labelPeptide_(const Feature& feature, FeatureMapSim& result) const
  {
    // modify with iTRAQ modification (needed for mass calculation and MS/MS signal)
    //site="Y" - low abundance
    //site="N-term"
    //site="K" - lysine
    String modification = (itraq_type_==ItraqConstants::FOURPLEX ? "iTRAQ4plex" : "iTRAQ8plex");
    vector<PeptideHit> pep_hits(feature.getPeptideIdentifications()[0].getHits());
    AASequence seq(pep_hits[0].getSequence());
    // N-term
    seq.setNTerminalModification(modification);
    // all "K":
    for (Size i=0;i<seq.size();++i)
    {
      if (seq[i]=='K' && !seq.isModified(i)) seq.setModification(i,modification);
    }
    result.resize(1);
    result[0] = feature;
    pep_hits[0].setSequence(seq);
    result[0].getPeptideIdentifications()[0].setHits(pep_hits);
    // some "Y":
    // for each "Y" create two new features, depending on labeling efficiency on "Y":
    if (y_labeling_efficiency_==0) return;

    for (Size i=0;i<seq.size();++i)
    {
      if ( seq[i]=='Y' && !seq.isModified(i))
      { 
        if (y_labeling_efficiency_==1)
        {
          addModificationToPeptideHit_(result.back(), modification, i);
        }
        else
        { // double number of features:
          Size f_count=result.size();
          for (Size f=0;f<f_count;++f)
          {
            // copy feature
            result.push_back(result[f]);
            // modify the copy
            addModificationToPeptideHit_(result.back(), modification, i);
            // adjust intensities:
            result.back().setIntensity(result.back().getIntensity() * y_labeling_efficiency_);
            result[f].setIntensity(result[f].getIntensity() * (1-y_labeling_efficiency_));
         }
        }
      }
    }
    

  }

  DoubleReal ITRAQLabeler::getRTProfileIntensity_(const Feature & f, const DoubleReal MS2_RT_time) const
  {
    // compute intensity correction factor for feature from RT profile
    const DoubleList& elution_bounds = f.getMetaValue("elution_profile_bounds");
    const DoubleList& elution_ints   = f.getMetaValue("elution_profile_intensities");

    // check that RT is within the elution bound:
    OPENMS_POSTCONDITION(f.getConvexHull().getBoundingBox().encloses(MS2_RT_time, f.getMZ()), "The MS2 spectrum has wrong parent features! The feature does not touch the spectrum's RT!")

    if (MS2_RT_time < elution_bounds[1] || elution_bounds[3] < MS2_RT_time)
  {
      LOG_WARN << "Warn: requesting MS2 RT for " << MS2_RT_time << ", but bounds are only from [" <<elution_bounds[1] << "," << elution_bounds[3] << "]\n";
      return 0;
    }

    // do linear interpolation
    DoubleReal width = elution_bounds[3] - elution_bounds[1];
    DoubleReal offset = MS2_RT_time - elution_bounds[1];
    Int index = floor(offset / (width / (elution_ints.size()-1)) + 0.5);

    OPENMS_POSTCONDITION(index < (Int)elution_ints.size(), "Wrong index computation! (Too large)")

    return elution_ints[index];
  }

  Matrix<SimIntensityType> ITRAQLabeler::getItraqIntensity_(const Feature & f, const DoubleReal MS2_RT_time) const
  {

    DoubleReal factor = getRTProfileIntensity_(f, MS2_RT_time);

    //std::cerr << "\n\nfactor is: " << factor << "\n";
		// fill map with values present (all missing ones remain 0)
		Matrix<SimIntensityType> m(ItraqConstants::CHANNEL_COUNT[itraq_type_], 1, 0);
    Size ch(0);
    Size ch_internal(0);
		for (ChannelMapType::ConstIterator it=channel_map_.begin();it!=channel_map_.end();++it) 
    {
      SimIntensityType intensity(0);
      if (it->second.active && f.metaValueExists(getChannelIntensityName(ch_internal)))
      { // peptide is present in this channel
        intensity = (DoubleReal) f.getMetaValue(getChannelIntensityName(ch_internal));
        ++ch_internal;
      }
      m.setValue(ch, 0, intensity * factor);
      ++ch;
    }

    return m;
  }

} // namespace OpenMS
