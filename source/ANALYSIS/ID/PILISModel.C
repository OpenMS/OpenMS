// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/SYSTEM/File.h>

#include <cmath>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <fstream>


#define PRECURSOR_DEBUG

#define TRAINING_DEBUG
#undef  TRAINING_DEBUG

#define SIM_DEBUG
#undef  SIM_DEBUG

using namespace std;

// TODO New QXP pathway handling
// TODO New loss modeling with P's in the sequence, multiple losses, ...
// TODO New XXD yk-2 enhancement
//
// TODO XPHXXX yk-2 enhancement
// TODO N-terminal Lysine NH3-loss
// 
// TODO if E N-term and charge directed state, abundant loss of water
// the bk-1/bk-2 pathway somtimes show bx+18 ions, especially when H is C-term to it
//
// TODO a-ion generation when suffix has high PA and prefix low PA
// TODO b Q1 loss

namespace OpenMS 
{

	PILISModel::PILISModel()
		: DefaultParamHandler("PILISModel"),
			valid_(false)
	{	
		defaults_.setValue("upper_mz", 2000.0, "Max m/z value of product ions in the simulated spectrum");
		defaults_.setValue("lower_mz", 200.0, "Lowest m/z value of product ions in the simulated spectrum");
		defaults_.setValue("charge_remote_threshold", 0.2, "If the probability for the proton at the N-terminus is lower than this value, enable charge-remote pathways");
		defaults_.setValue("charge_directed_threshold", 0.3, "Limit the proton availability at the N-terminus to at least this value for charge-directed pathways");
		
		//defaults_.setValue("side_chain_intensity_threshold", 0.0, "Side-chain pathways are active if this threshold is exceeded by the probability of N-terminal proton");
		defaults_.setValue("model_depth", 4, "The number of explicitly modeled backbone cleavages from N-terminus and C-terminus, would be 9 for the default value");
		defaults_.setValue("visible_model_depth", 30, "The maximal possible size of a peptide to be modeled");
		defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, used to identify the precursor peak and its loss peaks for training");
		defaults_.setValue("peak_mass_tolerance", 0.3, "Peak mass tolerance of the product ions, used to identify the ions for training");
		defaults_.setValue("fixed_modifications", "", "Fixed modifications in format '57.001@C'");

		defaults_.setValue("min_y_ion_intensity", 0.20, ".");
		defaults_.setValue("min_b_ion_intensity", 0.15, ".");
		defaults_.setValue("min_a_ion_intensity", 0.05, ".");
		defaults_.setValue("min_y_loss_intensity", 0.05, ".");
		defaults_.setValue("min_b_loss_intensity", 0.02, ".");

		defaults_.setValue("charge_loss_factor", 0.2, "Factor which accounts for the loss of a proton from a ++ ion, the higher the value the more ++ ions loose protons");
		defaults_.setValue("pseudo_counts", 1e-15, "Value which is added for every transition trained of the underlying hidden Markov model");
		defaults_.setValue("max_isotope", 2, "Maximal isotope peak which is reported by the model, 0 means all isotope peaks are reported");

		defaultsToParam_();
		initModels_();
	}

	PILISModel::~PILISModel()
	{
	}

	PILISModel::PILISModel(const PILISModel& model)
		: DefaultParamHandler(model),
			hmm_(model.hmm_),
			hmm_precursor_(model.hmm_precursor_),
			prot_dist_(model.prot_dist_),
			tsg_(model.tsg_),
			valid_(model.valid_),
			peaks_(model.peaks_),
			spectra_aligner_(model.spectra_aligner_)
	{
	}

	PILISModel& PILISModel::operator = (const PILISModel& model)
	{
		if (this != &model)
		{
			DefaultParamHandler::operator=(model);
			hmm_ = model.hmm_;
	    hmm_precursor_ = model.hmm_precursor_;
	    prot_dist_ = model.prot_dist_;
	    tsg_ = model.tsg_;
	    valid_ = model.valid_;
			peaks_ = model.peaks_;
			spectra_aligner_ = model.spectra_aligner_;
		}
		return *this;
	}
	
	void PILISModel::readFromFile(const String& filename)
	{
		// read the model
		String new_filename = File::find(filename);
    if (!File::readable(new_filename))
    {
     	throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, new_filename);
    }
    if (File::empty(new_filename))
    {
     	throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, new_filename);
    }

		TextFile file;
		file.load(new_filename, true);

		TextFile::Iterator it_begin(file.begin()), it_end(file.begin());
		it_begin = file.search(it_begin, "BASE_MODEL_BEGIN");
		it_end = file.search(it_begin, "BASE_MODEL_END");
		parseHMMModel_(++it_begin, it_end, hmm_);

/*
		// seek to next interval
		it_begin = file.search(it_end, "PRECURSOR_MODEL_BEGIN");
		it_end = file.search(it_begin, "PRECURSOR_MODEL_END");
		parseHMMModel_(++it_begin, it_end, hmm_pre_loss_);

		// loss models
		it_begin = file.search(it_end, "BION_LOSS_MODEL_BEGIN");
		it_end = file.search(it_begin, "BION_LOSS_MODEL_END");
		parseHMMLightModel_(++it_begin, it_end, hmms_losses_[BIon]);
	
		it_begin = file.search(it_end, "B2ION_LOSS_MODEL_BEGIN");
		it_end = file.search(it_begin, "B2ION_LOSS_MODEL_END");
		parseHMMLightModel_(++it_begin, it_end, hmms_losses_[B2Ion]);
		
		// y-ion loss model
		it_begin = file.search(it_end, "YION_LOSS_MODEL_BEGIN");
		it_end = file.search(it_begin, "YION_LOSS_MODEL_END");
		parseHMMLightModel_(++it_begin, it_end, hmms_losses_[YIon]);
*/
		valid_ = true;
		return;
	}

	void PILISModel::writeGraphMLFile(const String& filename)
	{
		hmm_.writeGraphMLFile(filename);
		return;
	}

	void PILISModel::writeToFile(const String& filename)
	{
		
		#ifdef SIM_DEBUG
		cerr << "writing to file '" << file_name << "'" << endl;
		#endif

		ofstream out(filename.c_str());
		out << "BASE_MODEL_BEGIN" << endl;
		hmm_.write(out);
		out << "BASE_MODEL_END" << endl;

		/*
		out << "PRECURSOR_MODEL_BEGIN" << endl;
		hmm_pre_loss_.write(out);
		out << "PRECURSOR_MODEL_END" << endl;

		out << "BION_LOSS_MODEL_BEGIN" << endl;
		hmms_losses_[BIon].write(out);
		out << "BION_LOSS_MODEL_END" << endl;

		out << "B2ION_LOSS_MODEL_BEGIN" << endl;
		hmms_losses_[B2Ion].write(out);
		out << "B2ION_LOSS_MODEL_END" << endl;
		
		out << "YION_LOSS_MODEL_BEGIN" << endl;
		hmms_losses_[YIon].write(out);
		out << "YION_LOSS_MODEL_END" << endl;
	*/	
		return;
	}

	void PILISModel::train(const RichPeakSpectrum& in_spec, const AASequence& peptide, UInt charge)
	{
		if (!valid_)
		{
			cerr << "PILISModel: cannot train, initialize model from file first, e.g. data/PILIS/PILIS_model_default.dat" << endl;
			return;
		}

		if (peptide.size() >= 40)
		{
			cerr << "PILISModel: cannot train peptide " << peptide << " of length " << peptide.size() << " (max is 39)" << endl;
			return;
		}


		RichPeakSpectrum train_spec = in_spec;
		train_spec.sortByPosition();
		
		#ifdef TRAINING_DEBUG
		cout << "peptide: " << peptide  << "(z=" << charge << ")" << endl;
		#endif
		double peptide_weight((peptide.getMonoWeight() + double(charge)) / double(charge));
		
		Map<String, double> pre_ints;
		getPrecursorIntensitiesFromSpectrum_(train_spec, pre_ints, peptide_weight, charge);

		// get the ions intensities, y and b ions and losses H2O, NH3 respectively
		IonPeaks_ ints_1, ints_2;
		double sum1 = getIntensitiesFromSpectrum_(train_spec, ints_1, peptide, 1);
			
		double sum2(0);
			
		if (charge > 1)
		{
			sum2 = getIntensitiesFromSpectrum_(train_spec, ints_2, peptide, 2);
		}

		// normalize the intensities
		double sum(sum1);

		for (Map<String, double>::ConstIterator it = pre_ints.begin(); it != pre_ints.end(); ++it)
		{
#ifdef PRECURSOR_DEBUG
			cerr << it->first << " " << it->second << endl;
#endif
			sum += it->second;
		}
		
		for (Map<String, double>::Iterator it = pre_ints.begin(); it != pre_ints.end(); ++it)
		{
			it->second /= sum;
		}

		for (Map<IonType_, vector<double> >::Iterator it1 = ints_1.ints.begin(); it1 != ints_1.ints.end(); ++it1)
		{
			for (vector<double>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				*it2 /= sum;
			}
		}

		if (charge > 1)
		{
			for (Map<IonType_, vector<double> >::Iterator it1 = ints_2.ints.begin(); it1 != ints_2.ints.end(); ++it1)
			{
				for (vector<double>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					*it2 /= sum;
				}
			}
		}
	
		if (sum == 0)
		{
			// does not make sense to proceed here
			cerr << "PILISModel::train warning: no peaks found which match the given sequence" << endl;
			return;
		}

		Map<UInt, double> bb_charge_full, sc_charge_full;
		prot_dist_.getProtonDistribution(bb_charge_full, sc_charge_full, peptide, charge, Residue::YIon);
		prot_dist_.setPeptideProtonDistribution(bb_charge_full, sc_charge_full);
	
		vector<double> bb_init, sc_init, cr_init;
		bool is_charge_remote = getInitialTransitionProbabilities_(bb_init, cr_init, sc_init, bb_charge_full, sc_charge_full, peptide);
		double bb_sum(0);
		for (Map<UInt, double>::ConstIterator iit = bb_charge_full.begin(); iit != bb_charge_full.end(); ++iit)
		{
			bb_sum += iit->second;
		}

		// activation of protons sitting at lysine and histidine side chains
		for (UInt i = 0; i != peptide.size(); ++i)
		{
			if (peptide[i].getOneLetterCode() == "H" || peptide[i].getOneLetterCode() == "K")
			{
				bb_sum += 0.5 * sc_charge_full[i]; // TODO
			}
		}

		if  (bb_sum > 1)
		{
			bb_sum = 1;		
			//bb_sum /= charge;
		}

		if (bb_sum < (double)param_.getValue("charge_directed_threshold") * charge)
		{
			bb_sum = (double)param_.getValue("charge_directed_threshold") * charge;
		}

		// clear the main Hidden Markov Model
		hmm_.clearInitialTransitionProbabilities();
		hmm_.clearTrainingEmissionProbabilities();
		
		double charge_sum(0);
		vector<AASequence> prefixes, suffixes;

		vector<double> bb_charges_all;
    for (UInt i = 0; i != bb_charge_full.size(); ++i)
    {
      bb_charges_all.push_back(bb_charge_full[i]);
    }
    sort(bb_charges_all.begin(), bb_charges_all.end());
    double bb_charges_median = bb_charges_all[UInt(bb_charges_all.size()/2.0)];

		
		
		// for each site: 1. set proton distribution, 2. initial training intensities, 3. train the model
		for (UInt i = 0; i != peptide.size() - 1; ++i)
		{
			String pos_name, y_name1, b_name1, a_name1, y_name, b_name;
			String y_name2, b_name2, a_name2, prefix_size(i + 1), suffix_size(peptide.size() - 1 - i);

			UInt suffix_pos(peptide.size() - i  - 1);
			
			if (i < floor(peptide.size()/2))
			{
				pos_name = prefix_size;
				y_name = "yk-"+prefix_size;
				b_name = "b"+prefix_size;
				y_name1 = "yk-"+prefix_size+"+";
				b_name1 = "b"+prefix_size+"+";
				a_name1 = "a"+prefix_size+"+";

        y_name2 = "yk-"+prefix_size+"++";
        b_name2 = "b"+prefix_size+"++";
        a_name2 = "a"+prefix_size+"++";

			}
			else
			{
				pos_name = "k-"+suffix_size;
				y_name = "y"+suffix_size;
				b_name = "bk-"+suffix_size;
				y_name1 = "y"+suffix_size+"+";
				b_name1 = "bk-"+suffix_size+"+";
				a_name1 = "ak-"+suffix_size+"+";

        y_name2 = "y"+suffix_size+"++";
        b_name2 = "bk-"+suffix_size+"++";
        a_name2 = "ak-"+suffix_size+"++";
			}
						
			AASequence prefix(peptide.getPrefix(i + 1)), suffix(peptide.getSuffix(peptide.size() - i - 1));
			//String aa1(peptide[i].getOneLetterCode()), aa2(peptide[i + 1].getOneLetterCode());
			AASequence aa1_seq, aa2_seq;
			aa1_seq += &peptide[i];
			aa2_seq += &peptide[i + 1];
			String aa1(aa1_seq.toString()), aa2(aa2_seq.toString());
			
			// calc PAs and get b/y ratios for bxyz pathway
			double bint1(0), yint1(0), bint2(0), yint2(0), b_sc_int1(0), b_sc_int2(0), y_sc_int1(0), y_sc_int2(0);
			prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, bint1, yint1, bint2, yint2, ProtonDistributionModel::ChargeDirected);

			double bb_enhance_factor(max(1.0, sqrt(bb_charge_full[i+1] / bb_charges_median)));
			//cerr << "BB_FACTOR=" << bb_enhance_factor << endl;
				hmm_.setInitialTransitionProbability("BB"+pos_name, /*bb_init[i]*/bb_sum /* + sc_charge_full[i]*/ * bb_enhance_factor /* + 1.0/(-log10(bb_charge_full[i+1]))*/);
				hmm_.setInitialTransitionProbability(aa1+aa2+"_bxyz"+pos_name, /*bb_init[i]*/ bb_sum /* + sc_charge_full[i]*/ * bb_enhance_factor /* + 1.0/(-log10(bb_charge_full[i+1]))*/);

				hmm_.setInitialTransitionProbability(aa1+aa2+"_axyz"+pos_name, bb_sum * bb_enhance_factor);
				//cerr << "CHARGES=" << aa1+aa2+"bxyz"+pos_name << " " << bb_sum + 1.0/(-log10(bb_charge_full[i+1])) << " " << bb_sum << " " << 1.0/(-log10(bb_charge_full[i+1])) << endl;

			//cerr << "ChargeStats: BB=" << BB_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;
			//charge_sum += BB_charges[i]; 

			prefixes.push_back(prefix);
			suffixes.push_back(suffix);
			
			double b_cr_int1(0), b_cr_int2(0), y_cr_int1(0), y_cr_int2(0);
			if ((aa1 == "D" || aa1 == "E" || pos_name == "k-1" || pos_name == "k-2") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("CR"+pos_name, /*cr_init[i]*/1 - bb_sum);
				charge_sum += cr_init[i];

				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,	Residue::BIon, b_cr_int1, y_cr_int1, b_cr_int2, y_cr_int2, ProtonDistributionModel::ChargeRemote);

				//cerr << "ChargeStats: CR=" << CR_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;
			}

			if ((aa1 == "K" || aa1 == "H" || aa1 == "R"))
			{
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_sc_int1, y_sc_int1, b_sc_int2, y_sc_int2, ProtonDistributionModel::SideChain);
				hmm_.setInitialTransitionProbability("SC"+pos_name, /*sc_init[i]*/sc_charge_full[i]);

				//cerr << "ChargeStats: SC=" << SC_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;	
				//charge_sum += SC_charges[i];
			}

			//cerr << "suffix_pos=" << suffix_pos << ", prefix_pos=" << i << ", peptide.size()=" << peptide.size() << endl;
			double y_sum1 = ints_1.ints[YIon][suffix_pos - 1] + ints_1.ints[YIon_H2O][suffix_pos - 1] + ints_1.ints[YIon_NH3][suffix_pos - 1];
			double b_sum1 = ints_1.ints[BIon][i] + ints_1.ints[BIon_H2O][i] + ints_1.ints[BIon_NH3][i] + ints_1.ints[AIon][i];
		
			double y_sum2(0), b_sum2(0);
			if (charge > 1)
			{
				y_sum2 = ints_2.ints[YIon][suffix_pos - 1] + ints_2.ints[YIon_H2O][suffix_pos - 1] + ints_2.ints[YIon_NH3][suffix_pos - 1];
				b_sum2 = ints_2.ints[BIon][i] + ints_2.ints[BIon_H2O][i] + ints_2.ints[BIon_NH3][i] + ints_2.ints[AIon][i];
			}

			//cerr << prefix << " - " << suffix << " " << b_sum1 << " " << y_sum1 << " " << i + 1 << " " << suffix_pos << " " << ints_1.ints[BIon][i] << endl;
			
			hmm_.setTrainingEmissionProbability(b_name1, b_sum1);
			hmm_.setTrainingEmissionProbability(a_name1, ints_1.ints[AIon][i]);
			hmm_.setTrainingEmissionProbability(y_name1, y_sum1);

			if (charge > 1)
			{
      	hmm_.setTrainingEmissionProbability(b_name2, b_sum2);
      	hmm_.setTrainingEmissionProbability(y_name2, y_sum2);
			}

			hmm_.setTrainingEmissionProbability("AA"+pos_name, b_sum1 + b_sum2 + y_sum1 + y_sum2);
			hmm_.setTrainingEmissionProbability("end"+pos_name, 0.5/(double(peptide.size() - 1)));

			//cerr << aa1 << aa2 << ": " << b_name1 << "=" << b_ints1[i] << ", " << y_name1 << "=" << y_ints1[y_ints1.size()-i-1] << 
			//		", y-H2O+ = " << y_H2O_ints1[y_H2O_ints1.size() - i - 1] << endl;
			//cerr << aa1 << aa2 << ": " << b_name2 << "=" << b_ints2[i] << ", " << y_name2 << "=" << y_ints2[y_ints2.size()-i-1] << 
			//		", y-H2O++ = " << y_H2O_ints2[y_H2O_ints2.size() - i - 1] << endl;
			//cerr << aa1 << aa2 << ": " << pos_name << " " << BB_charges[i] << " " << SC_charges[i] << endl;
	
			// now enable the states
				hmm_.enableTransition("BB"+pos_name, "AA"+pos_name);
				hmm_.enableTransition("BB"+pos_name, "end"+pos_name);

				hmm_.enableTransition(aa1+aa2+"_bxyz"+pos_name, "bxyz"+pos_name);
				hmm_.enableTransition(aa1+aa2+"_bxyz"+pos_name, "end"+pos_name);

				hmm_.setTransitionProbability("bxyz"+pos_name, b_name1, bint1);
				hmm_.setTransitionProbability("bxyz"+pos_name, y_name1, yint1);



				hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "axyz"+pos_name);
				hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "end"+pos_name);
				hmm_.setTransitionProbability("axyz"+pos_name, a_name1, bint1);
				hmm_.setTransitionProbability("axyz"+pos_name, y_name1, yint1);
				
				
			if (charge > 1 && bb_sum > 0)
			{
      	hmm_.setTransitionProbability("bxyz"+pos_name, b_name2, bint2);
      	hmm_.setTransitionProbability("bxyz"+pos_name, y_name2, yint2);

				hmm_.setTransitionProbability("axyz"+pos_name, b_name2, bint2);
				hmm_.setTransitionProbability("axyz"+pos_name, y_name2, yint2);
			}

			if (aa1 == "D" && is_charge_remote)
			{
				hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);

				hmm_.setTrainingEmissionProbability("A"+pos_name, b_sum1 * b_cr_int1 + b_sum2 * b_cr_int2 + y_sum1 * y_cr_int1 + y_sum2 * y_cr_int2);

				hmm_.setTransitionProbability("D"+pos_name, b_name1, b_cr_int1);
				hmm_.setTransitionProbability("D"+pos_name, y_name1, y_cr_int1);
				hmm_.setTransitionProbability("D"+pos_name, b_name2, b_cr_int2);
				hmm_.setTransitionProbability("D"+pos_name, y_name2, y_cr_int2);
				hmm_.setInitialTransitionProbability(aa2+"_D"+pos_name, /*cr_init[i]*/ 1 - bb_sum);
				hmm_.enableTransition(aa2+"_D"+pos_name, "D"+pos_name);
				hmm_.enableTransition(aa2+"_D"+pos_name, "end"+pos_name);
			}

			if (pos_name == "k-1" && is_charge_remote && aa1 != "D" && aa1 != "E" && (aa2 != "R" || peptide.has("K") || peptide.has("H")))
			{
				hmm_.enableTransition("CRk-1", "Ak-1");
				hmm_.enableTransition("CRk-1", "endk-1");

				hmm_.setTrainingEmissionProbability("Ak-1", b_sum1 + b_sum2 + y_sum1 + y_sum2);

				hmm_.setTransitionProbability("bk-1", b_name1, b_cr_int1);
				hmm_.setTransitionProbability("bk-1", y_name1, y_cr_int1);
				//hmm_.setTransitionProbability("bk-1", b_name2, b_cr_int2);
				//hmm_.setTransitionProbability("bk-1", y_name2, y_cr_int2);
				
				hmm_.setInitialTransitionProbability(aa1+"_bk-1", 1 - bb_sum);
				
				hmm_.enableTransition(aa1+"_bk-1", "bk-1");
				hmm_.enableTransition(aa1+"_bk-1", "endk-1");

				//cerr << aa1 << " " << b_sum1 + b_sum2 + y_sum1 + y_sum2 << endl;

				//cerr << aa1 << hmm_.getTransitionProbability(aa1 + "bk-1", "bk-1") << endl;
			}

			if (pos_name == "k-2" && is_charge_remote && aa1 != "D" && aa1 != "E" && (aa2 != "R" || peptide.has("K") || peptide.has("H")))
			{
				hmm_.enableTransition("CRk-2", "Ak-2");
				hmm_.enableTransition("CRk-2", "endk-2");
				
				hmm_.setTrainingEmissionProbability("Ak-2", b_sum1 + b_sum2 + y_sum1 + y_sum2);

				hmm_.setTransitionProbability("bk-2", b_name1, b_cr_int1);
				hmm_.setTransitionProbability("bk-2", y_name1, y_cr_int1);

				hmm_.setInitialTransitionProbability(aa1+"_bk-2", 1 - bb_sum);
				
				hmm_.enableTransition(aa1+"_bk-2", "bk-2");
        hmm_.enableTransition(aa1+"_bk-2", "endk-2");
			}
			
			
			if (aa1 == "E" && is_charge_remote)
      { 
        hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);
				hmm_.setTrainingEmissionProbability("A"+pos_name, b_sum1 * b_cr_int1 + b_sum2 * b_cr_int2 + y_sum1 * y_cr_int1 + y_sum2 * y_cr_int2);
        hmm_.setTransitionProbability("E"+pos_name, b_name1, b_cr_int1);
        hmm_.setTransitionProbability("E"+pos_name, y_name1, y_cr_int1);
				hmm_.setTransitionProbability("E"+pos_name, b_name2, b_cr_int2);
				hmm_.setTransitionProbability("E"+pos_name, y_name2, y_cr_int2);
				hmm_.setInitialTransitionProbability(aa2+"_E"+pos_name, /*cr_init[i]*/ 1 - bb_sum);
        hmm_.enableTransition(aa2+"_E"+pos_name, "E"+pos_name);
				hmm_.enableTransition(aa2+"_E"+pos_name, "end"+pos_name);
      }

			if (aa1 == "K" /*&& is_charge_remote*/)
			{
				hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
				hmm_.setTrainingEmissionProbability("ASC"+pos_name, b_sum1 * b_sc_int1 + b_sum2 * b_sc_int2 + y_sum1 * y_sc_int1 + y_sum2 * y_sc_int2);
				hmm_.setTransitionProbability("K"+pos_name, b_name1, b_sc_int1);
				hmm_.setTransitionProbability("K"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("K"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("K"+pos_name, y_name2, y_sc_int2);
				hmm_.setInitialTransitionProbability(aa2+"_K"+pos_name, /*sc_init[i]*/sc_charge_full[i]);
        hmm_.enableTransition(aa2+"_K"+pos_name, "K"+pos_name);
				hmm_.enableTransition(aa2+"_K"+pos_name, "end"+pos_name);
			}

			if (aa1 == "H" /*&& is_charge_remote*/)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
				hmm_.setTrainingEmissionProbability("ASC"+pos_name, b_sum1 * b_sc_int1 + b_sum2 * b_sc_int2 + y_sum1 * y_sc_int1 + y_sum2 * y_sc_int2);
        hmm_.setTransitionProbability("H"+pos_name, b_name1, b_sc_int1);
        hmm_.setTransitionProbability("H"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("H"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("H"+pos_name, y_name2, y_sc_int2);
				hmm_.setInitialTransitionProbability(aa2+"_H"+pos_name, /*sc_init[i]*/sc_charge_full[i]);
        hmm_.enableTransition(aa2+"_H"+pos_name, "H"+pos_name);
				hmm_.enableTransition(aa2+"_H"+pos_name, "end"+pos_name);
      }

      if (aa1 == "R" /*&& is_charge_remote*/)
      { 
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
				hmm_.setTrainingEmissionProbability("ASC"+pos_name, b_sum1 * b_sc_int1 + b_sum2 * b_sc_int2 + y_sum1 * y_sc_int1 + y_sum2 * y_sc_int2);
        hmm_.setTransitionProbability("R"+pos_name, b_name1, b_sc_int1);
        hmm_.setTransitionProbability("R"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("R"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("R"+pos_name, y_name2, y_sc_int2);
				hmm_.setInitialTransitionProbability(aa2+"_R"+pos_name, /*sc_init[i]*/sc_charge_full[i]);
        hmm_.enableTransition(aa2+"_R"+pos_name, "R"+pos_name);
				hmm_.enableTransition(aa2+"_R"+pos_name, "end"+pos_name);
      } 

			/*
			// losses
			double y_loss_sum_1 = ints_1.ints[YIon][suffix_pos - 1] + ints_1.ints[YIon_H2O][suffix_pos - 1] + ints_1.ints[YIon_NH3][suffix_pos - 1];
			if (y_loss_sum_1 > 0.0)
			{
				Map<String, double> y_intensities;
				y_intensities[LOSS_TYPE_H2O] = ints_1.ints[YIon_H2O][suffix_pos - 1];
				y_intensities[LOSS_TYPE_NH3] = ints_1.ints[YIon_NH3][suffix_pos - 1];
				double ion_intensity = ints_1.ints[YIon][suffix_pos - 1];
				trainNeutralLossesFromIon_(y_sum1, y_intensities, YIon, ion_intensity, suffix);
				//cerr << "YLOSS: " << suffix << " " << y_sum1 << " " << y_intensities[LOSS_TYPE_H2O] << " " << y_intensities[LOSS_TYPE_NH3] << " " << ion_intensity << endl;

				//cerr << "YLOSS: " << prefix << " " << suffix << " " << 1 << " " << y_sum1 << " y-H2O " << y_intensities[LOSS_TYPE_H2O] << " y-NH3 " << y_intensities[LOSS_TYPE_NH3] << " y " << ion_intensity << endl;
			}
		
			double y_loss_sum_2(0);
			if (charge > 1)
			{
				y_loss_sum_2 = ints_2.ints[YIon][suffix_pos - 1] + ints_2.ints[YIon_H2O][suffix_pos - 1] + ints_2.ints[YIon_NH3][suffix_pos - 1];
			}
			if (charge > 1 && y_loss_sum_2 > 0.0)
			{
				Map<String, double> y_intensities;
				y_intensities[LOSS_TYPE_H2O] = ints_2.ints[YIon_H2O][suffix_pos - 1];
				y_intensities[LOSS_TYPE_NH3] = ints_2.ints[YIon_NH3][suffix_pos - 1];
				double ion_intensity = ints_2.ints[YIon][suffix_pos - 1];

				trainNeutralLossesFromIon_(y_sum2, y_intensities, YIon, ion_intensity, suffix);

				//cerr << "YLOSS: " << prefix << " " << suffix << " " << 2 << " " << y_sum2 << " y-H2O " << y_intensities[LOSS_TYPE_H2O] << " y-NH3 " << y_intensities[LOSS_TYPE_NH3] << " y " << ion_intensity << endl;
			}
			

			double b_loss_sum_1 = ints_1.ints[BIon][i] + ints_1.ints[BIon_H2O][i] + ints_1.ints[BIon_NH3][i] + ints_1.ints[AIon][i];
			if (b_loss_sum_1 > 0.0)
			{
      	Map<String, double> b_intensities;
      	b_intensities[LOSS_TYPE_H2O] = ints_1.ints[BIon_H2O][i];
      	b_intensities[LOSS_TYPE_NH3] = ints_1.ints[BIon_NH3][i];
				b_intensities[LOSS_TYPE_CO] = ints_1.ints[AIon][i];
      	double ion_intensity = ints_1.ints[BIon][i];
				if (prefix.size() == 2)
				{
					trainNeutralLossesFromIon_(b_sum1, b_intensities, B2Ion, ion_intensity, prefix);
				}
				else
				{
					trainNeutralLossesFromIon_(b_sum1, b_intensities, BIon, ion_intensity, prefix);	
				}

				//cerr << "BLOSS: " << prefix << " " << suffix << " " << 1 << " " << b_sum1 << " b-H2O " << b_intensities[LOSS_TYPE_H2O] << " b-NH3 " << b_intensities[LOSS_TYPE_NH3] << " b-CO " << b_intensities[LOSS_TYPE_CO] << " b " << ion_intensity << endl;
			}

			double b_loss_sum_2(0);
			if (charge > 1)
			{
				b_loss_sum_2 = ints_2.ints[BIon][i] + ints_2.ints[BIon_H2O][i] + ints_2.ints[BIon_NH3][i] + ints_2.ints[AIon][i];
			}
			if (charge > 1 && b_loss_sum_2 > 0.0)
			{
				Map<String, double> b_intensities;
				b_intensities[LOSS_TYPE_H2O] = ints_2.ints[BIon_H2O][i];
				b_intensities[LOSS_TYPE_NH3] = ints_2.ints[BIon_NH3][i];
				b_intensities[LOSS_TYPE_CO] = ints_2.ints[AIon][i];
				double ion_intensity = ints_2.ints[BIon][i];
				//cerr << "BLOSS: " << prefix << " " << suffix << " " << 2 << " " << b_sum2 << " b-H2O " << b_intensities[LOSS_TYPE_H2O] << " b-NH3 " << b_intensities[LOSS_TYPE_NH3] << " b-CO " << b_intensities[LOSS_TYPE_CO] << " b " << ion_intensity << endl;
				if (prefix.size() == 2)
				{
					trainNeutralLossesFromIon_(b_sum2, b_intensities, B2Ion, ion_intensity, prefix);
				}
				else
				{
					trainNeutralLossesFromIon_(b_sum2, b_intensities, BIon, ion_intensity, prefix);
				}
			}
*/
		}

		//cerr << "PrecursorStats: " << pre_int1 << " " << pre_int2 << " " << pre_int_H2O_1 << " " << pre_int_H2O_2 << " " << pre_int_NH3_1 << " " << pre_int_NH3_2 << endl;

		// precursor handling
		double pre_sum(0);
		
		for (Map<String, double>::ConstIterator it = pre_ints.begin(); it != pre_ints.end(); ++it)
		{
			pre_sum += it->second;
		}
		
		/*
		if (peptide[0].getOneLetterCode() == "Q" && pre_ints.pre_H2O > 0.05)
		{
			trainPrecursorIons_(1, pre_ints, peptide);
		}*/
		if (is_charge_remote && charge < 3 && pre_sum > 0.05)
		{
			//cerr << "BB_SUM: " << bb_sum << " " << pre_sum << " " << peptide << endl;
			trainPrecursorIons_(max(0.01, 1.0 - bb_sum), pre_ints, peptide);
		}

		// now train the model with the data set
		hmm_.train();

		//stringstream ss;
		//ss << peptide;
		//hmm_.writetoYGFFile(String("stats/model_graph_train_"+ss.str()+"_"+String(charge)+".graphml").c_str());

		hmm_.disableTransitions();

		return;
	}

	void PILISModel::getPrecursorIntensitiesFromSpectrum_(const RichPeakSpectrum& train_spec, Map<String, double>& peak_ints, double peptide_weight, UInt charge)
	{
#ifdef PRECURSOR_DEBUG
		cerr << "PILISModel::getPrecursorIntensitiesFromSpectrum_(#peaks=" << train_spec.size() << ", weight=" << peptide_weight << ", charge=" << charge << ")" << endl;
#endif
		vector<String> precursor_losses;
		precursor_losses.push_back(EmpiricalFormula("H2O").getString());
		precursor_losses.push_back(EmpiricalFormula("NH3").getString());
		//precursor_losses.push_back(EmpiricalFormula("H2S").getString());
		//precursor_losses.push_back(EmpiricalFormula("SCH3").getString());
		
		//double pre_error = (double)param_.getValue("precursor_mass_tolerance");
		double pre_error = (double)param_.getValue("peak_mass_tolerance");

  	for (RichPeakSpectrum::ConstIterator it = train_spec.begin(); it != train_spec.end(); ++it)
    {
			// intact peptide
			if (fabs(it->getMZ() - peptide_weight / (double)charge) < pre_error)
			{
				if (peak_ints.has(""))
				{
					peak_ints[""] += it->getIntensity();
				}
				else
				{
					peak_ints[""] = it->getIntensity();
				}
			}

			// single losses
			for (vector<String>::const_iterator lit1 = precursor_losses.begin(); lit1 != precursor_losses.end(); ++lit1)
			{		
				if (fabs(it->getMZ() - (peptide_weight - EmpiricalFormula(*lit1).getMonoWeight()) / (double)charge) < pre_error)
				{
					if (peak_ints.has(*lit1))
					{
						peak_ints[*lit1] += it->getIntensity();
					}
					else
					{
						peak_ints[*lit1] = it->getIntensity();
					}
				}

				for (vector<String>::const_iterator lit2 = lit1; lit2 != precursor_losses.end(); ++lit2)
				{
					String name;
					if (*lit1 < *lit2)
					{
						name = *lit1 + "-" + *lit2;
					}
					else
					{
						name = *lit2 + "-" + *lit1;
					}

					if (fabs(it->getMZ() - (peptide_weight - EmpiricalFormula(*lit1).getMonoWeight() -  EmpiricalFormula(*lit2).getMonoWeight()) / (double)charge) < pre_error)
					{
						if (peak_ints.has(name))
						{
							peak_ints[name] += it->getIntensity();
						}
						else
						{
							peak_ints[name] = it->getIntensity();
						}
					}
				}
			}
			
		}
		return;
	}
	
	double PILISModel::getIntensitiesFromSpectrum_(const RichPeakSpectrum& train_spec, IonPeaks_& ion_ints, const AASequence& peptide, UInt z)
	{
		double sum(0);
    RichPeakSpectrum y_theo_spec, y_H2O_theo_spec, y_NH3_theo_spec;
		TheoreticalSpectrumGenerator tsg;
    tsg_.addPeaks(y_theo_spec, peptide, Residue::YIon, z);
    
		RichPeak1D p;
		p.setIntensity(1);
		static const double h2o_weight = EmpiricalFormula("H2O").getMonoWeight();
		static const double nh3_weight = EmpiricalFormula("NH3").getMonoWeight();
		static const double co_weight = EmpiricalFormula("CO").getMonoWeight();
    for (RichPeakSpectrum::ConstIterator it = y_theo_spec.begin(); it != y_theo_spec.end(); ++it)
    {
      p.setPosition(it->getPosition() - h2o_weight / (double)z);
      y_H2O_theo_spec.push_back(p);

      p.setPosition(it->getPosition() - nh3_weight / (double)z);
      y_NH3_theo_spec.push_back(p);
    }

		vector<double> y_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, y_theo_spec, y_ints);
		ion_ints.ints[YIon] = y_ints;

		vector<double> y_H2O_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, y_H2O_theo_spec, y_H2O_ints);
		ion_ints.ints[YIon_H2O] = y_H2O_ints;

		vector<double> y_NH3_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, y_NH3_theo_spec, y_NH3_ints);
		ion_ints.ints[YIon_NH3] = y_NH3_ints;

    RichPeakSpectrum a_theo_spec;
    tsg.addPeaks(a_theo_spec, peptide, Residue::AIon, z);
			
		vector<double> a_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, a_theo_spec, a_ints);
		ion_ints.ints[AIon] = a_ints;

    RichPeakSpectrum b_theo_spec, b_H2O_theo_spec, b_NH3_theo_spec, b_CO_theo_spec;
    tsg.addPeaks(b_theo_spec, peptide, Residue::BIon, z);
		
    for (RichPeakSpectrum::ConstIterator it = b_theo_spec.begin(); it != b_theo_spec.end(); ++it)
    {
      p.setPosition(it->getPosition() - h2o_weight / (double)z);
      b_H2O_theo_spec.push_back(p);
      p.setPosition(it->getPosition() - nh3_weight / (double)z);
      b_NH3_theo_spec.push_back(p);
			p.setPosition(it->getPosition() - co_weight / (double)z);
			b_CO_theo_spec.push_back(p);
    }

		vector<double> b_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, b_theo_spec, b_ints);
		ion_ints.ints[BIon] = b_ints;

		vector<double> b_H2O_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, b_H2O_theo_spec, b_H2O_ints);
		ion_ints.ints[BIon_H2O] = b_H2O_ints;

		vector<double> b_NH3_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, b_NH3_theo_spec, b_NH3_ints);
		ion_ints.ints[BIon_NH3] = b_NH3_ints;

		vector<double> b_CO_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, b_CO_theo_spec, b_CO_ints);
		ion_ints.ints[AIon] = b_CO_ints;
		
		return sum;
	}

	double PILISModel::getIntensitiesFromComparison_(const RichPeakSpectrum& train_spec, const RichPeakSpectrum& theo_spec, vector<double>& intensities)
	{
		double sum(0);
		vector<pair<UInt, UInt> > peak_map;
		spectra_aligner_.getSpectrumAlignment(peak_map, train_spec, theo_spec);

    for (vector<pair<UInt, UInt> >::const_iterator it = peak_map.begin(); it != peak_map.end(); ++it)
    {
			intensities[it->second] = train_spec[it->first].getIntensity();
			sum += train_spec[it->first].getIntensity();
    }
		return sum;
	}
	
	void PILISModel::trainPrecursorIons_(double initial_probability, const Map<String, double>& ints, const AASequence& peptide)
	{
#ifdef PRECURSOR_DEBUG
		cerr << "PILISModel::trainPrecursorIons_(" << initial_probability << ", " << ints.size() << ", " << peptide << ")" << endl;
#endif
		// clean up
		hmm_precursor_.clearInitialTransitionProbabilities();
		hmm_precursor_.clearTrainingEmissionProbabilities();

		// set start transition prob
		hmm_precursor_.setInitialTransitionProbability("start", initial_probability);
	
		// set emission probabilities from the precursor ions present in the spectrum
		for (Map<String, double>::ConstIterator it = ints.begin(); it != ints.end(); ++it)
		{
			if (it->first != "")
			{
				hmm_precursor_.setTrainingEmissionProbability("p-" + it->first, it->second);
			}
			else
			{
				hmm_precursor_.setTrainingEmissionProbability("p", it->second);
			}
		}
	
		// build the final model
		enablePrecursorIonStates_(peptide);
	
#ifdef PRECURSOR_DEBUG
		stringstream ss;
		ss << peptide;
		hmm_precursor_.writeGraphMLFile("graphs/model_graph_train_precursor_" + peptide.toUnmodifiedString() + ".graphml");
#endif
		
		// train
		hmm_precursor_.train();

		// reset the model
		hmm_precursor_.disableTransitions();

		return;
	}

	void PILISModel::trainNeutralLossesFromIon_(double /*initial_probability*/, 
																								const Map<String, double>& /*ints*/, 
																								IonType_ /*ion_type*/, 
																								double /*ion_intensity*/, 
																								const AASequence& /*ion*/)
	{
		//Map<String, double> intensities = ints;
		//cerr << "LOSS: " << initial_probability << " " << intensities[LOSS_TYPE_H2O] << " " << intensities[LOSS_TYPE_NH3] << " " << ion_intensity << " " << ion << endl;
		// b-ions loss
		/*
		HiddenMarkovModelLight* hmm = &hmms_losses_[ion_type];
		hmm->clearInitialTransitionProbabilities();
		hmm->clearTrainingEmissionProbabilities();
		if (ion_type == BIon || ion_type == B2Ion)
		{
			// normalize to sum = 1
			double sum = ion_intensity + intensities[LOSS_TYPE_H2O] + intensities[LOSS_TYPE_NH3] + intensities[LOSS_TYPE_CO];
			intensities[LOSS_TYPE_H2O] /= sum;
			intensities[LOSS_TYPE_NH3] /= sum;
			intensities[LOSS_TYPE_CO] /= sum;
			ion_intensity /= sum;	
			
			hmm->setInitialTransitionProbability(B_ION, initial_probability);

			// H2O loss
			if (intensities.has(LOSS_TYPE_H2O))
			{
				hmm->setTrainingEmissionProbability(B_H2O, intensities[LOSS_TYPE_H2O]);
			}
			
			// NH3 loss
			if (intensities.has(LOSS_TYPE_NH3))
			{
				hmm->setTrainingEmissionProbability(B_NH3, intensities[LOSS_TYPE_NH3]);
			}

			// a-ions (CO loss)
			if (intensities.has(LOSS_TYPE_CO))
			{
				hmm->setTrainingEmissionProbability(A_ION, intensities[LOSS_TYPE_CO]);
			}
			// end state
			hmm->setTrainingEmissionProbability(B_LOSS_END, ion_intensity);
			enableNeutralLossStates_(ion_type, ion);

			hmm->train();
			hmm->disableTransitions();

			return;
		}

		if (ion_type == YIon)
		{
			double sum = ion_intensity + intensities[LOSS_TYPE_H2O] + intensities[LOSS_TYPE_NH3];
      intensities[LOSS_TYPE_H2O] /= sum;
      intensities[LOSS_TYPE_NH3] /= sum;
      ion_intensity /= sum;

			hmm->setInitialTransitionProbability(Y_ION, initial_probability);
      if (intensities.has(LOSS_TYPE_H2O))
      {
        hmm->setTrainingEmissionProbability(Y_H2O, intensities[LOSS_TYPE_H2O]);
      }

      if (intensities.has(LOSS_TYPE_NH3))
      {
        hmm->setTrainingEmissionProbability(Y_NH3, intensities[LOSS_TYPE_NH3]);
      }
      hmm->setTrainingEmissionProbability(Y_LOSS_END, ion_intensity);

			enableNeutralLossStates_(ion_type, ion);

			hmm->train();
			hmm->disableTransitions();
			
      return;
		}
*/
		return;
	}

  void PILISModel::enablePrecursorIonStates_(const AASequence& peptide)
	{
#ifdef PRECURSOR_DEBUG
		cerr << "void PILISModel::enablePrecursorIonStates_(" << peptide << ")" << endl;
#endif
		UInt max_explicit(4);
		Map<String, UInt> double_losses;
		vector<String> nexts;

		for (AASequence::ConstIterator it1 = peptide.begin(); it1 != peptide.end(); ++it1)
		{
			AASequence aa1;
			aa1 += &*it1;

			String loss1 = it1->getLossFormula().getString();
			if (loss1 == "")
			{
				continue;
			}
		
			for (AASequence::ConstIterator it2 = it1 + 1; it2 != peptide.end(); ++it2)
			{
				AASequence aa2;
				aa2 += &*it2;

				String loss2 = it2->getLossFormula().getString();
				if (loss2 == "")
				{
					continue;
				}

				String name;
				if (aa1 < aa2)
				{
					name = aa1.toString() + aa2.toString();
				}
				else
				{
					name = aa2.toString() + aa1.toString();
				}

				String losses;
				if (loss1 < loss2)
				{
					losses = "-" + loss1 + "-" + loss2;
				}
				else
				{
					losses = "-" + loss2 + "-" + loss1;
				}

				String num;
				if (!double_losses.has(name))
				{
					double_losses[name] = 1;
					num = "1";
				}
				else
				{
					++double_losses[name];
					if (double_losses[name] > max_explicit)
					{
						break;
					}
					num = String(double_losses[name]);
				}

				name += losses + "_" + num;
			
				// enable transitions to emit state, and single emission states
				hmm_precursor_.enableTransition(name, "p" + losses);
				hmm_precursor_.enableTransition(name, aa1.toString() + "-" + loss1 + "_1");
				vector<String> split;
				name.split('_', split);
				hmm_precursor_.setTransitionProbability(name, aa1.toString() + "-" + loss1 + "_1", hmm_precursor_.getTransitionProbability(split[0], aa1.toString() + "-" + loss1));
				hmm_precursor_.enableTransition(name, aa2.toString() + "-" + loss2 + "_1");
				hmm_precursor_.setTransitionProbability(name, aa2.toString() + "-" + loss2 + "_1", hmm_precursor_.getTransitionProbability(split[0], aa2.toString() + "-" + loss2));
				
				// store transition to connect the next lines
				nexts.push_back(name);
			}
		}


		if (nexts.size() > 0)
		{
			hmm_precursor_.setTransitionProbability("start", *nexts.begin(), 1.0);
		}
		
		for (vector<String>::const_iterator it = nexts.begin(); it != nexts.end(); ++it)
		{
			if (it != (nexts.end() - 1))
			{
				vector<String>::const_iterator next_it = it + 1;
				hmm_precursor_.enableTransition(*it, *next_it);
				vector<String> split;
        it->split('_', split);
        hmm_precursor_.setTransitionProbability(*it, *next_it, hmm_precursor_.getTransitionProbability(split[0], split[0] + "-next"));

			}
		}
		
		Map<String, UInt> single_losses;
		vector<String> single_nexts;
		for (AASequence::ConstIterator it1 = peptide.begin(); it1 != peptide.end(); ++it1)
    {
      AASequence aa1;
      aa1 += &*it1;
			String name = aa1.toString();
			String loss = it1->getLossFormula().getString();
			String num;
			if (loss != "")
			{
				if (!single_losses.has(name))
				{
					single_losses[name] = 1;
					num = "1";
				}
				else
				{
					++single_losses[name];
					if (single_losses[name] > max_explicit)
					{
						break;
					}
					num = String(single_losses[name]);
				}

				name += "-" + loss + "_" + num;
				hmm_precursor_.enableTransition(name, "p-" + loss);
				single_nexts.push_back(name);
			}
		}

		for (vector<String>::const_iterator it = single_nexts.begin(); it != single_nexts.end(); ++it)
		{
			if (it != (single_nexts.end() - 1))
			{
				vector<String>::const_iterator next_it = it + 1;
				hmm_precursor_.enableTransition(*it, *next_it);
				vector<String> split;
				it->split('_', split);
				hmm_precursor_.setTransitionProbability(*it, *next_it, hmm_precursor_.getTransitionProbability(split[0], split[0] + "-next"));
			}
		}

		// last transition, no reaction took place
		if (single_nexts.size() > 0)
		{
			hmm_precursor_.enableTransition(single_nexts.back(), "p");
		}


		return;
	}

  void PILISModel::enableNeutralLossStates_(IonType_ /*ion_type*/, const AASequence& /*ion*/)
	{
		return;
	}
	

	void PILISModel::evaluate()
	{
		hmm_.evaluate();
		hmm_precursor_.evaluate();

		
		//hmm_.estimateUntrainedTransitions(); // TODO
	}

	void PILISModel::getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, UInt charge)
	{
		//cerr << "==============================================================================" << endl;
		//cerr << peptide << " " << charge << endl;
		if (!valid_)
		{
			cerr << "PILISModel: cannot simulate, initialize model from file first, e.g. data/PILIS/PILIS_model_default.dat" << endl;
			return;
		}

		if (peptide.size() > 39)
		{
			cerr << "PILISModel: cannot generate spectra of peptide '" << peptide << "' of length " << peptide.size() << " (max is 39)" << endl;
			return;
		}
		
		// calc proton distribution
		Map<UInt, double> sc_charge_full, bb_charge_full;
		prot_dist_.getProtonDistribution(bb_charge_full, sc_charge_full, peptide, charge, Residue::YIon);
		prot_dist_.setPeptideProtonDistribution(bb_charge_full, sc_charge_full);

		hmm_.clearInitialTransitionProbabilities();
		hmm_.clearTrainingEmissionProbabilities();
	
    // set charges
    vector<double> bb_init, sc_init, cr_init;
		bool is_charge_remote = getInitialTransitionProbabilities_(bb_init, cr_init, sc_init, bb_charge_full, sc_charge_full, peptide);
	
		double bb_sum(0), bb_sum_orig(0);
		for (Map<UInt, double>::ConstIterator iit = bb_charge_full.begin(); iit != bb_charge_full.end(); ++iit)
		{
			bb_sum += iit->second;
		}
	
    // activation of protons sitting at lysine and histidine side chains
    for (UInt i = 0; i != peptide.size(); ++i)
    {
      if (peptide[i].getOneLetterCode() == "H" || peptide[i].getOneLetterCode() == "K")
      {
        bb_sum += 0.5 * sc_charge_full[i]; // TODO
      }
    }
		
		if (bb_sum > 1)
		{
			bb_sum = 1;
		}
		bb_sum_orig = bb_sum;
		//bb_sum /= charge;
		if (bb_sum < (double)param_.getValue("charge_directed_threshold") * charge)
		{
			bb_sum = (double)param_.getValue("charge_directed_threshold") * charge;
		}

		double charge_sum(0);

		vector<AASequence> suffixes, prefixes;
		vector<String> suffix_names1, suffix_names2, prefix_names1, prefix_names2, a_names1, a_names2;
	
		vector<double> bb_charges_all;
		for (UInt i = 0; i != bb_charge_full.size(); ++i)
		{
			bb_charges_all.push_back(bb_charge_full[i]);
		}
		sort(bb_charges_all.begin(), bb_charges_all.end());
		double bb_charges_median = bb_charges_all[UInt(bb_charges_all.size()/2.0)];
		
		// get the paths
		for (UInt i = 0; i != peptide.size() - 1; ++i)
		{

      String pos_name, y_name1, b_name1, a_name1, y_name, b_name;
      String y_name2, b_name2, a_name2;
			String i_plus1(i+1), pep_size_i(peptide.size() - 1 - i);


      if (i < floor(peptide.size()/2))
      {
        pos_name = i_plus1;
				y_name = "yk-"+i_plus1;
				b_name = "b"+i_plus1;
        y_name1 = "yk-"+i_plus1+"+";
        b_name1 = "b"+i_plus1+"+";
        a_name1 = "a"+i_plus1+"+";
        
        y_name2 = "yk-"+i_plus1+"++";
        b_name2 = "b"+i_plus1+"++";
        a_name2 = "a"+i_plus1+"++";
      }
      else
      { 
        pos_name = "k-"+pep_size_i;
				y_name = "y"+pep_size_i;
				b_name = "bk-"+pep_size_i;
        y_name1 = "y"+pep_size_i+"+";
        b_name1 = "bk-"+pep_size_i+"+";
        a_name1 = "ak-"+pep_size_i+"+";
        
        y_name2 = "y"+pep_size_i+"++";
        b_name2 = "bk-"+pep_size_i+"++";
        a_name2 = "ak-"+pep_size_i+"++";
      }

			AASequence prefix(peptide.getPrefix(i + 1)), suffix(peptide.getSuffix(peptide.size() - 1 - i));
			suffixes.push_back(suffix);
			prefixes.push_back(prefix);
			suffix_names1.push_back(y_name1);
			suffix_names2.push_back(y_name2);
			prefix_names1.push_back(b_name1);
			prefix_names2.push_back(b_name2);
			a_names1.push_back(a_name1);
			a_names2.push_back(a_name2);

			//String aa1(peptide[i].getOneLetterCode()), aa2(peptide[i+1].getOneLetterCode());
			AASequence aa1_seq, aa2_seq;
			aa1_seq += &peptide[i];
			aa2_seq += &peptide[i + 1];
			String aa1(aa1_seq.toString()), aa2(aa2_seq.toString());

			double bb_enhance_factor(max(1.0, sqrt(bb_charge_full[i+1] / bb_charges_median)));
			hmm_.setInitialTransitionProbability("BB"+pos_name, /*bb_init[i]*/ bb_sum /* + sc_charge_full[i]*/ * bb_enhance_factor /* + 10 * 1.0/(-log10(bb_charge_full[i+1]))*/);
			charge_sum += bb_init[i];

			double charge_loss_factor = (double)param_.getValue("charge_loss_factor");
			double b_cr_int1(0), b_cr_int2(0), y_cr_int1(0), y_cr_int2(0);
			double bint1(0), yint1(0), bint2(0), yint2(0), aint1(0), aint2(0), ayint1(0), ayint2(0), b_sc_int1(0), b_sc_int2(0), y_sc_int1(0), y_sc_int2(0);
			if ((aa1 == "D" || aa1 == "E" || pos_name == "k-1" || pos_name == "k-2") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("CR"+pos_name, /*cr_init[i]*/ 1 - bb_sum);
				charge_sum += cr_init[i];
				prot_dist_.getChargeStateIntensities(	peptide, prefix, suffix, charge,
																							Residue::BIon, b_cr_int1, y_cr_int1, b_cr_int2, y_cr_int2, ProtonDistributionModel::ChargeRemote);

				//cerr << "charge state intensities (CR): " << peptide << " (" << prefix << "-" << suffix << ") " << b_cr_int1 << " " << y_cr_int1 << " " << b_cr_int2 << " " << y_cr_int2 << " " << charge << endl;
				
				// TODO
				//if (charge < 3)
				//{
					y_cr_int1 += charge_loss_factor * y_cr_int2;
					y_cr_int2 *= 1 - charge_loss_factor;
					b_cr_int1 += charge_loss_factor * b_cr_int2;
					b_cr_int2 *= 1 - charge_loss_factor;
				//}
			}

			if ((aa1 == "K" || aa1 == "H" || aa1 == "R")/* && is_charge_remote*/)
			{
				hmm_.setInitialTransitionProbability("SC"+pos_name, /*sc_init[i]*/sc_charge_full[i]);
				charge_sum += sc_init[i];
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,	Residue::BIon, b_sc_int1, y_sc_int1, b_sc_int2, y_sc_int2, ProtonDistributionModel::SideChain);

				// TODO
				//if (charge < 3)
				//{
					y_sc_int1 += charge_loss_factor * y_sc_int2;
        	y_sc_int2 *= 1 - charge_loss_factor;
        	b_sc_int1 += charge_loss_factor * b_sc_int2;
        	b_sc_int2 *= 1 - charge_loss_factor;
				//}

				//cerr << "charge state intensities (SC): " << peptide << " (" << prefix << "-" << suffix << ") " << b_sc_int1 << " " << y_sc_int1 << " " << b_sc_int2 << " " << y_sc_int2 << " " << charge << endl;
			}
		
      prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, bint1, yint1, bint2, yint2, ProtonDistributionModel::ChargeDirected);

			//cerr << "charge state intensities (CD): " << peptide << " (" << prefix << "-" << suffix << ") " << bint1 << " " << yint1 << " " << bint2 << " " << yint2 << " " << charge << endl;

      prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::AIon, aint1, ayint1, aint2, ayint2, ProtonDistributionModel::ChargeDirected);	
	
			// TODO correction of the doubly charged ions to singly charged ones, to account for proton loss (and not to need it model explicitly)
			//if (charge < 3)
			//{
				yint1 += charge_loss_factor * yint2;
				yint2 *= 1 - charge_loss_factor;
				bint1 += charge_loss_factor * bint2;
				bint2 *= 1 - charge_loss_factor;
			//}

      // now enable the states
			hmm_.enableTransition("BB"+pos_name, "AA"+pos_name);
			hmm_.enableTransition("BB"+pos_name, "end"+pos_name);
			
      hmm_.enableTransition("AA"+pos_name, aa1+aa2+"_bxyz"+pos_name);
      hmm_.enableTransition(aa1+aa2+"_bxyz"+pos_name, "bxyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"_bxyz"+pos_name, "end"+pos_name);

      hmm_.setTransitionProbability("bxyz"+pos_name, b_name1, bint1);
      hmm_.setTransitionProbability("bxyz"+pos_name, y_name1, yint1);

			hmm_.setTransitionProbability("bxyz"+pos_name, b_name2, bint2);
			hmm_.setTransitionProbability("bxyz"+pos_name, y_name2, yint2);



			// axyz
			hmm_.enableTransition("AA"+pos_name, aa1+aa2+"_axyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "axyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "end"+pos_name);

			hmm_.setTransitionProbability("axyz"+pos_name, a_name1, aint1);
			hmm_.setTransitionProbability("axyz"+pos_name, y_name1, yint1);

			hmm_.setTransitionProbability("axyz"+pos_name, a_name2, aint2);
			hmm_.setTransitionProbability("axyz"+pos_name, y_name2, yint2);
			
			//cerr << "neutral loss charge sums: " << H2O_b_charge_sum << " " << H2O_y_charge_sum << " " << NH3_b_charge_sum << " " << NH3_y_charge_sum << " " << charge_dir_tmp << endl;

      if (aa1 == "D" && is_charge_remote)
      {
        hmm_.setTransitionProbability("D"+pos_name, b_name1, b_cr_int1);
        hmm_.setTransitionProbability("D"+pos_name, y_name1, y_cr_int1);
				hmm_.setTransitionProbability("D"+pos_name, b_name2, b_cr_int2);
				hmm_.setTransitionProbability("D"+pos_name, y_name2, y_cr_int2);
				hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);
				hmm_.enableTransition("A"+pos_name, aa2+"_D"+pos_name);
        hmm_.enableTransition(aa2+"_D"+pos_name, "D"+pos_name);
				hmm_.enableTransition(aa2+"_D"+pos_name, "end"+pos_name);
			}

      if (pos_name == "k-1" && is_charge_remote && aa1 != "D" && ((aa2 != "R" && aa1 != "E") || peptide.has("K") || peptide.has("H")))
      {
				if (charge == 1)
				{
					hmm_.setTransitionProbability("bk-1", b_name1, 1.0);
					hmm_.setTransitionProbability("bk-1", y_name1, 0.0);
				}
				else
				{
					hmm_.setTransitionProbability("bk-1", b_name1, b_cr_int1);
					hmm_.setTransitionProbability("bk-1", y_name1, y_cr_int1);
					hmm_.setTransitionProbability("bk-1", b_name2, b_cr_int2);
					hmm_.setTransitionProbability("bk-1", y_name2, y_cr_int2);
				}

				hmm_.enableTransition("CRk-1", "Ak-1");
        hmm_.enableTransition("CRk-1", "endk-1");
				hmm_.enableTransition("Ak-1", aa1+"_bk-1");

				hmm_.enableTransition(aa1 + "_bk-1", "bk-1");
        hmm_.enableTransition(aa1 + "_bk-1", "end"+pos_name);
				//cerr << aa1 << " " << aa2 << hmm_.getTransitionProbability(aa1 + aa2 + "bk-1", "bk-1");
      }

			if (pos_name == "k-2" && is_charge_remote && aa1 != "D" && ((aa2 != "R" && aa1 != "E") || peptide.has("K") || peptide.has("H")))
			{
				if (charge == 1)
				{
					hmm_.setTransitionProbability("bk-2", b_name1, 1.0);
					hmm_.setTransitionProbability("bk-2", y_name1, 0.0);
				}
				else
				{
					hmm_.setTransitionProbability("bk-2", b_name1, b_cr_int1);
					hmm_.setTransitionProbability("bk-2", y_name1, y_cr_int1);
					hmm_.setTransitionProbability("bk-2", b_name2, b_cr_int2);
					hmm_.setTransitionProbability("bk-2", y_name2, y_cr_int2);
				}
				hmm_.enableTransition("CRk-2", "Ak-2");
				hmm_.enableTransition("CRk-2", "endk-2");
				hmm_.enableTransition("Ak-2", aa1+"_bk-2");

				hmm_.enableTransition(aa1 + "_bk-2", "bk-2");
				hmm_.enableTransition(aa1 + "_bk-2", "end"+pos_name);
			}

			if (aa1 == "E" && is_charge_remote)
      {
        hmm_.setTransitionProbability("E"+pos_name, b_name1, b_cr_int1);
        hmm_.setTransitionProbability("E"+pos_name, y_name1, y_cr_int1);
				hmm_.setTransitionProbability("E"+pos_name, b_name2, b_cr_int2);
				hmm_.setTransitionProbability("E"+pos_name, y_name2, y_cr_int2);
        hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);
        hmm_.enableTransition("A"+pos_name, aa2+"_E"+pos_name);
        hmm_.enableTransition(aa2+"_E"+pos_name, "E"+pos_name);
				hmm_.enableTransition(aa2+"_E"+pos_name, "end"+pos_name);
      }

			if (aa1 == "K"/* && is_charge_remote*/)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
        hmm_.setTransitionProbability("K"+pos_name, b_name1, b_sc_int1);
        hmm_.setTransitionProbability("K"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("K"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("K"+pos_name, y_name2, y_sc_int2);
        hmm_.enableTransition("ASC"+pos_name, aa2+"_K"+pos_name);
        hmm_.enableTransition(aa2+"_K"+pos_name, "K"+pos_name);
				hmm_.enableTransition(aa2+"_K"+pos_name, "end"+pos_name);
			}
				
			if (aa1 == "H"/* && is_charge_remote*/)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
        hmm_.setTransitionProbability("H"+pos_name, b_name1, b_sc_int1);
        hmm_.setTransitionProbability("H"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("H"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("H"+pos_name, y_name2, y_sc_int2);
        hmm_.enableTransition("ASC"+pos_name, aa2+"_H"+pos_name);
        hmm_.enableTransition(aa2+"_H"+pos_name, "H"+pos_name);
				hmm_.enableTransition(aa2+"_H"+pos_name, "end"+pos_name);
      }

      if (aa1 == "R"/* && is_charge_remote*/)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
        hmm_.setTransitionProbability("R"+pos_name, b_name1, b_sc_int1);
        hmm_.setTransitionProbability("R"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("R"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("R"+pos_name, y_name2, y_sc_int2);
        hmm_.enableTransition("ASC"+pos_name, aa2+"_R"+pos_name);
        hmm_.enableTransition(aa2+"_R"+pos_name, "R"+pos_name);
				hmm_.enableTransition(aa2+"_R"+pos_name, "end"+pos_name);
      }

		}

    Map<HMMState*, double> tmp;
		hmm_.calculateEmissionProbabilities(tmp);

		// clear peaks from last spectrum
		peaks_.clear();

		//stringstream peptide_ss;
		//peptide_ss << peptide;
		//hmm_.writeGraphMLFile(String("model_graph_train_"+peptide_ss.str()+"_"+String(charge)+".graphml").c_str());

		// losses
		vector<Map<String, double> > prefix_losses, suffix_losses;
		vector<double> prefix_ints1, prefix_ints2, suffix_ints1, suffix_ints2;
		double prefix_sum(0), suffix_sum(0);
		for (UInt i = 0; i != prefixes.size(); ++i)
		{
			Map<String, double> prefix_loss;
			double prefix_int1 = tmp[hmm_.getState(prefix_names1[i])];
			prefix_ints1.push_back(prefix_int1);
			double prefix_int2 = tmp[hmm_.getState(prefix_names2[i])];
			prefix_sum = prefix_int1 + prefix_int2;
			//cerr << "prefix loss ints, "<< prefix_names1[i] << " "  << prefix_int1 << " " << prefix_int2 << endl;
			prefix_ints2.push_back(prefix_int2);
			/*
			if (prefixes[i].size() == 2)
			{
				getNeutralLossesFromIon_(prefix_loss, prefix_int1 + prefix_int2, B2Ion, prefixes[i]);
			}
			else
			{
				getNeutralLossesFromIon_(prefix_loss, prefix_int1 + prefix_int2, BIon, prefixes[i]);
			}
			prefix_losses.push_back(prefix_loss);
			*/

			Map<String, double> suffix_loss;
			double suffix_int1 = tmp[hmm_.getState(suffix_names1[i])];
			suffix_ints1.push_back(suffix_int1);
			double suffix_int2 = tmp[hmm_.getState(suffix_names2[i])];
			suffix_sum = suffix_int1 + suffix_int2;
			suffix_ints2.push_back(suffix_int2);
			//cerr << "suffix loss ints, " << suffix_names1[i] << " "  << suffix_int1 << " " <<  suffix_int2 << endl;

			/*
			getNeutralLossesFromIon_(suffix_loss, suffix_int1 + suffix_int2, YIon, suffixes[i]);
			suffix_losses.push_back(suffix_loss);
			*/
		}

		//cerr << "PREFIX_SUM=" << prefix_sum << ", SUFFIX_SUM=" << suffix_sum << endl;

		hmm_.disableTransitions();

		// read the emission probs and put the peaks into a spec
		UInt max_isotope = (UInt)param_.getValue("max_isotope");
		IsotopeDistribution id(max_isotope);
	
		// register name
		RichPeak1D p;
		p.metaRegistry().registerName("IonName", "Name of the ion");
		
		for (UInt i = 0; i != prefixes.size(); ++i)
		{
			// prefix
			double weight = prefixes[i].getMonoWeight(Residue::BIon);
			id.estimateFromPeptideWeight(weight);
			
			//cerr << weight + 1.0 << " " << prefix_ints1[i] << " " << id[0].first << " ";
			//for (vector<pair<UInt, double> >::const_iterator it = id.begin(); it != id.end(); ++it)
			//{
			//	cerr << "(" << it->first << "|" << it->second << ") ";
			//}
			//cerr << endl;
			
			// first isotope peak
			addPeaks_(weight, 1, 0.0, prefix_ints1[i], spec, id, "b"+String(i+1) + "+");
			if (charge >= 2)
			{
				addPeaks_(weight, 2, 0.0, prefix_ints2[i], spec, id, "b"+String(i+1) + "++");

				// neutral losses
				// get fractions as the different charge states are treated together 
				//double loss_1_fraction = prefix_ints1[i] / (prefix_ints1[i] + prefix_ints2[i]);
				//double loss_2_fraction = prefix_ints2[i] / (prefix_ints1[i] + prefix_ints2[i]);
		
				/*
				if (prefix_losses[i].has(LOSS_TYPE_H2O))
				{
					addPeaks_(weight, 1, -18.0, prefix_losses[i][LOSS_TYPE_H2O] * loss_1_fraction, spec, id, "b" + String(i+1) + "-H2O+");

          // doubly charged
					addPeaks_(weight, 2, -18.0, prefix_losses[i][LOSS_TYPE_H2O] * loss_2_fraction, spec, id, "b" + String(i+1) + "-H2O++");
				}
				if (prefix_losses[i].has(LOSS_TYPE_NH3))
				{
					addPeaks_(weight, 1, -17.0, prefix_losses[i][LOSS_TYPE_NH3] * loss_1_fraction, spec, id, "b" + String(i+1) + "-NH3+");
          // doubly charged
					addPeaks_(weight, 2, -17.0, prefix_losses[i][LOSS_TYPE_NH3] * loss_2_fraction, spec, id, "b" + String(i+1) + "-NH3++");
				}

				if (prefix_losses[i].has(LOSS_TYPE_CO))
				{
					addPeaks_(weight, 1, -28.0, prefix_losses[i][LOSS_TYPE_CO] * loss_1_fraction, spec, id, "a" + String(i+1) + "+");
					addPeaks_(weight, 2, -28.0, prefix_losses[i][LOSS_TYPE_CO] * loss_2_fraction, spec, id, "a" + String(i+1) + "++");
				}
				*/
				
			}
			else
			{
				/*
				if (prefix_losses[i].has(LOSS_TYPE_H2O))
        {
					//cerr << weight << " " << prefix_losses[i][LOSS_TYPE_H2O] << "b" << String(i+1) << "-H2O" << endl;
					addPeaks_(weight, 1, -18.0, prefix_losses[i][LOSS_TYPE_H2O], spec, id, "b" + String(i+1) + "-H2O+");
				}
				if (prefix_losses[i].has(LOSS_TYPE_NH3))
        {
					//cerr << weight << " " << prefix_losses[i][LOSS_TYPE_NH3] << "b" << String(i+1) << "-NH3" << endl;
					addPeaks_(weight, 1, -17.0, prefix_losses[i][LOSS_TYPE_NH3], spec, id, "b" + String(i+1) + "-NH3+");
				}

				if (prefix_losses[i].has(LOSS_TYPE_CO))
				{
					addPeaks_(weight, 1, -28.0, prefix_losses[i][LOSS_TYPE_CO], spec, id, "a" + String(i+1) + "+");
				}
				*/
			}

			/*
			// a-ions
			//weight = prefixes[i].getMonoWeight(AIon);
			//id.estimateFromPeptideWeight(weight);
			double a_int1 = tmp[hmm_.getState(a_names1[i])];
			addPeaks_(weight - 28.0, 1, 0.0, a_int1, spec, id, "a" + String(i+1) + "+");

			if (charge >= 2)
			{
				double a_int2 = tmp[hmm_.getState(a_names2[i])];
				addPeaks_(weight - 28.0, 2, 0.0, a_int2, spec, id, "a" + String(i+1) + "++");
			}*/

			// suffix ions
			weight = suffixes[i].getMonoWeight(Residue::YIon);
      id.estimateFromPeptideWeight(weight);
			addPeaks_(weight, 1, 0.0, suffix_ints1[i], spec, id, suffix_names1[i]);
      if (charge >= 2)
      {
				addPeaks_(weight, 2, 0.0, suffix_ints2[i], spec, id, suffix_names2[i]);

				/*
        // neutral losses
        // get fractions as the different charge states are treated together
        double loss_1_fraction = suffix_ints1[i] / (suffix_ints1[i] + suffix_ints2[i]);
        double loss_2_fraction = suffix_ints2[i] / (suffix_ints1[i] + suffix_ints2[i]);

        if (suffix_losses[i].has(LOSS_TYPE_H2O))
	      {
					addPeaks_(weight, 1, -18.0, suffix_losses[i][LOSS_TYPE_H2O] * loss_1_fraction, spec, id, "y" + String(i + 1) + "-H2O+");
          // doubly charged
					addPeaks_(weight, 2, -18.0, suffix_losses[i][LOSS_TYPE_H2O] * loss_2_fraction, spec, id, "y" + String(i + 1) + "-H2O++");
        }

        if (suffix_losses[i].has(LOSS_TYPE_NH3))
        {
					addPeaks_(weight, 1, -17.0, suffix_losses[i][LOSS_TYPE_NH3] * loss_1_fraction, spec, id, "y" + String(i + 1) + "-NH3+");
          // doubly charged
					addPeaks_(weight, 2, -17.0, suffix_losses[i][LOSS_TYPE_NH3] * loss_2_fraction, spec, id, "y" + String(i + 1) + "-NH3++");
        }*/
      }
      else
      {
				/*
        if (suffix_losses[i].has(LOSS_TYPE_H2O))
        {
					//cerr << "H2O: " << suffix_losses[i][LOSS_TYPE_H2O] << " " << suffix_ints1[i] << endl;
					//cerr << "NH3: " << suffix_losses[i][LOSS_TYPE_H2O] << endl;
					addPeaks_(weight, 1, -18.0, suffix_losses[i][LOSS_TYPE_H2O], spec, id, "y" + String(i + 1) + "-H2O+");
        }
        if (suffix_losses[i].has(LOSS_TYPE_NH3))
        {
					addPeaks_(weight, 1, -17.0, suffix_losses[i][LOSS_TYPE_NH3], spec, id, "y" + String(i + 1) + "-NH3+");
        }*/
      }
		}

		if ((is_charge_remote && charge < 3 && !(peptide.has("D") && charge == 2)) || peptide[0].getOneLetterCode() == "Q")
		{
			Map<String, double> pre_ints;
			getPrecursorIons_(pre_ints, max(0.01, 1 - bb_sum - suffix_sum - prefix_sum), peptide);

			double weight = peptide.getMonoWeight();
			id.estimateFromPeptideWeight(weight);
			
			for (Map<String, double>::ConstIterator it = pre_ints.begin(); it != pre_ints.end(); ++it)
			{
				double loss_weight = 0;
				if (it->first == "p-O1H2-O1H2")
				{
					loss_weight = -36.0;
				}
				if (it->first == "p-O1H2")
				{
					loss_weight = -18.0;
				}
				if (it->first == "p-N1H3")
				{
					loss_weight = -17.0;
				}
				if (it->first == "p-N1H3-O1H2")
				{
					loss_weight = -35.0;
				}
				if (it->first == "p-N1H3-N1H3")
				{
					loss_weight = -34.0;
				}
				addPeaks_(weight, charge, loss_weight , it->second, spec, id, it->first);
			}
		}
		
		
		/*
		if ((is_charge_remote && charge < 3 && !(peptide.has("D") && charge == 2)) || peptide[0].getOneLetterCode() == "Q")
		{
			Map<String_, double> pre_ints;
			//cerr << "PRECURSOR_GET: " << peptide << " " <<  bb_sum << " " << 1 - bb_sum - suffix_sum - prefix_sum << " (" << peptide << ", " << charge << ")" << endl;

			if (peptide[0].getOneLetterCode() == "Q")
			{
				getPrecursorIons_(pre_ints, 1, peptide, false);
			}
			else
			{
				getPrecursorIons_(pre_ints, max(0.01, 1 - bb_sum - suffix_sum - prefix_sum), peptide, false);
			}

			double weight = peptide.getMonoWeight();
			id.estimateFromPeptideWeight(weight);
		
			
			if (pre_ints.has(LOSS_TYPE_NH2CHNH))
			{
				addPeaks_(weight, charge, -44.0, pre_ints[LOSS_TYPE_NH2CHNH], spec, id, "M-NH2-CH=NH");
			}
			
			if (pre_ints.has(LOSS_TYPE_H2O))
			{
				addPeaks_(weight, charge, -18.0, pre_ints[LOSS_TYPE_H2O], spec, id, "M-H2O"); 
				//cerr << "Int from Pre: " << pre_ints[LOSS_TYPE_H2O] << endl;
			}

			if (pre_ints.has(LOSS_TYPE_NH3))
			{
				addPeaks_(weight, charge, -17.0, pre_ints[LOSS_TYPE_NH3], spec, id, "M-NH3");
			}

			
			if (pre_ints.has(LOSS_TYPE_NONE))
			{
				addPeaks_(weight, charge, 0.0, pre_ints[LOSS_TYPE_NONE], spec, id, "M");
			}
			
			if (pre_ints.has(LOSS_TYPE_H2O_H2O))
			{
				addPeaks_(weight, charge, -36.0, pre_ints[LOSS_TYPE_H2O_H2O], spec, id, "M-H2O-H2O");
			}

			if (pre_ints.has(LOSS_TYPE_H2O_NH3))
			{
				addPeaks_(weight, charge, -35.0, pre_ints[LOSS_TYPE_H2O_NH3], spec, id, "M-H2O-NH3");
			}
		}
	*/
		// now build the spectrum with the peaks
		double intensity_max(0);
		for (Map<double, vector<RichPeak1D> >::ConstIterator it = peaks_.begin(); it != peaks_.end(); ++it)
		{
			if (it->second.size() == 1/* && it->second.begin()->getIntensity() != 0*/)
			{
				spec.push_back(*it->second.begin());
				if (intensity_max < spec.back().getIntensity())
				{
					intensity_max = spec.back().getIntensity();
				}
			}
			else
			{
				double int_sum(0);
				p = *it->second.begin();
				for (vector<RichPeak1D>::const_iterator pit = it->second.begin(); pit != it->second.end(); ++pit)
				{
					int_sum += pit->getIntensity();
					if (String(pit->getMetaValue("IonName")) != "")
					{
						p.setMetaValue("IonName", pit->getMetaValue("IonName"));
					}
				}

				p.setIntensity(int_sum);
				spec.push_back(p);
				if (intensity_max < int_sum)
				{
					intensity_max = int_sum;
				}
			}
		}

		spec.sortByPosition();

		double min_y_int((double)param_.getValue("min_y_ion_intensity"));
		double min_b_int((double)param_.getValue("min_b_ion_intensity"));
		double min_a_int((double)param_.getValue("min_a_ion_intensity"));
		double min_y_loss_int((double)param_.getValue("min_y_loss_intensity"));
		double min_b_loss_int((double)param_.getValue("min_b_loss_intensity"));

		// TODO switch to enable disable default
		// TODO consider ++ ions
		for (RichPeakSpectrum::Iterator it = spec.begin(); it != spec.end(); ++it)
		{
			it->setIntensity(it->getIntensity() / intensity_max);

			String ion_name(it->getMetaValue("IonName"));
			if (ion_name != "")
			{
				if (ion_name.hasSubstring("y") && (charge > 2 || !ion_name.hasSubstring("++")))
				{
					if (ion_name.hasSubstring("H2O") || ion_name.hasSubstring("NH3"))
					{
						if (it->getIntensity() < min_y_loss_int)
						{
							it->setIntensity(min_y_loss_int);
						}
					}
					else
					{
						if (it->getIntensity() < min_y_int)
						{
							it->setIntensity(min_y_int);
						}
					}
				}

				if (ion_name.hasSubstring("b") && (charge > 2 || !ion_name.hasSubstring("++")))
        {
          if (ion_name.hasSubstring("H2O") || ion_name.hasSubstring("NH3"))
          {
						if (it->getIntensity() < min_b_loss_int)
						{
	            it->setIntensity(min_b_loss_int);
						}
          }
          else
          {
						if (it->getIntensity() < min_b_int)
						{
	            it->setIntensity(min_b_int);
						}
          }
        }

				if (ion_name.hasSubstring("a") && (charge > 2 || !ion_name.hasSubstring("++")))
				{
					if (it->getIntensity() < min_a_int)
					{
						it->setIntensity(min_a_int);
					}
				}
			}
		}

		return;
	}

	void PILISModel::getPrecursorIons_(Map<String, double>& intensities, double initial_probability, const AASequence& precursor)
	{
		hmm_precursor_.setInitialTransitionProbability("start", initial_probability);

		enablePrecursorIonStates_(precursor);

		Map<HMMState*, double> tmp;
		hmm_precursor_.calculateEmissionProbabilities(tmp);

		for (Map<HMMState*, double>::ConstIterator it = tmp.begin(); it != tmp.end(); ++it)
		{
			intensities[it->first->getName()] = it->second;
		}
		
#ifdef PRECURSOR_DEBUG
		for (Map<HMMState*, double>::ConstIterator it = tmp.begin(); it != tmp.end(); ++it)
		{
			cerr << it->first->getName() << " -> " << it->second << endl;
		}		
		//stringstream ss;
		//ss << precursor;
		//hmm_pre_loss_.writeGraphMLFile(String("model_graph_train_"+ss.str()+"_precursor.graphml").c_str());
#endif
		
		hmm_precursor_.disableTransitions();
	}

  void PILISModel::getNeutralLossesFromIon_(Map<String, double>& /*intensities*/, double /*initial_probability*/, IonType_ /*ion_type*/, const AASequence& /*ion*/)
	{
		/*
		HiddenMarkovModelLight* hmm = &hmms_losses_[ion_type];
	
		if (ion_type == BIon || ion_type == B2Ion)
		{
			hmm->setInitialTransitionProbability(B_ION, initial_probability);
			enableNeutralLossStates_(ion_type, ion);

			Map<HMMStateLight*, double> tmp;

			hmm->calculateEmissionProbabilities(tmp);
	
			
			//stringstream ss;
			//ss << ion;
			//hmm->writeGraphMLFile(String("model_graph_train_bion_"+ss.str()+".graphml").c_str());
		
				
			
			if (tmp.has(hmm->getState(B_H2O)))
			{
				intensities[LOSS_TYPE_H2O] = tmp[hmm->getState(B_H2O)];
			}
			if (tmp.has(hmm->getState(B_NH3)))
			{
				intensities[LOSS_TYPE_NH3] = tmp[hmm->getState(B_NH3)];
			}
			if (tmp.has(hmm->getState(A_ION)))
			{
				intensities[LOSS_TYPE_CO] = tmp[hmm->getState(A_ION)];
			}

			hmm->disableTransitions();
			
			return;
		}

		if (ion_type == YIon)
    {
      hmm->setInitialTransitionProbability(Y_ION, initial_probability);
      enableNeutralLossStates_(ion_type, ion);

      Map<HMMStateLight*, double> tmp;
      hmm->calculateEmissionProbabilities(tmp);

			//stringstream ss;
			//ss << ion;
			//hmm->writeGraphMLFile(String("model_graph_train_yion_"+ss.str()+".graphml").c_str());
			
			if (tmp.has(hmm->getState(Y_H2O)))
			{
	      intensities[LOSS_TYPE_H2O] = tmp[hmm->getState(Y_H2O)];
			}
			if (tmp.has(hmm->getState(Y_NH3)))
			{
	      intensities[LOSS_TYPE_NH3] = tmp[hmm->getState(Y_NH3)];
			}

			hmm->disableTransitions();
			
      return;
		}
*/
		return;
	}

	bool PILISModel::getInitialTransitionProbabilities_(std::vector<double>& bb_init, 
																												std::vector<double>& cr_init, 
																												std::vector<double>& sc_init, 
																												const Map<UInt, double>& bb_charges, 
																												const Map<UInt, double>& sc_charges, 
																												const AASequence& peptide)
	{
		bool is_charge_remote(false);

		double charge_dir_tmp(bb_charges[0]);

    double bb_sum(0);
    for (UInt i = 0; i != bb_charges.size(); ++i)
    {
      bb_sum += bb_charges[i];
    }

    charge_dir_tmp = bb_sum;
    if (bb_sum < (double)param_.getValue("charge_directed_threshold"))
    {
      charge_dir_tmp = (double)param_.getValue("charge_directed_threshold");
			is_charge_remote = true;

			bb_sum = (double)param_.getValue("charge_directed_threshold");
    }

    double bb_avg = (bb_sum - bb_charges[0]) / (double)(bb_charges.size() - 1);

					
		for (UInt i = 0; i != peptide.size() - 1; ++i)
    {
      bb_init.push_back(bb_charges[i+1] * charge_dir_tmp/* * (double)param_.getValue("charge_directed_factor")*/);
      String aa(peptide[i].getOneLetterCode());
      if (sc_charges.has(i))
      {
        if ((aa == "K" || aa == "R" || aa == "H") /*&& bb_sum < (double)param_.getValue("charge_remote_threshold")*/)
        {
          sc_init.push_back(sc_charges[i] * bb_avg/* * (double)param_.getValue("side_chain_factor")*/);
        }
        else
        {
          sc_init.push_back(0.0);
        }
        if (bb_sum < (double)param_.getValue("charge_remote_threshold") && (aa == "D" || aa == "E"))
        {
          cr_init.push_back(((1 - bb_sum) * bb_avg /* sc_charge[i]*/)/* * (double)param_.getValue("charge_remote_factor")*/);
        }
        else
        {
          cr_init.push_back(0.0);
        }
      }
      else
      {
        sc_init.push_back(0.0);
        cr_init.push_back(0.0);
      }
    }

		// normalize the initial probability values
		double init_prob_sum(0);
		for (UInt i = 0; i != bb_init.size(); ++i)
		{
			init_prob_sum += bb_init[i] + sc_init[i] + cr_init[i];
		}
		// C-term bb positional acmino acid might have a proton associated with (is inactive, as now pathway uses it)
		if (sc_charges.has(peptide.size() - 1))
		{
			init_prob_sum += sc_charges[peptide.size() - 1];
		}
		for (UInt i = 0; i != bb_init.size(); ++i)
		{
			bb_init[i] /= init_prob_sum;
			sc_init[i] /= init_prob_sum;
			cr_init[i] /= init_prob_sum;
		}
		
		//return is_charge_remote;
		return true;
	}

	void PILISModel::addPeaks_(double mz, int charge, double offset, double intensity, RichPeakSpectrum& /*spectrum*/, const IsotopeDistribution& id, const String& name)
	{
		static RichPeak1D p;
		UInt i = 0;
		for (IsotopeDistribution::ConstIterator it = id.begin(); it != id.end(); ++it, ++i)
		{
			double pos = (mz + i + charge + offset) / (double)charge;
			p.setPosition(pos);
			if (it == id.begin())
			{
				p.setMetaValue("IonName", String(name.c_str()));
			}

			if (pos >= (double)param_.getValue("lower_mz") && pos <= (double)param_.getValue("upper_mz"))
			{
				p.setIntensity(intensity * it->second);
				peaks_[p.getMZ()].push_back(p);
			}

			if (it == id.begin())
			{
				p.setMetaValue("IonName", String(""));
			}
		}

		return;
	}
	
	
	void PILISModel::initModels_()
	{
		return;
	}

	void PILISModel::parseHMMModel_(const TextFile::ConstIterator& begin, const TextFile::ConstIterator& end, HiddenMarkovModel& hmm)
	{
		if (begin == end)
    {
      return;
    }
		
		//Size num_syn(0);
    for (TextFile::ConstIterator it = begin; it != end; ++it)
    {
      String line = *it;
      // comment?
			//cerr << line << endl;
      if (line[0] == '#')
      {
        continue;
      }

      vector<String> split;
      line.split(' ', split);

      if (split.size() > 0)
      {
        String id = split[0];

        if (id == "State")
        {
          bool hidden(true);
          if (split.size() != 2 && split[2] == "false")
          {
            hidden = false;
          }
          hmm.addNewState(new HMMState(split[1], hidden));
					//cerr << "added new state: '" << split[1] << "', " << hidden << endl;
          continue;
        }

        if (id == "Synonym")
        {
					//++num_syn;
					//hmm.addSynonymTransition(split[1], split[2], split[3], split[4]);
					hmm.addSynonymTransition(split[3], split[4], split[1], split[2]);
          continue;
        }

        if (id == "Transition")
        {
          hmm.setTransitionProbability(split[1], split[2], split[3].toFloat());
          continue;
        }
      }
    }
    //hmm_.disableTransitions();
		hmm.buildSynonyms();

		hmm.disableTransitions();
	
		//cerr << hmm_.getNumberOfStates() << endl;
		
		return;
	}

	void PILISModel::updateMembers_()
	{
		double pseudo_counts = (double)param_.getValue("pseudo_counts");
		hmm_.setPseudoCounts(pseudo_counts);

		
		hmm_precursor_.setPseudoCounts(pseudo_counts);

		/*
		hmm_pre_loss_.setPseudoCounts(pseudo_counts);
		for (Map<IonType_, HiddenMarkovModelLight>::Iterator it = hmms_losses_.begin(); it != hmms_losses_.end(); ++it)
		{
			it->second.setPseudoCounts(pseudo_counts);
		}
		*/
	}
} // namespace OpenMS


