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
#include <OpenMS/ANALYSIS/ID/PILISModelGenerator.h>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/SYSTEM/File.h>

#include <cmath>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <fstream>


#define TRAINING_DEBUG
#undef  TRAINING_DEBUG

#define SIM_DEBUG
#undef  SIM_DEBUG

using namespace std;

// new TODOS
// CR and CD distinction of loss models


// old ones
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
		defaults_.setValue("charge_remote_threshold", 0.15, "If the probability for the proton at the N-terminus is lower than this value, enable charge-remote pathways");
		defaults_.setValue("charge_directed_threshold", 0.15, "Limit the proton availability at the N-terminus to at least this value for charge-directed pathways");
		defaults_.setValue("model_depth", 4, "The number of explicitly modeled backbone cleavages from N-terminus and C-terminus, would be 9 for the default value");
		defaults_.setValue("visible_model_depth", 30, "The maximal possible size of a peptide to be modeled");
		defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, used to identify the precursor peak and its loss peaks for training");
		defaults_.setValue("fragment_mass_tolerance", 0.3, "Peak mass tolerance of the product ions, used to identify the ions for training");
		defaults_.setValue("variable_modifications", StringList::create("MOD:00719,MOD:09997"), "Modifications which should be included in the model, represented by PSI-MOD accessions.");
		defaults_.setValue("fixed_modifications", StringList::create(""), "Modifications which should replace the unmodified amino acid, represented by PSI-MOD accessions.");
						
		defaults_.setValue("min_y_ion_intensity", 0.0, ".");
		defaults_.setValue("min_b_ion_intensity", 0.0, ".");
		defaults_.setValue("min_a_ion_intensity", 0.0, ".");
		defaults_.setValue("min_y_loss_intensity", 0.0, ".");
		defaults_.setValue("min_b_loss_intensity", 0.0, ".");

		defaults_.setValue("charge_loss_factor", 0.2, "Factor which accounts for the loss of a proton from a ++ ion, the higher the value the more ++ ions loose protons");
		defaults_.setValue("side_chain_activation", 0.3, "Additional activation of proton sitting at side chain, especially important at lysin and histidine residues");
		defaults_.setValue("pseudo_counts", 1e-15, "Value which is added for every transition trained of the underlying hidden Markov model");
		defaults_.setValue("max_isotope", 2, "Maximal isotope peak which is reported by the model, 0 means all isotope peaks are reported");
		defaults_.setValue("max_fragment_charge", 2, "Maximal charge state allowed for fragment ions");
		defaults_.setValue("suppress_single_ions", "true", "If set to true, single ions are suppressed, e.g. y1 and b1");

		defaultsToParam_();
	}

	PILISModel::~PILISModel()
	{
	}

	PILISModel::PILISModel(const PILISModel& model)
		: DefaultParamHandler(model),
			hmm_(model.hmm_),
			prot_dist_(model.prot_dist_),
			tsg_(model.tsg_),
			valid_(model.valid_),
			peaks_(model.peaks_),
			spectra_aligner_(model.spectra_aligner_),
			precursor_model_cr_(model.precursor_model_cr_),
			precursor_model_cd_(model.precursor_model_cd_),
			
			a_ion_losses_cr_(model.a_ion_losses_cr_),
			a_ion_losses_cd_(model.a_ion_losses_cd_),
			
			b_ion_losses_cr_(model.b_ion_losses_cr_),
			b_ion_losses_cd_(model.b_ion_losses_cd_),
			
			b2_ion_losses_cr_(model.b2_ion_losses_cr_),
			b2_ion_losses_cd_(model.b2_ion_losses_cd_),
			
			y_ion_losses_cr_(model.y_ion_losses_cr_),
			y_ion_losses_cd_(model.y_ion_losses_cd_)
	{
	}

	PILISModel& PILISModel::operator = (const PILISModel& model)
	{
		if (this != &model)
		{
			DefaultParamHandler::operator=(model);
			hmm_ = model.hmm_;
	    prot_dist_ = model.prot_dist_;
	    tsg_ = model.tsg_;
	    valid_ = model.valid_;
			peaks_ = model.peaks_;
			spectra_aligner_ = model.spectra_aligner_;
			precursor_model_cr_ = model.precursor_model_cr_;
			precursor_model_cd_ = model.precursor_model_cd_;
			
			a_ion_losses_cr_ = model.a_ion_losses_cr_;
			a_ion_losses_cd_ = model.a_ion_losses_cd_;
			
			b_ion_losses_cr_ = model.b_ion_losses_cr_;
			b_ion_losses_cd_ = model.b_ion_losses_cd_;
      
			b2_ion_losses_cr_ = model.b2_ion_losses_cr_;
			b2_ion_losses_cd_ = model.b2_ion_losses_cd_;
			
			y_ion_losses_cr_ = model.y_ion_losses_cr_;
			y_ion_losses_cd_ = model.y_ion_losses_cd_;
		}
		return *this;
	}

	void PILISModel::setHMM(const HiddenMarkovModel& model)
	{
		hmm_ = model;
	}

	void PILISModel::init(bool generate_models)
	{
		if (generate_models)
		{
			PILISModelGenerator gen;
			Param gen_param(gen.getParameters());
			gen_param.setValue("variable_modifications", (StringList)param_.getValue("variable_modifications"));
			gen_param.setValue("fixed_modifications", (StringList)param_.getValue("fixed_modifications"));
			gen_param.setValue("model_depth", (UInt)param_.getValue("model_depth"));
			gen_param.setValue("visible_model_depth", (UInt)param_.getValue("visible_model_depth"));
			gen.setParameters(gen_param);
			gen.getModel(hmm_);
		}


		Param pre_param(precursor_model_cr_.getParameters());
		pre_param.setValue("fragment_mass_tolerance", (double)param_.getValue("fragment_mass_tolerance"));
		pre_param.setValue("variable_modifications", (StringList)param_.getValue("variable_modifications"));
		pre_param.setValue("fixed_modifications", (StringList)param_.getValue("fixed_modifications"));
		pre_param.setValue("pseudo_counts", (double)param_.getValue("pseudo_counts"));
		pre_param.setValue("C_term_H2O_loss", "true");
		pre_param.setValue("ion_name", "p");
		precursor_model_cr_.setParameters(pre_param);
		precursor_model_cd_.setParameters(pre_param);
	
		if (generate_models)
		{
			precursor_model_cr_.generateModel();
			precursor_model_cd_.generateModel();
		}

		pre_param.setValue("C_term_H2O_loss", "false");
		pre_param.setValue("enable_double_losses", "false");
		pre_param.setValue("ion_name", "b");
		b_ion_losses_cr_.setParameters(pre_param);
		b_ion_losses_cd_.setParameters(pre_param);
		if (generate_models)
		{
			b_ion_losses_cr_.generateModel();
			b_ion_losses_cd_.generateModel();
		}

    pre_param.setValue("C_term_H2O_loss", "false");
    pre_param.setValue("enable_double_losses", "false");
    pre_param.setValue("ion_name", "b2");
    b2_ion_losses_cr_.setParameters(pre_param);
		b2_ion_losses_cd_.setParameters(pre_param);
    if (generate_models)
    {
      b2_ion_losses_cr_.generateModel();
			b2_ion_losses_cd_.generateModel();
    }

    pre_param.setValue("C_term_H2O_loss", "false");
    pre_param.setValue("enable_double_losses", "false");
    pre_param.setValue("ion_name", "a");
    a_ion_losses_cr_.setParameters(pre_param);
		a_ion_losses_cd_.setParameters(pre_param);
		if (generate_models)
    {
      a_ion_losses_cr_.generateModel();
			a_ion_losses_cd_.generateModel();
    }
		
		pre_param.setValue("C_term_H2O_loss", "true");
		pre_param.setValue("ion_name", "y");
		y_ion_losses_cr_.setParameters(pre_param);
		y_ion_losses_cd_.setParameters(pre_param);
		if (generate_models)
		{
			y_ion_losses_cr_.generateModel();
			y_ion_losses_cd_.generateModel();
		}

		if (generate_models)
		{
			valid_ = true;
		}
	}
				
	
	void PILISModel::readFromFile(const String& filename)
	{
		// read the model
    if (!File::readable(filename))
    {
     	throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    if (File::empty(filename))
    {
     	throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

		init(false);
		
		TextFile file;
		file.load(filename, true);

		TextFile::Iterator it_begin(file.begin()), it_end(file.begin());
		it_begin = file.search(it_begin, "BASE_MODEL_BEGIN");
		it_end = file.search(it_begin, "BASE_MODEL_END");
		parseHMMModel_(++it_begin, it_end, hmm_);

		// seek to next interval
		it_begin = file.search(it_end, "PRECURSOR_MODEL_CR_BEGIN");
		it_end = file.search(it_begin, "PRECURSOR_MODEL_CR_END");
		HiddenMarkovModel precursor_model_cr_hmm;
		parseHMMModel_(++it_begin, it_end, precursor_model_cr_hmm);
		precursor_model_cr_.setHMM(precursor_model_cr_hmm);

    it_begin = file.search(it_end, "PRECURSOR_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "PRECURSOR_MODEL_CD_END");
		HiddenMarkovModel precursor_model_cd_hmm;
    parseHMMModel_(++it_begin, it_end, precursor_model_cd_hmm);
		precursor_model_cd_.setHMM(precursor_model_cd_hmm);
		
		// b loss models
    it_begin = file.search(it_end, "BION_LOSS_MODEL_CR_BEGIN");
    it_end = file.search(it_begin, "BION_LOSS_MODEL_CR_END");
    HiddenMarkovModel b_ion_loss_hmm_cr;
    parseHMMModel_(++it_begin, it_end, b_ion_loss_hmm_cr);
    b_ion_losses_cr_.setHMM(b_ion_loss_hmm_cr);

    it_begin = file.search(it_end, "BION_LOSS_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "BION_LOSS_MODEL_CD_END");
    HiddenMarkovModel b_ion_loss_hmm_cd;
    parseHMMModel_(++it_begin, it_end, b_ion_loss_hmm_cd);
    b_ion_losses_cd_.setHMM(b_ion_loss_hmm_cd);

		it_begin = file.search(it_end, "B2ION_LOSS_MODEL_CR_BEGIN");
		it_end = file.search(it_begin, "B2ION_LOSS_MODEL_CR_END");
		HiddenMarkovModel b2_ion_loss_hmm_cr;
		parseHMMModel_(++it_begin, it_end, b2_ion_loss_hmm_cr);
		b2_ion_losses_cr_.setHMM(b2_ion_loss_hmm_cr);

    it_begin = file.search(it_end, "B2ION_LOSS_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "B2ION_LOSS_MODEL_CD_END");
    HiddenMarkovModel b2_ion_loss_hmm_cd;
    parseHMMModel_(++it_begin, it_end, b2_ion_loss_hmm_cd);
    b2_ion_losses_cd_.setHMM(b2_ion_loss_hmm_cd);
	
		// a-ion loss model
		it_begin = file.search(it_end, "AION_LOSS_MODEL_CR_BEGIN");
		it_end = file.search(it_begin, "AION_LOSS_MODEL_CR_END");
		HiddenMarkovModel a_ion_loss_hmm_cr;
		parseHMMModel_(++it_begin, it_end, a_ion_loss_hmm_cr);
		a_ion_losses_cr_.setHMM(a_ion_loss_hmm_cr);

    it_begin = file.search(it_end, "AION_LOSS_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "AION_LOSS_MODEL_CD_END");
    HiddenMarkovModel a_ion_loss_hmm_cd;
    parseHMMModel_(++it_begin, it_end, a_ion_loss_hmm_cd);
    a_ion_losses_cd_.setHMM(a_ion_loss_hmm_cd);
		
		// y-ion loss model
		it_begin = file.search(it_end, "YION_LOSS_MODEL_CR_BEGIN");
		it_end = file.search(it_begin, "YION_LOSS_MODEL_CR_END");
		HiddenMarkovModel y_ion_loss_hmm_cr;
		parseHMMModel_(++it_begin, it_end, y_ion_loss_hmm_cr);
		y_ion_losses_cr_.setHMM(y_ion_loss_hmm_cr);

    it_begin = file.search(it_end, "YION_LOSS_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "YION_LOSS_MODEL_CD_END");
    HiddenMarkovModel y_ion_loss_hmm_cd;
    parseHMMModel_(++it_begin, it_end, y_ion_loss_hmm_cd);
    y_ion_losses_cd_.setHMM(y_ion_loss_hmm_cd);

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

		out << "PRECURSOR_MODEL_CR_BEGIN" << endl;
		precursor_model_cr_.getHMM().write(out);
		out << "PRECURSOR_MODEL_CR_END" << endl;
		
		out << "PRECURSOR_MODEL_CD_BEGIN" << endl;
		precursor_model_cd_.getHMM().write(out);
		out << "PRECURSOR_MODEL_CD_END" << endl;

		out << "BION_LOSS_MODEL_CR_BEGIN" << endl;
		b_ion_losses_cr_.getHMM().write(out);
		out << "BION_LOSS_MODEL_CR_END" << endl;

    out << "BION_LOSS_MODEL_CD_BEGIN" << endl;
    b_ion_losses_cd_.getHMM().write(out);
    out << "BION_LOSS_MODEL_CD_END" << endl;

		out << "B2ION_LOSS_MODEL_CR_BEGIN" << endl;
		b2_ion_losses_cr_.getHMM().write(out);
		out << "B2ION_LOSS_MODEL_CR_END" << endl;

    out << "B2ION_LOSS_MODEL_CD_BEGIN" << endl;
    b2_ion_losses_cd_.getHMM().write(out);
    out << "B2ION_LOSS_MODEL_CD_END" << endl;


    out << "AION_LOSS_MODEL_CR_BEGIN" << endl;
    a_ion_losses_cr_.getHMM().write(out);
    out << "AION_LOSS_MODEL_CR_END" << endl;

    out << "AION_LOSS_MODEL_CD_BEGIN" << endl;
    a_ion_losses_cd_.getHMM().write(out);
    out << "AION_LOSS_MODEL_CD_END" << endl;
		
		out << "YION_LOSS_MODEL_CR_BEGIN" << endl;
		y_ion_losses_cr_.getHMM().write(out);
		out << "YION_LOSS_MODEL_CR_END" << endl;

    out << "YION_LOSS_MODEL_CD_BEGIN" << endl;
    y_ion_losses_cd_.getHMM().write(out);
    out << "YION_LOSS_MODEL_CD_END" << endl;
	
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
		//double peptide_weight((peptide.getMonoWeight() + double(charge)) / double(charge));
		
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

		/*
		for (Map<String, double>::ConstIterator it = pre_ints.begin(); it != pre_ints.end(); ++it)
		{
			sum += it->second;
		}*/
		
		/*
		for (Map<String, double>::Iterator it = pre_ints.begin(); it != pre_ints.end(); ++it)
		{
			it->second /= sum;
		}*/

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
		double side_chain_activation(param_.getValue("side_chain_activation"));
		for (UInt i = 0; i != peptide.size(); ++i)
		{
			if (peptide[i].getOneLetterCode() != "R") // || peptide[i].getOneLetterCode() == "K")
			{
				bb_sum += side_chain_activation * sc_charge_full[i]; // TODO
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
			
			if (i < floor(peptide.size()/2.0))
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

			UInt max_fragment_charge = param_.getValue("max_fragment_charge");
			for (UInt z = 1; z <= max_fragment_charge; ++z)
			{
				double charge_remote_threshold((double)param_.getValue("charge_remote_threshold"));
				double avail_bb_sum_prefix = getAvailableBackboneCharge_(prefix, Residue::BIon, z);
				double avail_bb_sum_suffix = getAvailableBackboneCharge_(suffix, Residue::YIon, z);
			
				if (prefix.size() != 2)
				{
					double pre_weight = prefix.getMonoWeight(Residue::BIon) + 1.0;

					if (avail_bb_sum_prefix <= charge_remote_threshold)
					{
						//cerr << "Train b-ion losses, CR (avail=" << avail_bb_sum_prefix << ", z=" << z << ", prefix=" << prefix << endl;
						b_ion_losses_cr_.train(in_spec, prefix, pre_weight, z);
					}
					else
					{
						//cerr << "Train b-ion losses, CD (avail=" << avail_bb_sum_prefix << ", z=" << z << ", prefix=" << prefix << endl;
						b_ion_losses_cd_.train(in_spec, prefix, pre_weight, z);
					}
				}
				else
				{
					double pre_weight = prefix.getMonoWeight(Residue::BIon) + 1.0;

					if (avail_bb_sum_prefix <= charge_remote_threshold)
					{
						b2_ion_losses_cr_.train(in_spec, prefix, pre_weight, z);
					}
					else
					{
						b2_ion_losses_cd_.train(in_spec, prefix, pre_weight, z);
					}
				}

				double a_pre_weight = prefix.getMonoWeight(Residue::AIon) + 1.0;
				if (avail_bb_sum_prefix <= charge_remote_threshold)
				{
					a_ion_losses_cr_.train(in_spec, prefix, a_pre_weight, z);
				}
				else
				{
					a_ion_losses_cd_.train(in_spec, prefix, a_pre_weight, z);
				}
			
				double suf_weight = suffix.getMonoWeight(Residue::YIon) + 1.0;
				if (avail_bb_sum_suffix <= charge_remote_threshold)
				{
					//cerr << "Train y-ion losses, CR (avail=" << avail_bb_sum_suffix << ", z=" << z << ", suffix=" << suffix << endl;
					y_ion_losses_cr_.train(in_spec, suffix, suf_weight, z);
				}
				else
				{
					//cerr << "Train y-ion losses, CD (avail=" << avail_bb_sum_suffix << ", z=" << z << ", suffix=" << suffix << endl;
					y_ion_losses_cd_.train(in_spec, suffix, suf_weight, z);
				}
			}
		}
	
		double charge_remote_threshold((double)param_.getValue("charge_remote_threshold"));
		if (bb_sum <= charge_remote_threshold)
		{
			double peptide_weight((peptide.getMonoWeight() + charge) / (double)charge);
			precursor_model_cr_.train(in_spec, peptide, peptide_weight, charge);
		}
		else
		{
			double peptide_weight((peptide.getMonoWeight() + charge) / (double)charge);
			precursor_model_cd_.train(in_spec, peptide, peptide_weight, charge);
		}

		// now train the model with the data set
		hmm_.train();
		hmm_.disableTransitions();

		return;
	}

	// TODO Code cleanup, only b, y and a ions are needed
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
	
	void PILISModel::evaluate()
	{
		hmm_.evaluate();
		precursor_model_cr_.evaluate();
		precursor_model_cd_.evaluate();
		a_ion_losses_cr_.evaluate();
		a_ion_losses_cd_.evaluate();
		b_ion_losses_cr_.evaluate();
		b_ion_losses_cd_.evaluate();
		b2_ion_losses_cr_.evaluate();
		b2_ion_losses_cd_.evaluate();
		y_ion_losses_cr_.evaluate();
		y_ion_losses_cd_.evaluate();
		
		hmm_.estimateUntrainedTransitions(); // TODO reimplement in HiddenMarkovModel
	}

	void PILISModel::getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, UInt charge)
	{
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


      if (i < floor(peptide.size()/2.0))
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
		vector<double> prefix_ints1, prefix_ints2, suffix_ints1, suffix_ints2, a_ints1;
		double prefix_sum(0), suffix_sum(0);
		for (UInt i = 0; i != prefixes.size(); ++i)
		{
			double prefix_int1 = tmp[hmm_.getState(prefix_names1[i])];
			//cerr << prefix_names1[i] << " " << prefix_int1 << endl;
			prefix_ints1.push_back(prefix_int1);
			
			double prefix_int2 = tmp[hmm_.getState(prefix_names2[i])];
			//cerr << prefix_names2[i] << " " << prefix_int2 << endl;
			prefix_sum = prefix_int1 + prefix_int2;
			//cerr << "prefix loss ints, "<< prefix_names1[i] << " "  << prefix_int1 << " " << prefix_int2 << endl;
			prefix_ints2.push_back(prefix_int2);

			double a_int1 = tmp[hmm_.getState(a_names1[i])];
			a_ints1.push_back(a_int1);
			
			double suffix_int1 = tmp[hmm_.getState(suffix_names1[i])];
			//cerr << suffix_names1[i] << " " << suffix_int1 << endl;
			suffix_ints1.push_back(suffix_int1);
			
			double suffix_int2 = tmp[hmm_.getState(suffix_names2[i])];
			//cerr << suffix_names2[i] << " " << suffix_int2 << endl;
			suffix_sum = suffix_int1 + suffix_int2;
			suffix_ints2.push_back(suffix_int2);
			//cerr << "suffix loss ints, " << suffix_names1[i] << " "  << suffix_int1 << " " <<  suffix_int2 << endl;
		}

		//cerr << "PREFIX_SUM=" << prefix_sum << ", SUFFIX_SUM=" << suffix_sum << endl;

		hmm_.disableTransitions();

		// read the emission probs and put the peaks into a spec
		UInt max_isotope = (UInt)param_.getValue("max_isotope");
		IsotopeDistribution id(max_isotope);
	
		// register name
		RichPeak1D p;
		//p.metaRegistry().registerName("IonName", "Name of the ion");
		
		double charge_remote_threshold((double)param_.getValue("charge_remote_threshold"));
		bool suppress_single_ions(param_.getValue("suppress_single_ions").toBool());
		for (UInt i = 0; i != prefixes.size(); ++i)
		{
			// prefix
			double weight = prefixes[i].getMonoWeight(Residue::BIon);
			//id.estimateFromPeptideWeight(weight);
			id = prefixes[i].getFormula(Residue::BIon).getIsotopeDistribution(2);
			
			// first isotope peak
			//addPeaks_(weight, 1, 0.0, prefix_ints1[i], spec, id, "b"+String(i+1) + "+");

      vector<RichPeak1D> b_loss_peaks;

			double avail_bb_sum_prefix1 = getAvailableBackboneCharge_(prefixes[i], Residue::BIon, 1);
			if (avail_bb_sum_prefix1 <= charge_remote_threshold)
			{
				b_ion_losses_cr_.getIons(b_loss_peaks, prefixes[i], prefix_ints1[i]);
			}
			else
			{
				b_ion_losses_cd_.getIons(b_loss_peaks, prefixes[i], prefix_ints1[i]);
			}
      for (vector<RichPeak1D>::const_iterator it = b_loss_peaks.begin(); it != b_loss_peaks.end(); ++it)
      {
				String b_ion_name = it->getMetaValue("IonName");
				vector<String> split;
				b_ion_name.split('-', split);
				if (split.size() == 0)
				{
					b_ion_name += String(i + 1);
				}
				else
				{
					b_ion_name = split[0] + String(i + 1) + "-";
					for (UInt j = 1; j != split.size(); ++j)
					{
						b_ion_name += split[j];
					}
				}

				if (!suppress_single_ions || (suppress_single_ions && i + 1 > 1))
				{
					b_ion_name += "+";
        	addPeaks_(weight, 1, it->getMZ(), it->getIntensity(), spec, id, b_ion_name);
				}
      }

			if (charge > 1)
			{
				b_loss_peaks.clear();
				double avail_bb_sum_prefix2 = getAvailableBackboneCharge_(prefixes[i], Residue::BIon, 2);
				if (avail_bb_sum_prefix2 <= charge_remote_threshold)
				{
					b_ion_losses_cr_.getIons(b_loss_peaks, prefixes[i], prefix_ints2[i]);
				}
				else
				{
					b_ion_losses_cd_.getIons(b_loss_peaks, prefixes[i], prefix_ints2[i]);
				}
      	for (vector<RichPeak1D>::const_iterator it = b_loss_peaks.begin(); it != b_loss_peaks.end(); ++it)
      	{
        	String b_ion_name = it->getMetaValue("IonName");
        	vector<String> split;
        	b_ion_name.split('-', split);
        	if (split.size() == 0)
        	{
          	b_ion_name += String(i + 1);
        	}
        	else
       	 	{
          	b_ion_name = split[0] + String(i + 1) + "-";
          	for (UInt j = 1; j != split.size(); ++j)
          	{
            	b_ion_name += split[j];
          	}
        	}

					if (!suppress_single_ions || (suppress_single_ions && i + 1 > 1))
					{
						b_ion_name += "++";
        		addPeaks_(weight, 2, it->getMZ(), it->getIntensity(), spec, id, b_ion_name);
					}
      	}
			}

			weight = prefixes[i].getMonoWeight(Residue::AIon);
      vector<RichPeak1D> a_loss_peaks;
			if (avail_bb_sum_prefix1 <= charge_remote_threshold)
			{
				a_ion_losses_cr_.getIons(a_loss_peaks, prefixes[i], a_ints1[i]);
			}
			else
			{
				a_ion_losses_cd_.getIons(a_loss_peaks, prefixes[i], a_ints1[i]);
			}
      for (vector<RichPeak1D>::const_iterator it = a_loss_peaks.begin(); it != a_loss_peaks.end(); ++it)
      {
        String a_ion_name = it->getMetaValue("IonName");
        vector<String> split;
        a_ion_name.split('-', split);
        if (split.size() == 0)
        {
          a_ion_name += String(i + 1);
        }
        else
        {
          a_ion_name = split[0] + String(i + 1) + "-";
          for (UInt j = 1; j != split.size(); ++j)
          {
            a_ion_name += split[j];
          }
        }

				if (!suppress_single_ions || (suppress_single_ions && i + 1 > 1))
				{
        	a_ion_name += "+";
	        addPeaks_(weight, 1, it->getMZ(), it->getIntensity(), spec, id, a_ion_name);
				}
      }

			// suffix ions
			weight = suffixes[i].getMonoWeight(Residue::YIon);
      //id.estimateFromPeptideWeight(weight);
			id = suffixes[i].getFormula(Residue::YIon).getIsotopeDistribution(2);
			// neutral losses
      vector<RichPeak1D> y_loss_peaks;

			double avail_bb_sum_suffix1 = getAvailableBackboneCharge_(suffixes[i], Residue::YIon, 1);
			if (avail_bb_sum_suffix1 <= charge_remote_threshold)
			{
      	y_ion_losses_cr_.getIons(y_loss_peaks, suffixes[i], suffix_ints1[i]);
			}
			else
			{
				y_ion_losses_cd_.getIons(y_loss_peaks, suffixes[i], suffix_ints1[i]);
			}
			for (vector<RichPeak1D>::const_iterator it = y_loss_peaks.begin(); it != y_loss_peaks.end(); ++it)
			{
				String y_ion_name = it->getMetaValue("IonName");
        vector<String> split;
        y_ion_name.split('-', split);
				UInt num = suffixes.size() - i;
        if (split.size() == 0)
        {
          y_ion_name += String(num);
        }
        else
        {
          y_ion_name = split[0] + String(num) + "-";
          for (UInt j = 1; j != split.size(); ++j)
          {
            y_ion_name += split[j];
          }
        }
					
				if (!suppress_single_ions || (suppress_single_ions && num > 1))
				{							
					y_ion_name += "+";
					addPeaks_(weight, 1, it->getMZ(), it->getIntensity(), spec, id, y_ion_name);
				}
			}
		
			if (charge > 1)
			{
				y_loss_peaks.clear();
				double avail_bb_sum_suffix2 = getAvailableBackboneCharge_(suffixes[i], Residue::YIon, 2);
				if (avail_bb_sum_suffix2 <= charge_remote_threshold)
				{
					y_ion_losses_cr_.getIons(y_loss_peaks, suffixes[i], suffix_ints2[i]);
				}
				else
				{
					y_ion_losses_cd_.getIons(y_loss_peaks, suffixes[i], suffix_ints2[i]);
				}
      	for (vector<RichPeak1D>::const_iterator it = y_loss_peaks.begin(); it != y_loss_peaks.end(); ++it)
      	{
        	String y_ion_name = it->getMetaValue("IonName");
        	vector<String> split;
        	y_ion_name.split('-', split);
					UInt num = suffixes.size() - i;
        	if (split.size() == 0)
        	{
          	y_ion_name += String(num);
        	}
        	else
        	{
          	y_ion_name = split[0] + String(num) + "-";
          	for (UInt j = 1; j != split.size(); ++j)
         		{
            	y_ion_name += split[j];
          	}
        	}

					if (!suppress_single_ions || (suppress_single_ions && num > 1))
					{
						y_ion_name += "++";
	        	addPeaks_(weight, 2, it->getMZ(), it->getIntensity(), spec, id, y_ion_name);
					}
				}
			}
		}

		// precursor intensities
		vector<RichPeak1D> pre_peaks;
		//cerr << peptide << " ->getSpectrum: bb_sum=" << bb_sum << ", charge_remote=" << is_charge_remote<< endl;
		if (bb_sum <= charge_remote_threshold/*(is_charge_remote && charge < 3 && !(peptide.has("D") && charge == 2)) || peptide[0].getOneLetterCode() == "Q"*/)
		{
			precursor_model_cr_.getIons(pre_peaks, peptide, bb_sum);
		}
		else
		{
			precursor_model_cd_.getIons(pre_peaks, peptide, max(0.0, 1.0 - bb_sum));
		}
		
		double weight = peptide.getMonoWeight();
		id.estimateFromPeptideWeight(weight);	
		for (vector<RichPeak1D>::const_iterator it = pre_peaks.begin(); it != pre_peaks.end(); ++it)
		{
			addPeaks_(weight, charge, it->getMZ(), it->getIntensity(), spec, id, it->getMetaValue("IonName"));
		}
		
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
		//double min_a_int((double)param_.getValue("min_a_ion_intensity"));
		//double min_y_loss_int((double)param_.getValue("min_y_loss_intensity"));
		//double min_b_loss_int((double)param_.getValue("min_b_loss_intensity"));
		

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
								/*
					if (ion_name.hasSubstring("H2O") || ion_name.hasSubstring("NH3"))
					{
						if (it->getIntensity() < min_y_loss_int)
						{
							it->setIntensity(min_y_loss_int);
						}
					}
					else
					{*/
						if (it->getIntensity() < min_y_int)
						{
							it->setIntensity(min_y_int);
						}
					//}
				}

				if (ion_name.hasSubstring("b") && (charge > 2 || !ion_name.hasSubstring("++")))
        {
								/*
          if (ion_name.hasSubstring("H2O") || ion_name.hasSubstring("NH3"))
          {
						if (it->getIntensity() < min_b_loss_int)
						{
	            it->setIntensity(min_b_loss_int);
						}
          }
          else
          {*/
						if (it->getIntensity() < min_b_int)
						{
	            it->setIntensity(min_b_int);
						}
          //}
        }
/*
				if (ion_name.hasSubstring("a") && (charge > 2 || !ion_name.hasSubstring("++")))
				{
					if (it->getIntensity() < min_a_int)
					{
						it->setIntensity(min_a_int);
					}
				}*/
			}
			
		}

		return;
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


	double PILISModel::getAvailableBackboneCharge_(const AASequence& ion, Residue::ResidueType res_type, int charge)
	{
		double bb_sum(0);
		Map<UInt, double> bb_charges, sc_charges;
		prot_dist_.getProtonDistribution(bb_charges, sc_charges, ion, charge, res_type);

		for (Map<UInt, double>::ConstIterator it = bb_charges.begin(); it != bb_charges.end(); ++it)
		{
			bb_sum += it->second;
		}

    // activation of protons sitting at lysine and histidine side chains
    double side_chain_activation(param_.getValue("side_chain_activation"));
    for (UInt i = 0; i != ion.size(); ++i)
    {
      if (ion[i].getOneLetterCode() != "R") // || peptide[i].getOneLetterCode() == "K")
      {
        bb_sum += side_chain_activation * sc_charges[i]; // TODO
      }
    }

    if  (bb_sum > 1)
    {
      bb_sum = 1;
    }

    if (bb_sum < (double)param_.getValue("charge_directed_threshold") * charge)
    {
      bb_sum = (double)param_.getValue("charge_directed_threshold") * charge;
    }
		return bb_sum;
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
		
		return is_charge_remote;
		//return true;
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
		//hmm.buildSynonyms();
		hmm.disableTransitions();
	
		//cerr << hmm_.getNumberOfStates() << endl;
		
		return;
	}

	void PILISModel::updateMembers_()
	{
		double pseudo_counts = (double)param_.getValue("pseudo_counts");
		hmm_.setPseudoCounts(pseudo_counts);

		// pass parameters to precursor model
		Param precursor_param_cr(precursor_model_cr_.getParameters());
		precursor_param_cr.setValue("pseudo_counts", pseudo_counts);
		precursor_model_cr_.setParameters(precursor_param_cr);

		Param precursor_param_cd = precursor_model_cd_.getParameters();
		precursor_param_cd.setValue("pseudo_counts", pseudo_counts);
		precursor_model_cd_.setParameters(precursor_param_cd);

		Param b_ion_losses_param_cr = b_ion_losses_cr_.getParameters();
		b_ion_losses_param_cr.setValue("pseudo_counts", pseudo_counts);
		b_ion_losses_cr_.setParameters(b_ion_losses_param_cr);

    Param b_ion_losses_param_cd = b_ion_losses_cd_.getParameters();
    b_ion_losses_param_cd.setValue("pseudo_counts", pseudo_counts);
    b_ion_losses_cd_.setParameters(b_ion_losses_param_cd);

		
    Param b2_ion_losses_param_cr = b2_ion_losses_cr_.getParameters();
    b2_ion_losses_param_cr.setValue("pseudo_counts", pseudo_counts);
    b2_ion_losses_cr_.setParameters(b2_ion_losses_param_cr);

    Param b2_ion_losses_param_cd = b2_ion_losses_cd_.getParameters();
    b2_ion_losses_param_cd.setValue("pseudo_counts", pseudo_counts);
    b2_ion_losses_cd_.setParameters(b2_ion_losses_param_cd);

		
		Param y_ion_losses_param_cr = y_ion_losses_cr_.getParameters();
    y_ion_losses_param_cr.setValue("pseudo_counts", pseudo_counts);
    y_ion_losses_cr_.setParameters(y_ion_losses_param_cr);

    Param y_ion_losses_param_cd = y_ion_losses_cd_.getParameters();
    y_ion_losses_param_cd.setValue("pseudo_counts", pseudo_counts);
    y_ion_losses_cd_.setParameters(y_ion_losses_param_cd);

		
		Param a_ion_losses_param_cr = a_ion_losses_cr_.getParameters();
		a_ion_losses_param_cr.setValue("pseudo_counts", pseudo_counts);
		a_ion_losses_cr_.setParameters(a_ion_losses_param_cr);

    Param a_ion_losses_param_cd = a_ion_losses_cd_.getParameters();
    a_ion_losses_param_cd.setValue("pseudo_counts", pseudo_counts);
    a_ion_losses_cd_.setParameters(a_ion_losses_param_cd);
	}
} // namespace OpenMS


