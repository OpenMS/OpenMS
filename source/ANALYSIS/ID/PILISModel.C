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

#define MIN_DECIMAL_VALUE 1e-10

using namespace std;

// new TODOS
// - New proton activation function, e.g. which handles internal sc occupancy 
//   and multiple e.g. P in a peptide


// old ones
// New QXP pathway handling
// New XXD yk-2 enhancement
//
// XPHXXX yk-2 enhancement


namespace OpenMS 
{

	PILISModel::PILISModel()
		: DefaultParamHandler("PILISModel"),
			valid_(false)
	{	
		defaults_.setValue("upper_mz", 2000.0, "Max m/z value of product ions in the simulated spectrum");
		defaults_.setValue("lower_mz", 200.0, "Lowest m/z value of product ions in the simulated spectrum");
		defaults_.setValue("charge_remote_threshold", 0.05, "If the probability for the proton at the N-terminus is lower than this value, enable charge-remote pathways");
		defaults_.setValue("charge_directed_threshold", 0.05, "Limit the proton availability at the N-terminus to at least this value for charge-directed pathways");
		defaults_.setValue("model_depth", 6, "The number of explicitly modeled backbone cleavages from N-terminus and C-terminus, would be 9 for the default value");
		defaults_.setValue("visible_model_depth", 50, "The maximal possible size of a peptide to be modeled");
		defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, used to identify the precursor peak and its loss peaks for training");
		defaults_.setValue("fragment_mass_tolerance", 0.3, "Peak mass tolerance of the product ions, used to identify the ions for training");
		defaults_.setValue("variable_modifications", StringList::create("MOD:00719,MOD:09997"), "Modifications which should be included in the model, represented by PSI-MOD accessions.");
		defaults_.setValue("fixed_modifications", StringList::create(""), "Modifications which should replace the unmodified amino acid, represented by PSI-MOD accessions.");
						
		defaults_.setValue("min_y_ion_intensity", 0.0, ".");
		defaults_.setValue("min_b_ion_intensity", 0.0, ".");
		defaults_.setValue("min_a_ion_intensity", 0.0, ".");
		defaults_.setValue("min_y_loss_intensity", 0.0, ".");
		defaults_.setValue("min_b_loss_intensity", 0.0, ".");

		defaults_.setValue("side_chain_activation", 0.1, "Additional activation of proton sitting at side chain, especially important at lysin and histidine residues");
		defaults_.setValue("pseudo_counts", 1e-15, "Value which is added for every transition trained of the underlying hidden Markov model");
		defaults_.setValue("max_isotope", 2, "Maximal isotope peak which is reported by the model, 0 means all isotope peaks are reported");

		defaults_.setValue("max_fragment_charge_training", 1, "Maximal allowed charge states for ions to be considered for training");
		defaults_.setValue("max_fragment_charge", 4, "Maximal charge state allowed for fragment ions");
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

		if (peptide.size() >= (UInt)param_.getValue("visible_model_depth"))
		{
			cerr << "PILISModel: cannot train peptide " << peptide << " of length " << peptide.size() << " (max for this model is " << (UInt)param_.getValue("visible_model_depth") << ", as defined by parameter \"visible_model_depth\")" << endl;
			return;
		}

		RichPeakSpectrum train_spec = in_spec;
		train_spec.sortByPosition();
		
		#ifdef TRAINING_DEBUG
		cout << "peptide: " << peptide  << "(z=" << charge << ")" << endl;
		#endif

		// get proton distribution
		Map<UInt, double> bb_charge_full, sc_charge_full;
		prot_dist_.getProtonDistribution(bb_charge_full, sc_charge_full, peptide, charge, Residue::YIon);
		prot_dist_.setPeptideProtonDistribution(bb_charge_full, sc_charge_full);
	
		// get start probabilities
		vector<double> bb_init, sc_init, cr_init;
		bool is_charge_remote = getInitialTransitionProbabilities_(bb_init, cr_init, sc_init, bb_charge_full, sc_charge_full, peptide);
		double bb_sum(0);
		for (Map<UInt, double>::ConstIterator iit = bb_charge_full.begin(); iit != bb_charge_full.end(); ++iit)
		{
			bb_sum += iit->second;
		}

		// activation of protons sitting at lysine and histidine side chains
		double side_chain_activation(param_.getValue("side_chain_activation"));
		for (Size i = 0; i != peptide.size(); ++i)
		{
			if (peptide[i].getOneLetterCode() != "R")
			{
				bb_sum += side_chain_activation * sc_charge_full[i];
			}
		}

		// backbone proton availability thresholds
		if  (bb_sum > 1)
		{
			bb_sum = 1;		
		}
		if (bb_sum < (double)param_.getValue("charge_directed_threshold") * charge)
		{
			bb_sum = (double)param_.getValue("charge_directed_threshold") * charge;
		}

		// clear the main Hidden Markov Model
		hmm_.clearInitialTransitionProbabilities();
		hmm_.clearTrainingEmissionProbabilities();
		
		//double harge_sum(0);
		vector<AASequence> prefixes, suffixes;

		// get the median backbone proton occupancy, for normalization purposes
		vector<double> bb_charges_all;
    for (Size i = 0; i != bb_charge_full.size(); ++i)
    {
      bb_charges_all.push_back(bb_charge_full[i]);
    }
    sort(bb_charges_all.begin(), bb_charges_all.end());
    double bb_charges_median = bb_charges_all[UInt(bb_charges_all.size()/2.0)];

		// for each site: 
		// 1. set proton distribution, 
		// 2. initial training intensities, 
		// 3. train the model
		for (Size i = 0; i != peptide.size() - 1; ++i)
		{
			String pos_name, prefix_size(i + 1), suffix_size(peptide.size() - 1 - i);

			if (i < floor((peptide.size() - 1.0)/2.0))
			{
				pos_name = prefix_size;
			}
			else
			{
				pos_name = "k-"+suffix_size;
			}
						
			AASequence prefix(peptide.getPrefix(i + 1)), suffix(peptide.getSuffix(peptide.size() - i - 1));
			AASequence aa1_seq, aa2_seq;
			aa1_seq += &peptide[i];
			aa2_seq += &peptide[i + 1];
			String aa1(aa1_seq.toString()), aa2(aa2_seq.toString());
			
			// calc PAs and get b/y ratios for bxyz pathway
			vector<double> b_cr_ints(charge, 0), y_cr_ints(charge, 0), b_sc_ints(charge, 0), y_sc_ints(charge, 0), b_ints(charge, 0), y_ints(charge, 0), a_ints(charge, 0), ay_ints(charge, 0);
			prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_ints, y_ints, ProtonDistributionModel::ChargeDirected);

			prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::AIon, a_ints, ay_ints, ProtonDistributionModel::ChargeDirected);

			double bb_enhance_factor(max(1.0, sqrt(bb_charge_full[i+1] / bb_charges_median)));
			//cerr << "bb_sum=" << bb_sum << " " << bb_enhance_factor << endl;

			if (!is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("BB"+pos_name, bb_sum * bb_enhance_factor);
				hmm_.setInitialTransitionProbability(aa1+aa2+"_bxyz"+pos_name, bb_sum  * bb_enhance_factor);
				hmm_.setInitialTransitionProbability(aa1+aa2+"_axyz"+pos_name, bb_sum  * bb_enhance_factor);
			}

			prefixes.push_back(prefix);
			suffixes.push_back(suffix);
			
			if ((aa1 == "D" || aa1 == "E" || pos_name == "k-1" || pos_name == "k-2") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("CR"+pos_name, 1 - bb_sum);
				//charge_sum += cr_init[i];

				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,	Residue::BIon, b_cr_ints, y_cr_ints, ProtonDistributionModel::ChargeRemote);

				//cerr << "ChargeStats: CR=" << CR_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;
			}

			if ((aa1 == "K" || aa1 == "H" || aa1 == "R") && is_charge_remote)
			{
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_sc_ints, y_sc_ints, ProtonDistributionModel::SideChain);
				hmm_.setInitialTransitionProbability("SC"+pos_name, sc_charge_full[i]);

				//cerr << "ChargeStats: SC=" << SC_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;	
				//charge_sum += SC_charges[i];
			}

			
			UInt max_fragment_charge_training = param_.getValue("max_fragment_charge_training");
			Map<int, double> sum_a, sum_b, sum_y;
			double sum_a_ints(0.0), sum_b_ints(0.0), sum_y_ints(0.0);
      for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
      {
				sum_a[z] = 0;
				sum_b[z] = 0;
				sum_y[z] = 0;
        double charge_remote_threshold((double)param_.getValue("charge_remote_threshold"));
        double avail_bb_sum_prefix = getAvailableBackboneCharge_(prefix, Residue::BIon, z);
        double avail_bb_sum_suffix = getAvailableBackboneCharge_(suffix, Residue::YIon, z);

        if (prefix.size() != 2)
        {
          double pre_weight = prefix.getMonoWeight(Residue::BIon);

          if (avail_bb_sum_prefix <= charge_remote_threshold)
          {
            //cerr << "Train b-ion losses, CR (avail=" << avail_bb_sum_prefix << ", z=" << z << ", prefix=" << prefix << endl;
            sum_b[z] += b_ion_losses_cr_.train(in_spec, prefix, pre_weight, z);
          }
          else
          {
            //cerr << "Train b-ion losses, CD (avail=" << avail_bb_sum_prefix << ", z=" << z << ", prefix=" << prefix << endl;
            sum_b[z] += b_ion_losses_cd_.train(in_spec, prefix, pre_weight, z);
          }
        }
        else
        {
          double pre_weight = prefix.getMonoWeight(Residue::BIon);

          if (avail_bb_sum_prefix <= charge_remote_threshold)
          {
            sum_b[z] += b2_ion_losses_cr_.train(in_spec, prefix, pre_weight, z);
          }
          else
          {
            sum_b[z] += b2_ion_losses_cd_.train(in_spec, prefix, pre_weight, z);
          }
        }

        double a_pre_weight = prefix.getMonoWeight(Residue::AIon);
        if (avail_bb_sum_prefix <= charge_remote_threshold)
        {
          sum_a[z] += a_ion_losses_cr_.train(in_spec, prefix, a_pre_weight, z);
        }
        else
        {
					sum_a[z] += a_ion_losses_cd_.train(in_spec, prefix, a_pre_weight, z);
        }

        double suf_weight = suffix.getMonoWeight(Residue::YIon);
        if (avail_bb_sum_suffix <= charge_remote_threshold)
        {
          //cerr << "Train y-ion losses, CR (avail=" << avail_bb_sum_suffix << ", z=" << z << ", suffix=" << suffix << endl;
          sum_y[z] += y_ion_losses_cr_.train(in_spec, suffix, suf_weight, z);
        }
        else
        {
          //cerr << "Train y-ion losses, CD (avail=" << avail_bb_sum_suffix << ", z=" << z << ", suffix=" << suffix << endl;
          sum_y[z] += y_ion_losses_cd_.train(in_spec, suffix, suf_weight, z);
        }

				sum_a_ints += sum_a[z];
				sum_b_ints += sum_b[z];
				sum_y_ints += sum_y[z];
      }

			double emission_prob(0);

			// end state is always needed
			hmm_.setTrainingEmissionProbability("end"+pos_name, 0.5/(double(peptide.size() - 1)));
			if (!is_charge_remote)
			{
				emission_prob = 0;
				hmm_.setTrainingEmissionProbability("AA"+pos_name, sum_a_ints + sum_b_ints + sum_y_ints);

				// now enable the states
				hmm_.enableTransition("BB"+pos_name, "AA"+pos_name);
				hmm_.enableTransition("BB"+pos_name, "end"+pos_name);

				// bxyz path
				hmm_.enableTransition(aa1+aa2+"_bxyz" + pos_name, "bxyz"+pos_name);
				hmm_.enableTransition(aa1+aa2+"_bxyz" + pos_name, "end"+pos_name);
			
				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += b_ints[z - 1] * sum_b[z];
				 	emission_prob += y_ints[z - 1] * sum_y[z];
				}
				hmm_.setTransitionProbability("bxyz" + pos_name, "bxyz_" + pos_name + "-ions", 1.0);
				hmm_.setTrainingEmissionProbability("bxyz_" + pos_name + "-ions", emission_prob);

				// axyz path
				hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "axyz"+pos_name);
				hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "end"+pos_name);
			
				emission_prob = 0;
				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += a_ints[z - 1] * sum_a[z];
					//emission_prob += ay_ints[z] * sum_y[z]; // TODO test if this can be done, or ruins the a-ion intensities
				}
				hmm_.setTransitionProbability("axyz" + pos_name, "axyz_" + pos_name + "-ions", 1.0);
				hmm_.setTrainingEmissionProbability("axyz_" + pos_name + "-ions", emission_prob);
			}
			
			if (aa1 == "D" && is_charge_remote)
			{
				emission_prob = 0;
				hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);

				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += b_cr_ints[z - 1] * sum_b[z];
					emission_prob += y_cr_ints[z - 1] * sum_y[z];
				}

				hmm_.setTrainingEmissionProbability("A"+pos_name, emission_prob);
				hmm_.setTransitionProbability("D" + pos_name, "D_" + pos_name + "-ions", 1.0);
				hmm_.setTrainingEmissionProbability("D_" + pos_name + "-ions", emission_prob);

				hmm_.setInitialTransitionProbability(aa2+"_D"+pos_name, 1 - bb_sum);
				hmm_.enableTransition(aa2+"_D"+pos_name, "D"+pos_name);
				hmm_.enableTransition(aa2+"_D"+pos_name, "end"+pos_name);
			}

			// TODO test !has R rule
			if (pos_name == "k-1" && is_charge_remote && aa1 != "D" && aa1 != "E" && !peptide.has("R"))
			{
				emission_prob = 0;
				hmm_.enableTransition("CRk-1", "Ak-1");
				hmm_.enableTransition("CRk-1", "endk-1");

				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += b_cr_ints[z - 1] * sum_b[z];
					emission_prob += y_cr_ints[z - 1] * sum_y[z];
				}

				hmm_.setTrainingEmissionProbability("Ak-1", emission_prob);
				hmm_.setTransitionProbability("bk-1", "bk-1_-ions", 1.0);
				hmm_.setTrainingEmissionProbability("bk-1_-ions", emission_prob);

				hmm_.setInitialTransitionProbability(aa1+"_bk-1", 1 - bb_sum);
				hmm_.enableTransition(aa1+"_bk-1", "bk-1");
				hmm_.enableTransition(aa1+"_bk-1", "endk-1");
			}

			if (pos_name == "k-2" && is_charge_remote && aa1 != "D" && aa1 != "E" && !peptide.has("R"))
			{
				emission_prob = 0;
				hmm_.enableTransition("CRk-2", "Ak-2");
				hmm_.enableTransition("CRk-2", "endk-2");
				
				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += b_cr_ints[z - 1] * sum_b[z];
					emission_prob += y_cr_ints[z - 1] * sum_y[z];
				}

				hmm_.setTrainingEmissionProbability("Ak-2", emission_prob);
				hmm_.setTransitionProbability("bk-2", "bk-2_-ions", 1.0);
				hmm_.setTrainingEmissionProbability("bk-2_-ions", emission_prob);

				hmm_.setInitialTransitionProbability(aa1+"_bk-2", 1 - bb_sum);
				hmm_.enableTransition(aa1+"_bk-2", "bk-2");
        hmm_.enableTransition(aa1+"_bk-2", "endk-2");
			}
			
			if (aa1 == "E" && is_charge_remote)
      { 
				emission_prob = 0;
        hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);

				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += b_cr_ints[z - 1] * sum_b[z];
					emission_prob += y_cr_ints[z - 1] * sum_y[z];
				}

				hmm_.setTrainingEmissionProbability("A"+pos_name, emission_prob);
				hmm_.setTransitionProbability("E" + pos_name, "E_" + pos_name + "-ions", 1.0);
				hmm_.setTrainingEmissionProbability("E_" + pos_name + "-ions", emission_prob);

				hmm_.setInitialTransitionProbability(aa2+"_E"+pos_name, 1 - bb_sum);
        hmm_.enableTransition(aa2+"_E"+pos_name, "E"+pos_name);
				hmm_.enableTransition(aa2+"_E"+pos_name, "end"+pos_name);
      }

			if (aa1 == "K" && is_charge_remote)
			{
				emission_prob = 0;
				hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);

				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += b_sc_ints[z - 1] * sum_b[z];
					emission_prob += y_sc_ints[z - 1] * sum_y[z];
				}

				hmm_.setTrainingEmissionProbability("ASC"+pos_name, emission_prob);
				hmm_.setTransitionProbability("K" + pos_name, "K_" + pos_name + "-ions", 1.0);
				hmm_.setTrainingEmissionProbability("K_" + pos_name + "-ions", emission_prob);

				hmm_.setInitialTransitionProbability(aa2+"_K"+pos_name, sc_charge_full[i]);
        hmm_.enableTransition(aa2+"_K"+pos_name, "K"+pos_name);
				hmm_.enableTransition(aa2+"_K"+pos_name, "end"+pos_name);
			}

			if (aa1 == "H" && is_charge_remote)
      {
				emission_prob = 0;
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);

				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += b_sc_ints[z - 1] * sum_b[z];
					emission_prob += y_sc_ints[z - 1] * sum_y[z];
				}
				
				hmm_.setTrainingEmissionProbability("ASC"+pos_name, emission_prob);
				hmm_.setTransitionProbability("H" + pos_name, "H_" + pos_name + "-ions", 1.0);
				hmm_.setTrainingEmissionProbability("H_" + pos_name + "-ions", emission_prob);

				hmm_.setInitialTransitionProbability(aa2+"_H"+pos_name, sc_charge_full[i]);
        hmm_.enableTransition(aa2+"_H"+pos_name, "H"+pos_name);
				hmm_.enableTransition(aa2+"_H"+pos_name, "end"+pos_name);
      }

      if (aa1 == "R" && is_charge_remote)
      { 
				emission_prob = 0;
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);

				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					emission_prob += b_sc_ints[z - 1] * sum_b[z];
					emission_prob += y_sc_ints[z - 1] * sum_y[z];
				}
				
				hmm_.setTrainingEmissionProbability("ASC"+pos_name, emission_prob);
				hmm_.setTransitionProbability("R" + pos_name, "R_" + pos_name + "-ions", 1.0);
        hmm_.setTrainingEmissionProbability("R_" + pos_name + "-ions", emission_prob);

				hmm_.setInitialTransitionProbability(aa2+"_R"+pos_name, sc_charge_full[i]);
        hmm_.enableTransition(aa2+"_R"+pos_name, "R"+pos_name);
				hmm_.enableTransition(aa2+"_R"+pos_name, "end"+pos_name);
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
		
		hmm_.setVariableModifications((StringList)param_.getValue("variable_modifications"));
		hmm_.estimateUntrainedTransitions();
	}

	void PILISModel::getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, UInt charge)
	{
		if (!valid_)
		{
			cerr << "PILISModel: cannot simulate, initialize model from file first, e.g. data/PILIS/PILIS_model_default.dat" << endl;
			return;
		}

		if (peptide.size() > (UInt)param_.getValue("visible_model_depth"))
		{
			cerr << "PILISModel: cannot generate spectra of peptide '" << peptide << "' of length " << peptide.size() << " (max of this model is " << (UInt)param_.getValue("visible_model_depth") << ", as defined by \"visible_model_depth\")" << endl;
			return;
		}
		
		// calc proton distribution
		Map<UInt, double> sc_charge_full, bb_charge_full;
		//cerr << "getProtonDistribution: " << peptide << " " << charge << endl;
		prot_dist_.getProtonDistribution(bb_charge_full, sc_charge_full, peptide, charge, Residue::YIon);
		prot_dist_.setPeptideProtonDistribution(bb_charge_full, sc_charge_full);

		//for (Size i = 0; i != bb_charge_full.size() - 1; ++i)
		//{
		//	cerr << "i: bb=" << bb_charge_full[i] << ", sc=" << sc_charge_full[i] << endl;
		//}

		hmm_.clearInitialTransitionProbabilities();
		hmm_.clearTrainingEmissionProbabilities();
	
    // set charges
    vector<double> bb_init, sc_init, cr_init;
		bool is_charge_remote = getInitialTransitionProbabilities_(bb_init, cr_init, sc_init, bb_charge_full, sc_charge_full, peptide);
		//cerr << "is_charge_remote=" << is_charge_remote << endl;
		//for (Size i = 0; i != bb_init.size(); ++i)
		//{
		//	cerr << "bb=" << bb_init[i] << ", cr=" << cr_init[i] << ", sc=" << sc_init[i] << endl;
		//}
	
		double bb_sum(0), bb_sum_orig(0);
		for (Map<UInt, double>::ConstIterator iit = bb_charge_full.begin(); iit != bb_charge_full.end(); ++iit)
		{
			bb_sum += iit->second;
		}
	
    // activation of protons sitting at lysine and histidine side chains
		double side_chain_activation(param_.getValue("side_chain_activation"));
    for (Size i = 0; i != peptide.size(); ++i)
    {
      if (peptide[i].getOneLetterCode() == "H" || peptide[i].getOneLetterCode() == "K")
      {
        bb_sum += side_chain_activation * sc_charge_full[i];
      }
    }
		
		if (bb_sum > 1)
		{
			bb_sum = 1;
		}
		bb_sum_orig = bb_sum;
		if (bb_sum < (double)param_.getValue("charge_directed_threshold") * charge)
		{
			bb_sum = (double)param_.getValue("charge_directed_threshold") * charge;
		}

		vector<AASequence> suffixes, prefixes;
		vector<String> pos_names;
	
		vector<double> bb_charges_all;
		for (Size i = 0; i != bb_charge_full.size(); ++i)
		{
			bb_charges_all.push_back(bb_charge_full[i]);
		}
		sort(bb_charges_all.begin(), bb_charges_all.end());
		double bb_charges_median = bb_charges_all[UInt(bb_charges_all.size()/2.0)];
	
		vector<vector<double> > all_b_cr_ints, all_y_cr_ints, all_b_sc_ints, all_y_sc_ints, all_b_ints, all_y_ints, all_a_ints, all_ay_ints;
		// get the paths
		for (Size i = 0; i != peptide.size() - 1; ++i)
		{
      String pos_name, y_name1, b_name1, a_name1, y_name, b_name;
      String y_name2, b_name2, a_name2;
			String i_plus1(i+1), pep_size_i(peptide.size() - 1 - i);

      if (i < floor((peptide.size() - 1.0)/2.0))
      {
        pos_name = i_plus1;
      }
      else
      { 
        pos_name = "k-"+pep_size_i;
      }

			AASequence prefix(peptide.getPrefix(i + 1)), suffix(peptide.getSuffix(peptide.size() - 1 - i));
			suffixes.push_back(suffix);
			prefixes.push_back(prefix);
			pos_names.push_back(pos_name);

			AASequence aa1_seq, aa2_seq;
			aa1_seq += &peptide[i];
			aa2_seq += &peptide[i + 1];
			String aa1(aa1_seq.toString()), aa2(aa2_seq.toString());

			double bb_enhance_factor(max(1.0, sqrt(bb_charge_full[i+1] / bb_charges_median)));
			hmm_.setInitialTransitionProbability("BB"+pos_name, bb_sum * bb_enhance_factor);

			vector<double> b_cr_ints(charge, 0), y_cr_ints(charge, 0), b_sc_ints(charge, 0), y_sc_ints(charge, 0), b_ints(charge, 0), y_ints(charge, 0), a_ints(charge, 0), ay_ints(charge, 0);
			if ((aa1 == "D" || aa1 == "E" || pos_name == "k-1" || pos_name == "k-2") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("CR"+pos_name, 1 - bb_sum);
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_cr_ints, y_cr_ints, ProtonDistributionModel::ChargeRemote);
			}

			if ((aa1 == "K" || aa1 == "H" || aa1 == "R") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("SC"+pos_name, sc_charge_full[i]);
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,	Residue::BIon, b_sc_ints, y_sc_ints, ProtonDistributionModel::SideChain);
			}
		
      prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_ints, y_ints, ProtonDistributionModel::ChargeDirected);
      prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::AIon, a_ints, ay_ints, ProtonDistributionModel::ChargeDirected);	

			//if (charge == 3)
			//{
			//	cerr << prefix << " - " << suffix << " ";
			//	for (Size ions_charge = 0; ions_charge != b_ints.size(); ++ions_charge)
			//	{
			//		cerr << b_ints[ions_charge] << " - " << y_ints[ions_charge] << " | ";
			//	}
			//	cerr << endl;
			//}

			all_a_ints.push_back(a_ints);
			all_ay_ints.push_back(ay_ints);
			all_b_ints.push_back(b_ints);
			all_b_cr_ints.push_back(b_cr_ints);
			all_y_cr_ints.push_back(y_cr_ints);
			all_b_sc_ints.push_back(b_sc_ints);
			all_y_sc_ints.push_back(y_sc_ints);
			all_y_ints.push_back(y_ints);

      // now enable the states
			hmm_.enableTransition("BB"+pos_name, "AA"+pos_name);
			hmm_.enableTransition("BB"+pos_name, "end"+pos_name);
			
      hmm_.enableTransition("AA"+pos_name, aa1+aa2+"_bxyz"+pos_name);
      hmm_.enableTransition(aa1+aa2+"_bxyz"+pos_name, "bxyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"_bxyz"+pos_name, "end"+pos_name);

			hmm_.setTransitionProbability("bxyz" + pos_name, "bxyz_" + pos_name + "-ions", 1.0);

			// axyz
			hmm_.enableTransition("AA"+pos_name, aa1+aa2+"_axyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "axyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "end"+pos_name);
			hmm_.setTransitionProbability("axyz" + pos_name, "axyz_" + pos_name + "-ions", 1.0);

			//cerr << "neutral loss charge sums: " << H2O_b_charge_sum << " " << H2O_y_charge_sum << " " << NH3_b_charge_sum << " " << NH3_y_charge_sum << " " << charge_dir_tmp << endl;

      if (aa1 == "D" && is_charge_remote)
      {
				hmm_.setTransitionProbability("D" + pos_name, "D_" + pos_name + "-ions", 1.0);
				hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);
				hmm_.enableTransition("A"+pos_name, aa2+"_D"+pos_name);
        hmm_.enableTransition(aa2+"_D"+pos_name, "D"+pos_name);
				hmm_.enableTransition(aa2+"_D"+pos_name, "end"+pos_name);
			}

      if (pos_name == "k-1" && is_charge_remote && aa1 != "D" && aa1 != "E" && !peptide.has("R"))
      {
				hmm_.setTransitionProbability("bk-1", "bk-1_-ions", 1.0);

				hmm_.enableTransition("CRk-1", "Ak-1");
        hmm_.enableTransition("CRk-1", "endk-1");
				hmm_.enableTransition("Ak-1", aa1+"_bk-1");

				hmm_.enableTransition(aa1 + "_bk-1", "bk-1");
        hmm_.enableTransition(aa1 + "_bk-1", "end"+pos_name);
      }

			if (pos_name == "k-2" && is_charge_remote && aa1 != "D" && aa1 != "E" && !peptide.has("R"))
			{
				hmm_.setTransitionProbability("bk-2", "bk-2_-ions", 1.0);

				hmm_.enableTransition("CRk-2", "Ak-2");
				hmm_.enableTransition("CRk-2", "endk-2");
				hmm_.enableTransition("Ak-2", aa1+"_bk-2");

				hmm_.enableTransition(aa1 + "_bk-2", "bk-2");
				hmm_.enableTransition(aa1 + "_bk-2", "end"+pos_name);
			}

			if (aa1 == "E" && is_charge_remote)
      {
				hmm_.setTransitionProbability("E" + pos_name, "E_" + pos_name + "-ions", 1.0);
        hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);
        hmm_.enableTransition("A"+pos_name, aa2+"_E"+pos_name);
        hmm_.enableTransition(aa2+"_E"+pos_name, "E"+pos_name);
				hmm_.enableTransition(aa2+"_E"+pos_name, "end"+pos_name);
      }

			if (aa1 == "K" && is_charge_remote)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
				hmm_.setTransitionProbability("K" + pos_name, "K_" + pos_name + "-ions", 1.0);
        hmm_.enableTransition("ASC"+pos_name, aa2+"_K"+pos_name);
        hmm_.enableTransition(aa2+"_K"+pos_name, "K"+pos_name);
				hmm_.enableTransition(aa2+"_K"+pos_name, "end"+pos_name);
			}
				
			if (aa1 == "H" && is_charge_remote)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
				hmm_.setTransitionProbability("H" + pos_name, "H_" + pos_name + "-ions", 1.0);
        hmm_.enableTransition("ASC"+pos_name, aa2+"_H"+pos_name);
        hmm_.enableTransition(aa2+"_H"+pos_name, "H"+pos_name);
				hmm_.enableTransition(aa2+"_H"+pos_name, "end"+pos_name);
      }

      if (aa1 == "R" && is_charge_remote)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
				hmm_.setTransitionProbability("R" + pos_name, "R_" + pos_name + "-ions", 1.0);
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


		UInt max_isotope = (UInt)param_.getValue("max_isotope");
		UInt max_fragment_charge = (UInt)param_.getValue("max_fragment_charge");
		double charge_remote_threshold = (double)param_.getValue("charge_remote_threshold");
		bool suppress_single_ions = param_.getValue("suppress_single_ions").toBool();
		for (Size i = 0; i != prefixes.size(); ++i)
		{
			vector<double> path_intensities;
			path_intensities.push_back(tmp[hmm_.getState("bxyz_" + pos_names[i] + "-ions")]);
			//cerr << prefixes[i] << " - " << suffixes[i] <<  endl;
			//cerr << "bxyz-intensity: " << path_intensities.back() << endl;
      path_intensities.push_back(tmp[hmm_.getState("H_" + pos_names[i] + "-ions")]);

      path_intensities.push_back(tmp[hmm_.getState("K_" + pos_names[i] + "-ions")]);
      path_intensities.push_back(tmp[hmm_.getState("R_" + pos_names[i] + "-ions")]);
      path_intensities.push_back(tmp[hmm_.getState("D_" + pos_names[i] + "-ions")]);
			//cerr << "Dcr-intensity: " << path_intensities.back() << endl;
      path_intensities.push_back(tmp[hmm_.getState("E_" + pos_names[i] + "-ions")]);
			//cerr << "Ecr-intensity: " << path_intensities.back() << endl;
			
			/*
      cerr << prefixes[i] << " - " << suffixes[i] << ": bxyz=" << path_intensities[0] <<
                                                     ", H=" << path_intensities[1] <<
                                                     ", K=" << path_intensities[2] <<
                                                     ", R=" << path_intensities[3] <<
                                                     ", D=" << path_intensities[4] <<
                                                     ", E=" << path_intensities[5];
*/

			
			if (i == peptide.size() - 2)
			{
				path_intensities.push_back(tmp[hmm_.getState("bk-1_-ions")]);
				//cerr << ", bk-1=" << path_intensities.back();
			}
			if (i == peptide.size() - 3)
			{
				path_intensities.push_back(tmp[hmm_.getState("bk-2_-ions")]);
				//cerr << ", bk-2=" << path_intensities.back();
			}
			path_intensities.push_back(tmp[hmm_.getState("axyz_" + pos_names[i] + "-ions")]);
			//cerr << ", axzy=" << path_intensities.back() << endl;
	
			vector<vector<vector<double> > > prefix_intensities, suffix_intensities;
			prefix_intensities.push_back(all_b_ints);
			prefix_intensities.push_back(all_b_sc_ints);
			prefix_intensities.push_back(all_b_sc_ints);
			prefix_intensities.push_back(all_b_sc_ints);
			prefix_intensities.push_back(all_b_cr_ints);
			prefix_intensities.push_back(all_b_cr_ints);
			if (i == peptide.size() - 2 || i == peptide.size() - 3)
			{
				prefix_intensities.push_back(all_b_cr_ints);
			}
			prefix_intensities.push_back(all_a_ints);

      suffix_intensities.push_back(all_y_ints);
      suffix_intensities.push_back(all_y_sc_ints);
      suffix_intensities.push_back(all_y_sc_ints);
      suffix_intensities.push_back(all_y_sc_ints);
      suffix_intensities.push_back(all_y_cr_ints);
      suffix_intensities.push_back(all_y_cr_ints);
			if (i == peptide.size() - 2 || i == peptide.size() - 3)
			{
				suffix_intensities.push_back(all_y_cr_ints);
			}
			suffix_intensities.push_back(all_ay_ints);
		

			for (Size int_it = 0; int_it != path_intensities.size(); ++int_it)
			{
				if (path_intensities[int_it] < MIN_DECIMAL_VALUE)
				{
					continue;
				}
				double prefix_weight = prefixes[i].getMonoWeight(Residue::BIon);
				double suffix_weight = suffixes[i].getMonoWeight(Residue::YIon);
				IsotopeDistribution prefix_id;
				
				if (int_it != path_intensities.size() - 1)
				{
					prefix_id = prefixes[i].getFormula(Residue::BIon).getIsotopeDistribution(max_isotope);
				}
				else
				{
					prefix_id = prefixes[i].getFormula(Residue::AIon).getIsotopeDistribution(max_isotope);
				}
				IsotopeDistribution suffix_id = suffixes[i].getFormula(Residue::YIon).getIsotopeDistribution(max_isotope);

				for (UInt z = 1; z <= charge && z <= max_fragment_charge; ++z)
				{
					if (prefix_intensities[int_it][i][z - 1] != 0)
					{
						vector<RichPeak1D> b_loss_peaks;
						double avail_bb_sum_prefix(0);
						
						if (int_it != path_intensities.size() - 1)
						{
							avail_bb_sum_prefix = getAvailableBackboneCharge_(prefixes[i], Residue::BIon, z);
							if (avail_bb_sum_prefix <= charge_remote_threshold)
            	{
								if (i == 1) // second BB position
								{
									b2_ion_losses_cr_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * path_intensities[int_it]);
								}
								else
								{
              		b_ion_losses_cr_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * path_intensities[int_it]);
								}
            	}
            	else
            	{
								if (i == 1)
								{
									b2_ion_losses_cd_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * path_intensities[int_it]);
								}
              	else
								{
									b_ion_losses_cd_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * path_intensities[int_it]);
								}
            	}
						}
						else
						{
							avail_bb_sum_prefix = getAvailableBackboneCharge_(prefixes[i], Residue::AIon, z);
							if (avail_bb_sum_prefix <= charge_remote_threshold)
              {
                a_ion_losses_cr_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * path_intensities[int_it]);
              }
              else
              {
                a_ion_losses_cd_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * path_intensities[int_it]);
              }
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
          			for (Size j = 1; j != split.size(); ++j)
          			{
            			b_ion_name += split[j];
          			}
        			}

        			//if (!suppress_single_ions || (suppress_single_ions && i + 1 > 1))
        			//{
								b_ion_name += String(z, '+');

								if (int_it != path_intensities.size() - 1)
								{
          				addPeaks_(prefix_weight, z, it->getMZ(), it->getIntensity(), spec, prefix_id, b_ion_name);
								}
								else
								{
									addPeaks_(prefix_weight - 28.0, z, it->getMZ(), it->getIntensity(), spec, prefix_id, b_ion_name);
								}
        			//}
						}
					}

					if (suffix_intensities[int_it][i][z - 1] > MIN_DECIMAL_VALUE && int_it != path_intensities.size() - 1) // no ay 
					{
						vector<RichPeak1D> y_loss_peaks;
						double avail_bb_sum_suffix = getAvailableBackboneCharge_(suffixes[i], Residue::YIon, z);
						if (avail_bb_sum_suffix <= charge_remote_threshold)
						{
							y_ion_losses_cr_.getIons(y_loss_peaks, suffixes[i], suffix_intensities[int_it][i][z - 1] * path_intensities[int_it]);
						}
						else
						{
							y_ion_losses_cd_.getIons(y_loss_peaks, suffixes[i], suffix_intensities[int_it][i][z - 1] * path_intensities[int_it]);
						}

						for (vector<RichPeak1D>::const_iterator it = y_loss_peaks.begin(); it != y_loss_peaks.end(); ++it)
        		{
          		String y_ion_name = it->getMetaValue("IonName");
          		vector<String> split;
          		y_ion_name.split('-', split);
          		if (split.size() == 0)
          		{
            		y_ion_name += String(peptide.size() - i - 1);
          		}
          		else
          		{
            		y_ion_name = split[0] + String(peptide.size() - i - 1) + "-";
            		for (Size j = 1; j != split.size(); ++j)
            		{
              		y_ion_name += split[j];
            		}
          		}

          		//if (!suppress_single_ions || (suppress_single_ions && i + 1 > 1))
          		//{
								y_ion_name += String(z, '+');
            		addPeaks_(suffix_weight, z, it->getMZ(), it->getIntensity(), spec, suffix_id, y_ion_name);
          		//}
						}
        	}
				}
			}
		}

		hmm_.disableTransitions();

		// precursor intensities
		vector<RichPeak1D> pre_peaks;
		//cerr << peptide << " ->getSpectrum: bb_sum=" << bb_sum << ", charge_remote=" << is_charge_remote<< endl;
		if (bb_sum <= charge_remote_threshold)
		{
			precursor_model_cr_.getIons(pre_peaks, peptide, bb_sum);
		}
		else
		{
			precursor_model_cd_.getIons(pre_peaks, peptide, max(0.0, 1.0 - bb_sum));
		}
		
		double weight = peptide.getMonoWeight();
		IsotopeDistribution id;
		id.estimateFromPeptideWeight(weight);	
		for (vector<RichPeak1D>::const_iterator it = pre_peaks.begin(); it != pre_peaks.end(); ++it)
		{
			addPeaks_(weight, charge, it->getMZ(), it->getIntensity(), spec, id, it->getMetaValue("IonName"));
		}
	
	
		// now build the spectrum with the peaks
		double intensity_max(0);
		for (Map<double, vector<RichPeak1D> >::ConstIterator it = peaks_.begin(); it != peaks_.end(); ++it)
		{
			if (it->second.size() == 1 && it->second.begin()->getIntensity() != 0)
			{
				spec.push_back(*it->second.begin());
				if (intensity_max < spec.back().getIntensity())
				{
					intensity_max = spec.back().getIntensity();
				}
			}
			else
			{
				RichPeak1D p;
				double int_sum(0);
				p = *it->second.begin();
				set<String> names;
				for (vector<RichPeak1D>::const_iterator pit = it->second.begin(); pit != it->second.end(); ++pit)
				{
					int_sum += pit->getIntensity();
					if (String(pit->getMetaValue("IonName")) != "")
					{
						names.insert(pit->getMetaValue("IonName"));
					}
				}

				String name;
				for (set<String>::const_iterator nit = names.begin(); nit != names.end(); ++nit)
				{
					name += *nit + "/";
				}
				p.setMetaValue("IonName", name);
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
    for (Size i = 0; i != ion.size(); ++i)
    {
      if (ion[i].getOneLetterCode() != "R")
      {
        bb_sum += side_chain_activation * sc_charges[i];
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

		//double charge_dir_tmp(bb_charges[0]); TODO obsolete????

    double bb_sum(0);
    for (Size i = 0; i != bb_charges.size(); ++i)
    {
			//cerr << i << " " << bb_charges[i] << endl;
      bb_sum += bb_charges[i];
    }
		//cerr << "bb_sum=" << bb_sum << endl;

    //charge_dir_tmp = bb_sum;
    if (bb_sum < (double)param_.getValue("charge_directed_threshold"))
    {
      //charge_dir_tmp = (double)param_.getValue("charge_directed_threshold");
			is_charge_remote = true;

			bb_sum = (double)param_.getValue("charge_directed_threshold");
    }

    double bb_avg = (bb_sum - bb_charges[0]) / (double)(bb_charges.size() - 1);

					
		for (Size i = 0; i != peptide.size() - 1; ++i)
    {
      bb_init.push_back(bb_charges[i+1] * bb_sum /* charge_dir_tmp*/);
      String aa(peptide[i].getOneLetterCode());
      if (sc_charges.has(i))
      {
        if ((aa == "K" || aa == "R" || aa == "H"))
        {
          sc_init.push_back(sc_charges[i] * bb_avg);
        }
        else
        {
          sc_init.push_back(0.0);
        }
        if (bb_sum <= (double)param_.getValue("charge_remote_threshold") && (aa == "D" || aa == "E"))
        {
          cr_init.push_back(((1 - bb_sum) * bb_avg));
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
		for (Size i = 0; i != bb_init.size(); ++i)
		{
			init_prob_sum += bb_init[i] + sc_init[i] + cr_init[i];
		}

		/*
		// C-term bb positional acmino acid might have a proton associated with (is inactive, as no pathway uses it)
		if (sc_charges.has(peptide.size() - 1))
		{
			init_prob_sum += sc_charges[peptide.size() - 1];
		}*/
		for (Size i = 0; i != bb_init.size(); ++i)
		{
			bb_init[i] /= init_prob_sum;
			sc_init[i] /= init_prob_sum;
			cr_init[i] /= init_prob_sum;
		}
		
		return is_charge_remote;
	}

	void PILISModel::addPeaks_(double mz, int charge, double offset, double intensity, RichPeakSpectrum& /*spectrum*/, const IsotopeDistribution& id, const String& name)
	{
		if (intensity < MIN_DECIMAL_VALUE)
		{
			return;
		}
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


