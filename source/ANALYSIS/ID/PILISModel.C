// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISModelGenerator.h>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
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

#define INIT_CHARGE_DEBUG
#undef  INIT_CHARGE_DEBUG

#define MIN_DECIMAL_VALUE 1e-8

using namespace std;

// new TODOS
// - New proton activation function, e.g. which handles internal sc occupancy 
//   and multiple e.g. P in a peptide
// - Pathway length scaling through root function, e.g. for loss models ????


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
		defaults_.setValue("charge_remote_threshold", 0.3, "If the probability for the proton at the N-terminus is lower than this value, enable charge-remote pathways");
		defaults_.setValue("charge_directed_threshold", 0.3, "Limit the proton availability at the N-terminus to at least this value for charge-directed pathways");
		defaults_.setValue("model_depth", 10, "The number of explicitly modeled backbone cleavages from N-terminus and C-terminus, would be 9 for the default value", StringList::create("advanced"));
		defaults_.setValue("visible_model_depth", 50, "The maximal possible size of a peptide to be modeled", StringList::create("advanced"));

		// tolerances
		defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, used to identify the precursor peak and its loss peaks for training");
		defaults_.setValue("fragment_mass_tolerance", 0.3, "Peak mass tolerance of the product ions, used to identify the ions for training");
		
		// modification parameters
		vector<String> all_mods;
		ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
		defaults_.setValue("variable_modifications", StringList::create("Oxidation (M),Carbamyl (K)"), "Modifications which should be included in the model, represented by PSI-MOD accessions.");
		defaults_.setValidStrings("variable_modifications", all_mods);
		defaults_.setValue("fixed_modifications", StringList::create(""), "Modifications which should replace the unmodified amino acid, represented by PSI-MOD accessions.");
		defaults_.setValidStrings("fixed_modifications", all_mods);


		defaults_.setValue("min_enhancement_factor", 2.0, "Minimal factor for bxyz backbone cleavages.", StringList::create("advanced"));
						
		defaults_.setValue("min_y_ion_intensity", 0.0, "Minimal intensity for y ions.", StringList::create("advanced"));
		defaults_.setValue("min_b_ion_intensity", 0.0, "Minimal intensity for b ions.", StringList::create("advanced"));
		defaults_.setValue("min_a_ion_intensity", 0.0, "Minimal intensity for a ions.", StringList::create("advanced"));
		defaults_.setValue("min_y_loss_intensity", 0.0, "Minimal intensity for y ions with neutral loss.", StringList::create("advanced"));
		defaults_.setValue("min_b_loss_intensity", 0.0, "Minimal intensity for b ions with neutral loss.", StringList::create("advanced"));

		defaults_.setValue("side_chain_activation", 0.2, "Additional activation of proton sitting at side chain, especially important at lysine and histidine residues");
		defaults_.setValue("pseudo_counts", 1e-15, "Value which is added for every transition trained of the underlying hidden Markov model", StringList::create("advanced"));
		defaults_.setValue("max_isotope", 2, "Maximal isotope peak which is reported by the model, 0 means all isotope peaks are reported");

		defaults_.setValue("max_fragment_charge_training", 2, "Maximal allowed charge states for ions to be considered for training");
		defaults_.setValue("max_fragment_charge", 4, "Maximal charge state allowed for fragment ions");

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

	void PILISModel::init(bool generate_models)
	{
		if (generate_models)
		{
			PILISModelGenerator gen;
			Param gen_param(gen.getParameters());
			gen_param.setValue("variable_modifications", (StringList)param_.getValue("variable_modifications"));
			gen_param.setValue("fixed_modifications", (StringList)param_.getValue("fixed_modifications"));
			gen_param.setValue("model_depth", param_.getValue("model_depth"));
			gen_param.setValue("visible_model_depth", param_.getValue("visible_model_depth"));
			gen.setParameters(gen_param);
			gen.getModel(hmm_);
		}


		Param pre_param(precursor_model_cr_.getParameters());
		pre_param.setValue("fragment_mass_tolerance", (DoubleReal)param_.getValue("fragment_mass_tolerance"));
		pre_param.setValue("variable_modifications", (StringList)param_.getValue("variable_modifications"));
		pre_param.setValue("fixed_modifications", (StringList)param_.getValue("fixed_modifications"));
		pre_param.setValue("pseudo_counts", (DoubleReal)param_.getValue("pseudo_counts"));
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
		parseHMMModel_(++it_begin, it_end, hmm_, param_);

		// seek to next interval
		it_begin = file.search(it_end, "PRECURSOR_MODEL_CR_BEGIN");
		it_end = file.search(it_begin, "PRECURSOR_MODEL_CR_END");
		Param p;
		HiddenMarkovModel precursor_model_cr_hmm;
		parseHMMModel_(++it_begin, it_end, precursor_model_cr_hmm, p);
		precursor_model_cr_.setHMM(precursor_model_cr_hmm);
		precursor_model_cr_.setParameters(p);

    it_begin = file.search(it_end, "PRECURSOR_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "PRECURSOR_MODEL_CD_END");
		p.clear();
		HiddenMarkovModel precursor_model_cd_hmm;
    parseHMMModel_(++it_begin, it_end, precursor_model_cd_hmm, p);
		precursor_model_cd_.setHMM(precursor_model_cd_hmm);
		precursor_model_cd_.setParameters(p);

		
		// b loss models
    it_begin = file.search(it_end, "BION_LOSS_MODEL_CR_BEGIN");
    it_end = file.search(it_begin, "BION_LOSS_MODEL_CR_END");
		p.clear();
    HiddenMarkovModel b_ion_loss_hmm_cr;
    parseHMMModel_(++it_begin, it_end, b_ion_loss_hmm_cr, p);
    b_ion_losses_cr_.setHMM(b_ion_loss_hmm_cr);
    b_ion_losses_cr_.setParameters(p);

    it_begin = file.search(it_end, "BION_LOSS_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "BION_LOSS_MODEL_CD_END");
		p.clear();
    HiddenMarkovModel b_ion_loss_hmm_cd;
    parseHMMModel_(++it_begin, it_end, b_ion_loss_hmm_cd, p);
    b_ion_losses_cd_.setHMM(b_ion_loss_hmm_cd);
    b_ion_losses_cd_.setParameters(p);

		it_begin = file.search(it_end, "B2ION_LOSS_MODEL_CR_BEGIN");
		it_end = file.search(it_begin, "B2ION_LOSS_MODEL_CR_END");
		p.clear();
		HiddenMarkovModel b2_ion_loss_hmm_cr;
		parseHMMModel_(++it_begin, it_end, b2_ion_loss_hmm_cr, p);
		b2_ion_losses_cr_.setHMM(b2_ion_loss_hmm_cr);
		b2_ion_losses_cr_.setParameters(p);

    it_begin = file.search(it_end, "B2ION_LOSS_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "B2ION_LOSS_MODEL_CD_END");
		p.clear();
    HiddenMarkovModel b2_ion_loss_hmm_cd;
    parseHMMModel_(++it_begin, it_end, b2_ion_loss_hmm_cd, p);
    b2_ion_losses_cd_.setHMM(b2_ion_loss_hmm_cd);
    b2_ion_losses_cd_.setParameters(p);
	
		// a-ion loss model
		it_begin = file.search(it_end, "AION_LOSS_MODEL_CR_BEGIN");
		it_end = file.search(it_begin, "AION_LOSS_MODEL_CR_END");
		p.clear();
		HiddenMarkovModel a_ion_loss_hmm_cr;
		parseHMMModel_(++it_begin, it_end, a_ion_loss_hmm_cr, p);
		a_ion_losses_cr_.setHMM(a_ion_loss_hmm_cr);
		a_ion_losses_cr_.setParameters(p);

    it_begin = file.search(it_end, "AION_LOSS_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "AION_LOSS_MODEL_CD_END");
		p.clear();
    HiddenMarkovModel a_ion_loss_hmm_cd;
    parseHMMModel_(++it_begin, it_end, a_ion_loss_hmm_cd, p);
    a_ion_losses_cd_.setHMM(a_ion_loss_hmm_cd);
    a_ion_losses_cd_.setParameters(p);
		
		// y-ion loss model
		it_begin = file.search(it_end, "YION_LOSS_MODEL_CR_BEGIN");
		it_end = file.search(it_begin, "YION_LOSS_MODEL_CR_END");
		p.clear();
		HiddenMarkovModel y_ion_loss_hmm_cr;
		parseHMMModel_(++it_begin, it_end, y_ion_loss_hmm_cr, p);
		y_ion_losses_cr_.setHMM(y_ion_loss_hmm_cr);
		y_ion_losses_cr_.setParameters(p);

    it_begin = file.search(it_end, "YION_LOSS_MODEL_CD_BEGIN");
    it_end = file.search(it_begin, "YION_LOSS_MODEL_CD_END");
		p.clear();
    HiddenMarkovModel y_ion_loss_hmm_cd;
    parseHMMModel_(++it_begin, it_end, y_ion_loss_hmm_cd, p);
    y_ion_losses_cd_.setHMM(y_ion_loss_hmm_cd);
    y_ion_losses_cd_.setParameters(p);

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
		cerr << "writing to file '" << filename << "'" << endl;
		#endif

		ofstream out(filename.c_str());
		out << "BASE_MODEL_BEGIN" << endl;
		writeParameters_(out, param_);
		hmm_.write(out);
		out << "BASE_MODEL_END" << endl;

		out << "PRECURSOR_MODEL_CR_BEGIN" << endl;
		writeParameters_(out, precursor_model_cr_.getParameters());
		precursor_model_cr_.getHMM().write(out);
		out << "PRECURSOR_MODEL_CR_END" << endl;
		
		out << "PRECURSOR_MODEL_CD_BEGIN" << endl;
		writeParameters_(out, precursor_model_cd_.getParameters());
		precursor_model_cd_.getHMM().write(out);
		out << "PRECURSOR_MODEL_CD_END" << endl;

		out << "BION_LOSS_MODEL_CR_BEGIN" << endl;
		writeParameters_(out, b_ion_losses_cr_.getParameters());
		b_ion_losses_cr_.getHMM().write(out);
		out << "BION_LOSS_MODEL_CR_END" << endl;

    out << "BION_LOSS_MODEL_CD_BEGIN" << endl;
		writeParameters_(out, b_ion_losses_cd_.getParameters());
    b_ion_losses_cd_.getHMM().write(out);
    out << "BION_LOSS_MODEL_CD_END" << endl;

		out << "B2ION_LOSS_MODEL_CR_BEGIN" << endl;
		writeParameters_(out, b2_ion_losses_cr_.getParameters());
		b2_ion_losses_cr_.getHMM().write(out);
		out << "B2ION_LOSS_MODEL_CR_END" << endl;

    out << "B2ION_LOSS_MODEL_CD_BEGIN" << endl;
		writeParameters_(out, b2_ion_losses_cd_.getParameters());
    b2_ion_losses_cd_.getHMM().write(out);
    out << "B2ION_LOSS_MODEL_CD_END" << endl;

    out << "AION_LOSS_MODEL_CR_BEGIN" << endl;
		writeParameters_(out, a_ion_losses_cr_.getParameters());
    a_ion_losses_cr_.getHMM().write(out);
    out << "AION_LOSS_MODEL_CR_END" << endl;

    out << "AION_LOSS_MODEL_CD_BEGIN" << endl;
		writeParameters_(out, a_ion_losses_cd_.getParameters());
    a_ion_losses_cd_.getHMM().write(out);
    out << "AION_LOSS_MODEL_CD_END" << endl;
		
		out << "YION_LOSS_MODEL_CR_BEGIN" << endl;
		writeParameters_(out, y_ion_losses_cr_.getParameters());
		y_ion_losses_cr_.getHMM().write(out);
		out << "YION_LOSS_MODEL_CR_END" << endl;

    out << "YION_LOSS_MODEL_CD_BEGIN" << endl;
		writeParameters_(out, y_ion_losses_cd_.getParameters());
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

		if (peptide.size() >= (Size)param_.getValue("visible_model_depth"))
		{
			cerr << "PILISModel: cannot train peptide " << peptide << " of length " << peptide.size() << " (max for this model is " << param_.getValue("visible_model_depth") << ", as defined by parameter \"visible_model_depth\")" << endl;
			return;
		}

		RichPeakSpectrum train_spec = in_spec;
		train_spec.sortByPosition();
		
		#ifdef TRAINING_DEBUG
		cout << "peptide: " << peptide  << "(z=" << charge << ")" << endl;
		#endif

		// get proton distribution
		vector<DoubleReal> bb_charge_full, sc_charge_full;
		prot_dist_.getProtonDistribution(bb_charge_full, sc_charge_full, peptide, charge, Residue::YIon);
		prot_dist_.setPeptideProtonDistribution(bb_charge_full, sc_charge_full);

		// get start probabilities
		vector<DoubleReal> bb_init, sc_init, cr_init;
		DoubleReal precursor_init(0);
		bool is_charge_remote = getInitialTransitionProbabilities_(bb_init, cr_init, sc_init, precursor_init, bb_charge_full, sc_charge_full, peptide);

		// clear the main Hidden Markov Model
		hmm_.clearInitialTransitionProbabilities();
		hmm_.clearTrainingEmissionProbabilities();
		
		//DoubleReal harge_sum(0);
		vector<AASequence> prefixes, suffixes;

		DoubleReal pep_weight(peptide.getMonoWeight());
		DoubleReal peptide_mz((pep_weight + charge)/(DoubleReal)charge);

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
			vector<DoubleReal> b_cr_ints(charge, 0), y_cr_ints(charge, 0), b_sc_ints(charge, 0), y_sc_ints(charge, 0), b_ints(charge, 0), y_ints(charge, 0), a_ints(charge, 0), ay_ints(charge, 0);
			prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_ints, y_ints, ProtonDistributionModel::ChargeDirected);

			prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::AIon, a_ints, ay_ints, ProtonDistributionModel::ChargeDirected);

			if (!is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("BB"+pos_name, bb_init[i]);
				hmm_.setInitialTransitionProbability(aa1+aa2+"_bxyz"+pos_name, bb_init[i]);
				hmm_.setInitialTransitionProbability(aa1+aa2+"_axyz"+pos_name, bb_init[i]);
			}

			prefixes.push_back(prefix);
			suffixes.push_back(suffix);
			
			if ((aa1 == "D" || aa1 == "E" || pos_name == "k-1" || pos_name == "k-2") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("CR"+pos_name, cr_init[i]);
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,	Residue::BIon, b_cr_ints, y_cr_ints, ProtonDistributionModel::ChargeRemote);
				//cerr << "ChargeStats: CR=" << CR_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;
			}

			if ((aa1 == "K" || aa1 == "H" || aa1 == "R") && is_charge_remote)
			{
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_sc_ints, y_sc_ints, ProtonDistributionModel::SideChain);
				hmm_.setInitialTransitionProbability("SC"+pos_name, /*sc_charge_full[i]*/ sc_init[i]);
				//cerr << "ChargeStats: SC=" << SC_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;	
			}

			
			UInt max_fragment_charge_training = param_.getValue("max_fragment_charge_training");
			Map<int, DoubleReal> sum_a, sum_b, sum_y;
			DoubleReal sum_a_ints(0.0), sum_b_ints(0.0), sum_y_ints(0.0);
      for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
      {
				sum_a[z] = 0;
				sum_b[z] = 0;
				sum_y[z] = 0;
        DoubleReal charge_remote_threshold((DoubleReal)param_.getValue("charge_remote_threshold"));
        DoubleReal avail_bb_sum_prefix = getAvailableBackboneCharge_(prefix, Residue::BIon, z);
        DoubleReal avail_bb_sum_suffix = getAvailableBackboneCharge_(suffix, Residue::YIon, z);

        if (prefix.size() != 2)
        {
          DoubleReal pre_weight = prefix.getMonoWeight(Residue::BIon);

          if (avail_bb_sum_prefix <= charge_remote_threshold)
          {
            //cerr << "Train b-ion losses, CR (avail=" << avail_bb_sum_prefix << ", z=" << z << ", prefix=" << prefix << endl;
            sum_b[z] += b_ion_losses_cr_.train(in_spec, prefix, pre_weight, z, pep_weight);
          }
          else
          {
            //cerr << "Train b-ion losses, CD (avail=" << avail_bb_sum_prefix << ", z=" << z << ", prefix=" << prefix << endl;
            sum_b[z] += b_ion_losses_cd_.train(in_spec, prefix, pre_weight, z, pep_weight);
          }
        }
        else
        {
          DoubleReal pre_weight = prefix.getMonoWeight(Residue::BIon);

          if (avail_bb_sum_prefix <= charge_remote_threshold)
          {
            sum_b[z] += b2_ion_losses_cr_.train(in_spec, prefix, pre_weight, z, pep_weight);
          }
          else
          {
            sum_b[z] += b2_ion_losses_cd_.train(in_spec, prefix, pre_weight, z, pep_weight);
          }
        }

        DoubleReal a_pre_weight = prefix.getMonoWeight(Residue::AIon);
        if (avail_bb_sum_prefix <= charge_remote_threshold)
        {
          sum_a[z] += a_ion_losses_cr_.train(in_spec, prefix, a_pre_weight, z, pep_weight);
        }
        else
        {
					sum_a[z] += a_ion_losses_cd_.train(in_spec, prefix, a_pre_weight, z, pep_weight);
        }

        DoubleReal suf_weight = suffix.getMonoWeight(Residue::YIon);
        if (avail_bb_sum_suffix <= charge_remote_threshold)
        {
          //cerr << "Train y-ion losses, CR (avail=" << avail_bb_sum_suffix << ", z=" << z << ", suffix=" << suffix << endl;
          sum_y[z] += y_ion_losses_cr_.train(in_spec, suffix, suf_weight, z, pep_weight);
        }
        else
        {
          //cerr << "Train y-ion losses, CD (avail=" << avail_bb_sum_suffix << ", z=" << z << ", suffix=" << suffix << endl;
          sum_y[z] += y_ion_losses_cd_.train(in_spec, suffix, suf_weight, z, pep_weight);
        }

				sum_a_ints += sum_a[z];
				sum_b_ints += sum_b[z];
				sum_y_ints += sum_y[z];
      }

			// end state is always needed
			hmm_.setTrainingEmissionProbability("end"+pos_name, 0.5/(DoubleReal(peptide.size() - 1)));
			if (!is_charge_remote)
			{
				//hmm_.setTrainingEmissionProbability("AA"+pos_name, sum_a_ints + sum_b_ints + sum_y_ints);
				hmm_.enableTransition("AA" + pos_name, aa1 + aa2 + "_bxyz" + pos_name);
				hmm_.enableTransition("AA" + pos_name, aa1 + aa2 + "_axyz" + pos_name);

				// now enable the states
				hmm_.enableTransition("BB"+pos_name, "AA"+pos_name);
				hmm_.enableTransition("BB"+pos_name, "end"+pos_name);

				// bxyz path
				hmm_.enableTransition(aa1+aa2+"_bxyz" + pos_name, "bxyz"+pos_name);
				hmm_.enableTransition(aa1+aa2+"_bxyz" + pos_name, "end"+pos_name);
			
				DoubleReal pre_emission_prob(0), suf_emission_prob(0);
				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					pre_emission_prob += b_ints[z - 1] * sum_b[z];
				 	suf_emission_prob += y_ints[z - 1] * sum_y[z];
				}
				hmm_.setTransitionProbability("bxyz" + pos_name, "bxyz_" + pos_name + "-prefix", 1.0);
				hmm_.setTransitionProbability("bxyz" + pos_name, "bxyz_" + pos_name + "-suffix", 1.0);
				hmm_.enableTransition("bxyz_" + pos_name + "-prefix", "bxyz_" + pos_name + "-prefix-ions");
				hmm_.enableTransition("bxyz_" + pos_name + "-prefix", "end" + pos_name);
				hmm_.enableTransition("bxyz_" + pos_name + "-suffix", "bxyz_" + pos_name + "-suffix-ions");
				hmm_.enableTransition("bxyz_" + pos_name + "-suffix", "end" + pos_name);

				hmm_.setTrainingEmissionProbability("bxyz_" + pos_name + "-prefix-ions", pre_emission_prob);
				hmm_.setTrainingEmissionProbability("bxyz_" + pos_name + "-suffix-ions", suf_emission_prob);

				// axyz path
				hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "axyz"+pos_name);
				hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "end"+pos_name);
			
				pre_emission_prob = 0;
				suf_emission_prob = 0;
				for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
				{
					pre_emission_prob += a_ints[z - 1] * sum_a[z];
					//suf_emission_prob += ay_ints[z] * sum_y[z]; // TODO test if this can be done, or ruins the a-ion intensities
				}
				hmm_.setTransitionProbability("axyz" + pos_name, "axyz_" + pos_name + "-prefix", 1.0);
				hmm_.enableTransition("axyz_" + pos_name + "-prefix", "axyz_" + pos_name + "-prefix-ions");
				hmm_.enableTransition("axyz_" + pos_name + "-prefix", "end" + pos_name);
				hmm_.setTrainingEmissionProbability("axyz_" + pos_name + "-prefix-ions", pre_emission_prob);

			}
			
			const String cr_aa("DE");
			for (String::ConstIterator it = cr_aa.begin(); it != cr_aa.end(); ++it)
			{
				if (is_charge_remote && aa1 == String(*it))
				{
					hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
					hmm_.enableTransition("CR"+pos_name, "end"+pos_name);

					DoubleReal pre_emission_prob(0), suf_emission_prob(0);
					for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
					{
						pre_emission_prob += b_cr_ints[z - 1] * sum_b[z];
						suf_emission_prob += y_cr_ints[z - 1] * sum_y[z];
					}

					//hmm_.setTrainingEmissionProbability("A"+pos_name, pre_emission_prob +  suf_emission_prob);
					hmm_.enableTransition("A" + pos_name, aa2 + "_" + String(*it) + pos_name);

					hmm_.setTransitionProbability(String(*it) + pos_name, String(*it) + "_" + pos_name + "-prefix", 1.0);
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-prefix", String(*it) + "_" + pos_name + "-prefix-ions");
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-prefix", "end" + pos_name);
					hmm_.setTrainingEmissionProbability(String(*it) + "_" + pos_name + "-prefix-ions", pre_emission_prob);

					hmm_.setTransitionProbability(String(*it) + pos_name, String(*it) + "_" + pos_name + "-suffix", 1.0);
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-suffix", String(*it) + "_" + pos_name + "-suffix-ions");
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-suffix", "end" + pos_name);
					hmm_.setTrainingEmissionProbability(String(*it) + "_" + pos_name + "-suffix-ions", suf_emission_prob);

					hmm_.setInitialTransitionProbability(aa2+"_" + String(*it) + pos_name, cr_init[i]);
					hmm_.enableTransition(aa2+"_" + String(*it) + pos_name, String(*it) + pos_name);
					hmm_.enableTransition(aa2+"_" + String(*it) + pos_name, "end" + pos_name);
				}
			}

			if (is_charge_remote && !peptide.has("D") && !peptide.has("E"))
			{

				if (pos_name == "k-1")
				{
					hmm_.enableTransition("CRk-1", "Ak-1");
					hmm_.enableTransition("CRk-1", "endk-1");

					DoubleReal pre_emission_prob(0), suf_emission_prob(0);
					for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
					{
						pre_emission_prob += b_cr_ints[z - 1] * sum_b[z];
						suf_emission_prob += y_cr_ints[z - 1] * sum_y[z];
					}

					//hmm_.setTrainingEmissionProbability("Ak-1", pre_emission_prob + suf_emission_prob);
					hmm_.enableTransition("Ak-1", aa1 + "_bk-1");


					hmm_.setTransitionProbability("bk-1", "bk-1_-prefix", 1.0);
					hmm_.enableTransition("bk-1_-prefix", "bk-1_-prefix-ions");
					hmm_.enableTransition("bk-1_-prefix", "endk-1");
					hmm_.setTrainingEmissionProbability("bk-1_-prefix-ions", pre_emission_prob);

					hmm_.setTransitionProbability("bk-1", "bk-1_-suffix", 1.0);
          hmm_.enableTransition("bk-1_-suffix", "bk-1_-suffix-ions");
          hmm_.enableTransition("bk-1_-suffix", "endk-1");
          hmm_.setTrainingEmissionProbability("bk-1_-suffix-ions", suf_emission_prob);

					hmm_.setTransitionProbability("bk-1", "bk-1_-ions", 1.0);

					hmm_.setInitialTransitionProbability(aa1+"_bk-1", cr_init[i]);
					hmm_.enableTransition(aa1+"_bk-1", "bk-1");
					hmm_.enableTransition(aa1+"_bk-1", "endk-1");
				}

				if (pos_name == "k-2")
				{
					hmm_.enableTransition("CRk-2", "Ak-2");
          hmm_.enableTransition("CRk-2", "endk-2");

          DoubleReal pre_emission_prob(0), suf_emission_prob(0);
          for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
          {
            pre_emission_prob += b_cr_ints[z - 1] * sum_b[z];
            suf_emission_prob += y_cr_ints[z - 1] * sum_y[z];
          }

          //hmm_.setTrainingEmissionProbability("Ak-2", pre_emission_prob + suf_emission_prob);
					hmm_.enableTransition("Ak-2", aa1 + "_bk-2");


          hmm_.setTransitionProbability("bk-2", "bk-2_-prefix", 1.0);
          hmm_.enableTransition("bk-2_-prefix", "bk-2_-prefix-ions");
          hmm_.enableTransition("bk-2_-prefix", "endk-2");
          hmm_.setTrainingEmissionProbability("bk-2_-prefix-ions", pre_emission_prob);

          hmm_.setTransitionProbability("bk-2", "bk-2_-suffix", 1.0);
          hmm_.enableTransition("bk-2_-suffix", "bk-2_-suffix-ions");
          hmm_.enableTransition("bk-2_-suffix", "endk-2");
          hmm_.setTrainingEmissionProbability("bk-2_-suffix-ions", suf_emission_prob);

          hmm_.setTransitionProbability("bk-2", "bk-2_-ions", 1.0);

          hmm_.setInitialTransitionProbability(aa1+"_bk-2", cr_init[i]);
          hmm_.enableTransition(aa1+"_bk-2", "bk-2");
          hmm_.enableTransition(aa1+"_bk-2", "endk-2");
				}
			}
			
			const String sc_aa("KRH");
			for (String::ConstIterator it = sc_aa.begin(); it != sc_aa.end(); ++it)
			{
				if (is_charge_remote && aa1 == String(*it))
				{

					hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
					hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
				
					DoubleReal pre_emission_prob(0), suf_emission_prob(0);
					for (UInt z = 1; z <= max_fragment_charge_training && z <= charge; ++z)
					{
						pre_emission_prob += b_sc_ints[z - 1] * sum_b[z];
						suf_emission_prob += y_sc_ints[z - 1] * sum_y[z];
					}

					//hmm_.setTrainingEmissionProbability("ASC"+pos_name, pre_emission_prob + suf_emission_prob);
					hmm_.enableTransition("ASC" + pos_name, aa2 + "_" + String(*it) + pos_name);

					hmm_.setTransitionProbability(String(*it) + pos_name, String(*it) + "_" + pos_name + "-prefix", 1.0);
					hmm_.setTransitionProbability(String(*it) + pos_name, String(*it) + "_" + pos_name + "-suffix", 1.0);
					
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-prefix", String(*it) + "_" + pos_name + "-prefix-ions");
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-prefix", "end" + pos_name);
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-suffix", String(*it) + "_" + pos_name + "-suffix-ions");
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-suffix", "end" + pos_name);

					hmm_.setTrainingEmissionProbability(String(*it) + "_" + pos_name + "-prefix-ions", pre_emission_prob);
					hmm_.setTrainingEmissionProbability(String(*it) + "_" + pos_name + "-suffix-ions", suf_emission_prob);


					hmm_.setInitialTransitionProbability(aa2 + "_" + String(*it) + pos_name, sc_init[i]);
        	hmm_.enableTransition(aa2+"_" + String(*it) + pos_name, String(*it) + pos_name);
					hmm_.enableTransition(aa2+"_" + String(*it) + pos_name, "end"+pos_name);

				}
			}
		}
	
		if (is_charge_remote)
		{
			precursor_model_cr_.train(in_spec, peptide, peptide_mz, charge, pep_weight);
		}
		else
		{
			precursor_model_cd_.train(in_spec, peptide, peptide_mz, charge, pep_weight);
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

		if (peptide.size() > (Size)param_.getValue("visible_model_depth"))
		{
			cerr << "PILISModel: cannot generate spectra of peptide '" << peptide << "' of length " << peptide.size() << " (max of this model is " << param_.getValue("visible_model_depth") << ", as defined by \"visible_model_depth\")" << endl;
			return;
		}
		
		// calc proton distribution
		vector<DoubleReal> sc_charge_full, bb_charge_full;
		//cerr << "getProtonDistribution: " << peptide << " " << charge << endl;
		prot_dist_.getProtonDistribution(bb_charge_full, sc_charge_full, peptide, charge, Residue::YIon);
		prot_dist_.setPeptideProtonDistribution(bb_charge_full, sc_charge_full);
		
		#ifdef INIT_CHARGE_DEBUG
		for (Size i = 0; i != bb_charge_full.size() - 1; ++i)
		{
			cerr << "i: bb=" << bb_charge_full[i] << ", sc=" << sc_charge_full[i] << endl;
		}
		#endif

		hmm_.clearInitialTransitionProbabilities();
		hmm_.clearTrainingEmissionProbabilities();
	
    // set charges
    vector<DoubleReal> bb_init, sc_init, cr_init;
		DoubleReal precursor_init(0);
		bool is_charge_remote = getInitialTransitionProbabilities_(bb_init, cr_init, sc_init, precursor_init, bb_charge_full, sc_charge_full, peptide);
		#ifdef INIT_CHARGE_DEBUG
		cerr << "is_charge_remote=" << is_charge_remote << endl;
		for (Size i = 0; i != bb_init.size(); ++i)
		{
			cerr << "bb=" << bb_init[i] << ", cr=" << cr_init[i] << ", sc=" << sc_init[i] << endl;
		}
		cerr << "precursor_init=" << precursor_init << endl;
		#endif
	
		vector<AASequence> suffixes, prefixes;
		vector<String> pos_names;
	
		vector<vector<DoubleReal> > all_b_cr_ints, all_y_cr_ints, all_b_sc_ints, all_y_sc_ints, all_b_ints, all_y_ints, all_a_ints, all_ay_ints;
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

			hmm_.setInitialTransitionProbability("BB"+pos_name, bb_init[i]);

			vector<DoubleReal> b_cr_ints(charge, 0), y_cr_ints(charge, 0), b_sc_ints(charge, 0), y_sc_ints(charge, 0), b_ints(charge, 0), y_ints(charge, 0), a_ints(charge, 0), ay_ints(charge, 0);
			if ((aa1 == "D" || aa1 == "E" || pos_name == "k-1" || pos_name == "k-2") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("CR"+pos_name, cr_init[i]);
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_cr_ints, y_cr_ints, ProtonDistributionModel::ChargeRemote);
			}

			if ((aa1 == "K" || aa1 == "H" || aa1 == "R") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("SC"+pos_name, sc_init[i]);
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,	Residue::BIon, b_sc_ints, y_sc_ints, ProtonDistributionModel::SideChain);
			}
		
      prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_ints, y_ints, ProtonDistributionModel::ChargeDirected);
      prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::AIon, a_ints, ay_ints, ProtonDistributionModel::ChargeDirected);	

			//if (charge == 3)
			//{
			#ifdef INIT_CHARGE_DEBUG
				cerr << prefix << " - " << suffix << " ";
				for (Size ions_charge = 0; ions_charge != b_ints.size(); ++ions_charge)
				{
					cerr << b_ints[ions_charge] << " - " << y_ints[ions_charge] << " | ";
				}
				cerr << endl;
			#endif
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

			hmm_.setTransitionProbability("bxyz" + pos_name, "bxyz_" + pos_name + "-prefix", 1.0);
			hmm_.setTransitionProbability("bxyz" + pos_name, "bxyz_" + pos_name + "-suffix", 1.0);

			hmm_.enableTransition("bxyz_" + pos_name + "-prefix", "bxyz_" + pos_name + "-prefix-ions");
			hmm_.enableTransition("bxyz_" + pos_name + "-prefix", "end" + pos_name);
			hmm_.enableTransition("bxyz_" + pos_name + "-suffix", "bxyz_" + pos_name + "-suffix-ions");
			hmm_.enableTransition("bxyz_" + pos_name + "-suffix", "end" + pos_name);


			// axyz
			hmm_.enableTransition("AA"+pos_name, aa1+aa2+"_axyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "axyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"_axyz"+pos_name, "end"+pos_name);
			hmm_.setTransitionProbability("axyz" + pos_name, "axyz_" + pos_name + "-prefix", 1.0);
			hmm_.enableTransition("axyz_" + pos_name + "-prefix", "axyz_" + pos_name + "-prefix-ions");
			hmm_.enableTransition("axyz_" + pos_name + "-prefix", "end" + pos_name);

			const String cr_aa("DE");
			for (String::ConstIterator it = cr_aa.begin(); it != cr_aa.end(); ++it)
			{
      	if (is_charge_remote && aa1 == String(*it))
				{      	
					hmm_.setTransitionProbability(String(*it) + pos_name, String(*it) + "_" + pos_name + "-prefix", 1.0);
					hmm_.setTransitionProbability(String(*it) + pos_name, String(*it) + "_" + pos_name + "-suffix", 1.0);

					hmm_.enableTransition(String(*it) + "_" + pos_name + "-prefix", String(*it) + "_" + pos_name + "-prefix-ions");
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-prefix", "end" + pos_name);
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-suffix", String(*it) + "_" + pos_name + "-suffix-ions");
					hmm_.enableTransition(String(*it) + "_" + pos_name + "-suffix", "end" + pos_name);

					hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
					hmm_.enableTransition("CR"+pos_name, "end"+pos_name);
					hmm_.enableTransition("A"+pos_name, aa2+"_" + String(*it) + pos_name);
        	hmm_.enableTransition(aa2+"_" + String(*it) + pos_name, String(*it) + pos_name);
					hmm_.enableTransition(aa2+"_" + String(*it) + pos_name, "end"+pos_name);
				}
			}

			if (!peptide.has("D") && !peptide.has("E") && is_charge_remote)
			{
      	if (pos_name == "k-1")
      	{
					hmm_.setTransitionProbability("bk-1", "bk-1_-prefix", 1.0);
					hmm_.setTransitionProbability("bk-1", "bk-1_-suffix", 1.0);

					hmm_.enableTransition("bk-1_-prefix", "bk-1_-prefix-ions");
					hmm_.enableTransition("bk-1_-prefix", "endk-1" );
					hmm_.enableTransition("bk-1_-suffix", "bk-1_-suffix-ions");
					hmm_.enableTransition("bk-1_-suffix", "endk-1");

					hmm_.enableTransition("CRk-1", "Ak-1");
        	hmm_.enableTransition("CRk-1", "endk-1");
					hmm_.enableTransition("Ak-1", aa1+"_bk-1");

					hmm_.enableTransition(aa1 + "_bk-1", "bk-1");
        	hmm_.enableTransition(aa1 + "_bk-1", "endk-1");
      	}

				if (pos_name == "k-2")
				{
          hmm_.setTransitionProbability("bk-2", "bk-2_-prefix", 1.0);
          hmm_.setTransitionProbability("bk-2", "bk-2_-suffix", 1.0);

          hmm_.enableTransition("bk-2_-prefix", "bk-2_-prefix-ions");
          hmm_.enableTransition("bk-2_-prefix", "endk-2" );
          hmm_.enableTransition("bk-2_-suffix", "bk-2_-suffix-ions");
          hmm_.enableTransition("bk-2_-suffix", "endk-2");

          hmm_.enableTransition("CRk-2", "Ak-2");
          hmm_.enableTransition("CRk-2", "endk-2");
          hmm_.enableTransition("Ak-2", aa1+"_bk-2");

          hmm_.enableTransition(aa1 + "_bk-2", "bk-2");
          hmm_.enableTransition(aa1 + "_bk-2", "endk-2");
				}
			}

			const String sc_aa("HKR");
			for (String::ConstIterator it = sc_aa.begin(); it != sc_aa.end(); ++it)
			{
				if (is_charge_remote && aa1 == String(*it))
				{
					hmm_.setTransitionProbability(String(*it) + pos_name, String(*it) + "_" + pos_name + "-prefix", 1.0);
          hmm_.setTransitionProbability(String(*it) + pos_name, String(*it) + "_" + pos_name + "-suffix", 1.0);

          hmm_.enableTransition(String(*it) + "_" + pos_name + "-prefix", String(*it) + "_" + pos_name + "-prefix-ions");
          hmm_.enableTransition(String(*it) + "_" + pos_name + "-prefix", "end" + pos_name);
          hmm_.enableTransition(String(*it) + "_" + pos_name + "-suffix", String(*it) + "_" + pos_name + "-suffix-ions");
          hmm_.enableTransition(String(*it) + "_" + pos_name + "-suffix", "end" + pos_name);



        	hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
					hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
        	hmm_.enableTransition("ASC"+pos_name, aa2+"_" + String(*it) + pos_name);
        	hmm_.enableTransition(aa2+"_" + String(*it) + pos_name, String(*it) + pos_name);
					hmm_.enableTransition(aa2+"_" + String(*it) + pos_name, "end"+pos_name);
				}
			}
		}

		//cerr << "Collecting peaks..." << endl;

    Map<HMMState*, DoubleReal> tmp;
		hmm_.calculateEmissionProbabilities(tmp);

		// clear peaks from last spectrum
		peaks_.clear();
		
		//stringstream peptide_ss;
		//peptide_ss << peptide;
		//hmm_.writeGraphMLFile(String("model_graph_train_"+peptide_ss.str()+"_"+String(charge)+".graphml").c_str());

		UInt max_isotope = (UInt)param_.getValue("max_isotope");
		UInt max_fragment_charge = (UInt)param_.getValue("max_fragment_charge");
		DoubleReal charge_remote_threshold = (DoubleReal)param_.getValue("charge_remote_threshold");
		for (Size i = 0; i != prefixes.size(); ++i)
		{
			String aa;
			AASequence aa_seq;
			aa_seq += peptide[i].getOneLetterCode();
			aa = aa_seq.toString();

			vector<DoubleReal> pre_path_intensities, suf_path_intensities;
			pre_path_intensities.push_back(tmp[hmm_.getState("bxyz_" + pos_names[i] + "-prefix-ions")]);
			suf_path_intensities.push_back(tmp[hmm_.getState("bxyz_" + pos_names[i] + "-suffix-ions")]);
			#ifdef INIT_CHARGE_DEBUG
			cerr << prefixes[i] << " - " << suffixes[i] << ": bxyz=" << pre_path_intensities.back() << "|" << suf_path_intensities.back() << " ";
			#endif

			const String cr_and_sc("DEHKR");
			for (String::ConstIterator it = cr_and_sc.begin(); it != cr_and_sc.end(); ++it)
			{
				if (aa == String(*it))
				{
					pre_path_intensities.push_back(tmp[hmm_.getState(String(*it) + "_" + pos_names[i] + "-prefix-ions")]);
					suf_path_intensities.push_back(tmp[hmm_.getState(String(*it) + "_" + pos_names[i] + "-suffix-ions")]);
					#ifdef INIT_CHARGE_DEBUG
					cerr << ": " << aa << "=" << pre_path_intensities.back() << "|" << suf_path_intensities.back() << " ";
					#endif
				}
			}
		
			if (i == peptide.size() - 2)
			{
				pre_path_intensities.push_back(tmp[hmm_.getState("bk-1_-prefix-ions")]);
				suf_path_intensities.push_back(tmp[hmm_.getState("bk-1_-suffix-ions")]);
			}
			if (i == peptide.size() - 3)
			{
				pre_path_intensities.push_back(tmp[hmm_.getState("bk-2_-prefix-ions")]);
				suf_path_intensities.push_back(tmp[hmm_.getState("bk-2_-suffix-ions")]);
			}

			pre_path_intensities.push_back(tmp[hmm_.getState("axyz_" + pos_names[i] + "-prefix-ions")]);
			suf_path_intensities.push_back(tmp[hmm_.getState("axyz_" + pos_names[i] + "-suffix-ions")]);

			vector<vector<vector<DoubleReal> > > prefix_intensities, suffix_intensities;
			prefix_intensities.push_back(all_b_ints);
			suffix_intensities.push_back(all_y_ints);
			const String sc_aa("HKR");
			if (sc_aa.hasSubstring(aa))
			{
				prefix_intensities.push_back(all_b_sc_ints);
				suffix_intensities.push_back(all_y_sc_ints);
			}
			const String cr_aa("DE");
			if (cr_aa.hasSubstring(aa))
			{
				prefix_intensities.push_back(all_b_cr_ints);
				suffix_intensities.push_back(all_y_cr_ints);	
			}

			if (i == peptide.size() - 2 || i == peptide.size() - 3)
			{
				prefix_intensities.push_back(all_b_cr_ints);
				suffix_intensities.push_back(all_y_cr_ints);
			}

			prefix_intensities.push_back(all_a_ints);
			suffix_intensities.push_back(all_ay_ints);

			for (Size int_it = 0; int_it != pre_path_intensities.size(); ++int_it)
			{
				if (pre_path_intensities[int_it] < MIN_DECIMAL_VALUE && suf_path_intensities[int_it] < MIN_DECIMAL_VALUE)
				{
					continue;
				}
				DoubleReal prefix_weight = prefixes[i].getMonoWeight(Residue::BIon);
				DoubleReal suffix_weight = suffixes[i].getMonoWeight(Residue::YIon);
				IsotopeDistribution prefix_id;
				
				if (int_it != pre_path_intensities.size() - 1)
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
					if (prefix_intensities[int_it][i][z - 1] > MIN_DECIMAL_VALUE)
					{
						vector<RichPeak1D> b_loss_peaks;
						DoubleReal avail_bb_sum_prefix(0);
						
						if (int_it != pre_path_intensities.size() - 1)
						{
							avail_bb_sum_prefix = getAvailableBackboneCharge_(prefixes[i], Residue::BIon, z);
							if (avail_bb_sum_prefix <= charge_remote_threshold)
            	{
								if (i == 1) // second BB position
								{
									b2_ion_losses_cr_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * pre_path_intensities[int_it]);
								}
								else
								{
              		b_ion_losses_cr_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * pre_path_intensities[int_it]);
								}
            	}
            	else
            	{
								if (i == 1)
								{
									b2_ion_losses_cd_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * pre_path_intensities[int_it]);
								}
              	else
								{
									b_ion_losses_cd_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * pre_path_intensities[int_it]);
								}
            	}
						}
						else
						{
							avail_bb_sum_prefix = getAvailableBackboneCharge_(prefixes[i], Residue::AIon, z);
							if (avail_bb_sum_prefix <= charge_remote_threshold)
              {
                a_ion_losses_cr_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * pre_path_intensities[int_it]);
              }
              else
              {
                a_ion_losses_cd_.getIons(b_loss_peaks, prefixes[i], prefix_intensities[int_it][i][z - 1] * pre_path_intensities[int_it]);
              }
						}

						for (vector<RichPeak1D>::const_iterator it = b_loss_peaks.begin(); it != b_loss_peaks.end(); ++it)
						{
        			String b_ion_name = it->getMetaValue("IonName");
							#ifdef INIT_CHARGE_DEBUG
							cerr << b_ion_name << " " << it->getMZ() << " " << it->getIntensity() << endl;
							#endif
        			vector<String> split;
        			b_ion_name.split('-', split);
        			if (split.empty())
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

							b_ion_name += String(z, '+');
							if (int_it != pre_path_intensities.size() - 1)
							{
          			addPeaks_(prefix_weight, z, it->getMZ(), it->getIntensity(), spec, prefix_id, b_ion_name);
							}
							else
							{
								addPeaks_(prefix_weight - 28.0, z, it->getMZ(), it->getIntensity(), spec, prefix_id, b_ion_name);
							}
						}
					}

					if (suffix_intensities[int_it][i][z - 1] > MIN_DECIMAL_VALUE && int_it != suf_path_intensities.size() - 1) // no ay 
					{
						vector<RichPeak1D> y_loss_peaks;
						DoubleReal avail_bb_sum_suffix = getAvailableBackboneCharge_(suffixes[i], Residue::YIon, z);
						if (avail_bb_sum_suffix <= charge_remote_threshold)
						{
							y_ion_losses_cr_.getIons(y_loss_peaks, suffixes[i], suffix_intensities[int_it][i][z - 1] * suf_path_intensities[int_it]);
						}
						else
						{
							y_ion_losses_cd_.getIons(y_loss_peaks, suffixes[i], suffix_intensities[int_it][i][z - 1] * suf_path_intensities[int_it]);
						}

						for (vector<RichPeak1D>::const_iterator it = y_loss_peaks.begin(); it != y_loss_peaks.end(); ++it)
        		{
          		String y_ion_name = it->getMetaValue("IonName");
							#ifdef INIT_CHARGE_DEBUG
							cerr << y_ion_name << " " << it->getMZ() << " " << it->getIntensity() << endl;
							#endif
          		vector<String> split;
          		y_ion_name.split('-', split);
          		if (split.empty())
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

							y_ion_name += String(z, '+');
            	addPeaks_(suffix_weight, z, it->getMZ(), it->getIntensity(), spec, suffix_id, y_ion_name);
						}
        	}
				}
			}
		}

		hmm_.disableTransitions();

		// precursor intensities
		vector<RichPeak1D> pre_peaks;
		if (is_charge_remote)
		{
			precursor_model_cr_.getIons(pre_peaks, peptide, precursor_init);
		}
		else
		{
			precursor_model_cd_.getIons(pre_peaks, peptide, precursor_init);
		}
		
		DoubleReal weight = peptide.getMonoWeight();
		IsotopeDistribution id(max_isotope);
		id.estimateFromPeptideWeight(weight);	
		for (vector<RichPeak1D>::const_iterator it = pre_peaks.begin(); it != pre_peaks.end(); ++it)
		{
			addPeaks_(weight, charge, it->getMZ(), it->getIntensity(), spec, id, it->getMetaValue("IonName"));
		}
	
		// now build the spectrum with the peaks
		DoubleReal intensity_max(0);
		for (Map<DoubleReal, vector<RichPeak1D> >::ConstIterator it = peaks_.begin(); it != peaks_.end(); ++it)
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
				DoubleReal int_sum(0);
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

		RichPeakSpectrum new_spec;
		// add up peak intensities which are close together
		vector<RichPeak1D> close_peaks;
		for (RichPeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
		{
			// empty
			if (close_peaks.empty())
			{
				close_peaks.push_back(*it);
				continue;
			}
			
			// include peak in cluster
			if (it->getPosition() - close_peaks.begin()->getMZ() < 0.1)
			{
				close_peaks.push_back(*it);
				continue;
			}
			else
			{
				// move cluster to 
				if (close_peaks.size() > 1)
				{
					DoubleReal mz(0), intensity(0);
					String name;
					for (vector<RichPeak1D>::const_iterator pit = close_peaks.begin(); pit != close_peaks.end(); ++pit)
					{
						mz += pit->getMZ();
						intensity += pit->getIntensity();
						String p_name = pit->getMetaValue("IonName");
						if (p_name != "")
						{
							name += p_name;
							if (pit != close_peaks.end() - 1)
							{
								name += "/";
							}
						}
					}
					RichPeak1D peak;
					peak.setPosition(mz / (DoubleReal)close_peaks.size());
					peak.setIntensity(intensity);
					peak.setMetaValue("IonName", name);
					new_spec.push_back(peak);

					// clear the actual cluster and actual peak
					close_peaks.clear();
					close_peaks.push_back(*it);
					continue;
				}
				// cluster has only one peak, add it to the spec and start new cluster
				new_spec.push_back(*close_peaks.begin());
				close_peaks.clear();
				close_peaks.push_back(*it);
			}
		}
		spec = new_spec;

		
		DoubleReal min_y_int((DoubleReal)param_.getValue("min_y_ion_intensity"));
		DoubleReal min_b_int((DoubleReal)param_.getValue("min_b_ion_intensity"));
		DoubleReal min_a_int((DoubleReal)param_.getValue("min_a_ion_intensity"));
		DoubleReal min_y_loss_int((DoubleReal)param_.getValue("min_y_loss_intensity"));
		DoubleReal min_b_loss_int((DoubleReal)param_.getValue("min_b_loss_intensity"));
		

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

	DoubleReal PILISModel::getAvailableBackboneCharge_(const AASequence& ion, Residue::ResidueType res_type, int charge)
	{
		DoubleReal bb_sum(0);
		vector<DoubleReal> bb_charges, sc_charges;
		prot_dist_.getProtonDistribution(bb_charges, sc_charges, ion, charge, res_type);

		for (vector<DoubleReal>::const_iterator it = bb_charges.begin(); it != bb_charges.end(); ++it)
		{
			bb_sum += *it;
		}

    // activation of protons sitting at lysine and histidine side chains
    DoubleReal side_chain_activation(param_.getValue("side_chain_activation"));
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

    if (bb_sum < (DoubleReal)param_.getValue("charge_directed_threshold") * charge)
    {
      bb_sum = (DoubleReal)param_.getValue("charge_directed_threshold") * charge;
    }
		return bb_sum;
	}
	
	
	bool PILISModel::getInitialTransitionProbabilities_(std::vector<DoubleReal>& bb_init, 
																											std::vector<DoubleReal>& cr_init, 
																											std::vector<DoubleReal>& sc_init, 
																											DoubleReal& precursor_init,
																											const vector<DoubleReal>& bb_charges, 
																											const vector<DoubleReal>& sc_charges, 
																											const AASequence& peptide)
	{
		bool is_charge_remote(false);

    DoubleReal bb_sum(0), bb_sum_orig(0);
    for (vector<DoubleReal>::const_iterator it = bb_charges.begin(); it != bb_charges.end(); ++it)
    {
      bb_sum += *it;
    }

		if (bb_sum > 1)
    {
      bb_sum = 1;
    }

		#ifdef INIT_CHARGE_DEBUG
		cerr << "bb_sum=" << bb_sum << endl;
		#endif
		bb_sum_orig = bb_sum;

		if (bb_sum < (DoubleReal)param_.getValue("charge_remote_threshold"))
    {
			is_charge_remote = true;
		}

		if (bb_sum < (DoubleReal)param_.getValue("charge_directed_threshold"))
		{
			bb_sum = (DoubleReal)param_.getValue("charge_directed_threshold");
		}
    
    // side-chain activiation
		DoubleReal side_chain_activation(param_.getValue("side_chain_activation"));
		for (Size i = 0; i != peptide.size(); ++i)
		{
			if (peptide[i].getOneLetterCode() != "R")
			{
				bb_sum += side_chain_activation * sc_charges[i];
			}
		}

		if (bb_sum > 1)
		{
			bb_sum = 1;
		}
	
		#ifdef INIT_CHARGE_DEBUG
		cerr << "bb_sum after side-chain-activiation=" << bb_sum << endl;
		#endif

		vector<DoubleReal> bb_charges_all = bb_charges;
		sort(bb_charges_all.begin(), bb_charges_all.end());
		DoubleReal bb_charges_median = bb_charges_all[(Size)(bb_charges_all.size()/2.0)];
		#ifdef INIT_CHARGE_DEBUG
		cerr << "bb_charges_median=" << bb_charges_median << endl;
		#endif

		DoubleReal min_enhancement_factor = param_.getValue("min_enhancement_factor");
		DoubleReal blocker_sum(1.0);
		for (Size i = 0; i != peptide.size() - 1; ++i)
    {
			DoubleReal bb_enhance_factor(max(min_enhancement_factor, sqrt(bb_charges[i+1] / bb_charges_median)));
			#ifdef INIT_CHARGE_DEBUG
			cerr << "bb_enhance_factor=" << bb_enhance_factor << endl;
			#endif
			
			if (sc_charges[i] != 0)
			{
				blocker_sum += 10.0 * sc_charges[i];
			}

      bb_init.push_back(/*bb_charges[i+1]  **/ bb_sum * bb_enhance_factor / blocker_sum);
      String aa(peptide[i].getOneLetterCode());
      if ((aa == "K" || aa == "R" || aa == "H"))
      {
        sc_init.push_back(sc_charges[i] /* * bb_charges_median*/);
      }
      else
      {
        sc_init.push_back(0.0);
      }

			if (is_charge_remote && (aa == "D" || aa == "E" || i == peptide.size() - 2 || i == peptide.size() - 3))
      {
        cr_init.push_back(((1 - bb_sum_orig) /** bb_charges_median*/));
      }
      else
      {
        cr_init.push_back(0.0);
      }
    }
		precursor_init = (1 - bb_sum) /** bb_charges_median*//*bb_avg*/ / 10.0;

		// normalize the initial probability values
		DoubleReal init_prob_sum(0);
		for (Size i = 0; i != bb_init.size(); ++i)
		{
			init_prob_sum += bb_init[i] + sc_init[i] + cr_init[i];
		}
		init_prob_sum += precursor_init;

		for (Size i = 0; i != bb_init.size(); ++i)
		{
			bb_init[i] /= init_prob_sum;
			sc_init[i] /= init_prob_sum;
			cr_init[i] /= init_prob_sum;
		}
		precursor_init /= init_prob_sum;
		
		return is_charge_remote;
	}

	void PILISModel::addPeaks_(DoubleReal mz, int charge, DoubleReal offset, DoubleReal intensity, RichPeakSpectrum& /*spectrum*/, const IsotopeDistribution& id, const String& name)
	{
		if (intensity < MIN_DECIMAL_VALUE)
		{
			return;
		}
		static RichPeak1D p;
		UInt i = 0;
		for (IsotopeDistribution::ConstIterator it = id.begin(); it != id.end(); ++it, ++i)
		{
			DoubleReal pos = (mz + i + charge + offset) / (DoubleReal)charge;
			p.setPosition(pos);
			if (it == id.begin())
			{
				p.setMetaValue("IonName", String(name.c_str()));
			}

			if (pos >= (DoubleReal)param_.getValue("lower_mz") && pos <= (DoubleReal)param_.getValue("upper_mz"))
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
	
	void PILISModel::parseHMMModel_(const TextFile::ConstIterator& begin, const TextFile::ConstIterator& end, HiddenMarkovModel& hmm, Param& param)
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
      line.split(' ', split, true);

      if ( !split.empty() )
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

				if (id == "Parameter")
				{
					// Parameter charge_remote_threshold 0.3 float
					if (split[split.size() - 1] == "float")
					{
						param.setValue(split[1], split[2].toDouble());
					}
					else if (split[split.size() - 1] == "int")
					{
						param.setValue(split[1], split[2].toInt());
					}
					else if (split[split.size() - 1] == "string_list")
					{
						String tmp_list;
						for (Size i = 2; i < split.size() - 1; ++i)
						{
							tmp_list += split[i];
						}
						param.setValue(split[1], StringList::create(tmp_list));
					}
					else if (split[split.size() - 1] == "string")
					{
						param.setValue(split[1], split[2]);
					}
					else 
					{
						throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, line);
					}
				}
      }
    }
    //hmm_.disableTransitions();
		//hmm.buildSynonyms();
		hmm.disableTransitions();
	
		//cerr << hmm_.getNumberOfStates() << endl;
		
		return;
	}

	void PILISModel::writeParameters_(ostream& os, const Param& param)
	{
		for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
		{
			os << "Parameter ";
			if (it->value.valueType() == DataValue::DOUBLE_VALUE)
			{
				os << it->name << " \"" << it->value << "\" float" << endl;
				continue;
			}
			if (it->value.valueType() == DataValue::INT_VALUE)
			{
				os << it->name << " \"" << it->value << "\" int" << endl;
				continue;
			}
			if (it->value.valueType() == DataValue::STRING_LIST)
			{
				StringList value = (StringList)it->value;
				String tmp;
				tmp.concatenate(value.begin(), value.end(), ",");
				os << it->name << " \"" << tmp << "\" string_list" << endl;
				continue;
			}
			if (it->value.valueType() == DataValue::STRING_VALUE)
			{
				os << it->name << " \"" << it->value << "\" string" << endl;
				continue;
			}
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Value-type of parameter " + it->name + " not implemented!");
		}
		return;
	}

	void PILISModel::updateMembers_()
	{
		DoubleReal pseudo_counts = (DoubleReal)param_.getValue("pseudo_counts");
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


