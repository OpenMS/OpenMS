// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#define TRAINING_DEBUG
#undef  TRAINING_DEBUG

#define SIM_DEBUG
#undef  SIM_DEBUG

#define ALIGMENT_DEBUG
#undef  ALIGMENT_DEBUG

// TODO parameter for upper and lower m/z boundaries
// TODO rewrite of neutral losses (depend on which charges?), do it like Pre?
// TODO charge directed/charge remote more fading, any more needed? Or ok?
// TODO position dependance of pathways (model explicitly?); hm, values that makes sense only for charge directed

using namespace std;

namespace OpenMS 
{

	ResidueDB PILISModel::res_db_;

	PILISModel::PILISModel()
		: DefaultParamHandler("PILISModel"),
			valid_(false)
	{	
		defaults_.setValue("upper_mz", 2000.0);
		defaults_.setValue("lower_mz", 200.0);
		defaults_.setValue("charge_remote_threshold", 0.2);
		defaults_.setValue("charge_remote_factor", 1.0);
		defaults_.setValue("charge_directed_threshold", 0.3);
		defaults_.setValue("charge_directed_factor", 1.0);
		defaults_.setValue("side_chain_intensity_threshold", 0.0);
		defaults_.setValue("side_chain_factor", 1.0);
		defaults_.setValue("model_depth", 4);
		defaults_.setValue("visible_model_depth", 30);
		defaults_.setValue("precursor_error", 3.0);
		
		//param_ = default_;
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
			hmms_losses_(model.hmms_losses_),
			prot_dist_(model.prot_dist_),
			tsg_(model.tsg_),
			name_to_enum_(model.name_to_enum_),
			enum_to_name_(model.enum_to_name_),
			valid_(model.valid_)
	{
		initModels_();
	}

	PILISModel& PILISModel::operator = (const PILISModel& model)
	{
		if (this != &model)
		{
			DefaultParamHandler::operator=(model);
			hmm_ = model.hmm_;
	    hmm_precursor_ = model.hmm_precursor_;
  	  hmms_losses_ = model.hmms_losses_;
	    prot_dist_ = model.prot_dist_;
	    tsg_ = model.tsg_;
	    name_to_enum_ = model.name_to_enum_;
	    enum_to_name_ = model.enum_to_name_;
	    valid_ = model.valid_;
		}
		return *this;
	}
	
	void PILISModel::readFromFile(const String& filename)
	{

		// read the model
		if (!File::exists(filename))
 	  {
 	   	throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    if (!File::readable(filename))
    {
     	throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    if (File::empty(filename))
    {
     	throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

		TextFile file;
		file.load(filename, true);

		TextFile::Iterator it_begin(file.begin()), it_end(file.begin());
		it_begin = file.search(it_begin, "BASE_MODEL_BEGIN");
		it_end = file.search(it_begin, "BASE_MODEL_END");
		parseBaseModel_(++it_begin, it_end);

		// seek to next interval
		it_begin = file.search(it_end, "PRECURSOR_MODEL_BEGIN");
		it_end = file.search(it_begin, "PRECURSOR_MODEL_END");
		parseHMMLightModel_(++it_begin, it_end, hmm_precursor_);

		// loss models
		it_begin = file.search(it_end, "BION_LOSS_MODEL_BEGIN");
		it_end = file.search(it_begin, "BION_LOSS_MODEL_END");
		parseHMMLightModel_(++it_begin, it_end, hmms_losses_[Residue::BIon]);
		
		// y-ion loss model
		it_begin = file.search(it_end, "YION_LOSS_MODEL_BEGIN");
		it_end = file.search(it_begin, "YION_LOSS_MODEL_END");
		parseHMMLightModel_(++it_begin, it_end, hmms_losses_[Residue::YIon]);
		
		valid_ = true;
		return;
	}

	void PILISModel::writetoYGFFile(const String& filename)
	{
		hmm_.writetoYGFFile(filename);
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

		out << "PRECURSOR_MODEL_BEGIN" << endl;
		hmm_precursor_.write(out);
		out << "PRECURSOR_MODEL_END" << endl;

		out << "BION_LOSS_MODEL_BEGIN" << endl;
		hmms_losses_[Residue::BIon].write(out);
		out << "BION_LOSS_MODEL_END" << endl;

		out << "YION_LOSS_MODEL_BEGIN" << endl;
		hmms_losses_[Residue::YIon].write(out);
		out << "YION_LOSS_MODEL_END" << endl;
		
		return;
	}

	void PILISModel::train(const PeakSpectrum& in_spec, const AASequence& peptide, UnsignedInt charge)
	{
		if (!valid_)
		{
			cerr << "PILISModel: cannot train, initialize model from file first, e.g. data/PILIS/PILIS_model_default.dat" << endl;
			return;
		}

		PeakSpectrum train_spec = in_spec;
		train_spec.getContainer().sortByPosition();
		
		#ifdef TRAINING_DEBUG
		cout << "peptide: " << peptide  << "(z=" << charge << ")" << endl;
		#endif
		double peptide_weight((peptide.getMonoWeight() + double(charge)) / double(charge));
		//double pre_int1(0), pre_int_H2O_1(0), pre_int_NH3_1(0), pre_int2(0), pre_int_H2O_2(0), pre_int_NH3_2(0);
		PrecursorPeaks_ pre_ints_1, pre_ints_2;
		getPrecursorIntensitiesFromSpectrum_(train_spec, pre_ints_1, peptide_weight, 1);
		if (charge > 1)
		{
			getPrecursorIntensitiesFromSpectrum_(train_spec, pre_ints_2, peptide_weight, 2);
		}
		// get the ions intensities, y and b ions and losses H2O, NH3 respectively
		
		IonPeaks_ ints_1, ints_2;
		//for (Size z = 1; z <= charge; ++z)
		//{
			// TODO
			double sum1 = getIntensitiesFromSpectrum_(train_spec, ints_1, peptide, 1);
			
			double sum2(0);
			
			if (charge > 1)
			{
				sum2 = getIntensitiesFromSpectrum_(train_spec, ints_2, peptide, 2);
			}
		//}

		// normalize the intensities
		//cerr << "pre_ints="	<< pre_ints_1.pre << " " <<  pre_ints_2.pre << " " << pre_ints_1.pre_NH3 << " " <<  pre_ints_2.pre_NH3 << " " <<  pre_ints_1.pre_H2O << " " <<  pre_ints_2.pre_H2O << endl;
		double sum(0);
		sum += sum1;
		sum += pre_ints_1.pre + pre_ints_1.pre_NH3 + pre_ints_1.pre_H2O;
		if (charge > 1)
		{
			sum += sum2 + pre_ints_2.pre + pre_ints_2.pre_NH3 + pre_ints_2.pre_H2O;
		}
		pre_ints_1.pre /= sum;
		pre_ints_1.pre_NH3 /= sum;
		pre_ints_1.pre_H2O /= sum;

		if (charge > 1)
		{
			pre_ints_2.pre  /= sum;
			pre_ints_2.pre_NH3 /= sum;
			pre_ints_2.pre_H2O /= sum;
		}

		for (HashMap<IonType_, vector<double> >::Iterator it1 = ints_1.ints.begin(); it1 != ints_1.ints.end(); ++it1)
		{
			for (vector<double>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				*it2 /= sum;
			}
		}

		if (charge > 1)
		{
			for (HashMap<IonType_, vector<double> >::Iterator it1 = ints_2.ints.begin(); it1 != ints_2.ints.end(); ++it1)
			{
				for (vector<double>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					*it2 /= sum;
				}
			}
		}
	
		//cerr << "sum=" << sum << endl;
		
		if (sum == 0)
		{
			// does not make sense to proceed here
			cerr << "warning: no peaks found which match the given sequence" << endl;
			return;
		}

		HashMap<Size, double> bb_charge_full;
		HashMap<Size, double> sc_charge_full;
		prot_dist_.getProtonDistribution(bb_charge_full, sc_charge_full, peptide, charge, Residue::YIon);
		prot_dist_.setPeptideProtonDistribution(bb_charge_full, sc_charge_full);
	
		vector<double> bb_init, sc_init, cr_init;
		bool is_charge_remote = getInitialTransitionProbabilities_(bb_init, cr_init, sc_init, bb_charge_full, sc_charge_full, peptide);

		// clear the main Hidden Markov Model
		hmm_.clearInitialTransitionProbabilities();
		hmm_.clearTrainingEmissionProbabilities();
		
		double charge_sum(0);

		vector<AASequence> prefixes, suffixes;
		
		// for each site: 1. set proton distribution, 2. initial training intensities, 3. train the model
		for (Size i = 0; i != peptide.size() - 1; ++i)
		{
			//cerr << i << " " << peptide << endl;
			String pos_name, y_name1, b_name1, a_name1, y_name, b_name;
			String y_name2, b_name2, a_name2, prefix_size(i + 1), suffix_size(peptide.size() - 1 - i);

			Size suffix_pos(peptide.size() - i  - 1);
			
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
						
			//cerr << "pos_name=" << pos_name << ", b_name=" << b_name1 << ", y_name=" << y_name1 << endl;
	
			AASequence prefix(peptide.getPrefix(i + 1)), suffix(peptide.getSuffix(peptide.size() - i - 1));
			String aa1(peptide[i]->getOneLetterCode()), aa2(peptide[i + 1]->getOneLetterCode());
			
			// calc PAs and get b/y ratios for bxyz pathway
			double bint1(0), yint1(0), bint2(0), yint2(0), b_sc_int1(0), b_sc_int2(0), y_sc_int1(0), y_sc_int2(0);
			prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, bint1, yint1, bint2, yint2, ProtonDistributionModel::ChargeDirected);

			hmm_.setInitialTransitionProbability("BB"+pos_name, bb_init[i]);
			// TODO disable this if necessary
			hmm_.setInitialTransitionProbability(aa1+aa2+"bxyz"+pos_name, bb_init[i]);
			hmm_.setInitialTransitionProbability(aa1+aa2+"axyz"+pos_name, bb_init[i]);

			//cerr << "ChargeStats: BB=" << BB_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;
			
			//charge_sum += BB_charges[i]; 

			prefixes.push_back(prefix);
			suffixes.push_back(suffix);
			
			double b_cr_int1(0), b_cr_int2(0), y_cr_int1(0), y_cr_int2(0);
			if ((aa1 == "D" || aa1 == "E") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("CR"+pos_name, cr_init[i]);
				charge_sum += cr_init[i];

				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,	Residue::BIon, b_cr_int1, y_cr_int1, b_cr_int2, y_cr_int2, ProtonDistributionModel::ChargeRemote);

				//cerr << "ChargeStats: CR=" << CR_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;
			}

			if (aa1 == "K" || aa1 == "H" || aa1 == "R")
			{
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, b_sc_int1, y_sc_int1, b_sc_int2, y_sc_int2, ProtonDistributionModel::SideChain);
				hmm_.setInitialTransitionProbability("SC"+pos_name, sc_init[i]);

				//cerr << "ChargeStats: SC=" << SC_charges[i] << ", " << peptide.getPrefix(i+1) << "-" << peptide.getSuffix(peptide.size() - 1 - i) << ", " << charge << endl;
				
				//charge_sum += SC_charges[i];
			}

			double y_sum1 = ints_1.ints[YIon][suffix_pos - 1] + ints_1.ints[YIon_H2O][suffix_pos - 1] + ints_1.ints[YIon_NH3][suffix_pos - 1];
			double b_sum1 = ints_1.ints[BIon][i] + ints_1.ints[BIon_H2O][i] + ints_1.ints[BIon_NH3][i]; 
		
			double y_sum2(0), b_sum2(0);
			if (charge > 1)
			{
				y_sum2 = ints_2.ints[YIon][suffix_pos - 1] + ints_2.ints[YIon_H2O][suffix_pos - 1] + ints_2.ints[YIon_NH3][suffix_pos - 1];
				b_sum2 = ints_2.ints[BIon][i] + ints_2.ints[BIon_H2O][i] + ints_2.ints[BIon_NH3][i];
			}

			//cerr << prefix << " - " << suffix << " " << b_sum1 << " " << y_sum1 << " " << i + 1 << " " << suffix_pos << endl;
			
			hmm_.setTrainingEmissionProbability(b_name1, b_sum1);
			hmm_.setTrainingEmissionProbability(a_name1, ints_1.ints[AIon][i]);
			hmm_.setTrainingEmissionProbability(y_name1, y_sum1);

			if (charge > 1)
			{
      	hmm_.setTrainingEmissionProbability(b_name2, b_sum2);
      	hmm_.setTrainingEmissionProbability(y_name2, y_sum2);
			}

			hmm_.setTrainingEmissionProbability("AA"+pos_name, b_sum1 + b_sum2 + y_sum1 + y_sum2);

			hmm_.setTrainingEmissionProbability("end"+pos_name, 0.5/(double(peptide.size()-1)));

			//cerr << aa1 << aa2 << ": " << b_name1 << "=" << b_ints1[i] << ", " << y_name1 << "=" << y_ints1[y_ints1.size()-i-1] << 
			//		", y-H2O+ = " << y_H2O_ints1[y_H2O_ints1.size() - i - 1] << endl;
			//cerr << aa1 << aa2 << ": " << b_name2 << "=" << b_ints2[i] << ", " << y_name2 << "=" << y_ints2[y_ints2.size()-i-1] << 
			//		", y-H2O++ = " << y_H2O_ints2[y_H2O_ints2.size() - i - 1] << endl;
			//cerr << aa1 << aa2 << ": " << pos_name << " " << BB_charges[i] << " " << SC_charges[i] << endl;
	
			// now enable the states
			hmm_.enableTransition("BB"+pos_name, "AA"+pos_name);
			hmm_.enableTransition("BB"+pos_name, "end"+pos_name);

			hmm_.enableTransition(aa1+aa2+"bxyz"+pos_name, "bxyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"bxyz"+pos_name, "end"+pos_name);

			hmm_.enableTransition(aa1+aa2+"axyz"+pos_name, "axyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"axyz"+pos_name, "end"+pos_name);
			
			hmm_.setTransitionProbability("bxyz"+pos_name, b_name1, bint1);
			hmm_.setTransitionProbability("bxyz"+pos_name, y_name1, yint1);

			hmm_.setTransitionProbability("axyz"+pos_name, a_name1, 1); // simple the one and only way to go

			if (charge > 1)
			{
      	hmm_.setTransitionProbability("bxyz"+pos_name, b_name2, bint2);
      	hmm_.setTransitionProbability("bxyz"+pos_name, y_name2, yint2);
			}

      //hmm_.setTransitionProbability("axyz"+pos_name, a_name2, aint2);
		

			if (aa1 == "D" && is_charge_remote)
			{
				hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);

				hmm_.setTrainingEmissionProbability("A"+pos_name, b_sum1 * b_cr_int1 + b_sum2 * b_cr_int2 + y_sum1 * y_cr_int1 + y_sum2 * y_cr_int2);

				hmm_.setTransitionProbability("D"+pos_name, b_name1, b_cr_int1);
				hmm_.setTransitionProbability("D"+pos_name, y_name1, y_cr_int1);
				hmm_.setTransitionProbability("D"+pos_name, b_name2, b_cr_int2);
				hmm_.setTransitionProbability("D"+pos_name, y_name2, y_cr_int2);
				hmm_.setInitialTransitionProbability(aa2+"D"+pos_name, cr_init[i]);
				hmm_.enableTransition(aa2+"D"+pos_name, "D"+pos_name);
				hmm_.enableTransition(aa2+"D"+pos_name, "end"+pos_name);
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
				hmm_.setInitialTransitionProbability(aa2+"E"+pos_name, cr_init[i]);
        hmm_.enableTransition(aa2+"E"+pos_name, "E"+pos_name);
				hmm_.enableTransition(aa2+"E"+pos_name, "end"+pos_name);
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
				hmm_.setInitialTransitionProbability(aa2+"K"+pos_name, sc_init[i]);
        hmm_.enableTransition(aa2+"K"+pos_name, "K"+pos_name);
				hmm_.enableTransition(aa2+"K"+pos_name, "end"+pos_name);
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
				hmm_.setInitialTransitionProbability(aa2+"H"+pos_name, sc_init[i]);
        hmm_.enableTransition(aa2+"H"+pos_name, "H"+pos_name);
				hmm_.enableTransition(aa2+"H"+pos_name, "end"+pos_name);
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
				hmm_.setInitialTransitionProbability(aa2+"RSC"+pos_name, sc_init[i]);
        hmm_.enableTransition(aa2+"RSC"+pos_name, "R"+pos_name);
				hmm_.enableTransition(aa2+"RSC"+pos_name, "end"+pos_name);
      } 

			// losses
			if (ints_1.ints[YIon][suffix_pos - 1] > 0.01)
			{
				HashMap<NeutralLossType_, double> y_intensities;
				y_intensities[LOSS_TYPE_H2O] = ints_1.ints[YIon_H2O][suffix_pos - 1];
				y_intensities[LOSS_TYPE_NH3] = ints_1.ints[YIon_NH3][suffix_pos - 1];
				double ion_intensity = ints_1.ints[YIon][suffix_pos - 1];

				if (charge > 1)
				{
					y_intensities[LOSS_TYPE_H2O] += ints_2.ints[YIon_H2O][suffix_pos - 1];
					y_intensities[LOSS_TYPE_NH3] += ints_2.ints[YIon_NH3][suffix_pos - 1];
					ion_intensity += ints_2.ints[YIon][suffix_pos - 1];
				}
				//cerr << "YLOSS: " << y_sum1 + y_sum2 << " " << y_intensities[LOSS_TYPE_H2O] << " " << y_intensities[LOSS_TYPE_NH3] << " " << ion_intensity << endl;
				trainNeutralLossesFromIon_(y_sum1 + y_sum2, y_intensities, Residue::YIon, ion_intensity, suffix);
			}
			

			if (ints_1.ints[BIon][i] > 0.01)
			{
      	HashMap<NeutralLossType_, double> b_intensities;
      	b_intensities[LOSS_TYPE_H2O] = ints_1.ints[BIon_H2O][i];
      	b_intensities[LOSS_TYPE_NH3] = ints_1.ints[BIon_NH3][i];
      	double ion_intensity = ints_1.ints[BIon][i];

				if (charge > 1)
				{
					b_intensities[LOSS_TYPE_H2O] += ints_2.ints[BIon_H2O][i];
					b_intensities[LOSS_TYPE_NH3] += ints_2.ints[BIon_NH3][i];
					ion_intensity += ints_2.ints[BIon][i];
				}
				//cerr << "BLOSS: " << b_sum1 + b_sum2 << " " << b_intensities[LOSS_TYPE_H2O] << " " << b_intensities[LOSS_TYPE_NH3] << " " << ion_intensity << endl;
				trainNeutralLossesFromIon_(b_sum1 + b_sum2, b_intensities, Residue::BIon, ion_intensity, prefix);
			}

		}

		//cerr << "PrecursorStats: " << pre_int1 << " " << pre_int2 << " " << pre_int_H2O_1 << " " << pre_int_H2O_2 << " " << pre_int_NH3_1 << " " << pre_int_NH3_2 << endl;

		// precursor handling
		if (is_charge_remote || peptide[0]->getOneLetterCode() == "Q")
		{
			double bb_avg(0);
			for (Size i = 1; i != bb_charge_full.size(); ++i)
			{
				bb_avg += bb_charge_full[i];
			}
			bb_avg /= (peptide.size() - 1);
			if (charge == 1)
			{
				trainPrecursorIons_(bb_avg, pre_ints_1.pre, pre_ints_1.pre_NH3, pre_ints_1.pre_H2O, peptide);
			}
			else
			{
				trainPrecursorIons_(bb_avg, pre_ints_2.pre, pre_ints_2.pre_NH3, pre_ints_2.pre_H2O, peptide);
			}
		}
						
		// now train the model with the data set
		hmm_.train();

		//stringstream ss;
		//ss << peptide;
		//hmm_.writetoYGFFile(String("stats/model_graph_train_"+ss.str()+"_"+String(charge)+".graphml").c_str());

		hmm_.disableTransitions();

		return;
	}

	void PILISModel::getPrecursorIntensitiesFromSpectrum_(const PeakSpectrum& train_spec, PrecursorPeaks_& peak_ints, double peptide_weight, UnsignedInt charge)
	{
		static double H2O_weight = Formulas::H2O.getMonoWeight();
		static double NH3_weight = Formulas::NH3.getMonoWeight();
		double pre_error = (double)param_.getValue("precursor_error");
		peak_ints.pre = 0;
		peak_ints.pre_H2O = 0;
		peak_ints.pre_NH3 = 0;
  	for (PeakSpectrum::ConstIterator it = train_spec.begin(); it != train_spec.end(); ++it)
    {
	    if (fabs(it->getPosition()[0] - peptide_weight) < pre_error)
	    {
	      peak_ints.pre += it->getIntensity();
  	  }

      if (fabs(it->getPosition()[0] - (peptide_weight - H2O_weight / double(charge))) < pre_error)
      {
      	peak_ints.pre_H2O += it->getIntensity();
      }

      if (fabs(it->getPosition()[0] - (peptide_weight - NH3_weight / double(charge))) < pre_error)
      {
       	peak_ints.pre_NH3 += it->getIntensity();
      }
		}
	}
	
	double PILISModel::getIntensitiesFromSpectrum_(const PeakSpectrum& train_spec, IonPeaks_& ion_ints, const AASequence& peptide, UnsignedInt z)
	{
		double sum(0);
    PeakSpectrum y_theo_spec, y_H2O_theo_spec, y_NH3_theo_spec;
		TheoreticalSpectrumGenerator tsg;
    tsg_.addPeaks(y_theo_spec, peptide, Residue::YIon, z);
    
		DPeak<1> p;
		p.setIntensity(1);
    for (PeakSpectrum::ConstIterator it = y_theo_spec.begin(); it != y_theo_spec.end(); ++it)
    {
      p.setPosition(it->getPosition() - Formulas::H2O.getMonoWeight() / (double)z);
      y_H2O_theo_spec.getContainer().push_back(p);

      p.setPosition(it->getPosition() - Formulas::NH3.getMonoWeight() / (double)z);
      y_NH3_theo_spec.getContainer().push_back(p);
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

    PeakSpectrum a_theo_spec;
    tsg.addPeaks(a_theo_spec, peptide, Residue::AIon, z);
			
		vector<double> a_ints(peptide.size() - 1, 0.0);
		sum += getIntensitiesFromComparison_(train_spec, a_theo_spec, a_ints);
		ion_ints.ints[AIon] = a_ints;

    PeakSpectrum b_theo_spec, b_H2O_theo_spec, b_NH3_theo_spec;
    tsg.addPeaks(b_theo_spec, peptide, Residue::BIon, z);
    for (PeakSpectrum::ConstIterator it = b_theo_spec.begin(); it != b_theo_spec.end(); ++it)
    {
      p.setPosition(it->getPosition() - Formulas::H2O.getMonoWeight() / (double)z);
      b_H2O_theo_spec.getContainer().push_back(p);
      p.setPosition(it->getPosition() - Formulas::NH3.getMonoWeight() / (double)z);
      b_NH3_theo_spec.getContainer().push_back(p);
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

		return sum;
	}

	double PILISModel::getIntensitiesFromComparison_(const PeakSpectrum& train_spec, const PeakSpectrum& theo_spec, vector<double>& intensities)
	{
		double sum(0);
		vector<pair<Size, Size> > peak_map;
		spectra_aligner_.getSpectrumAlignment(peak_map, train_spec, theo_spec);

    for (vector<pair<Size, Size> >::const_iterator it = peak_map.begin(); it != peak_map.end(); ++it)
    {
			intensities[it->second] = train_spec.getContainer()[it->first].getIntensity();
			sum += train_spec.getContainer()[it->first].getIntensity();
    }
		return sum;
	}
	
	void PILISModel::trainPrecursorIons_(double initial_probability, double intensity, double intensity_NH3, double intensity_H2O, const AASequence& peptide)
	{
		hmm_precursor_.clearInitialTransitionProbabilities();
		hmm_precursor_.clearTrainingEmissionProbabilities();
		hmm_precursor_.setTrainingEmissionProbability(PRE_MH, intensity);
		hmm_precursor_.setTrainingEmissionProbability(PRE_MH_H2O, intensity_H2O);
		hmm_precursor_.setTrainingEmissionProbability(PRE_MH_NH3, intensity_NH3);
		hmm_precursor_.setInitialTransitionProbability(PRE_ION, initial_probability);
		hmm_precursor_.setTrainingEmissionProbability(PRE_END, intensity + intensity_H2O + intensity_NH3);

		enablePrecursorIonStates_(peptide);

		// train the hmm on the above enabled states
		hmm_precursor_.train();

		hmm_precursor_.disableTransitions();
	}

	void PILISModel::trainNeutralLossesFromIon_(double initial_probability, 
																								const HashMap<NeutralLossType_, double>& intensities, 
																								Residue::ResidueType ion_type, 
																								double ion_intensity, 
																								const AASequence& ion)
	{
		//cerr << "LOSS: " << initial_probability << " " << intensities[LOSS_TYPE_H2O] << " " << intensities[LOSS_TYPE_NH3] << " " << ion_intensity << endl;
		// H2O loss
		HiddenMarkovModelLight* hmm = &hmms_losses_[ion_type];
		hmm->clearInitialTransitionProbabilities();
		hmm->clearTrainingEmissionProbabilities();
		if (ion_type == Residue::BIon)
		{
			hmm->setInitialTransitionProbability(B_ION, initial_probability);
			if (intensities.has(PILISModel::LOSS_TYPE_H2O))
			{
				hmm->setTrainingEmissionProbability(B_H2O, intensities[PILISModel::LOSS_TYPE_H2O]);
			}
			
			if (intensities.has(PILISModel::LOSS_TYPE_NH3))
			{
				hmm->setTrainingEmissionProbability(B_NH3, intensities[PILISModel::LOSS_TYPE_NH3]);
			}
			hmm->setTrainingEmissionProbability(B_LOSS_END, ion_intensity);
			enableNeutralLossStates_(ion_type, ion);

			hmm->train();
			hmm->disableTransitions();

			return;
		}

		if (ion_type == Residue::YIon)
		{
			hmm->setInitialTransitionProbability(Y_ION, initial_probability);
      if (intensities.has(PILISModel::LOSS_TYPE_H2O))
      {
        hmm->setTrainingEmissionProbability(Y_H2O, intensities[PILISModel::LOSS_TYPE_H2O]);
      }

      if (intensities.has(PILISModel::LOSS_TYPE_NH3))
      {
        hmm->setTrainingEmissionProbability(Y_NH3, intensities[PILISModel::LOSS_TYPE_NH3]);
      }
      hmm->setTrainingEmissionProbability(Y_LOSS_END, ion_intensity);

			enableNeutralLossStates_(ion_type, ion);

			hmm->train();
			hmm->disableTransitions();
			
      return;
		}

		return;
	}

  void PILISModel::enablePrecursorIonStates_(const AASequence& peptide)
	{
   	if (peptide.has("S"))
    {
      hmm_precursor_.enableTransition(PRE_ION, PRE_H2O_S);
      hmm_precursor_.enableTransition(PRE_H2O_S, PRE_MH_H2O);
      hmm_precursor_.enableTransition(PRE_H2O_S, PRE_END);
      hmm_precursor_.enableTransition(PRE_H2O_S, PRE_MH);
		}

    if (peptide.has("T"))
    {
			hmm_precursor_.enableTransition(PRE_ION, PRE_H2O_T);
      hmm_precursor_.enableTransition(PRE_H2O_T, PRE_MH_H2O);
      hmm_precursor_.enableTransition(PRE_H2O_T, PRE_END);
      hmm_precursor_.enableTransition(PRE_H2O_T, PRE_MH);
    }

    if (peptide.has("E"))
    {
			hmm_precursor_.enableTransition(PRE_ION, PRE_H2O_E);
			hmm_precursor_.enableTransition(PRE_H2O_E, PRE_MH_H2O);
			hmm_precursor_.enableTransition(PRE_H2O_E, PRE_END);
			hmm_precursor_.enableTransition(PRE_H2O_E, PRE_MH);
    }

    if (peptide.has("D"))
    {
			hmm_precursor_.enableTransition(PRE_ION, PRE_H2O_D);
      hmm_precursor_.enableTransition(PRE_H2O_D, PRE_MH_H2O);
      hmm_precursor_.enableTransition(PRE_H2O_D, PRE_END);
      hmm_precursor_.enableTransition(PRE_H2O_D, PRE_MH);
    }

    // NH3 loss pathways
    if (peptide.has("K"))
    {
      hmm_precursor_.enableTransition(PRE_ION, PRE_NH3_K);
      hmm_precursor_.enableTransition(PRE_NH3_K, PRE_MH_NH3);
      hmm_precursor_.enableTransition(PRE_NH3_K, PRE_END);
      hmm_precursor_.enableTransition(PRE_NH3_K, PRE_MH);
    }

    if (peptide.has("R"))
    {
      hmm_precursor_.enableTransition(PRE_ION, PRE_NH3_R);
      hmm_precursor_.enableTransition(PRE_NH3_R, PRE_MH_NH3);
      hmm_precursor_.enableTransition(PRE_NH3_R, PRE_END);
      hmm_precursor_.enableTransition(PRE_NH3_R, PRE_MH);
    }

    // pyroglutamic acid formation
    if (peptide[0]->getOneLetterCode() == "Q")
    {
      hmm_precursor_.enableTransition(PRE_ION, PRE_H2O_Q1);
      hmm_precursor_.enableTransition(PRE_H2O_Q1, PRE_MH_H2O);
      hmm_precursor_.enableTransition(PRE_H2O_Q1, PRE_END);
      hmm_precursor_.enableTransition(PRE_H2O_Q1, PRE_MH);
    }

    // common loss from C-terminus
    hmm_precursor_.enableTransition(PRE_ION, PRE_H2O_CTERM);
    hmm_precursor_.enableTransition(PRE_H2O_CTERM, PRE_MH_H2O);
    hmm_precursor_.enableTransition(PRE_H2O_CTERM, PRE_END);
    hmm_precursor_.enableTransition(PRE_H2O_CTERM, PRE_MH);

	}

  void PILISModel::enableNeutralLossStates_(Residue::ResidueType ion_type, const AASequence& ion)
	{
		
		if (ion_type == Residue::BIon)
		{
			HiddenMarkovModelLight* hmm = &hmms_losses_[ion_type];
			if (ion.has("S"))
      {
        hmm->enableTransition(B_ION, B_H2O_S);
        hmm->enableTransition(B_H2O_S, B_H2O);
        hmm->enableTransition(B_H2O_S, B_LOSS_END);
      }

      if (ion.has("T"))
      {
        hmm->enableTransition(B_ION, B_H2O_T);
        hmm->enableTransition(B_H2O_T, B_H2O);
        hmm->enableTransition(B_H2O_T, B_LOSS_END);
      }

      if (ion.has("E"))
      {
        hmm->enableTransition(B_ION, B_H2O_E);
        hmm->enableTransition(B_H2O_E, B_H2O);
        hmm->enableTransition(B_H2O_E, B_LOSS_END);
      }

      if (ion.has("D"))
      {
        hmm->enableTransition(B_ION, B_H2O_D);
        hmm->enableTransition(B_H2O_D, B_H2O);
        hmm->enableTransition(B_H2O_D, B_LOSS_END);
      }

      if (ion.has("R"))
      {
        hmm->enableTransition(B_ION, B_NH3_R);
        hmm->enableTransition(B_NH3_R, B_NH3);
        hmm->enableTransition(B_NH3_R, B_LOSS_END);
      }

      if (ion.has("K"))
      {
        hmm->enableTransition(B_ION, B_NH3_K);
        hmm->enableTransition(B_NH3_K, B_NH3);
        hmm->enableTransition(B_NH3_K, B_LOSS_END);
      }
      return;
		}

   	if (ion_type == Residue::YIon)
    {
			HiddenMarkovModelLight* hmm = &hmms_losses_[ion_type];
      if (ion.has("S"))
      {
        hmm->enableTransition(Y_ION, Y_H2O_S);
        hmm->enableTransition(Y_H2O_S, Y_H2O);
        hmm->enableTransition(Y_H2O_S, Y_LOSS_END);
      }
			
      if (ion.has("T"))
      {
        hmm->enableTransition(Y_ION, Y_H2O_T);
        hmm->enableTransition(Y_H2O_T, Y_H2O);
        hmm->enableTransition(Y_H2O_T, Y_LOSS_END);
      }

      if (ion.has("E"))
      {
        hmm->enableTransition(Y_ION, Y_H2O_E);
        hmm->enableTransition(Y_H2O_E, Y_H2O);
        hmm->enableTransition(Y_H2O_E, Y_LOSS_END);
      }

      if (ion.has("D"))
      {
        hmm->enableTransition(Y_ION, Y_H2O_D);
        hmm->enableTransition(Y_H2O_D, Y_H2O);
        hmm->enableTransition(Y_H2O_D, Y_LOSS_END);
      }

      if (ion[0]->getOneLetterCode() == "Q")
      {
        hmm->enableTransition(Y_ION, Y_H2O_Q1);
        hmm->enableTransition(Y_H2O_Q1, Y_H2O);
        hmm->enableTransition(Y_H2O_Q1, Y_LOSS_END);
      }

      hmm->enableTransition(Y_ION, Y_H2O_CTERM);
      hmm->enableTransition(Y_H2O_CTERM, Y_H2O);
      hmm->enableTransition(Y_H2O_CTERM, Y_LOSS_END);

      if (ion.has("R"))
      {
        hmm->enableTransition(Y_ION, Y_NH3_R);
        hmm->enableTransition(Y_NH3_R, Y_NH3);
        hmm->enableTransition(Y_NH3_R, Y_LOSS_END);
      }

      if (ion.has("K"))
      {
        hmm->enableTransition(Y_ION, Y_NH3_K);
        hmm->enableTransition(Y_NH3_K, Y_NH3);
        hmm->enableTransition(Y_NH3_K, Y_LOSS_END);
      }
      return;
    }

		return;
	}
	

	void PILISModel::evaluate()
	{
		hmm_.evaluate();
		hmm_.estimateUntrainedTransitions();

		hmm_precursor_.evaluate();
		//hmm_precursor_.dump();

		for (HashMap<Residue::ResidueType, HiddenMarkovModelLight>::Iterator it = hmms_losses_.begin(); it != hmms_losses_.end(); ++it)
		{
			it->second.evaluate();
			//it->second.dump();
		}
	}

/*
	void PILISModel::getSpectrumAlignment(HashMap<Size, Size>& peak_map, const PeakSpectrum& spec1, const PeakSpectrum& spec2)
	{
		//cerr << "HashMap<Size, Size> HMMModel2::getSpectrumAlignment(const PeakSpectrum& spec1, const PeakSpectrum& spec2)" << endl;
		// spec1 should be the fixed (theo spec)
		// "threaded" algorithm
	
		//const float rel_error(0.004);
		const float rel_error(0.004);
		
		const PeakSpectrum::ContainerType& a1 = spec1.getContainer();
		const PeakSpectrum::ContainerType& a2 = spec2.getContainer();
		
		Size last_j(0);
		for (Size i = 0; i != a1.size(); ++i)
		{
			float pos1 = a1[i].getPosition()[0];
			float error = pos1 * rel_error;
			//float error = 0.3;
			float diff = numeric_limits<float>::max();
			for (Size j = last_j; j != a2.size(); ++j)
			{
				float pos2 = a2[j].getPosition()[0];
				if (abs(pos1 - pos2) < error && abs(pos1 - pos2) < diff)
				{
					diff = abs(pos1-pos2);
					last_j = j;
				}
				
				if (pos2 > pos1 && abs(pos1 - pos2) > error)
				{
					break;
				}
			}
			if (diff != numeric_limits<float>::max())
			{
				peak_map[i] = last_j;
			}
		}
		#ifdef ALIGMENT_DEBUG
		for (HashMap<Size, Size>::ConstIterator it = peak_map.begin(); it != peak_map.end(); ++it)
		{
			cerr << it->first << " " << a1[it->first].getPosition()[0] << " <> " << a2[it->second].getPosition()[0] << " " << a2[it->second].getIntensity() << endl;
		}
		#endif
		//cerr << peak_map.size() << endl;
	}

	*/

	void PILISModel::getSpectrum(PeakSpectrum& spec, const AASequence& peptide, UnsignedInt charge)
	{
		if (!valid_)
		{
			cerr << "PILISModel: cannot train, initialize model from file first, e.g. data/PILIS/PILIS_model_default.dat" << endl;
		}
		
		// calc proton distribution
		HashMap<Size, double> sc_charge_full, bb_charge_full;
		prot_dist_.getProtonDistribution(bb_charge_full, sc_charge_full, peptide, charge, Residue::YIon);
		prot_dist_.setPeptideProtonDistribution(bb_charge_full, sc_charge_full);

		hmm_.clearInitialTransitionProbabilities();
		hmm_.clearTrainingEmissionProbabilities();
	
    // set charges
    vector<double> bb_init, sc_init, cr_init;
		bool is_charge_remote = getInitialTransitionProbabilities_(bb_init, cr_init, sc_init, bb_charge_full, sc_charge_full, peptide);

		double charge_sum(0);

		vector<AASequence> suffixes, prefixes;
		vector<String> suffix_names1, suffix_names2, prefix_names1, prefix_names2, a_names1, a_names2;
		
		// get the paths
		for (Size i = 0; i != peptide.size() - 1; ++i)
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

			String aa1(peptide[i]->getOneLetterCode()), aa2(peptide[i+1]->getOneLetterCode());

			hmm_.setInitialTransitionProbability("BB"+pos_name, bb_init[i]);
			charge_sum += bb_init[i];

			double b_cr_int1(0), b_cr_int2(0), y_cr_int1(0), y_cr_int2(0);
			double bint1(0), yint1(0), bint2(0), yint2(0), aint1(0), aint2(0), ayint1(0), ayint2(0), b_sc_int1(0), b_sc_int2(0), y_sc_int1(0), y_sc_int2(0);
			if ((aa1 == "D" || aa1 == "E") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("CR"+pos_name, cr_init[i]);
				charge_sum += cr_init[i];
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,
																		Residue::BIon, b_cr_int1, y_cr_int1, b_cr_int2, y_cr_int2, ProtonDistributionModel::ChargeRemote);

				// TODO
				y_cr_int1 += 0.5 * y_cr_int2;
				y_cr_int2 *= 0.5;
				b_cr_int1 += 0.5 * b_cr_int2;
				b_cr_int2 *= 0.5;
			}

			if ((aa1 == "K" || aa1 == "H" || aa1 == "R") && is_charge_remote)
			{
				hmm_.setInitialTransitionProbability("SC"+pos_name, sc_init[i]);
				charge_sum += sc_init[i];
				prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge,	Residue::BIon, b_sc_int1, y_sc_int1, b_sc_int2, y_sc_int2, ProtonDistributionModel::SideChain);

				// TODO
				y_sc_int1 += 0.5 * y_sc_int2;
        y_sc_int2 *= 0.5;
        b_sc_int1 += 0.5 * b_sc_int2;
        b_sc_int2 *= 0.5;
			}
		
      prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::BIon, bint1, yint1, bint2, yint2, ProtonDistributionModel::ChargeDirected);

			//cerr << "charge state intensities: " << peptide << " (" << prefix << "-" << suffix << ") " << bint1 << " " << yint1 << " " << bint2 << " " << yint2 << " " << charge << endl;

      prot_dist_.getChargeStateIntensities(peptide, prefix, suffix, charge, Residue::AIon, aint1, ayint1, aint2, ayint2, ProtonDistributionModel::ChargeDirected);	
	
			// TODO correction of the doubly charged ions to singly charged ones, to account for proton loss (and not to need it model explicitly)
			yint1 += 0.5 * yint2;
			yint2 *= 0.5;
			bint1 += 0.5 * bint2;
			bint2 *= 0.5;

      // now enable the states
			hmm_.enableTransition("BB"+pos_name, "AA"+pos_name);
			hmm_.enableTransition("BB"+pos_name, "end"+pos_name);
			
      hmm_.enableTransition("AA"+pos_name, aa1+aa2+"bxyz"+pos_name);
			hmm_.enableTransition("AA"+pos_name, aa1+aa2+"axyz"+pos_name);

      hmm_.enableTransition(aa1+aa2+"bxyz"+pos_name, "bxyz"+pos_name);
			hmm_.enableTransition(aa1+aa2+"axyz"+pos_name, "axyz"+pos_name);

			hmm_.enableTransition(aa1+aa2+"bxyz"+pos_name, "end"+pos_name);
			hmm_.enableTransition(aa1+aa2+"axyz"+pos_name, "end"+pos_name);

      hmm_.setTransitionProbability("bxyz"+pos_name, b_name1, bint1);
      hmm_.setTransitionProbability("bxyz"+pos_name, y_name1, yint1);

			hmm_.setTransitionProbability("bxyz"+pos_name, b_name2, bint2);
			hmm_.setTransitionProbability("bxyz"+pos_name, y_name2, yint2);

			// TODO reactivate a-ions (train only a-ions not in together with y_ions
			hmm_.setTransitionProbability("axyz"+pos_name, a_name1, aint1);
			hmm_.setTransitionProbability("axyz"+pos_name, a_name2, aint2);

			//cerr << "neutral loss charge sums: " << H2O_b_charge_sum << " " << H2O_y_charge_sum << " " << NH3_b_charge_sum << " " << NH3_y_charge_sum << " " << charge_dir_tmp << endl;

      if (aa1 == "D" && is_charge_remote)
      {
        hmm_.setTransitionProbability("D"+pos_name, b_name1, b_cr_int1);
        hmm_.setTransitionProbability("D"+pos_name, y_name1, y_cr_int1);
				hmm_.setTransitionProbability("D"+pos_name, b_name2, b_cr_int2);
				hmm_.setTransitionProbability("D"+pos_name, y_name2, y_cr_int2);
				hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);
				hmm_.enableTransition("A"+pos_name, aa2+"D"+pos_name);
        hmm_.enableTransition(aa2+"D"+pos_name, "D"+pos_name);
				hmm_.enableTransition(aa2+"D"+pos_name, "end"+pos_name);
			}

			if (aa1 == "E" && is_charge_remote)
      {
        hmm_.setTransitionProbability("E"+pos_name, b_name1, b_cr_int1);
        hmm_.setTransitionProbability("E"+pos_name, y_name1, y_cr_int1);
				hmm_.setTransitionProbability("E"+pos_name, b_name2, b_cr_int2);
				hmm_.setTransitionProbability("E"+pos_name, y_name2, y_cr_int2);
        hmm_.enableTransition("CR"+pos_name, "A"+pos_name);
				hmm_.enableTransition("CR"+pos_name, "end"+pos_name);
        hmm_.enableTransition("A"+pos_name, aa2+"E"+pos_name);
        hmm_.enableTransition(aa2+"E"+pos_name, "E"+pos_name);
				hmm_.enableTransition(aa2+"E"+pos_name, "end"+pos_name);
      }

			if (aa1 == "K"/* && is_charge_remote*/)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
        hmm_.setTransitionProbability("K"+pos_name, b_name1, b_sc_int1);
        hmm_.setTransitionProbability("K"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("K"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("K"+pos_name, y_name2, y_sc_int2);
        hmm_.enableTransition("ASC"+pos_name, aa2+"K"+pos_name);
        hmm_.enableTransition(aa2+"K"+pos_name, "K"+pos_name);
				hmm_.enableTransition(aa2+"K"+pos_name, "end"+pos_name);
			}
				
			if (aa1 == "H"/* && is_charge_remote*/)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
        hmm_.setTransitionProbability("H"+pos_name, b_name1, b_sc_int1);
        hmm_.setTransitionProbability("H"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("H"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("H"+pos_name, y_name2, y_sc_int2);
        hmm_.enableTransition("ASC"+pos_name, aa2+"H"+pos_name);
        hmm_.enableTransition(aa2+"H"+pos_name, "H"+pos_name);
				hmm_.enableTransition(aa2+"H"+pos_name, "end"+pos_name);
      }

      if (aa1 == "R"/* && is_charge_remote*/)
      {
        hmm_.enableTransition("SC"+pos_name, "ASC"+pos_name);
				hmm_.enableTransition("SC"+pos_name, "end"+pos_name);
        hmm_.setTransitionProbability("R"+pos_name, b_name1, b_sc_int1);
        hmm_.setTransitionProbability("R"+pos_name, y_name1, y_sc_int1);
				hmm_.setTransitionProbability("R"+pos_name, b_name2, b_sc_int2);
				hmm_.setTransitionProbability("R"+pos_name, y_name2, y_sc_int2);
        hmm_.enableTransition("ASC"+pos_name, aa2+"RSC"+pos_name);
        hmm_.enableTransition(aa2+"RSC"+pos_name, "R"+pos_name);
				hmm_.enableTransition(aa2+"RSC"+pos_name, "end"+pos_name);
      }

		}

    HashMap<HMMState*, double> tmp;
		hmm_.calculateEmissionProbabilities(tmp);

		// clear peaks from last spectrum
		peaks_.clear();

		//stringstream peptide_ss;
		//peptide_ss << peptide;
		//hmm_.writetoYGFFile(String("stats/model_graph_train_"+peptide_ss.str()+"_"+String(charge)+".graphml").c_str());

		// losses
		vector<HashMap<NeutralLossType_, double> > prefix_losses, suffix_losses;
		vector<double> prefix_ints1, prefix_ints2, suffix_ints1, suffix_ints2;
		for (Size i = 0; i != prefixes.size(); ++i)
		{
			HashMap<NeutralLossType_, double> prefix_loss;
			double prefix_int1 = tmp[hmm_.getState(prefix_names1[i])];
			prefix_ints1.push_back(prefix_int1);
			double prefix_int2 = tmp[hmm_.getState(prefix_names2[i])];
			prefix_ints2.push_back(prefix_int2);
			getNeutralLossesFromIon_(prefix_loss, prefix_int1 + prefix_int2, Residue::BIon, prefixes[i]);
			prefix_losses.push_back(prefix_loss);

			HashMap<NeutralLossType_, double> suffix_loss;
			double suffix_int1 = tmp[hmm_.getState(suffix_names1[i])];
			suffix_ints1.push_back(suffix_int1);
			double suffix_int2 = tmp[hmm_.getState(suffix_names2[i])];
			suffix_ints2.push_back(suffix_int2);
			//cerr << "getting loss ints: " << suffix_int1 + suffix_int2 << endl;
			getNeutralLossesFromIon_(suffix_loss, suffix_int1 + suffix_int2, Residue::YIon, suffixes[i]);
			suffix_losses.push_back(suffix_loss);
		}

		hmm_.disableTransitions();

		// read the emission probs and put the peaks into a spec
		IsotopeDistribution id(2);
	
		// register name
		DPeak<1> p;
		p.metaRegistry().registerName("IonName", "Name of the ion");
		
		for (Size i = 0; i != prefixes.size(); ++i)
		{
			// prefix
			double weight = prefixes[i].getMonoWeight(Residue::BIon);
			id.estimateFromPeptideWeight(weight);
			
			//cerr << weight + 1.0 << " " << prefix_ints1[i] << " " << id.getContainer()[0].first << " ";
			//for (vector<pair<UnsignedInt, double> >::const_iterator it = id.getContainer().begin(); it != id.getContainer().end(); ++it)
			//{
			//	cerr << "(" << it->first << "|" << it->second << ") ";
			//}
			//cerr << endl;
			
			// first isotope peak
			addPeaks_(weight, 1, 0.0, prefix_ints1[i], spec, id, "b"+String(i+1));
			if (charge == 2)
			{
				addPeaks_(weight, 2, 0.0, prefix_ints2[i], spec, id, "b"+String(i+1));

				// neutral losses
				// get fractions as the different charge states are treated together 
				double loss_1_fraction = prefix_ints1[i] / (prefix_ints1[i] + prefix_ints2[i]);
				double loss_2_fraction = prefix_ints2[i] / (prefix_ints1[i] + prefix_ints2[i]);
		
				if (prefix_losses[i].has(LOSS_TYPE_H2O))
				{
					addPeaks_(weight, 1, -18.0, prefix_losses[i][LOSS_TYPE_H2O] * loss_1_fraction, spec, id, "b" + String(i+1) + "-H2O");

          // doubly charged
					addPeaks_(weight, 2, -18.0, prefix_losses[i][LOSS_TYPE_H2O] * loss_2_fraction, spec, id, "b" + String(i+1) + "-H2O++");
				}
				if (prefix_losses[i].has(LOSS_TYPE_NH3))
				{
					addPeaks_(weight, 1, -17.0, prefix_losses[i][LOSS_TYPE_NH3] * loss_1_fraction, spec, id, "b" + String(i+1) + "-NH3");
          // doubly charged
					addPeaks_(weight, 2, -17.0, prefix_losses[i][LOSS_TYPE_NH3] * loss_2_fraction, spec, id, "b" + String(i+1) + "-NH3++");
				}
			}
			else
			{
				if (prefix_losses[i].has(LOSS_TYPE_H2O))
        {
					addPeaks_(weight, 1, -18.0, prefix_losses[i][LOSS_TYPE_H2O], spec, id, "b" + String(i+1) + "-H2O");
				}
				if (prefix_losses[i].has(LOSS_TYPE_NH3))
        {
					addPeaks_(weight, 1, -17.0, prefix_losses[i][LOSS_TYPE_NH3], spec, id, "b" + String(i+1) + "-NH3");
				}
			}

			// a-ions
			//weight = prefixes[i].getMonoWeight(Residue::AIon);
			//id.estimateFromPeptideWeight(weight);
			double a_int1 = tmp[hmm_.getState(a_names1[i])];
			addPeaks_(weight, 1, 0.0, a_int1, spec, id, "a" + String(i+1));

			if (charge == 2)
			{
				double a_int2 = tmp[hmm_.getState(a_names2[i])];
				addPeaks_(weight, 2, 0.0, a_int2, spec, id, "a" + String(i+1) + "++");
			}

			// suffix ions
			weight = suffixes[i].getMonoWeight(Residue::YIon);
      id.estimateFromPeptideWeight(weight);
			addPeaks_(weight, 1, 0.0, suffix_ints1[i], spec, id, suffix_names1[i]);
      if (charge == 2)
      {
				addPeaks_(weight, 2, 0.0, suffix_ints2[i], spec, id, suffix_names2[i] + "++");

        // neutral losses
        // get fractions as the different charge states are treated together
        double loss_1_fraction = suffix_ints1[i] / (suffix_ints1[i] + suffix_ints2[i]);
        double loss_2_fraction = suffix_ints2[i] / (suffix_ints1[i] + suffix_ints2[i]);

        if (suffix_losses[i].has(LOSS_TYPE_H2O))
	      {
					addPeaks_(weight, 1, -18.0, suffix_losses[i][LOSS_TYPE_H2O] * loss_1_fraction, spec, id, "y" + String(i + 1) + "-H2O");
          // doubly charged
					addPeaks_(weight, 2, -18.0, suffix_losses[i][LOSS_TYPE_H2O] * loss_2_fraction, spec, id, "y" + String(i + 1) + "-H2O++");
        }

        if (suffix_losses[i].has(LOSS_TYPE_NH3))
        {
					addPeaks_(weight, 1, -17.0, suffix_losses[i][LOSS_TYPE_NH3] * loss_1_fraction, spec, id, "y" + String(i + 1) + "-NH3");
          // doubly charged
					addPeaks_(weight, 2, -17.0, suffix_losses[i][LOSS_TYPE_NH3] * loss_2_fraction, spec, id, "y" + String(i + 1) + "-NH3++");
        }
      }
      else
      {
        if (suffix_losses[i].has(LOSS_TYPE_H2O))
        {
					//cerr << "H2O: " << suffix_losses[i][LOSS_TYPE_H2O] << " " << suffix_ints1[i] << endl;
					//cerr << "NH3: " << suffix_losses[i][LOSS_TYPE_H2O] << endl;
					addPeaks_(weight, 1, -18.0, suffix_losses[i][LOSS_TYPE_H2O], spec, id, "y" + String(i + 1) + "-H2O");
        }
        if (suffix_losses[i].has(LOSS_TYPE_NH3))
        {
					addPeaks_(weight, 1, -17.0, suffix_losses[i][LOSS_TYPE_NH3], spec, id, "y" + String(i + 1) + "-NH3");
        }
      }
		}

		// precursor
		if (is_charge_remote || peptide[0]->getOneLetterCode() == "Q")
		{
			// get bb_avg
			double bb_avg(0);
      for (Size i = 1; i != bb_charge_full.size(); ++i)
      {
        bb_avg += bb_charge_full[i];
      }
      bb_avg /= (peptide.size() - 1);

			HashMap<NeutralLossType_, double> pre_ints;
			// TODO initial transition prob
			getPrecursorIons_(pre_ints, bb_avg, peptide);
	
			double weight = peptide.getMonoWeight();
			id.estimateFromPeptideWeight(weight);
			
			if (pre_ints.has(LOSS_TYPE_H2O))
			{
				addPeaks_(weight, charge, -18.0, pre_ints[LOSS_TYPE_H2O], spec, id, "M-H2O"); 
			}

			if (pre_ints.has(LOSS_TYPE_NH3))
			{
				addPeaks_(weight, charge, -17.0, pre_ints[LOSS_TYPE_NH3], spec, id, "M-NH3");
			}

			if (pre_ints.has(LOSS_TYPE_NONE))
			{
				addPeaks_(weight, charge, 0.0, pre_ints[LOSS_TYPE_NONE], spec, id, "M");
			}
		}
	
		// now build the spectrum with the peaks
		//Peak p;
		for (HashMap<double, vector<Peak> >::ConstIterator it = peaks_.begin(); it != peaks_.end(); ++it)
		{
			if (it->second.size() == 1/* && it->second.begin()->getIntensity() != 0*/)
			{
				spec.getContainer().push_back(*it->second.begin());
			}
			else
			{
				double int_sum(0);
				p = *it->second.begin();
				for (vector<Peak>::const_iterator pit = it->second.begin(); pit != it->second.end(); ++pit)
				{
					int_sum += pit->getIntensity();
					if (pit->getMetaValue("IonName") != "")
					{
						p.setMetaValue("IonName", pit->getMetaValue("IonName"));
					}
				}

				/*if (int_sum != 0)
				{*/
					//p = *it->second.begin();
					p.setIntensity(int_sum);
					//p.setPosition(pit->first);
				//}
				spec.getContainer().push_back(p);
			}
		}

		spec.getContainer().sortByPosition();

		#ifdef SIM_DEBUG
		cerr << "now mapping the intensities to positions" << endl;
		for (HashMap<HMMState*, double>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
		{
			cerr << it->first->getName() << "\t" << it->second << endl;
		}
		#endif


		return;
	}

	void PILISModel::getPrecursorIons_(HashMap<NeutralLossType_, double>& intensities, double initial_probability, const AASequence& precursor)
	{
		hmm_precursor_.setInitialTransitionProbability(PRE_ION, initial_probability);
		
		enablePrecursorIonStates_(precursor);

		HashMap<HMMStateLight*, double> tmp;
		hmm_precursor_.calculateEmissionProbabilities(tmp);

		intensities[PILISModel::LOSS_TYPE_H2O] = tmp[hmm_precursor_.getState(PRE_MH_H2O)];
		intensities[PILISModel::LOSS_TYPE_NH3] = tmp[hmm_precursor_.getState(PRE_MH_NH3)];
		intensities[PILISModel::LOSS_TYPE_NONE] =	tmp[hmm_precursor_.getState(PRE_MH)];

		hmm_precursor_.disableTransitions();
	}

  void PILISModel::getNeutralLossesFromIon_(HashMap<NeutralLossType_, double>& intensities, double initial_probability, Residue::ResidueType ion_type, const AASequence& ion)
	{
		HiddenMarkovModelLight* hmm = &hmms_losses_[ion_type];
		
		if (ion_type == Residue::BIon)
		{
			hmm->setInitialTransitionProbability(B_ION, initial_probability);
			enableNeutralLossStates_(ion_type, ion);

			HashMap<HMMStateLight*, double> tmp;
			hmm->calculateEmissionProbabilities(tmp);
			
			intensities[PILISModel::LOSS_TYPE_H2O] = tmp[hmm->getState(B_H2O)];
			intensities[PILISModel::LOSS_TYPE_NH3] = tmp[hmm->getState(B_NH3)];

			hmm->disableTransitions();
			
			return;
		}

		if (ion_type == Residue::YIon)
    {
      hmm->setInitialTransitionProbability(Y_ION, initial_probability);
      enableNeutralLossStates_(ion_type, ion);

      HashMap<HMMStateLight*, double> tmp;
      hmm->calculateEmissionProbabilities(tmp);

      intensities[PILISModel::LOSS_TYPE_H2O] = tmp[hmm->getState(Y_H2O)];
      intensities[PILISModel::LOSS_TYPE_NH3] = tmp[hmm->getState(Y_NH3)];

			hmm->disableTransitions();
			
      return;
    }
		
		return;
	}

	bool PILISModel::getInitialTransitionProbabilities_(std::vector<double>& bb_init, 
																												std::vector<double>& cr_init, 
																												std::vector<double>& sc_init, 
																												const HashMap<Size, double>& bb_charges, 
																												const HashMap<Size, double>& sc_charges, 
																												const AASequence& peptide)
	{
		bool is_charge_remote(false);

		double charge_dir_tmp(bb_charges[0]);

    double bb_sum(0);
    for (Size i = 0; i != bb_charges.size(); ++i)
    {
      bb_sum += bb_charges[i];
    }

    charge_dir_tmp = bb_sum;
    if (bb_sum < (double)param_.getValue("charge_directed_threshold")/*CHARGE_DIRECTED_THRESHOLD*/)
    {
      charge_dir_tmp = (double)param_.getValue("charge_directed_threshold")/*CHARGE_DIRECTED_THRESHOLD*/;
			is_charge_remote = true;
    }

    double bb_avg = (bb_sum - bb_charges[0]) / (double)(bb_charges.size() - 1);

					
		for (Size i = 0; i != peptide.size() - 1; ++i)
    {
      bb_init.push_back(bb_charges[i+1] * charge_dir_tmp * (double)param_.getValue("charge_directed_factor")/*CHARGE_DIRECTED_FACTOR*/);
      String aa(peptide[i]->getOneLetterCode());
      if (sc_charges.has(i))
      {
        if ((aa == "K" || aa == "R" || aa == "H") /*&& bb_sum < (double)param_.getValue("charge_remote_threshold")*/)
        {
          sc_init.push_back(sc_charges[i] * bb_avg * (double)param_.getValue("side_chain_factor") /*SIDE_CHAIN_FACTOR*/);
        }
        else
        {
          sc_init.push_back(0.0);
        }
        if (bb_sum < (double)param_.getValue("charge_remote_threshold") /*CHARGE_REMOTE_THRESHOLD*/ && (aa == "D" || aa == "E"))
        {
          cr_init.push_back(((1 - bb_sum) * bb_avg /** sc_charge[i]*/) * (double)param_.getValue("charge_remote_factor")/*CHARGE_REMOTE_FACTOR*/);
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
		// C-term bb positional acmino acid might have a proton associated with (is inactive, as now pathway uses it)
		if (sc_charges.has(peptide.size() - 1))
		{
			init_prob_sum += sc_charges[peptide.size() - 1];
		}
		for (Size i = 0; i != bb_init.size(); ++i)
		{
			bb_init[i] /= init_prob_sum;
			sc_init[i] /= init_prob_sum;
			cr_init[i] /= init_prob_sum;
		}
		
		return is_charge_remote;
	}

	void PILISModel::addPeaks_(double mz, int charge, double offset, double intensity, PeakSpectrum& /*spectrum*/, const IsotopeDistribution& id, const String& name)
	{
		static DPeak<1> p;
		Size i = 0;
		for (IsotopeDistribution::ConstIterator it = id.begin(); it != id.end(); ++it, ++i)
		{
			double pos = (mz + i + charge + offset) / (double)charge;
			p.setPosition(pos);
			if (it == id.begin())
			{
				p.setMetaValue("IonName", string(name.c_str()));
			}

			if (pos >= (double)param_.getValue("lower_mz") && pos <= (double)param_.getValue("upper_mz"))
			{
				p.setIntensity(intensity * it->second);
				//spectrum.getContainer().push_back(p);
				peaks_[p.getPosition()[0]].push_back(p);
			}

			if (it == id.begin())
			{
				p.setMetaValue("IonName", string(""));
			}
		}

		return;
	}
	
	/*
	void PILISModel::initModel_()
	{
		cerr << ">> Init model " << endl;

		hmm_.addNewState(new HMMState(String("endcenter"), false));
		hmm_.addNewState(new HMMState(String("end"), false));
		
		hmm_.addNewState(new HMMState(String("BBcenter")));
		hmm_.addNewState(new HMMState(String("AAcenter")));
		hmm_.addNewState(new HMMState(String("CRcenter")));
		hmm_.addNewState(new HMMState(String("Acenter")));
		hmm_.addNewState(new HMMState(String("SCcenter")));
		hmm_.addNewState(new HMMState(String("ASCcenter")));
		
		hmm_.addNewState(new HMMState(String("bxyz")));
		hmm_.addNewState(new HMMState(String("axyz")));
		hmm_.addNewState(new HMMState(String("D")));
		hmm_.addNewState(new HMMState(String("E")));
		
		hmm_.addNewState(new HMMState(String("AABase1")));
		hmm_.addNewState(new HMMState(String("AABase2")));
		
		hmm_.addNewState(new HMMState(String("K")));
		hmm_.addNewState(new HMMState(String("H")));
		hmm_.addNewState(new HMMState(String("R")));
	
		const String residues("ACDEFGHIKMNPQRSTVWY");
	
		// create the residue states states
		for (Size ii = 0; ii != residues.size(); ++ii)
    {
      String first;
      first += residues[ii];
      hmm_.addNewState(new HMMState(first+"D"));
			hmm_.addNewState(new HMMState(first+"E"));

			hmm_.addNewState(new HMMState(first+"K"));
			hmm_.addNewState(new HMMState(first+"H"));
			hmm_.addNewState(new HMMState(first+"RSC"));
			
      for (Size j = 0; j != residues.size(); ++j)
      {
        String second;
        second += residues[j];
        hmm_.addNewState(new HMMState(first+second+"bxyz"));
				hmm_.addNewState(new HMMState(first+second+"axyz"));
      }
    }

		for (Size i = 1; i <= VISIBLE_MODEL_DEPTH; ++i)
		{
			// these states are really created
			// charge states
			hmm_.addNewState(new HMMState("BB"+String(i)));
			hmm_.addNewState(new HMMState("BBk-"+String(i)));
			hmm_.addNewState(new HMMState("CR"+String(i))); 
			hmm_.addNewState(new HMMState("CRk-"+String(i)));

			hmm_.addNewState(new HMMState("SC"+String(i)));
			hmm_.addNewState(new HMMState("SCk-"+String(i)));

			// states for trans mapping 
			hmm_.addNewState(new HMMState("AA"+String(i)));
			hmm_.addNewState(new HMMState("AAk-"+String(i)));

			hmm_.addNewState(new HMMState("A"+String(i)));
			hmm_.addNewState(new HMMState("Ak-"+String(i)));

			hmm_.addNewState(new HMMState("ASC"+String(i)));
			hmm_.addNewState(new HMMState("ASCk-"+String(i)));

			// emitting ion states
			hmm_.addNewState(new HMMState("b"+String(i)+"+", false));
			hmm_.addNewState(new HMMState("bk-"+String(i)+"+", false));
			hmm_.addNewState(new HMMState("y"+String(i)+"+", false));
			hmm_.addNewState(new HMMState("yk-"+String(i)+"+", false));
			hmm_.addNewState(new HMMState("a"+String(i)+"+", false));
			hmm_.addNewState(new HMMState("ak-"+String(i)+"+", false));

			hmm_.addNewState(new HMMState("b"+String(i)+"++", false));
      hmm_.addNewState(new HMMState("bk-"+String(i)+"++", false));
      hmm_.addNewState(new HMMState("y"+String(i)+"++", false));
      hmm_.addNewState(new HMMState("yk-"+String(i)+"++", false));
      hmm_.addNewState(new HMMState("a"+String(i)+"++", false));
      hmm_.addNewState(new HMMState("ak-"+String(i)+"++", false));

			hmm_.addNewState(new HMMState("end"+String(i), false));
			hmm_.addNewState(new HMMState("endk-"+String(i), false));

			// emitting neutral loss states
		
			// post AA collector states
			hmm_.addNewState(new HMMState("bxyz"+String(i)));
			hmm_.addNewState(new HMMState("bxyzk-"+String(i)));
			hmm_.addNewState(new HMMState("axyz"+String(i)));
			hmm_.addNewState(new HMMState("axyzk-"+String(i)));
			hmm_.addNewState(new HMMState("D"+String(i)));
			hmm_.addNewState(new HMMState("Dk-"+String(i)));
			hmm_.addNewState(new HMMState("E"+String(i)));
			hmm_.addNewState(new HMMState("Ek-"+String(i)));

			hmm_.addNewState(new HMMState("K"+String(i)));
			hmm_.addNewState(new HMMState("Kk-"+String(i)));
			hmm_.addNewState(new HMMState("H"+String(i)));
			hmm_.addNewState(new HMMState("Hk-"+String(i)));
			hmm_.addNewState(new HMMState("R"+String(i)));
			hmm_.addNewState(new HMMState("Rk-"+String(i)));

			// map the residue states
			for (Size ii = 0; ii != residues.size(); ++ii)
			{
				String first;
				first += residues[ii];
				hmm_.addNewState(new HMMState(first+"D"+String(i)));
				hmm_.addNewState(new HMMState(first+"Dk-"+String(i)));

        hmm_.addNewState(new HMMState(first+"E"+String(i)));
        hmm_.addNewState(new HMMState(first+"Ek-"+String(i)));

				hmm_.addNewState(new HMMState(first+"K"+String(i)));
				hmm_.addNewState(new HMMState(first+"Kk-"+String(i)));

				hmm_.addNewState(new HMMState(first+"H"+String(i)));
        hmm_.addNewState(new HMMState(first+"Hk-"+String(i)));

				hmm_.addNewState(new HMMState(first+"RSC"+String(i)));
        hmm_.addNewState(new HMMState(first+"RSCk-"+String(i)));
				
				for (Size j = 0; j != residues.size(); ++j)
				{
					String second;
					second += residues[j];
					hmm_.addNewState(new HMMState(first+second+"bxyz"+String(i)));
					hmm_.addNewState(new HMMState(first+second+"bxyzk-"+String(i)));
					hmm_.addNewState(new HMMState(first+second+"axyz"+String(i)));
					hmm_.addNewState(new HMMState(first+second+"axyzk-"+String(i)));
				}
			}
		}

		cerr << "setting transitions" << endl;


		hmm_.setTransitionProbability("AABase1", "AABase2", 1);

		// set the initial transitions
		for (Size i = 1; i <= VISIBLE_MODEL_DEPTH; ++i)
		{
			if (i <= MODEL_DEPTH)
			{
		
				hmm_.setTransitionProbability("BB"+String(i), "end"+String(i), 0.5);
				hmm_.setTransitionProbability("BBk-"+String(i), "endk-"+String(i), 0.5);
			
				hmm_.setTransitionProbability("SC"+String(i), "end"+String(i), 0.5);
				hmm_.setTransitionProbability("SCk-"+String(i), "endk-"+String(i), 0.5);

				hmm_.setTransitionProbability("CR"+String(i), "end"+String(i), 0.5);
				hmm_.setTransitionProbability("CRk-"+String(i), "endk-"+String(i), 0.5);
	
				hmm_.setTransitionProbability("BB"+String(i), "AA"+String(i), 0.5);
				hmm_.setTransitionProbability("BBk-"+String(i), "AAk-"+String(i), 0.5);

				hmm_.setTransitionProbability("CR"+String(i), "A"+String(i), 0.5);
				hmm_.setTransitionProbability("CRk-"+String(i), "Ak-"+String(i), 0.5);

        hmm_.setTransitionProbability("SC"+String(i), "ASC"+String(i), 0.5);
        hmm_.setTransitionProbability("SCk-"+String(i), "ASCk-"+String(i), 0.5);
			}
			else
			{
				hmm_.addSynonymTransition("BBcenter", "endcenter", "BB"+String(i), "end"+String(i));
				hmm_.addSynonymTransition("BBcenter", "endcenter", "BBk-"+String(i), "endk-"+String(i));
				hmm_.setTransitionProbability("BBcenter", "endcenter", 0.5);
			
				hmm_.addSynonymTransition("CRcenter", "endcenter", "CR"+String(i), "end"+String(i));
				hmm_.addSynonymTransition("CRcenter", "endcenter", "CRk-"+String(i), "endk-"+String(i));
				hmm_.setTransitionProbability("CRcenter", "endcenter", 0.5);

				hmm_.addSynonymTransition("SCcenter", "endcenter", "SC"+String(i), "end"+String(i));
				hmm_.addSynonymTransition("SCcenter", "endcenter", "SCk-"+String(i), "endk-"+String(i));
				hmm_.setTransitionProbability("SCcenter", "endcenter", 0.5);

				hmm_.addSynonymTransition("BBcenter", "AAcenter", "BB"+String(i), "AA"+String(i));
				hmm_.addSynonymTransition("BBcenter", "AAcenter", "BBk-"+String(i), "AAk-"+String(i));
				hmm_.setTransitionProbability("BBcenter", "AAcenter", 0.5);

				hmm_.addSynonymTransition("CRcenter", "Acenter", "CR"+String(i), "A"+String(i));
				hmm_.addSynonymTransition("CRcenter", "Acenter", "CRk-"+String(i), "Ak-"+String(i));
				hmm_.setTransitionProbability("CRcenter", "Acenter", 0.5);

        hmm_.addSynonymTransition("SCcenter", "ASCcenter", "SC"+String(i), "ASC"+String(i));
        hmm_.addSynonymTransition("SCcenter", "ASCcenter", "SCk-"+String(i), "ASCk-"+String(i));
        hmm_.setTransitionProbability("SCcenter", "ASCcenter", 0.5);
			}
			
			for (Size ii = 0; ii != residues.size(); ++ii)
			{
				String first;
				first += residues[ii];

				// CR D
				hmm_.addSynonymTransition("AABase1", "AABase2", "A"+String(i), first+"D"+String(i));
				hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-"+String(i), first+"Dk-"+String(i));
			
				hmm_.addSynonymTransition(first+"D", "D", first+"D"+String(i), "D"+String(i));
				hmm_.addSynonymTransition(first+"D", "end", first+"D"+String(i), "end"+String(i));
				hmm_.addSynonymTransition(first+"D", "D", first+"Dk-"+String(i), "Dk-"+String(i));
				hmm_.addSynonymTransition(first+"D", "end", first+"Dk-"+String(i), "endk-"+String(i));

				hmm_.setTransitionProbability(first+"D", "D", 0.5);
				hmm_.setTransitionProbability(first+"D", "end", 0.5);
		
				// CR E
				hmm_.addSynonymTransition("AABase1", "AABase2", "A"+String(i), first+"E"+String(i));
        hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-"+String(i), first+"Ek-"+String(i));
        
        hmm_.addSynonymTransition(first+"E", "E", first+"E"+String(i), "E"+String(i));
				hmm_.addSynonymTransition(first+"E", "end", first+"E"+String(i), "end"+String(i));
        hmm_.addSynonymTransition(first+"E", "E", first+"Ek-"+String(i), "Ek-"+String(i));
				hmm_.addSynonymTransition(first+"E", "end", first+"Ek-"+String(i), "endk-"+String(i));
        
        hmm_.setTransitionProbability(first+"E", "E", 0.5);
       	hmm_.setTransitionProbability(first+"E", "end", 0.5);
		
				// SC K
				hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+String(i), first+"K"+String(i));
        hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+String(i), first+"Kk-"+String(i));

        hmm_.addSynonymTransition(first+"K", "K", first+"K"+String(i), "K"+String(i));
				hmm_.addSynonymTransition(first+"K", "end", first+"K"+String(i), "end"+String(i));
        hmm_.addSynonymTransition(first+"K", "K", first+"Kk-"+String(i), "Kk-"+String(i));
				hmm_.addSynonymTransition(first+"K", "end", first+"Kk-"+String(i), "endk-"+String(i));

        hmm_.setTransitionProbability(first+"K", "K", 0.5);
        hmm_.setTransitionProbability(first+"K", "end", 0.5);
			
				// SC H
				hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+String(i), first+"H"+String(i));
        hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+String(i), first+"Hk-"+String(i));

        hmm_.addSynonymTransition(first+"H", "H", first+"H"+String(i), "H"+String(i));
				hmm_.addSynonymTransition(first+"H", "end", first+"H"+String(i), "end"+String(i));
        hmm_.addSynonymTransition(first+"H", "H", first+"Hk-"+String(i), "Hk-"+String(i));
				hmm_.addSynonymTransition(first+"H", "end", first+"Hk-"+String(i), "endk-"+String(i));

        hmm_.setTransitionProbability(first+"H", "H", 0.5);
        hmm_.setTransitionProbability(first+"H", "end", 0.5);

				// SC R	
				hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+String(i), first+"RSC"+String(i));
        hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+String(i), first+"RSCk-"+String(i));

        hmm_.addSynonymTransition(first+"RSC", "R", first+"RSC"+String(i), "R"+String(i));
				hmm_.addSynonymTransition(first+"RSC", "end", first+"RSC"+String(i), "end"+String(i));
        hmm_.addSynonymTransition(first+"RSC", "R", first+"RSCk-"+String(i), "Rk-"+String(i));
				hmm_.addSynonymTransition(first+"RSC", "end", first+"RSCk-"+String(i), "endk-"+String(i));

        hmm_.setTransitionProbability(first+"RSC", "R", 0.5);
        hmm_.setTransitionProbability(first+"RSC", "end", 0.5);
				

				for (Size j = 0; j != residues.size(); ++j)
				{
					String second;
					second += residues[j];

					hmm_.addSynonymTransition("AABase1", "AABase2", "AA"+String(i), first+second+"bxyz"+String(i));
					hmm_.addSynonymTransition("AABase1", "AABase2", "AA"+String(i), first+second+"axyz"+String(i));

					hmm_.addSynonymTransition("AABase1", "AABase2", "AAk-"+String(i), first+second+"bxyzk-"+String(i));
					hmm_.addSynonymTransition("AABase1", "AABase2", "AAk-"+String(i), first+second+"axyzk-"+String(i));
				
					if ((second == "P")&& i <= 2)
					{
            hmm_.setTransitionProbability(first+second+"bxyz"+String(i), "bxyz"+String(i), 0.01);
						hmm_.addSynonymTransition(first+second+"bxyz"+String(i), "end", first+second+"bxyz"+String(i), "end"+String(i));
						hmm_.setTransitionProbability(first+second+"bxyz"+String(i), "end", 0.99);

            hmm_.setTransitionProbability(first+second+"axyz"+String(i), "axyz"+String(i), 0.01);
						hmm_.addSynonymTransition(first+second+"axyz"+String(i), "end", first+second+"axyz"+String(i), "end"+String(i));
						hmm_.setTransitionProbability(first+second+"axyz"+String(i), "end", 0.99);
					}
					else
					{
					
						hmm_.addSynonymTransition(first+second+"bxyz", "bxyz", first+second+"bxyz"+String(i), "bxyz"+String(i));
						hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyz"+String(i), "end"+String(i));

						hmm_.addSynonymTransition(first+second+"axyz", "axyz", first+second+"axyz"+String(i), "axyz"+String(i));
						hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyz"+String(i), "end"+String(i)); 

						hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyz"+String(i), "end"+String(i));
						hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyz"+String(i), "end"+String(i));
					}
					
					hmm_.addSynonymTransition(first+second+"bxyz", "bxyz", first+second+"bxyzk-"+String(i), "bxyzk-"+String(i));
					hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyzk-"+String(i), "endk-"+String(i));
					
					hmm_.addSynonymTransition(first+second+"axyz", "axyz", first+second+"axyzk-"+String(i), "axyzk-"+String(i));
					hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyzk-"+String(i), "endk-"+String(i));
					
					hmm_.setTransitionProbability(first+second+"bxyz", "bxyz", 0.5);
					hmm_.setTransitionProbability(first+second+"axyz", "axyz", 0.5);

					hmm_.setTransitionProbability(first+second+"bxyz", "end", 0.5);
					hmm_.setTransitionProbability(first+second+"axyz", "end", 0.5);
					hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyzk-"+String(i), "end"+String(i));
          hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyzk-"+String(i), "end"+String(i));
				}
			}
		}

		hmm_.disableTransitions();
		hmm_.buildSynonyms();
		cerr << "#states=" << hmm_.getNumberOfStates() << endl;
	}*/

	void PILISModel::initModels_()
	{
		name_to_enum_["b-H2O"] = B_H2O;
    name_to_enum_["b-NH3"] = B_NH3;
    name_to_enum_["b_loss_end"] = B_LOSS_END;
    name_to_enum_["bion"] = B_ION;
    name_to_enum_["b_Base1"] = B_BASE1;
    name_to_enum_["b_Base2"] = B_BASE2;
    name_to_enum_["b_H2O_S"] = B_H2O_S;
    name_to_enum_["b_H2O_T"] = B_H2O_T;
    name_to_enum_["b_H2O_E"] = B_H2O_E;
    name_to_enum_["b_H2O_D"] = B_H2O_D;
    name_to_enum_["b_NH3_K"] = B_NH3_K;
    name_to_enum_["b_NH3_R"] = B_NH3_R;
	
		enum_to_name_[B_H2O] = "b-H2O";
		enum_to_name_[B_NH3] = "b-NH3";
		enum_to_name_[B_LOSS_END] = "b_loss_end";
		enum_to_name_[B_ION] = "bion";
		enum_to_name_[B_BASE1] = "b_Base1";
		enum_to_name_[B_BASE2] = "b_Base2";
		enum_to_name_[B_H2O_S] = "b_H2O_S";
		enum_to_name_[B_H2O_T] = "b_H2O_T";
		enum_to_name_[B_H2O_E] = "b_H2O_E";
		enum_to_name_[B_H2O_D] = "b_H2O_D";
		enum_to_name_[B_NH3_K] = "b_NH3_K";
		enum_to_name_[B_NH3_R] = "b_NH3_R";
		
    name_to_enum_["y-H2O"] = Y_H2O;
    name_to_enum_["y-NH3"] = Y_NH3;
    name_to_enum_["y_loss_end"] = Y_LOSS_END;
    name_to_enum_["yion"] = Y_ION;
    name_to_enum_["y_Base1"] = Y_BASE1;
    name_to_enum_["y_Base2"] = Y_BASE2;
    name_to_enum_["y_H2O_S"] = Y_H2O_S;
    name_to_enum_["y_H2O_T"] = Y_H2O_T;
    name_to_enum_["y_H2O_E"] = Y_H2O_E;
    name_to_enum_["y_H2O_D"] = Y_H2O_D;
    name_to_enum_["y_NH3_K"] = Y_NH3_K;
    name_to_enum_["y_NH3_R"] = Y_NH3_R;
		name_to_enum_["y_H2O_Q1"] = Y_H2O_Q1;
		name_to_enum_["y_H2O_Cterm"] = Y_H2O_CTERM;

		enum_to_name_[Y_H2O] = "y-H2O";
		enum_to_name_[Y_NH3] = "y-NH3";
		enum_to_name_[Y_LOSS_END] = "y_loss_end";
		enum_to_name_[Y_ION] = "yion";
		enum_to_name_[Y_BASE1] = "y_Base1";
		enum_to_name_[Y_BASE2] = "y_Base2";
		enum_to_name_[Y_H2O_S] = "y_H2O_S";
		enum_to_name_[Y_H2O_T] = "y_H2O_T";
		enum_to_name_[Y_H2O_E] = "y_H2O_E";
		enum_to_name_[Y_H2O_D] = "y_H2O_D";
		enum_to_name_[Y_NH3_K] = "y_NH3_K";
		enum_to_name_[Y_NH3_R] = "y_NH3_R";
		enum_to_name_[Y_H2O_Q1] = "y_H2O_Q1";
		enum_to_name_[Y_H2O_CTERM] = "y_H2O_Cterm";

    name_to_enum_["[M+H]-H2O"] = PRE_MH_H2O;
    name_to_enum_["[M+H]-NH3"] = PRE_MH_NH3;
		name_to_enum_["[M+H]"] = PRE_MH;
    name_to_enum_["Preend"] = PRE_END;
    name_to_enum_["Pre"] = PRE_ION;
    name_to_enum_["Pre_Base1"] = PRE_BASE1;
    name_to_enum_["Pre_Base2"] = PRE_BASE2;
    name_to_enum_["Pre_H2O_S"] = PRE_H2O_S;
    name_to_enum_["Pre_H2O_T"] = PRE_H2O_T;
    name_to_enum_["Pre_H2O_E"] = PRE_H2O_E;
    name_to_enum_["Pre_H2O_D"] = PRE_H2O_D;
    name_to_enum_["Pre_NH3_K"] = PRE_NH3_K;
    name_to_enum_["Pre_NH3_R"] = PRE_NH3_R;
    name_to_enum_["Pre_H2O_Q1"] = PRE_H2O_Q1;
    name_to_enum_["Pre_H2O_Cterm"] = PRE_H2O_CTERM;

		enum_to_name_[PRE_MH_H2O] = "[M+H]-H2O";
		enum_to_name_[PRE_MH_NH3] = "[M+H]-NH3";
		enum_to_name_[PRE_MH] = "[M+H]";
		enum_to_name_[PRE_END] = "Preend";
		enum_to_name_[PRE_ION] = "Pre";
		enum_to_name_[PRE_BASE1] = "Pre_Base1";
		enum_to_name_[PRE_BASE2] = "Pre_Base2";
		enum_to_name_[PRE_H2O_S] = "Pre_H2O_S";
		enum_to_name_[PRE_H2O_T] = "Pre_H2O_T";
		enum_to_name_[PRE_H2O_E] = "Pre_H2O_E";
		enum_to_name_[PRE_H2O_D] = "Pre_H2O_D";
		enum_to_name_[PRE_NH3_K] = "Pre_NH3_K";
		enum_to_name_[PRE_NH3_R] = "Pre_NH3_R";
		enum_to_name_[PRE_H2O_Q1] = "Pre_H2O_Q1";
		enum_to_name_[PRE_H2O_CTERM] = "Pre_H2O_Cterm";
	
    hmms_losses_[Residue::BIon] = HiddenMarkovModelLight();
    //parseModelFile(String("model_losses_bions.dat"), &hmms_losses_[Residue::BIon]);
    hmms_losses_[Residue::YIon] = HiddenMarkovModelLight();
    //parseModelFile(String("model_losses_yions.dat"), &hmms_losses_[Residue::YIon]);

		for (HashMap<States_, String>::ConstIterator it = enum_to_name_.begin(); it != enum_to_name_.end(); ++it)
		{
			hmm_precursor_.addIdToName(it->first, it->second);
			hmms_losses_[Residue::BIon].addIdToName(it->first, it->second);
			hmms_losses_[Residue::YIon].addIdToName(it->first, it->second);
		}
		
	}

	void PILISModel::parseBaseModel_(const TextFile::ConstIterator& begin, const TextFile::ConstIterator& end)
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
          hmm_.addNewState(new HMMState(split[1], hidden));
					//cerr << "added new state: '" << split[1] << "', " << hidden << endl;
          continue;
        }

        if (id == "Synonym")
        {
					//++num_syn;
					//hmm_.addSynonymTransition(split[1], split[2], split[3], split[4]);
					hmm_.addSynonymTransition(split[3], split[4], split[1], split[2]);
          continue;
        }

        if (id == "Transition")
        {
          hmm_.setTransitionProbability(split[1], split[2], split[3].toFloat());
          continue;
        }
      }
    }
    //hmm_.disableTransitions();
		hmm_.buildSynonyms();

		hmm_.disableTransitions();
	
		//cerr << hmm_.getNumberOfStates() << endl;
		
		return;
	}
	
	void PILISModel::parseHMMLightModel_(const TextFile::ConstIterator& begin, const TextFile::ConstIterator& end, HiddenMarkovModelLight& model)
	{
		if (begin == end)
		{
			return;
		}

		for (TextFile::ConstIterator it = begin; it != end; ++it)
		{	
			String line = *it;
			// comment?
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
					if (name_to_enum_.has(split[1]))
					{
						model.addNewState(new HMMStateLight(name_to_enum_[split[1]], hidden));
					}
					else
					{
						throw Exception::ParseError(__FILE__, __LINE__, "", "'" + split[1] + "'", " state not found");
					}
																								
					continue;
				}

				if (id == "Synonym")
				{
					if (!name_to_enum_.has(split[1]))
          {
            throw Exception::ParseError(__FILE__, __LINE__, "", "'" + split[1] + "'", " state not found");
          }
          if (!name_to_enum_.has(split[2]))
          {
            throw Exception::ParseError(__FILE__, __LINE__, "", "'" + split[2] + "'", " state not found");
          }
          if (!name_to_enum_.has(split[3]))
          {
            throw Exception::ParseError(__FILE__, __LINE__, "", "'" + split[3] + "'", " state not found");
          }
          if (!name_to_enum_.has(split[4]))
          {
            throw Exception::ParseError(__FILE__, __LINE__, "", "'" + split[4] + "'", " state not found");
          }

					model.addSynonymTransition(name_to_enum_[split[3]], name_to_enum_[split[4]], name_to_enum_[split[1]], name_to_enum_[split[2]]);
					continue;
				}
			
				if (id == "Transition")
				{
					if (!name_to_enum_.has(split[1]))
					{
						throw Exception::ParseError(__FILE__, __LINE__, "", "'" + split[1] + "'", " state not found");
					}
					if (!name_to_enum_.has(split[2]))
					{
						throw Exception::ParseError(__FILE__, __LINE__, "", "'" + split[2] + "'", " state not found");
					}
					model.setTransitionProbability(name_to_enum_[split[1]], name_to_enum_[split[2]], split[3].toFloat());
					continue;
				}
			}
		}
		model.disableTransitions();
		model.buildSynonyms();
	}

/*
	Param& PILISModel::getParam()
	{
		return param_;
	}

	const Param& PILISModel::getParam() const
	{
		return param_;
	}

	void PILISModel::setParam(const Param& param)
	{
		param_ = param;
	}

	void PILISModel::resetToDefaultParam()
	{
		param_ = default_;
	}
*/
} // namespace OpenMS

