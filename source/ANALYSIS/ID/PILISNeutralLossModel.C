// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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


#include <OpenMS/ANALYSIS/ID/PILISNeutralLossModel.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

// TODOs
//
// - Fix damn losses!
//
// N-terminal Lysine NH3-loss
// 
// if E N-term and charge directed state, abundant loss of water
// the bk-1/bk-2 pathway somtimes show bx+18 ions, especially when H is C-term to it
//

#define NEUTRAL_LOSS_MODEL_DEBUG
#undef NEUTRAL_LOSS_MODEL_DEBUG

#define NEUTRAL_LOSS_MODEL_FILE_DEBUG
#undef NEUTRAL_LOSS_MODEL_FILE_DEBUG

using namespace std;

namespace OpenMS 
{

	PILISNeutralLossModel::PILISNeutralLossModel()
		: DefaultParamHandler("PILISNeutralLossModel"),
			num_explicit_(0)
	{	
		defaults_.setValue("fragment_mass_tolerance", 0.4, "Peak mass tolerance of the product ions, used to identify the ions for training");
		
		defaults_.setValue("fixed_modifications", StringList::create(""), "Fixed modifications");
		defaults_.setValue("variable_modifications", StringList::create(""), "Variable modifications");

		defaults_.setValue("pseudo_counts", 1e-15, "Value which is added for every transition trained of the underlying hidden Markov model");
		defaults_.setValue("num_explicit", 2, "Number of explicitely modeled losses from the same kind of amino acid or combinations thereof");

		defaults_.setValue("min_int_to_train", 0.1, "Minimal intensity a ion and its losses must have to be considered for training.");
		
		defaults_.setValue("C_term_H2O_loss", "true", "enable water loss of the C-terminus");
		defaults_.setValidStrings("C_term_H2O_loss", StringList::create("true,false"));
		
		defaults_.setValue("ion_name", "p", "Ion base names used to set in meta values");
		defaults_.setValidStrings("ion_name", StringList::create("p,a,b,b2,y"));

		defaults_.setValue("enable_double_losses", "true", "if true, two different losses can occur at the same time, e.g. -H2O and -NH3 forming loss of -35Da");
		defaults_.setValidStrings("enable_double_losses", StringList::create("true,false"));

		defaultsToParam_();
	}

	PILISNeutralLossModel::~PILISNeutralLossModel()
	{
	}

	PILISNeutralLossModel::PILISNeutralLossModel(const PILISNeutralLossModel& model)
		: DefaultParamHandler(model),
			hmm_precursor_(model.hmm_precursor_),
			num_explicit_(model.num_explicit_)
	{
	}

	PILISNeutralLossModel& PILISNeutralLossModel::operator = (const PILISNeutralLossModel& model)
	{
		if (this != &model)
		{
			DefaultParamHandler::operator=(model);
	    hmm_precursor_ = model.hmm_precursor_;
			num_explicit_ = model.num_explicit_;
		}
		return *this;
	}

	void PILISNeutralLossModel::setHMM(const HiddenMarkovModel& model)
	{
		hmm_precursor_ = model;
	}

	const HiddenMarkovModel& PILISNeutralLossModel::getHMM() const
	{
		return hmm_precursor_;
	}
	
	double PILISNeutralLossModel::train(const RichPeakSpectrum& spec, const AASequence& peptide, double ion_weight, UInt charge, double peptide_weight)
	{
	#ifdef NEUTRAL_LOSS_MODEL_DEBUG
		cerr << "PILISNeutralLossModel::train(#spec.size()=" << spec.size() << ", peptide=" << peptide << ", ion_weight=" << ion_weight << ", charge=" << charge << ", peptide_weight=" << peptide_weight << ")" << endl;
	#endif
		Map<String, double> peak_ints;
		double intensity_sum = getIntensitiesFromSpectrum_(spec, peak_ints, ion_weight, peptide, charge);

		String ion_name((String)param_.getValue("ion_name"));
		double min_int_to_train((double)param_.getValue("min_int_to_train"));
		#ifdef NEUTRAL_LOSS_MODEL_DEBUG
		cerr << ion_name << " intensity_sum=" << intensity_sum << ", min_int_to_train=" << min_int_to_train << endl;
		#endif
		if (intensity_sum < min_int_to_train || (ion_name != "p" && peak_ints[ion_name] == 0) || (ion_weight < peptide_weight / 2.0))
		{
			return intensity_sum;
		}
		
		double max_int(0);
		for (Map<String, double>::ConstIterator it = peak_ints.begin(); it != peak_ints.end(); ++it)
		{
			//it->second /= intensity_sum;
			if (it->second > max_int)
			{
				max_int = it->second;
			}
#ifdef NEUTRAL_LOSS_MODEL_DEBUG
			cerr << "LOSS: " << it->first << " " << it->second << " " << peptide << ", sum=" << intensity_sum <<endl;
#endif
		}

		for (Map<String, double>::Iterator it = peak_ints.begin(); it != peak_ints.end(); ++it)
		{
			it->second /= max_int;
		}
		
		trainIons_(1.0, peak_ints, peptide);
		
		return intensity_sum;
	}

	void PILISNeutralLossModel::getIons(vector<RichPeak1D>& peaks, const AASequence& peptide, double initial_prob)
	{
		Map<String, double> pre_ints;
		getIons_(pre_ints, initial_prob, peptide);
		for (Map<String, double>::ConstIterator it = pre_ints.begin(); it != pre_ints.end(); ++it)
		{
			RichPeak1D p;
			p.setIntensity(it->second);
			p.setMetaValue("IonName", it->first);
			// get the position of the peak
			vector<String> split;
			it->first.split('-', split);
			if (split.size() == 0)
			{
				p.setPosition(0);
			}
			else
			{
				if (split.size() == 2)
				{
					p.setPosition(- EmpiricalFormula(split[1]).getMonoWeight());
				}
				else
				{
					if (split.size() == 3)
					{
						p.setPosition(- (EmpiricalFormula(split[1]).getMonoWeight() + EmpiricalFormula(split[2]).getMonoWeight()));
					}
				}
			}
			peaks.push_back(p);
		}
	}
	
	double PILISNeutralLossModel::getIntensitiesFromSpectrum_(const RichPeakSpectrum& train_spec, Map<String, double>& peak_ints, double ion_weight, const AASequence& peptide, UInt charge)
	{
#ifdef NEUTRAL_LOSS_MODEL_DEBUG
		cerr << "PILISNeutralLossModel::getIntensitiesFromSpectrum_(#peaks=" << train_spec.size() << ", weight=" << ion_weight << ", peptide=" << peptide  <<  ", charge=" << charge << ")" << endl;
#endif
		set<String> precursor_losses;
		bool enable_double_losses(param_.getValue("enable_double_losses").toBool());
		for (AASequence::ConstIterator it = peptide.begin(); it != peptide.end(); ++it)
		{
			vector<EmpiricalFormula> losses = it->getLossFormulas();
			for (vector<EmpiricalFormula>::const_iterator lit = losses.begin(); lit != losses.end(); ++lit)
			{
				//cerr << "AA-LOSS: " << it->getOneLetterCode() << " " << lit->getString() << endl;
				precursor_losses.insert(lit->getString());
			}
		}

		bool enable_COOH(param_.getValue("C_term_H2O_loss").toBool());
		if (enable_COOH)
		{
			String h2o(EmpiricalFormula("H2O").getString());
			precursor_losses.insert(h2o);
		}

		//vector<EmpiricalFormula> pre_losses;
		vector<double> pre_loss_weights;
		vector<String> pre_loss_names;

		for (set<String>::const_iterator it = precursor_losses.begin(); it != precursor_losses.end(); ++it)
		{
			//pre_losses.push_back(*it);
			pre_loss_weights.push_back(EmpiricalFormula(*it).getMonoWeight());
			pre_loss_names.push_back(*it);
		}

		// init all possible losses with zero intensity
		String ion_name = param_.getValue("ion_name");
		Map<String, double> loss_weights;
		for (UInt i = 0; i != pre_loss_names.size(); ++i)
		{
			String name1(pre_loss_names[i]);
			peak_ints[ion_name + "-" + name1] = 0;
			loss_weights[ion_name + "-" + name1] = pre_loss_weights[i];

			if (enable_double_losses)
			{
				for (UInt j = 0; j != pre_loss_names.size(); ++j)
				{
					String name2(pre_loss_names[j]);
        	String name;
        	if (name1 < name2)
        	{
          	name = ion_name + "-" + name1 + "-" + name2;
        	}
        	else
        	{
          	name = ion_name + "-" + name2 + "-" + name1;
        	}
					peak_ints[name] = 0;
					loss_weights[name] = pre_loss_weights[i] + pre_loss_weights[j];
				}
			}
		}
	
		loss_weights[ion_name] = 0;
		peak_ints[ion_name] = 0;

	
		// now match peaks to the different losses and combinations	
		double pre_error = (double)param_.getValue("fragment_mass_tolerance");
		double intensity_sum(0);
  	for (RichPeakSpectrum::ConstIterator it = train_spec.begin(); it != train_spec.end(); ++it)
    {
			for (Map<String, double>::ConstIterator loss_it = loss_weights.begin(); loss_it != loss_weights.end(); ++loss_it)
			{
				double diff = fabs(it->getMZ() - (ion_weight - loss_it->second + (double)charge) / (double)charge);
				if (diff < pre_error)
				{
					double factor = (pre_error - diff) / pre_error;
					#ifdef NEUTRAL_LOSS_MODEL_DEBUG
					cerr << "Found: ";
					cerr << "factor=" << factor
							 << " intensity=" << it->getIntensity()
							 << " @m/z=" << it->getMZ() 
							 << " ion_weight=" << ion_weight 
							 << " loss_weight=" << loss_it->second 
							 << " theo_ion_pos=" <<  (ion_weight - loss_it->second + (double)charge) / (double)charge
							 << " diff=" << diff << endl;
					#endif
					peak_ints[loss_it->first] += it->getIntensity() * factor;
					intensity_sum += it->getIntensity() * factor;
				}
			}
		}
		return intensity_sum;
	}
	
	void PILISNeutralLossModel::trainIons_(double initial_probability, const Map<String, double>& ints, const AASequence& peptide)
	{
#ifdef NEUTRAL_LOSS_MODEL_DEBUG
		cerr << "PILISNeutralLossModel::trainIons_(" << initial_probability << ", " << ints.size() << ", " << peptide << ")" << endl;
#endif
		// clean up
		hmm_precursor_.clearInitialTransitionProbabilities();
		hmm_precursor_.clearTrainingEmissionProbabilities();

		// set start transition prob
		hmm_precursor_.setInitialTransitionProbability("start", initial_probability);
	
		// set emission probabilities from the precursor ions present in the spectrum
		for (Map<String, double>::ConstIterator it = ints.begin(); it != ints.end(); ++it)
		{
			hmm_precursor_.setTrainingEmissionProbability(it->first, it->second);
		}
	
		// build the final model
		enableIonStates_(peptide);
	
#ifdef NEUTRAL_LOSS_MODEL_FILE_DEBUG
		stringstream ss;
		ss << peptide;
		hmm_precursor_.writeGraphMLFile("graphs/model_graph_train_" + (String)param_.getValue("ion_name") + "_" + peptide.toUnmodifiedString() + ".graphml");
#endif
		
		// train
		hmm_precursor_.train();

		//hmm_precursor_.dump();
		
		// reset the model
		hmm_precursor_.disableTransitions();

		return;
	}

  void PILISNeutralLossModel::enableIonStates_(const AASequence& peptide)
	{
#ifdef NEUTRAL_LOSS_MODEL_DEBUG
		cerr << "void PILISNeutralLossModel::enableIonStates_(" << peptide << ")" << endl;
#endif
	
		String ion_name = param_.getValue("ion_name");
		bool enable_COOH(param_.getValue("C_term_H2O_loss").toBool());
		String h2o;
		if (enable_COOH)
		{
			h2o = EmpiricalFormula("H2O").getString();
		}
		
		Map<String, UInt> double_losses;
		vector<String> nexts;

		bool enable_double_losses(param_.getValue("enable_double_losses").toBool());
		if (enable_double_losses)
		{
		for (AASequence::ConstIterator it1 = peptide.begin(); it1 != peptide.end(); ++it1)
		{
			AASequence aa1;
			aa1 += &*it1;

			//bool has_NTerm_loss(false);
			//has_NTerm_loss = it1->hasNTermNeutralLosses();

			vector<EmpiricalFormula> losses1 = it1->getLossFormulas();

			//if (has_NTerm_loss)
			//{
			//	losses1 = it1->getNTermLossFormulas();
			//}

			for (vector<EmpiricalFormula>::const_iterator loss_it1 = losses1.begin(); loss_it1 != losses1.end(); ++loss_it1)
			{
				String loss1 = loss_it1->getString();
				if (loss1 == "")
				{
					continue;
				}
		
				for (AASequence::ConstIterator it2 = it1 + 1; it2 != peptide.end(); ++it2)
				{
					AASequence aa2;
					aa2 += &*it2;


					vector<EmpiricalFormula> losses2 = it2->getLossFormulas();
					for (vector<EmpiricalFormula>::const_iterator loss_it2 = losses2.begin(); loss_it2 != losses2.end(); ++loss_it2)
					{
						String loss2 = loss_it2->getString();
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
							if (double_losses[name] > num_explicit_)
							{
								continue;
							}
							num = String(double_losses[name]);
						}
	
						name += losses + "_" + num;
				
						// enable transitions to emit state, and single emission states
						hmm_precursor_.enableTransition(name, ion_name + losses);

						//if (has_NTerm_loss)
						//{
						//	cerr << "Enabling NTermNeutralLoss: " << name << endl;
						//	hmm_precursor_.enableTransition(name, aa1.toString() + "-NTerm-" + loss1);
						//}
						//else
						//{
							hmm_precursor_.enableTransition(name, aa1.toString() + "-" + loss1 + "_1");
						//}
						vector<String> split;
						name.split('_', split);
						//if (!has_NTerm_loss)
						//{
							hmm_precursor_.addSynonymTransition(split[0], aa1.toString() + "-" + loss1, name, aa1.toString() + "-" + loss1 + "_1");
							hmm_precursor_.enableTransition(name, aa2.toString() + "-" + loss2 + "_1");
							hmm_precursor_.addSynonymTransition(split[0], aa2.toString() + "-" + loss2, name, aa2.toString() + "-" + loss2 + "_1");
						//}
				
						// store transition to connect the next lines
						nexts.push_back(name);
					}
				}
		
				String losses;
				if (h2o < loss1)
				{
					losses = "-" + h2o + "-" + loss1;
				}
				else
				{
					losses = "-" + loss1 + "-" + h2o;
				}
			
				if (enable_COOH)
				{
					// add loss - COOH loss
					String name = aa1.toString() + "COOH";
					// get num
 	       	String num;
 	       	if (!double_losses.has(name))
 	       	{
 	        	double_losses[name] = 1;
          	num = "1";
        	}
        	else
        	{
          	++double_losses[name];
          	if (double_losses[name] > num_explicit_)
          	{
            	continue;
          	}
          	num = String(double_losses[name]);
        	}
        	name += losses + "_" + num;

					hmm_precursor_.enableTransition(name, ion_name + losses);
					hmm_precursor_.enableTransition(name, aa1.toString() + "-" + loss1 + "_1");
					vector<String> split;
					name.split('_', split);
					hmm_precursor_.addSynonymTransition(split[0], aa1.toString() + "-" + loss1, name, aa1.toString() + "-" + loss1 + "_1");
				
					hmm_precursor_.enableTransition(name, "COOH-" + h2o + "_1");
					hmm_precursor_.addSynonymTransition(split[0], "COOH-" + h2o, name, "COOH-" + h2o + "_1");

					nexts.push_back(name);
				}
			}
		}
		}


		if (enable_double_losses)
		{
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
        	//hmm_precursor_.addSynonymTransition(*it, *next_it, hmm_precursor_.getTransitionProbability(split[0], split[0] + "-next"));
					hmm_precursor_.addSynonymTransition(split[0], split[0] + "-next", *it, *next_it);
				}
			}
		}
		
		Map<String, UInt> single_losses;
		vector<String> single_nexts;
		for (AASequence::ConstIterator it1 = peptide.begin(); it1 != peptide.end(); ++it1)
    {
      AASequence aa1;
      aa1 += &*it1;
			String name = aa1.toString();

			if (it1 == peptide.begin() && it1->hasNTermNeutralLosses())
			{
				vector<EmpiricalFormula> losses = it1->getNTermLossFormulas();

				for (vector<EmpiricalFormula>::const_iterator loss_it = losses.begin(); loss_it != losses.end(); ++loss_it)
				{
					String loss = loss_it->getString();

					if (loss == "")
					{
						continue;
					}
					String new_name = name + "-NTerm-" + loss;
					new_name = name + "-NTerm-" + loss;
					//cerr << "Enabling NTermNeutralLoss: " << new_name << " at peptide: " << peptide << endl;
					hmm_precursor_.enableTransition(name + "-NTerm-" + loss, ion_name + "-" + loss);
					single_nexts.push_back(new_name);
				}
			}
			else
			{
				vector<EmpiricalFormula> losses = it1->getLossFormulas();
				for (vector<EmpiricalFormula>::const_iterator loss_it = losses.begin(); loss_it != losses.end(); ++loss_it)
				{
			
					String loss = loss_it->getString();
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
							if (single_losses[name] > num_explicit_)
							{
								continue;
							}
							num = String(single_losses[name]);
						}	
						String new_name;

						new_name = name + "-" + loss + "_" + num;
						hmm_precursor_.enableTransition(new_name, ion_name + "-" + loss);
						single_nexts.push_back(new_name);
					}
				}
			}
		}

		// C-terminal water loss
		if (enable_COOH)
		{
			hmm_precursor_.enableTransition("COOH-" + h2o + "_1", ion_name + "-" + h2o);
			single_nexts.push_back("COOH-" + h2o + "_1");
		}

		if (!enable_double_losses && single_nexts.size() > 0)
		{
			hmm_precursor_.setTransitionProbability("start", *single_nexts.begin(), 1.0);
		}
		
		for (vector<String>::const_iterator it = single_nexts.begin(); it != single_nexts.end(); ++it)
		{
			if (it != (single_nexts.end() - 1))
			{
				vector<String>::const_iterator next_it = it + 1;
				hmm_precursor_.enableTransition(*it, *next_it);

				if (it->has('_'))
				{
					vector<String> split;
					it->split('_', split);
					hmm_precursor_.addSynonymTransition(split[0], split[0] + "-next", *it, *next_it);
				}
				else
				{
					hmm_precursor_.addSynonymTransition(*it, *it + "-next", *it, *next_it);
				}
			}
		}

		// last transition, no reaction took place
		if (single_nexts.size() > 0)
		{
			hmm_precursor_.enableTransition(single_nexts.back(), ion_name);
		}
		else
		{
			hmm_precursor_.setTransitionProbability("start", ion_name, 1.0);
		}


		return;
	}

	void PILISNeutralLossModel::evaluate()
	{
		hmm_precursor_.evaluate();
	}

	void PILISNeutralLossModel::getIons_(Map<String, double>& intensities, double initial_probability, const AASequence& precursor)
	{
		//cerr << "getIons_: " << initial_probability << " " << precursor << endl;
		hmm_precursor_.setInitialTransitionProbability("start", 1.0);

		enableIonStates_(precursor);

		Map<HMMState*, double> tmp;
		hmm_precursor_.calculateEmissionProbabilities(tmp);

		double max_prob(0);
		for (Map<HMMState*, double>::ConstIterator it = tmp.begin(); it != tmp.end(); ++it)
		{
			intensities[it->first->getName()] = it->second;
			if (it->second > max_prob)
			{
				max_prob = it->second;
			}
		}

		for (Map<String, double>::Iterator it = intensities.begin(); it != intensities.end(); ++it)
		{
			it->second = it->second / max_prob * initial_probability;
		}
		
#ifdef NEUTRAL_LOSS_MODEL_DEBUG
		for (Map<HMMState*, double>::ConstIterator it = tmp.begin(); it != tmp.end(); ++it)
		{
			cerr << it->first->getName() << " -> " << it->second << endl;
		}		
#endif

		

#ifdef NEUTRAL_LOSS_MODEL_FILE_DEBUG
		String pep_seq(precursor.toString());
		pep_seq.substitute('(', '_');
		pep_seq.substitute(')', '_');
		pep_seq.substitute(':', '_');
		hmm_precursor_.writeGraphMLFile("graphs/model_graph_predict_precursor_" + (String)param_.getValue("ion_name") + pep_seq + ".graphml");
#endif
		
		hmm_precursor_.disableTransitions();
	}

  void PILISNeutralLossModel::generateModel()
  {
    set<const Residue*> residues(ResidueDB::getInstance()->getResidues("Natural20"));
    StringList variable_modifications = param_.getValue("variable_modifications");

		// add variable modified residues
    for (StringList::const_iterator it = variable_modifications.begin(); it != variable_modifications.end(); ++it)
    {
      residues.insert(ResidueDB::getInstance()->getModifiedResidue(*it));
#ifdef NEUTRAL_LOSS_MODEL_DEBUG
      AASequence aa;
      aa += ResidueDB::getInstance()->getModifiedResidue(*it);
      cerr << "AddingModifiedResidue: " << aa.toString() << endl;
#endif
    }

		// TODO add also NTerm loss formulas
		// collect the different possible neutral loss molecules from the residues
    set<String> losses;
    for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
    {
      vector<EmpiricalFormula> res_losses = (*it)->getLossFormulas();
      for (vector<EmpiricalFormula>::const_iterator loss_it = res_losses.begin(); loss_it != res_losses.end(); ++loss_it)
      {
        String loss = loss_it->getString();
        losses.insert(loss);
#ifdef NEUTRAL_LOSS_MODEL_DEBUG
        cerr << "Loss: " << loss << ", of residue: " << (*it)->getName() << endl;
#endif
      }
    }

#ifdef NEUTRAL_LOSS_MODEL_DEBUG
    cerr << "Adding states..." << endl;
#endif

    // precursor is not fragmented
		String ion_name = param_.getValue("ion_name");
    hmm_precursor_.addNewState(new HMMState(ion_name, false));

    // emitting nodes for single losses
    for (set<String>::const_iterator it = losses.begin(); it != losses.end(); ++it)
    {
      hmm_precursor_.addNewState(new HMMState(ion_name + "-" + *it, false));
    }

		// emitting nodes for double losses
		bool enable_double_losses(param_.getValue("enable_double_losses").toBool());
    set<String> double_losses;
		if (enable_double_losses)
		{
    	for (set<String>::const_iterator it1 = losses.begin(); it1 != losses.end(); ++it1)
    	{
      	set<String>::const_iterator it2 = it1;
      	++it2;
      	for (; it2 != losses.end(); ++it2)
      	{
        	double_losses.insert(*it1 + "-" + *it2);
      	}
    	}
    	for (set<String>::const_iterator it1 = losses.begin(); it1 != losses.end(); ++it1)
    	{
      	for (set<String>::const_iterator it2 = it1; it2 != losses.end(); ++it2)
      	{
        	hmm_precursor_.addNewState(new HMMState(ion_name + "-" + *it1 + "-" + *it2, false));
      	}
    	}
		}

		// the next section creates the different reaction pathways
		bool enable_COOH(param_.getValue("C_term_H2O_loss").toBool());
		String h2o;
		if (enable_COOH)
		{
    	// H2O loss from the C-terminus
    	h2o = EmpiricalFormula("H2O").getString();
    	hmm_precursor_.addNewState(new HMMState("COOH-" + h2o, true));
		}
		hmm_precursor_.addNewState(new HMMState("start", true));

    // add double loss states
    // add edges from double loss states to single loss states and emitting states
#ifdef NEUTRAL_LOSS_MODEL_DEBUG
    cerr << "Adding double loss states" << endl;
#endif

		if (enable_COOH)
		{
    	for (UInt i = 0; i != num_explicit_; ++i)
    	{
      	hmm_precursor_.addNewState(new HMMState("COOH-" + h2o + "_" + String(i + 1)));
    	}
		}

    for (set<const Residue*>::const_iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
    {
			AASequence aa1;
			aa1 += *it1;
			
			// aa is located at the N-terminus, this can lead to a different fragmentation behaviour
			vector<EmpiricalFormula> NTerm_res_losses = (*it1)->getNTermLossFormulas();
			for (vector<EmpiricalFormula>::const_iterator loss_it1 = NTerm_res_losses.begin(); loss_it1 != NTerm_res_losses.end(); ++loss_it1)
			{
				String loss1 = loss_it1->getString();
				if (loss1 == "")
				{
					continue;
				}
				//cerr << "Adding NTermNeutralLoss: " << aa1.toString() << "-NTerm-" << loss1 << endl;
				hmm_precursor_.addNewState(new HMMState(aa1.toString() + "-NTerm-" + loss1));
				hmm_precursor_.addNewState(new HMMState(aa1.toString() + "-NTerm-" + loss1 + "-next"));
			}
						
      vector<EmpiricalFormula> res_losses1 = (*it1)->getLossFormulas();
      for (vector<EmpiricalFormula>::const_iterator loss_it1 = res_losses1.begin(); loss_it1 != res_losses1.end(); ++loss_it1)
      {
        String loss1 = loss_it1->getString();
        if (loss1 == "")
        {
          continue;
        }

				// aa is located somewhere in the peptide
        hmm_precursor_.addNewState(new HMMState(aa1.toString() + "-" + loss1));
        hmm_precursor_.addNewState(new HMMState(aa1.toString() + "-" + loss1 + "-next"));

        for (UInt i = 0; i != num_explicit_; ++i)
        {
          hmm_precursor_.addNewState(new HMMState(aa1.toString() + "-" + loss1 + "_" + String(i + 1)));
        }

				if (enable_COOH && enable_double_losses)
				{
        	String losses;
        	if (h2o < loss1)
        	{
          	losses = "-" + h2o  + "-" + loss1;
        	}
        	else
        	{
          	losses = "-" + loss1 + "-" + h2o;
        	}
        	String cooh_name = aa1.toString() + "COOH" + losses;
        	hmm_precursor_.addNewState(new HMMState(cooh_name));
        	hmm_precursor_.addNewState(new HMMState(cooh_name + "-next"));
        	hmm_precursor_.setTransitionProbability(cooh_name, ion_name + losses, 0.1);
        	hmm_precursor_.setTransitionProbability(cooh_name, aa1.toString() + "-" + loss1, 0.1);
        	hmm_precursor_.setTransitionProbability(cooh_name, "COOH-" + h2o, 0.1);
        	hmm_precursor_.setTransitionProbability(cooh_name, cooh_name + "-next", 0.7);

        	for (UInt i = 0; i != num_explicit_; ++i)
        	{
          	String cooh_name_num = cooh_name + "_" + String(i + 1);
          	hmm_precursor_.addNewState(new HMMState(cooh_name_num));
          	hmm_precursor_.addSynonymTransition(cooh_name, ion_name + losses, cooh_name_num, ion_name + losses);
        	}
				}
      }
    }

		if (enable_double_losses)
		{
     	for (set<const Residue*>::const_iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
      {
        AASequence aa1;
        aa1 += *it1;

        vector<EmpiricalFormula> res_losses1 = (*it1)->getLossFormulas();

        for (vector<EmpiricalFormula>::const_iterator loss_it1 = res_losses1.begin(); loss_it1 != res_losses1.end(); ++loss_it1)
        {
          String loss1 = loss_it1->getString();

          if (loss1 == "")
          {
            continue;
          }

          for (set<const Residue*>::const_iterator it2 = it1; it2 != residues.end(); ++it2)
          {
            AASequence aa2;
            aa2 += *it2;

            vector<EmpiricalFormula> res_losses2 = (*it2)->getLossFormulas();
            for (vector<EmpiricalFormula>::const_iterator loss_it2 = res_losses2.begin(); loss_it2 != res_losses2.end(); ++loss_it2)
            {
              String loss2 = loss_it2->getString();
              if (loss2 != "")
              {
                String name;
                String losses;
                if (loss1 < loss2)
                {
                  losses = "-" + loss1 + "-" + loss2;
                }
                else
                {
                  losses = "-" + loss2 + "-" + loss1;
                }
                if (aa1 < aa2)
                {
                  name = aa1.toString() + aa2.toString();
                }
                else
                {
                  name = aa2.toString() + aa1.toString();
                }
                hmm_precursor_.addNewState(new HMMState(name + losses, true));
                hmm_precursor_.addNewState(new HMMState(name + losses + "-next", true));

								/*
								hmm_precursor_.addNewState(new HMMState(name + losses + "-Nterm", true));
								hmm_precursor_.addNewState(new HMMState(name + losses + "-NTerm-next", true));
								*/
								
               	for (UInt i = 0; i != num_explicit_; ++i)
                {
                  hmm_precursor_.addNewState(new HMMState(name + losses + "_" + String(i + 1), true));
                }

                hmm_precursor_.setTransitionProbability(name + losses, ion_name + losses, 0.1);
                hmm_precursor_.setTransitionProbability(name + losses, aa1.toString() + "-" + loss1, 0.1);
                hmm_precursor_.setTransitionProbability(name + losses, aa2.toString() + "-" + loss2, 0.1);
                hmm_precursor_.setTransitionProbability(name + losses, name + losses + "-next", 0.7);

                for (UInt i = 0; i != num_explicit_; ++i)
                {
                  String state_name_num = name + losses + "_" + String(i + 1);
                  String state_name = name + losses;
                  hmm_precursor_.addSynonymTransition(state_name, ion_name + losses,  state_name_num, ion_name + losses);
                  hmm_precursor_.addSynonymTransition(state_name, aa1.toString() + "-" + loss1, state_name_num, aa1.toString() + "-" + loss1 + "_" + String(i + 1));
                  hmm_precursor_.addSynonymTransition(state_name, aa2.toString() + "-" + loss2, state_name_num, aa2.toString() + "-" + loss2 + "_" + String(i + 1));
                }
              }
            }
          }
        }
    	}
		}

#ifdef NEUTRAL_LOSS_MODEL_DEBUG
    cerr << "Adding single loss states" << endl;
#endif
		if (enable_COOH)
		{
    	String cooh_name = "COOH-" + h2o;
    	hmm_precursor_.addNewState(new HMMState(cooh_name + "-next"));

			if (!enable_double_losses)
			{
    		hmm_precursor_.setTransitionProbability(cooh_name, ion_name + "-" + h2o, 0.001);
    		hmm_precursor_.setTransitionProbability(cooh_name, ion_name, 0.199);
    		hmm_precursor_.setTransitionProbability(cooh_name, "COOH-" + h2o + "-next", 0.8);
			}
			else
			{
				hmm_precursor_.setTransitionProbability(cooh_name, ion_name + "-" + h2o, 0.25);
        hmm_precursor_.setTransitionProbability(cooh_name, ion_name, 0.25);
        hmm_precursor_.setTransitionProbability(cooh_name, "COOH-" + h2o + "-next", 0.5);
			}

    	hmm_precursor_.addSynonymTransition(cooh_name, ion_name, cooh_name + "_1", ion_name);
    	hmm_precursor_.addSynonymTransition(cooh_name, ion_name + "-" + h2o, cooh_name + "_1", ion_name + "-" + h2o);
		}
		
    for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
    {
			AASequence aa;
      aa += *it;
			vector<EmpiricalFormula> NTerm_res_losses = (*it)->getNTermLossFormulas();
			for (vector<EmpiricalFormula>::const_iterator loss_it = NTerm_res_losses.begin(); loss_it != NTerm_res_losses.end(); ++loss_it)
			{
				String loss = loss_it->getString();
				//cerr << "Enabling NTermNeutralLoss: " << aa.toString() << "-NTerm-" << loss << endl;
	
				if (loss == "")
				{
					continue;
				}
				if (!enable_double_losses)
				{
          hmm_precursor_.setTransitionProbability(aa.toString() + "-NTerm-" + loss, ion_name + "-" + loss, 0.25);
          hmm_precursor_.setTransitionProbability(aa.toString() + "-NTerm-" + loss, ion_name, 0.25);
          hmm_precursor_.setTransitionProbability(aa.toString() + "-NTerm-" + loss, aa.toString() + "-NTerm-" + loss + "-next", 0.5);
				}
				else
				{
          hmm_precursor_.setTransitionProbability(aa.toString() + "-NTerm-" + loss, ion_name + "-" + loss, 0.25);
          hmm_precursor_.setTransitionProbability(aa.toString() + "-NTerm-" + loss, ion_name, 0.25);
          hmm_precursor_.setTransitionProbability(aa.toString() + "-NTerm-" + loss, aa.toString() + "-NTerm-" + loss + "-next", 0.5);
				}
			}
			
						
      vector<EmpiricalFormula> res_losses = (*it)->getLossFormulas();
      for (vector<EmpiricalFormula>::const_iterator loss_it = res_losses.begin(); loss_it != res_losses.end(); ++loss_it)
      {
        String loss = loss_it->getString();

        if (loss == "")
        {
          continue;
        }

				if (!enable_double_losses)
				{
        	hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, ion_name + "-" + loss, 0.001);
        	hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, ion_name, 0.199);
        	hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, aa.toString() + "-" + loss + "-next", 0.8);
				}
				else
				{
					hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, ion_name + "-" + loss, 0.25);
          hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, ion_name, 0.25);
          hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, aa.toString() + "-" + loss + "-next", 0.5);
				}

        for (UInt i = 0; i != num_explicit_; ++i)
        {
          String name_num = aa.toString() + "-" + loss + "_" + String(i + 1);
          String name = aa.toString() + "-" + loss;
          hmm_precursor_.addSynonymTransition(name, ion_name + "-" + loss, name_num, ion_name + "-" + loss);
          hmm_precursor_.addSynonymTransition(name, ion_name, name_num, ion_name);
        }
      }
    }

#ifdef NEUTRAL_LOSS_MODEL_DEBUG
    cerr << "Finalizing HMM" << endl;
#endif
    hmm_precursor_.disableTransitions();
    //hmm_precursor_.buildSynonyms();


#ifdef NEUTRAL_LOSS_MODEL_FILE_DEBUG
    cerr << "#States: " << hmm_precursor_.getNumberOfStates() << endl;
    hmm_precursor_.writeGraphMLFile("hmm_" + (String)param_.getValue("ion_name") + ".graphML");
#endif

    return;
  }

	
	void PILISNeutralLossModel::updateMembers_()
	{
		double pseudo_counts = (double)param_.getValue("pseudo_counts");
		hmm_precursor_.setPseudoCounts(pseudo_counts);
		num_explicit_ = (UInt)param_.getValue("num_explicit");
	}
} // namespace OpenMS


