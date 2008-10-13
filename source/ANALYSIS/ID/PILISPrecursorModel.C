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


#include <OpenMS/ANALYSIS/ID/PILISPrecursorModel.h>
//#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
//#include <OpenMS/CHEMISTRY/AASequence.h>
//#include <OpenMS/SYSTEM/File.h>

//#include <cmath>
//#include <sstream>
//#include <algorithm>
//#include <numeric>
//#include <fstream>


#define PRECURSOR_DEBUG

using namespace std;

namespace OpenMS 
{

	PILISPrecursorModel::PILISPrecursorModel()
		: DefaultParamHandler("PILISPrecursorModel")
	{	
		defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, used to identify the precursor peak and its loss peaks for training");
		defaults_.setValue("peak_mass_tolerance", 0.3, "Peak mass tolerance of the product ions, used to identify the ions for training");
		defaults_.setValue("fixed_modifications", "", "Fixed modifications in format '57.001@C'");

		defaults_.setValue("pseudo_counts", 1e-15, "Value which is added for every transition trained of the underlying hidden Markov model");

		defaultsToParam_();
	}

	PILISPrecursorModel::~PILISPrecursorModel()
	{
	}

	PILISPrecursorModel::PILISPrecursorModel(const PILISPrecursorModel& model)
		: DefaultParamHandler(model),
			hmm_precursor_(model.hmm_precursor_)
	{
	}

	PILISPrecursorModel& PILISPrecursorModel::operator = (const PILISPrecursorModel& model)
	{
		if (this != &model)
		{
			DefaultParamHandler::operator=(model);
	    hmm_precursor_ = model.hmm_precursor_;
		}
		return *this;
	}
	

	/*
	void PILISPrecursorModel::writeGraphMLFile(const String& filename)
	{
		hmm_precursor_.writeGraphMLFile(filename);
		return;
	}
	*/

	void PILISPrecursorModel::getPrecursorIntensitiesFromSpectrum_(const RichPeakSpectrum& train_spec, Map<String, double>& peak_ints, double peptide_weight, UInt charge)
	{
#ifdef PRECURSOR_DEBUG
		cerr << "PILISPrecursorModel::getPrecursorIntensitiesFromSpectrum_(#peaks=" << train_spec.size() << ", weight=" << peptide_weight << ", charge=" << charge << ")" << endl;
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
	
	void PILISPrecursorModel::trainPrecursorIons_(double initial_probability, const Map<String, double>& ints, const AASequence& peptide)
	{
#ifdef PRECURSOR_DEBUG
		cerr << "PILISPrecursorModel::trainPrecursorIons_(" << initial_probability << ", " << ints.size() << ", " << peptide << ")" << endl;
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

  void PILISPrecursorModel::enablePrecursorIonStates_(const AASequence& peptide)
	{
#ifdef PRECURSOR_DEBUG
		cerr << "void PILISPrecursorModel::enablePrecursorIonStates_(" << peptide << ")" << endl;
#endif
		UInt max_explicit(4);
		Map<String, UInt> double_losses;
		vector<String> nexts;

		for (AASequence::ConstIterator it1 = peptide.begin(); it1 != peptide.end(); ++it1)
		{
			AASequence aa1;
			aa1 += &*it1;

			vector<EmpiricalFormula> losses1 = it1->getLossFormulas();

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

			vector<EmpiricalFormula> losses;
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

	void PILISPrecursorModel::evaluate()
	{
		hmm_precursor_.evaluate();
	}

	void PILISPrecursorModel::getPrecursorIons_(Map<String, double>& intensities, double initial_probability, const AASequence& precursor)
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

	void PILISPrecursorModel::updateMembers_()
	{
		double pseudo_counts = (double)param_.getValue("pseudo_counts");
		hmm_precursor_.setPseudoCounts(pseudo_counts);
	}
} // namespace OpenMS


