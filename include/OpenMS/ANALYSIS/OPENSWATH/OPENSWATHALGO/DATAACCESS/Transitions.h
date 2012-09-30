/*
 * Transitions.h
 *
 *  Created on: Jul 30, 2012
 *      Author: witold
 */

#ifndef TRANSITIONS_H_
#define TRANSITIONS_H_

namespace OpenSwath
{



	struct Peptide
	{
			double rt;
			int charge;
			std::string sequence;
			std::string id;
			int getChargeState() const
			{
				return charge;
			}
			std::vector<LightModification> modifications;
			std::vector<LightTransition> transitions;
	};

	struct Protein
	{
			std::string id;
			std::string sequence;
			std::vector<Peptide> peptides;
	};

	///
	struct TargetedExperiment
	{
			std::vector<Protein> proteins;
	};

}

#endif /* TRANSITIONS_H_ */
