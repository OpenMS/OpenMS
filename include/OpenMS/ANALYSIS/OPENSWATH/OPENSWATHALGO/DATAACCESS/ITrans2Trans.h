/*
 * ITrans2Trans.h
 *
 *  Created on: Jul 30, 2012
 *      Author: witold
 */

#ifndef ITRANS2TRANS_H_
#define ITRANS2TRANS_H_

#include <map>
#include "TransitionExperiment.h"
#include "Transitions.h"

namespace OpenSwath
{

	inline void convert(LightTargetedExperiment & lte,
			std::map<std::string, std::vector<OpenSwath::LightTransition> > & transmap)
	{

		typedef std::pair<std::string, std::vector<OpenSwath::LightTransition> > Mpair;
		typedef std::map<std::string, std::vector<OpenSwath::LightTransition> > Mmap;
		std::vector<LightTransition> ltrans = lte.getTransitions();

		/*iterate over transitions*/
		std::vector<LightTransition>::iterator ltrit = ltrans.begin();
		std::vector<LightTransition>::iterator ltrend = ltrans.end();
		for (; ltrit != ltrend; ++ltrit) {
			std::string pepref = ltrit->getPeptideRef();

			Mmap::iterator it = transmap.find(pepref);
			if (it == transmap.end()) {
				std::vector<LightTransition> ltv;
				ltv.push_back(*ltrit);
				transmap.insert(Mpair(pepref, ltv));
			} else {
				it->second.push_back(*ltrit);
			}
		}
	} //end convert


	// spiegel
	inline bool findPeptide(const LightTargetedExperiment & lte, const std::string & peptideRef,
			LightPeptide & pep)
	{
		std::vector<LightPeptide>::const_iterator beg =lte.peptides.begin();
		std::vector<LightPeptide>::const_iterator end =lte.peptides.end();
		for (; beg != end; ++beg) {
			//std::cout << beg->id << " " << peptideRef << std::endl;
			if (beg->id.compare(peptideRef) == 0) {
				pep = *beg;
				return true;
			}
		}
		return false;
	}

} //end namespace

#endif /* ITRANS2TRANS_H_ */
