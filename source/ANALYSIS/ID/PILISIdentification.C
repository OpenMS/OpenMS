// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
//

#include <OpenMS/ANALYSIS/ID/PILISIdentification.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

using namespace std;

namespace OpenMS
{
	
	PILISIdentification::PILISIdentification()
		:	sequence_db_(0),
			hmm_model_(0)
	{
		param_.setValue("prcr_m_tol", double(3.0));
		param_.setValue("max_candidates", 200);
		param_.setValue("score_name", "ZhangSimilarityScore");
	}

	PILISIdentification::PILISIdentification(const PILISIdentification& /*PILIS_id*/)
	{
		// TODO
	}
	
	PILISIdentification::~PILISIdentification()
	{
	}

	void PILISIdentification::setSequenceDB(PILISSequenceDB* sequence_db)
	{
		sequence_db_ = sequence_db;
	}

	void PILISIdentification::setModel(PILISModel* hmm_model)
	{
		hmm_model_ = hmm_model;
	}

	void PILISIdentification::getIdentifications(vector<Identification>& ids, const PeakMap& exp)
	{
		// get the parameters
		double pre_tol = (double)param_.getValue("prcr_m_tol");
		//double peak_tol = (double)param_.getValue("pk_m_tol");
		unsigned int max_candidates = (unsigned int)param_.getValue("max_candidates");
		//unsigned int hits = (unsigned int)param_.getValue("hits");
		String score_name = param_.getValue("score_name");

		// scoring
		CompareFunctor* scorer = Factory<CompareFunctor>::create(score_name);
			
		for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it)
		{
			double pre_pos = it->getPrecursorPeak().getPosition()[0];
			vector<PILISSequenceDB::PepStruct> cand_peptides;
			sequence_db_->getPeptides(cand_peptides, pre_pos - pre_tol, pre_pos + pre_tol);
			cerr << "#cand peptides: " << cand_peptides.size() << endl;
	
			HashMap<String, Size> sequence_to_charge;

			// get simple spectra for pre-eliminate most of the candidates
			TheoreticalSpectrumGenerator tsg;
			Identification pre_id;
			for (vector<PILISSequenceDB::PepStruct>::const_iterator it1 = cand_peptides.begin(); it1 != cand_peptides.end(); ++it1)
			{
				AASequence peptide_sequence(it1->peptide);
				// TODO parameter settings
				PeakSpectrum spec;
				tsg.getSpectrum(spec, peptide_sequence, it1->charge);
				double score = (*scorer)(*it, spec);
				PeptideHit peptide_hit(score, "Zhang", 0, it1->peptide);
				pre_id.insertPeptideHit(peptide_hit);

				sequence_to_charge[it1->peptide] = it1->charge;
			}

			pre_id.assignRanks();


			Identification id;
			for (Size i = 0; i < pre_id.getPeptideHits().size() && i < max_candidates; ++i)
			{
				String sequence = pre_id.getPeptideHits()[i].getSequence();
				AASequence peptide_sequence(sequence);
				PeakSpectrum spec;
				hmm_model_->getSpectrum(spec, peptide_sequence, sequence_to_charge[sequence]);
				
				// normalize the spectra and add intensity to too small peaks
				// TODO remove cheating
    		double max(0);
  			for (PeakSpectrum::ConstIterator it1 = spec.begin(); it1 != spec.end(); ++it1)
  			{
    			if (max < it1->getIntensity())
    			{
      			max = it1->getIntensity();
    			}
  			}

  			for (PeakSpectrum::Iterator it1 = spec.begin(); it1 != spec.end(); ++it1)
  			{
    			it1->setIntensity(it1->getIntensity()/max);
  			}

    		for (PeakSpectrum::Iterator it1 = spec.begin(); it1 != spec.end(); ++it1)
    		{
      		if (it1->getIntensity() < 0.01 && it1->getMetaValue("IonName") != "")
      		{
        		it1->setIntensity(0.01);
      		}
    		}

				double score = (*scorer)(*it, spec);
				PeptideHit peptide_hit(score, "Zhang", 0, sequence);
				id.insertPeptideHit(peptide_hit);
				//cerr << peptide_sequence << " " << score << endl;
			}

			id.assignRanks();
			ids.push_back(id);
		}

		return;
	}

	Param& PILISIdentification::getParam()
	{
		return param_;
	}

	void PILISIdentification::setParam(const Param& param)
	{
		param_ = param;
	}
}

