// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_SIMTYPES_H
#define OPENMS_SIMULATION_SIMTYPES_H

#include <vector>
#include <utility>
#include <map>
#include <utility>

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace OpenMS 
{
  /// Coordinate type in mz and rt dimension
  typedef Peak2D::CoordinateType SimCoordinateType;
  
	/// Abundance of proteins/peptides
	typedef Peak2D::IntensityType SimIntensityType;
	
  /// Charge of a peptide
  typedef Feature::ChargeType SimChargeType;
  
  /// Raw data point
	typedef RichPeak1D SimPointType;
	
	/// stores abundance information supported by the simulator
	typedef Map<String, SimIntensityType> FASTAEntryEnhanced;

	/// Container for FASTAEntry & addtional sim specific information
	typedef std::vector< std::pair<FASTAFile::FASTAEntry, FASTAEntryEnhanced> > SampleProteins;

	/// Sim FeatureMap
	typedef FeatureMap<> FeatureMapSim;

  /// Sim MSExperiment type
  typedef MSExperiment< SimPointType > MSSimExperiment;
  
	/// Probability of a modification to occur
	typedef Real ProbabilityType;

  /// A posttranslational modification
  struct PTM
  {
		/// Residue (modified)
    Residue residue;
    
    /// Probability to occur (0-1)
    ProbabilityType probability;

		/// Constructor
		PTM(const Residue p_residue, const ProbabilityType p_probability)
			: residue(p_residue),
				probability(p_probability)
		{}

  };  
 
	/// stores statistical properties of modifications on a single AA
	/// (all @p candidates will have the same AA residue)
	struct PTMRow
	{
		/// possible modifications (all affecting the same AA)
		std::vector<PTM> candidates;
		/// p that no modification will take place
		ProbabilityType probability_none;
		/// sum of all modifications' probabilites (can be >1)
		ProbabilityType probability_sum;
		/// number of AA affected
		Size rnd_amount;

		PTMRow()
			:candidates(),
			 probability_none(1),
			 probability_sum(0),
			 rnd_amount(0)
		{}

		/// add a candidate and update members
		void add(const PTM& ptm)
		{
			candidates.push_back( ptm );
			probability_none *= 1-ptm.probability;
			probability_sum += ptm.probability;
		}

		/// draw one of the candidates according to their p's (sum normalized to one)
		Residue draw(const gsl_rng* rnd_gen)
		{
			double p_random = gsl_ran_flat(rnd_gen,0.0,1.0);
			// normalize to one
			p_random *= probability_sum;
			// the candidate that covers the rnd number wins
			double left=0,right=0;
			for (Size i=0;i<candidates.size();++i)
			{
				left=right;
				right+=candidates[i].probability;
				if (left<=p_random && p_random<=right) return candidates[i].residue;
			}
			return candidates.back().residue;
		}

	};

	/// mapping from AA-one letter code to list of possible PTM's
	typedef Map<String, PTMRow> PTMTable;

}

#endif
