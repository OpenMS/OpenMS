// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_PILISCROSSVALIDATION_H
#define OPENMS_ANALYSIS_ID_PILISCROSSVALIDATION_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <iostream>
#include <vector>

namespace OpenMS
{
	class PeakSpectrumCompareFunctor;
	class PILISModel;

	class OPENMS_DLLAPI PILISCrossValidation : public DefaultParamHandler
	{

		public:

		/** @brief this struct represents a peptide spectrum pair
		*/
		struct Peptide
		{
			Peptide()
				: charge(0)
			{

			}

			Peptide(const Peptide& rhs)
				:	sequence(rhs.sequence),
					charge(rhs.charge),
					spec(rhs.spec),
					hits(rhs.hits)
			{
			}

			virtual ~Peptide()
			{
			}

			Peptide& operator = (const Peptide& rhs)
			{
				if (&rhs != this)
				{
					sequence = rhs.sequence;
					charge = rhs.charge;
					spec = rhs.spec;
					hits = rhs.hits;
				}
				return *this;
			}

  		AASequence sequence;
  		Int charge;
  		RichPeakSpectrum spec;

  		std::vector<PeptideHit> hits;

  		bool operator < (const Peptide& peptide) const
  		{
    		if (sequence < peptide.sequence)
    		{
      		return true;
    		}
    		else
    		{
      		if (sequence == peptide.sequence)
      		{
        		return charge < peptide.charge;
      		}
    		}
    		return false;
  		}

		};


		struct Option
		{
  		enum Type
  		{
    		INT = 0,
    		DOUBLE = 1,
    		BOOL = 2,
    		STRINGLIST = 3
  		};

		  Option()
    		: type(INT),
      		int_min(0),
      		int_max(0),
      		int_stepsize(0),
      		dbl_min(0),
      		dbl_max(0),
      		dbl_stepsize(0)
  		{
  		}

			Option(const Option& rhs)
				: type(rhs.type),
					int_min(rhs.int_min),
					int_max(rhs.int_max),
					int_stepsize(rhs.int_stepsize),
					dbl_min(rhs.dbl_min),
					dbl_max(rhs.dbl_max),
					dbl_stepsize(rhs.dbl_stepsize)
			{
			}

		  Option(Type t, double min, double max, double stepsize)
  		{
    		type = t;
    		if (type == INT)
    		{
      		int_min = (int)min;
      		int_max = (int)max;
      	int_stepsize = (int)stepsize;
   			}
    		else
    		{
      		if (type == DOUBLE)
      		{
        		dbl_min = min;
        		dbl_max = max;
        		dbl_stepsize = stepsize;
      		}
      		else
     			{
        		std::cerr << "Type: " << t << " is not known!" << std::endl;
      		}
    		}
  		}

  		Option& operator = (const Option& rhs)
 			{
    		if (&rhs != this)
    		{
      		type = rhs.type;
      		int_min = rhs.int_min;
      		int_max = rhs.int_max;
      		int_stepsize = rhs.int_stepsize;
      		dbl_min = rhs.dbl_min;
      		dbl_max = rhs.dbl_max;
      		dbl_stepsize = rhs.dbl_stepsize;
    		}
    		return *this;
  		}

  		Type type;
  		int int_min;
  		int int_max;
  		int int_stepsize;
  		double dbl_min;
  		double dbl_max;
  		double dbl_stepsize;
		};



		/** @brief Implementation of a cross valdidation training for the PILIS model

				This class serves as an implementation of a cross validation training for
				the PILIS model. It includes a range of parameters which can be set to 
				perform a GridSearch additionally.
		*/
		PILISCrossValidation();

		PILISCrossValidation(const PILISCrossValidation& rhs);

		virtual ~PILISCrossValidation();

		PILISCrossValidation& operator = (const PILISCrossValidation& rhs);

		void setOptions(const Map<String, Option>& rhs)
		{
			cv_options_ = rhs;
		}

		void setOption(const String& name, const Option& option)
		{
			cv_options_[name] = option;
		}

		void apply(Param& PILIS_param, const PILISModel& base_model, const std::vector<Peptide>& peptides);

		double scoreHits(const std::vector<std::vector<std::vector<RichPeakSpectrum> > >& sim_spectra, const std::vector<std::vector<RichPeakSpectrum> >& exp_spectra);

		protected:

		double scoreSpectra_(const RichPeakSpectrum& spec1, const RichPeakSpectrum& spec2);

		void partition_(std::vector<std::vector<Peptide> >& parts, const std::vector<Peptide>& source);

		void generateParameters_(const Param& param, const Map<String, Option>& options, std::vector<Param>& parameters);

		Map<String, Option> cv_options_;

		void updateMembers_();

		PeakSpectrumCompareFunctor* pscf_;
	};
}

#endif

