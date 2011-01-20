// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H
#define OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <string>
#include <vector>

namespace OpenMS
{
	/**
		 @brief Class for the enzymatic digestion of proteins
		 
     Digestion can be performed using simple regular expressions, 
     e.g. [KR] | [^P]
     for trypsin. Also missed cleavages can be modelled, i.e. adjacent peptides are not cleaved
     due to enzyme malfunction/access restrictions. If @em n missed cleavages are given, all possible resulting
     peptides (cleaved and uncleaved) with up to @em n missed cleavages are returned.
     Thus @b no random selection of just @em n specific missed cleavage sites is performed.

     An alternative model is also available, where the protein is cleaved only at positions where a cleavage model
     trained on real data, exceeds a certain threshold. The model is published in 
     Siepen et al. (2007), "Prediction of missed cleavage sites in tryptic peptides aids protein identification in proteomics.", doi: 10.1021/pr060507u
     The model is only available for trypsin and ignores the missed cleavage setting. You should however use setLogThreshold()
     to adjust FP vs FN rates. A higher threshold increases the number of cleavages predicted.

		 @ingroup Chemistry
	*/
	class OPENMS_DLLAPI EnzymaticDigestion
	{
		public: 
			/// Possible enzymes for the digestion (adapt NamesOfEnzymes & getEnzymeByName() & nextCleavageSite_() if you add more enzymes here)
			enum Enzyme
			{
				TRYPSIN,
				SIZE_OF_ENZYMES
			};
			
			/// Names of the Enzymes
			static const std::string NamesOfEnzymes[SIZE_OF_ENZYMES];
			
			/// Default constructor
			EnzymaticDigestion();
			
			/// Returns the number of missed cleavages for the digestion
			SignedSize getMissedCleavages() const;

			/// Sets the number of missed cleavages for the digestion (default is 0). This setting is ignored when log model is used.
			void setMissedCleavages(SignedSize missed_cleavages);	

			/// Returns the enzyme for the digestion
			Enzyme getEnzyme() const;
			
			/// Sets the enzyme for the digestion (default is TRYPSIN).
			void setEnzyme(Enzyme enzyme);		

			/// convert enzyme string name to enum
			/// returns SIZE_OF_ENZYMES if @p name is not valid
			Enzyme getEnzymeByName(const String& name);

			/// Performs the enzymatic digestion of a protein.
			void digest(const AASequence& protein, std::vector<AASequence>& output);
			
			/// Returns the number of peptides a digestion of @p protein would yield.
			Size peptideCount(const AASequence& protein);

      /// use trained model when digesting?
			bool isLogModelEnabled() const;

      /// enables/disabled the trained model
			void setLogModelEnabled(bool enabled);
		
      /// Returns the threshold which needs to be exceeded to call a cleavage (only for the trained cleavage model on real data)
			DoubleReal getLogThreshold() const;
			
			/// Sets the threshold which needs to be exceeded to call a cleavage (only for the trained cleavage model on real data)
      /// Default is 0.25
			void setLogThreshold(DoubleReal threshold);	    

		protected:
			/// Number of missed cleavages
			SignedSize missed_cleavages_;
			/// Used enzyme
			Enzyme enzyme_;
      /// use the log model or naive digestion (with missed cleavages)
      bool use_log_model_;
      /// Threshold to decide if position is cleaved or missed (only for the model)
      DoubleReal log_model_threshold_;

      // define a binding site by position and AA
      struct BindingSite
      {
        Size position;
        String AAname;
        
        BindingSite ()
          : position(), AAname() {}

        BindingSite (const Size& p, const String& name)
          : position(p), AAname(name) {}

        bool operator < (const BindingSite& rhs) const
        {
          return (position < rhs.position) || ((position == rhs.position) && (AAname < rhs.AAname));
        }

        bool operator == (const BindingSite& rhs) const
        {
          return position==rhs.position && AAname==rhs.AAname;
        }

      };

      // define the log likelihood for missed and cleavage model
      struct CleavageModel
      {
        DoubleReal p_cleave;
        DoubleReal p_miss;

        CleavageModel()
          : p_cleave(0), p_miss(0) {}
        CleavageModel (const DoubleReal& p_c, const DoubleReal& p_m)
          : p_cleave(p_c), p_miss(p_m) {}
      };

      /// Holds the cleavage model
      Map<BindingSite, CleavageModel> model_data_;
			
			///moves the iterator @p it after the next cleavage site of the @p sequence
			void nextCleavageSite_(const AASequence& sequence, AASequence::ConstIterator& iterator);
	};

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H
