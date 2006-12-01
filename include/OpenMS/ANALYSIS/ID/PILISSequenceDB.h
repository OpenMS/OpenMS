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


#ifndef OPENMS_ANALYSIS_ID_PILISSEQUENCEDB_H
#define OPENMS_ANALYSIS_ID_PILISSEQUENCEDB_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>

#include <vector>
#include <limits>

namespace OpenMS
{
	class PILISSequenceDB
	{

		public:

			struct PepStruct
			{
			  String peptide;
			  float weight;
			  short charge;
			};

			PILISSequenceDB();
			
			PILISSequenceDB(const PILISSequenceDB&);
			
			virtual ~PILISSequenceDB();
			
			const PILISSequenceDB& operator = (const PILISSequenceDB&);

			/// I/O
			void addFASTAFile(const String& filename);

			// 
			void addPeptidesFromFile(const String& filename);

			void digestSequencesTryptic(Size missed_cleavages = 0);

			void clearProteins();

			void clearPeptides();

			// range query (m/z range)
			void getPeptides(std::vector<PepStruct>& peptides, double range_start = 0, double range_stop = std::numeric_limits<double>::max());

			// granularity factor of the peptides
			void setFactor(double factor);

			double getFactor() const;

			void setReplaceXandL(bool replace = true);

			bool isReplaceXandL() const;

			unsigned int countPeptides() const;

			bool has(const String& peptide);

		protected:

			void digest_(const String& sequence, Size missed_cleavages);

			std::vector<std::pair<String, String> > proteins_;

			HashMap<Size, std::vector<PepStruct> > peptides_;

			double factor_;

			bool replace_X_and_L_;
	};
}

#endif
