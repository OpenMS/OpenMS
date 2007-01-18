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


#ifndef OPENMS_ANALYSIS_ID_PILISSEQUENCEDB_H
#define OPENMS_ANALYSIS_ID_PILISSEQUENCEDB_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>

#include <vector>
#include <limits>

namespace OpenMS
{
	/**
	  @brief Simple sequence data base class
	*/
	class PILISSequenceDB
	{

		public:

			/// contains a peptide sequence, the charge and the weight of the peptide
			struct PepStruct
			{
			  String peptide;
			  float weight;
			  short charge;
			};

			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			PILISSequenceDB();
			
			/// copy constructor
			PILISSequenceDB(const PILISSequenceDB&);
			
			/// destructor
			virtual ~PILISSequenceDB();
			//@}
			
			/// assignment operator
			PILISSequenceDB& operator = (const PILISSequenceDB& rhs);

			/** Accessors
			*/
			//@{
			/// add peptides from a fasta file
			void addFASTAFile(const String& filename);

			/** @brief adds peptides from a file

					The file should have the following format:
					Three columns filled with <code> sequence weight charge </code>
			*/
			void addPeptidesFromFile(const String& filename);

			///
			void digestProteinsTryptic(Size missed_cleavages = 0);

			/// deletes all proteins
			void clearProteins();

			/// deletes all peptides
			void clearPeptides();

			// range query (m/z range); all peptides in the range [range_start, range_stop] are copied into peptides
			void getPeptides(std::vector<PepStruct>& peptides, double range_start = 0, double range_stop = std::numeric_limits<double>::max());

			/// granularity factor of the peptides (binsize)
			void setFactor(double factor);

			/// returns the granularity factor
			double getFactor() const;

			/// X (which stand for leucine or isoleucine) and L (leucine) are replaced by and I (isoleucine)
			void setReplaceXandL(bool replace = true);

			/// return true if X and L are replaced
			bool isReplaceXandL() const;

			/// returns the number of peptides in the database
			unsigned int countPeptides() const;

			/// returns the number of proteins in the database
			unsigned int countProteins() const;
			
			/// returns true if the peptide is contained in the database
			bool has(const String& peptide) const;
			//@}

		protected:

			/// performs a tryptic digest
			void digestTryptic_(const String& seq, std::vector<String>& peptides, Size missedcleavages);

			/// adds a peptide into the database if not already contained
			void addPeptide_(const String& peptide);
			
			std::vector<std::pair<String, String> > proteins_;

			HashMap<Size, std::vector<PepStruct> > peptides_;

			double factor_;

			bool replace_X_and_L_;
	};
}

#endif

