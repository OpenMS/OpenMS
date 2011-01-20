// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_ANALYSIS_TARGETED_INCLUSIONEXCLUSIONLIST_H
#define OPENMS_ANALYSIS_TARGETED_INCLUSIONEXCLUSIONLIST_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  /**
     @brief


		 @TODO allow modifications (fixed?)

     
  */
  class OPENMS_DLLAPI InclusionExclusionList 
  {
  public:
    /** @name Constructors and destructors
     */
    //@{
    /// default constructor
    InclusionExclusionList();
    
    
    //@}

//     void loadTargets(FeatureMap<>& map, std::vector<IncludeExcludeTarget>& targets,TargetedExperiment& exp);

//     void loadTargets(std::vector<FASTAFile::FASTAEntry>& fasta_entries, std::vector<IncludeExcludeTarget>& targets,
//                      TargetedExperiment& exp, Size missed_cleavages = 0);


		/**
			 @brief Writes inclusion or exclusion list of tryptic peptides of the given proteins (tab-delimited).

			 @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
		 */
    void writeTargets(std::vector<FASTAFile::FASTAEntry>& fasta_entries, String& out_path,IntList& charges,
											String rt_model_path,DoubleReal rel_rt_window_size,bool rt_in_seconds,Size missed_cleavages = 0);

		/**
			 @brief Writes inclusion or exclusion list of given feature map.

			 @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
		 */
    void writeTargets(FeatureMap<>& map,String& out_path,DoubleReal rel_rt_window_size,bool rt_in_seconds);
		
		/**
			 @brief Writes inclusion or exclusion list of given peptide ids (tab-delimited).

			 @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
			 @exception Exception::InvalidSize is thrown if a peptide id contains more than one hit
			 @exception Exception::MissingInformation is thrown if a peptide id contains no RT information
		 */
		void writeTargets(std::vector<PeptideIdentification>& pep_ids,String& out_path,DoubleReal rel_rt_window_size,IntList& charges,bool rt_in_seconds);
  };

}

#endif
