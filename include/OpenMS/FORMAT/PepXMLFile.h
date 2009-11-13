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
// $Maintainer: Chris Bielow, Hendrik Weisser $
// $Authors: Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEPXMLFILE_H
#define OPENMS_FORMAT_PEPXMLFILE_H

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>
#include <map>
#include <set>


namespace OpenMS 
{
  /**
    @brief Used to load and store PepXML files
    
    This class is used to load and store documents that implement the schema of PepXML files.
		
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI PepXMLFile
  	: protected Internal::XMLHandler,
  		public Internal::XMLFile
  {
		public:
		
			/// Constructor
			PepXMLFile();

			/// Destructor
			virtual ~PepXMLFile();
			
			/**
				@brief Loads peptide sequences with modifications out of a PepXML file
				
				@param filename PepXML file to load
				@param proteins Protein identification output
				@param peptides Peptide identification output
				@param experiment_name Experiment file name, which is used to extract the corresponding search results from the PepXML file. 
				@param experiment MS run to extract the retention times from (PepXML contains only scan numbers). If the experiment is empty, it is read from @a experiment_name.
				
				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename, std::vector<ProteinIdentification>& proteins, std::vector<PeptideIdentification>& peptides, const String& experiment_name, MSExperiment<>& experiment);
									
			/**
				@brief @a load function with empty defaults for some parameters (see above)
				
				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename, std::vector<ProteinIdentification>& proteins, std::vector<PeptideIdentification>& peptides, const String& experiment_name = "");

			/**
				@brief Stores idXML as PepXML file

				@exception Exception::UnableToCreateFile is thrown if the file could not be opened for writing
			*/
			void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids);

  	protected:
		
			// Docu in base class
			virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
			virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

		private:
			
		  void matchModification_(DoubleReal mass, String& modification_description, const String& origin);
	
			struct AminoAcidModification
			{
				String aminoacid;
				String massdiff;
				DoubleReal mass;
				bool variable;
				String description;
				String terminus;

				AminoAcidModification()
					: mass(0),
						variable(false)
				{
				}

				AminoAcidModification(const AminoAcidModification& rhs)
					:	aminoacid(rhs.aminoacid),
						massdiff(rhs.massdiff),
						mass(rhs.mass),
						variable(rhs.variable),
						description(rhs.description),
						terminus(rhs.terminus)
				{
				}

				virtual ~AminoAcidModification()
				{
				}

				AminoAcidModification& operator = (const AminoAcidModification& rhs)
				{
					if (this != &rhs)
					{
						aminoacid = rhs.aminoacid;
						massdiff = rhs.massdiff;
						mass = rhs.mass;
						variable = rhs.variable;
						description = rhs.description;
						terminus = rhs.terminus;
					}
					return *this;
				}

			};
	
			/// Pointer to the list of identified proteins
			std::vector<ProteinIdentification>* proteins_;
						
			/// Pointer to the list of identified peptides
			std::vector<PeptideIdentification>* peptides_;
			
			/// Pointer to the experiment from which the pepXML file was generated
			const MSExperiment<>* experiment_;

			/// Name of the associated experiment (filename of the data file, extension will be removed)
			String exp_name_;	

			/// Pointer to the mapping between scan number in the pepXML file and index in the corresponding MSExperiment
			std::map<Size, Size>* scan_map_;

			/// Retention time and mass-to-charge tolerance
			DoubleReal rt_tol_, mz_tol_;

			/// Hydrogen data (for mass types)
			Element hydrogen_;

			/// Do current entries belong to the experiment of interest (for pepXML files that bundle results from different experiments)?
			bool wrong_experiment_;

			/// References to currently active ProteinIdentifications
			std::vector<std::vector<ProteinIdentification>::iterator> current_proteins_;
						
			/// Search parameters of the current identification run
			ProteinIdentification::SearchParameters params_;

			/// Enyzme associated with the current identification run
			ProteinIdentification::DigestionEnzyme enzyme_;

			/// PeptideIdentification instance currently being processed
			PeptideIdentification current_peptide_;

			/// PeptideHit instance currently being processed
			PeptideHit peptide_hit_;

			/// Sequence of the current peptide hit				
			String current_sequence_;
			
			/// RT and m/z of current PeptideIdentification
			DoubleReal rt_, mz_;
						
			/// Precursor ion charge
			Int charge_;
	
			/// ID of current search result
			Size search_id_;

			/// Identifier linking PeptideIdentifications and ProteinIdentifications
			String prot_id_;

			/// Date the pepXML file was generated
			DateTime date_;
		
			/// Mass of a hydrogen atom (monoisotopic/average depending on case)
			DoubleReal hydrogen_mass_;

			/// The modifications of the current peptide hit (position is 1-based)
			std::vector<std::pair<String, Size> > current_modifications_;
			
			/// Fixed aminoacid modifications
			std::vector<AminoAcidModification> fixed_modifications_;	
			
			/// Variable aminoacid modifications
			std::vector<AminoAcidModification> variable_modifications_;
		
			//@}
									
	};
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_PEPXMLFILE_H
