// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_MODIFICATIONSDB_H
#define OPENMS_CHEMISTRY_MODIFICATIONSDB_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <set>

namespace OpenMS
{
	// forward declarations
	class ResidueModification;
	class Residue;

	/** @ingroup Chemistry
	
			@brief database which holds all residue modifications from UniMod
			
			This singleton class serves as a storage of the available modifications
			represented by UniMod (www.unimod.org). The modifications are identified
			by there name and possibly other ids from UniMod or the PSI-MOD ontology.
			Modifications can have different specificities, e.g. they can occur only
			at the termini, anywhere or only at specific amino acids.

			The modifications are defined in share/OpenMS/CHEMISTRY/unimod.xml and 
			in share/OpenMS/CHEMISTRY/PSI-MOD.obo. The unimod file can be directly
			downloaded from unimod.org and replaced if the modifications change.

			To add a new modification, not contained in UniMod, one should follow
			the way described at the unimod.org website and download the file then
			from unimod.org. The same can be done to add support for the modifications
			to search engines, e.g. Mascot.
	*/
	class OPENMS_DLLAPI ModificationsDB
	{					
		public:

			inline static ModificationsDB* getInstance()
      {
        static ModificationsDB* db_ = 0;
        if (db_ == 0)
        {
          db_ = new ModificationsDB;
        }
        return db_;
      }

			/// returns the number of modifications read from the unimod.xml file
			Size getNumberOfModifications() const;

			/// returns the modification with the given index 
			const ResidueModification& getModification(Size index) const;

			/// returns all modifications which have the given name as synonym
			void searchTerminalModifications(std::set<const ResidueModification*>& mods, const String& name, ResidueModification::Term_Specificity term_spec) const;

			/// returns all modification which have the given name as synonym and the given origin
			void searchModifications(std::set<const ResidueModification*>& mods, const String& orgin, const String& mod_name, ResidueModification::Term_Specificity term_spec) const;
			
      /// returns all modification which have the given name as synonym
      void searchModifications(std::set<const ResidueModification*>& mods, const String& mod_name, ResidueModification::Term_Specificity term_spec) const;

			/** @brief returns the modifications of the given name

					This can either be the PSI-MOD identifier or every other unique
					identifier which can be found in the PSI-MOD definitions file.

					To search for more than one modification searchModifications() can be used!
					
					@exception ElementNotFound is thrown if no or more than one element is found
			*/
			const ResidueModification& getTerminalModification(const String& name, ResidueModification::Term_Specificity term_spec) const;

			/// returns the modification with the given name and given residue 
			const ResidueModification& getModification(const String& residue_name, const String& mod_name, ResidueModification::Term_Specificity term_spec) const;

			const ResidueModification& getModification(const String& modification) const;

			/// returns the index of the modification in the mods_ vector; a unique name must be given
			Size findModificationIndex(const String& mod_name) const;

			/// query the modifications DB to get the terminal modifications with mass
			void getTerminalModificationsByDiffMonoMass(std::vector<String>& mods, DoubleReal mass, DoubleReal error, ResidueModification::Term_Specificity term_spec);

			/// query the modifications DB to get the modifications with mass, without any specific origin
			void getModificationsByDiffMonoMass(std::vector<String>& mods, DoubleReal mass, DoubleReal error = 0.0);
			
			/// query the modifications DB to get modifications with the given mass at the given residue
			void getModificationsByDiffMonoMass(std::vector<String>& mods, const String& residue, DoubleReal mass, DoubleReal error = 0.0);

			/// adds modifications from a given file in OBO format
			void readFromOBOFile(const String& filename);

			/// adds modifications from a given file in Unimod XML format
			void readFromUnimodXMLFile(const String& filename);
			
			/// get all modifications that can be used for identification searches
			void getAllSearchModifications(std::vector<String>& modifications);

      /// print information stored in modification db
       void print();
		protected:

			/// stores the modifications
			std::vector<ResidueModification*> mods_;

			/// stores the mappings of (unique) names to the modifications
			Map<String, std::set<const ResidueModification*> > modification_names_;
			
			
		private:

			/** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      ModificationsDB();

      ///copy constructor
      ModificationsDB(const ModificationsDB& residue_db);

      /// destructor
      virtual ~ModificationsDB();
      //@}

      /** @name Assignment
      */
      //@{
      /// assignment operator
      ModificationsDB& operator = (const ModificationsDB& aa);
      //@}		
	};
}
#endif
