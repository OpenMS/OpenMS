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
// $Maintainer: Sandro Andreotti $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEPNOVOINFILE_H
#define OPENMS_FORMAT_PEPNOVOINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <map>


namespace OpenMS
{
	/**
		@brief PepNovo input file adapter.
		
		Creates a PepNovo_PTMs.txt file for PepNovo search.
  	  	
  	@ingroup FileIO
	*/
  class OPENMS_DLLAPI PepNovoInfile
  {
		public:
			/// default constructor
			PepNovoInfile();

			/// copy constructor
			PepNovoInfile(const PepNovoInfile& pepnovo_infile);

			/// destructor
			virtual ~PepNovoInfile();

			/// assignment operator
			PepNovoInfile& operator=(const PepNovoInfile& pepnovo_infile);

			/// equality operator
			bool operator==(const PepNovoInfile& pepnovo_infile) const;

			/** stores the experiment data in a PepNovo input file that can be used as input for PepNovo shell execution
					
				@param filename the file which the input file is stored into
				@throw Exception::UnableToCreateFile is thrown if the given file could not be created
			*/
			void store(const String& filename);

			/** @brief generates the PepNovo Infile for given fixed and variable modifications			 *
			 *
			 * @param fixed_mods StringList of fixed modifications unique identifiers
			 * @param variable_mods StringList of variable modifications unique identifiers
			 */
			void setModifications(const StringList &fixed_mods, const StringList &variable_mods);

			/** @brief return the modifications.
			 *
			 *  the modification unique identifiers are mapped to the keys used
			 *  in the PepNovo Infile (origin+rounded monoisotopic mass of modification ).
			 *  (e.g. modification_key_map["K+16"]=="Oxidation (K)" )
			 */
			void getModifications(std::map<String,String>& modification_key_map) const;

		private:
			ModificationDefinitionsSet mods_;
			std::map<String,String>mods_and_keys_;
			TextFile ptm_file_;


		 /** retrieves the name of modification, and generates the corresponding line for the
			 PepNovo infile.
			 @param modification the modification
			 @param variable should be set to true if it variable
		*/
		  String handlePTMs_(const String &modification, const bool variable);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEPNOVOINFILE_H
