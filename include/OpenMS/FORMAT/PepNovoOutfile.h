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

#ifndef OPENMS_FORMAT_PEPNOVOOUTFILE_H
#define OPENMS_FORMAT_PEPNOVOOUTFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <map>

namespace OpenMS
{
	class ProteinIdentification;

	/**
		@brief Representation of a PepNovo output file
	
		This class serves to read in a PepNovo outfile. The information can be
		retrieved via the load function.

		@ingroup FileIO
	*/
	class OPENMS_DLLAPI PepNovoOutfile
	{
		public:
			/// Constructor
			PepNovoOutfile();

			/// copy constructor
			PepNovoOutfile(const PepNovoOutfile& pepnovo_outfile);

			/// destructor
			virtual ~PepNovoOutfile();

			/// assignment operator
			PepNovoOutfile& operator=(const PepNovoOutfile& pepnovo_outfile);

			/// equality operator
			bool operator==(const PepNovoOutfile& pepnovo_outfile) const;

			 /**
				@brief loads data from a PepNovo outfile

				@param result_filename the file to be loaded
				@param peptide_identifications the peptide identifications
				@param protein_identification the protein identification
				@param score_threshold cutoff threshold for the PepNovo score (PnvScr)
				@param id_rt_mz map the spectrum identifiers returned by PepNovo
				to the rt and mz values of the spectrum (used to map the identifications back to the spectra). key= &lt;PepNovo Id&gt;, value= &lt;pair&lt;rt,mz&gt; &gt;.
				For spectra not present in this map identifications cannot be mapped back.
				@param mod_id_map map the OpenMS id for modifications (FullId) to the ids returned by PepNovo key= &lt;PepNovo_key&gt;, value= &lt;OpenMS FullId&gt;
			*/
			void load(const std::string& result_filename, std::vector< PeptideIdentification >& peptide_identifications,
			    ProteinIdentification& protein_identification,
			    const DoubleReal& score_threshold,
			    const std::map< String, std::pair<DoubleReal,DoubleReal> >& id_rt_mz,
			    const std::map<String, String> &mod_id_map);

			/** @brief get the search engine version and search parameters from a PepNovo output file
			 *
			 * search parameters (precursor tolerance, peak mass tolerance, allowed modifications)are stored in the protein_identification.

				@param pepnovo_output_without_parameters_filename
				@param protein_identification
			*/
			void getSearchEngineAndVersion(const String& pepnovo_output_without_parameters_filename, ProteinIdentification& protein_identification);

	};

} //namespace OpenMS

#endif // OPENMS_FORMAT_PEPNOVOOUTFILE_H

