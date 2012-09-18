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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_SemanticValidator SemanticValidator
	
	@brief SemanticValidator for analysisXML and mzML files.
	
	This util is able to validate analysisXML and mzML files
	using an instance document and a mapping file.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_SemanticValidator.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude UTILS_SemanticValidator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSemanticValidator
	: public TOPPBase
{
 public:
	TOPPSemanticValidator()
		: TOPPBase("SemanticValidator","SemanticValidator for analysisXML and mzML files.",false)
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "Input file, either analysisXML or mzML.");
		registerInputFile_("mapping_file", "<file>", "", "Mapping file which is used to semantically validate the given XML file against this mapping file (see 'share/OpenMS/MAPPING' for templates).");
	}	
	
	ExitCodes main_(int , const char**)
	{
		String in_file = getStringOption_("in");
		String mapping_file = getStringOption_("mapping_file");

		CVMappings mappings;
		CVMappingFile().load(mapping_file, mappings, false);

		ControlledVocabulary cv;
		cv.loadFromOBO("PSI-MOD", File::find("/CHEMISTRY/PSI-MOD.obo"));
	  cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
		cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
		cv.loadFromOBO("brenda",File::find("/CV/brenda.obo"));
		cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));
		cv.loadFromOBO("UNIMOD",File::find("/CV/unimod.obo"));
		cv.loadFromOBO("PSI-MS",File::find("/CV/psi-ms.obo"));

		// check cv params
		Internal::SemanticValidator semantic_validator(mappings, cv);
		semantic_validator.setCheckTermValueTypes(true);
		semantic_validator.setCheckUnits(true);

		StringList errors, warnings;

		/*bool valid =*/ semantic_validator.validate(in_file, errors, warnings);
    for (Size i=0; i<warnings.size(); ++i)
    {
    	cout << "Warning: " << warnings[i] << endl;
    }
    for (Size i=0; i<errors.size(); ++i)
    {
      cout << "Error: " << errors[i] << endl;
    }

		if (warnings.empty() && errors.empty())
		{
			cout << "Congratulations, the file is valid!" << endl;
		}

		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPSemanticValidator tool;
	return tool.main(argc,argv);
}

/// @endcond



