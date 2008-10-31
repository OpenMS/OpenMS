// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/CVMappings.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page CVInspector CVInspector
	
	@todo Docu (Andreas)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPCVInspector
	: public TOPPBase
{
 public:
	TOPPCVInspector()
		: TOPPBase("CVInspector","TODO")
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("OBO_file","<file>","","Ontology file in OBO format where the CV terms are stored");
		registerInputFile_("mapping_file","<file>","","Mapping file in CVMapping (XML) format, where the mappings are defined", false);
	}	
	
	ExitCodes main_(int , const char**)
	{
		String OBO_file = getStringOption_("OBO_file");
		String mapping_file = getStringOption_("mapping_file");
	
		// load cv terms
		ControlledVocabulary cv;
		cv.loadFromOBO("PSI-PI", OBO_file);
//		cv.loadFromOBO("PATO", "quality.obo");
//		cv.loadFromOBO("UO", "unit.obo");
//		cv.loadFromOBO("brenda", "brenda.obo");
//		cv.loadFromOBO("GO", "goslim_goa.obo");
		Map<String, ControlledVocabulary::CVTerm> terms = cv.getTerms();
		cerr << "Loaded " << terms.size() << " terms from file" << endl;
	
		// load mappings from mapping file
		CVMappings mappings;
		CVMappingFile().load(mapping_file, mappings);
		cerr << "Loaded " << mappings.getMappingRules().size() << " mapping rules" << endl;

		// iterator over all mapping rules and store the mentioned terms
		set<String> used_terms;
		for (vector<CVMappings::CVMappingRule>::const_iterator it = mappings.getMappingRules().begin(); it != mappings.getMappingRules().end(); ++it)
		{
			set<String> allowed_terms;
			// iterate over all allowed terms
			for (vector<CVMappings::CVTerm>::const_iterator tit = it->getCVTerms().begin(); tit != it->getCVTerms().end(); ++tit)
			{
				// check whether the term itself it allowed, or only its children
				if (tit->getUseTerm())
				{
					allowed_terms.insert(tit->getAccession());
				}
							
				// check whether we need the whole tree, or just the term itself
				if (tit->getAllowChildren())
				{
					cv.getAllChildTerms(allowed_terms, tit->getAccession());

					// also add the term itself to the used_terms, because all the children are allowed
					used_terms.insert(tit->getAccession());
				}
			}

			// print the allowed terms for the rule
			cout << "MappingRule: id=" << it->getIdentifier() << ", elementPath=" << it->getElementPath() << ", #terms=" << it->getCVTerms().size() << endl;
			for (set<String>::const_iterator ait = allowed_terms.begin(); ait != allowed_terms.end(); ++ait)
			{
				cout << *ait << " " << terms[*ait].name << endl;
			}
			used_terms.insert(allowed_terms.begin(), allowed_terms.end());
		}

		// find unused terms, which CANNOT be used in the XML due to the mapping file
		set<String> unused_terms;
		for (Map<String, ControlledVocabulary::CVTerm>::ConstIterator it = terms.begin(); it != terms.end(); ++it)
		{
			if (used_terms.find(it->first) == used_terms.end())
			{
				unused_terms.insert(it->first);
			}
		}

		cout << "\n\nCVTerms which are unused in the mapping file and therefore MUST NOT be used in an instance document" << endl;
		for (set<String>::const_iterator it = unused_terms.begin(); it != unused_terms.end(); ++it)
		{
			cout << *it << " " << terms[*it].name;

			// print also parent names
			for (set<String>::const_iterator pit = terms[*it].parents.begin(); pit != terms[*it].parents.end(); ++pit)
			{
				cout << " " << terms[*pit].id << " " << terms[*pit].name;
			}
			cout << endl;
		}
		
		
		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPCVInspector tool;
	return tool.main(argc,argv);
}

/// @endcond



