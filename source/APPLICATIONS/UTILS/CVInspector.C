// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_CVInspector CVInspector
	
	@brief A tool for visualization and validation of PSI mapping and CV files.
	
	This tool is used to validate the correct use of mapping files and CV files.
	
	It can also generate a HTML representation of mapping file and CV.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_CVInspector.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPCVInspector
	: public TOPPBase
{
 public:
	TOPPCVInspector()
		: TOPPBase("CVInspector","A tool for visualization and validation of PSI mapping and CV files.",false)
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFileList_("cv_files", "<files>", StringList(), "List of ontology files in OBO format.");
		registerStringList_("cv_names", "<names>", StringList(), "List of identifiers (one for each ontology file).");
		registerInputFile_("mapping_file", "<file>", "", "Mapping file in CVMapping (XML) format.", false);
		registerStringList_("ignore_cv", "<list>", StringList::create("UO,PATO,BTO"), "A list of CV identifiers which should be ignored.", false);
		registerOutputFile_("html","<file>", "", "Writes an HTML version of the mapping file with annotated CV terms", false);
	}
	
	void writeTermTree_(const String& accession, const ControlledVocabulary& cv, TextFile& file, UInt indent)
	{
		const ControlledVocabulary::CVTerm& term = cv.getTerm(accession);
		for (set<String>::const_iterator it=term.children.begin(); it!=term.children.end(); ++it)
		{
			const ControlledVocabulary::CVTerm& child_term = cv.getTerm(*it);
			String subterm_line;
			for (Size i=0; i<4*indent; ++i) subterm_line += "&nbsp;";
			String description = child_term.description;
			if (child_term.synonyms.size()!=0)
			{
				description += String(" -- Synonyms: '") + child_term.synonyms.concatenate(", ") + "'";
			}
			subterm_line += "- <span title=\"" + description + "\">" + child_term.id + " ! " + child_term.name + "</span>";
			StringList tags;
			if (child_term.obsolete)
			{
				tags.push_back("<font color=darkred>obsolete</font>");
			}
			if (child_term.xref_type!=ControlledVocabulary::CVTerm::NONE)
			{
				tags.push_back("value-type=" + ControlledVocabulary::CVTerm::getXRefTypeName(child_term.xref_type));
			}
			if (child_term.units.size()>0)
			{
				StringList units;
				for (set<String>::const_iterator u_it=child_term.units.begin(); u_it!=child_term.units.end(); ++u_it)
				{
					units.push_back(*u_it + "!" + cv.getTerm(*u_it).name);
				}
				tags.push_back(String("units=") + units.concatenate(","));
			}
			if (child_term.xref_binary.size()>0)
			{
				StringList types;
				for (StringList::const_iterator u_it=child_term.xref_binary.begin(); u_it!=child_term.xref_binary.end(); ++u_it)
				{
					types.push_back(*u_it + "!" + cv.getTerm(*u_it).name);
				}
				tags.push_back(String("binary-array-types=") + types.concatenate(","));
			}
			if (tags.size()!=0)
			{
				subterm_line += String("<FONT color=\"grey\"> (") + tags.concatenate(", ") + ")</FONT>";
			}
			file.push_back(subterm_line + "<BR>");
			writeTermTree_(child_term.id, cv, file, indent+1);
		}
	}
	
	ExitCodes main_(int , const char**)
	{
		StringList cv_files = getStringList_("cv_files");
		StringList cv_names = getStringList_("cv_names");
		if (cv_files.size()!=cv_names.size())
		{
			cerr << "Error: You have to specify an identifier for each CV file. Aborting!" << endl;
			return ILLEGAL_PARAMETERS;
		}
	
		// load cv terms
		ControlledVocabulary cv;
		for (Size i=0; i<cv_files.size(); ++i)
		{
			cv.loadFromOBO(cv_names[i], cv_files[i]);
		}
		Map<String, ControlledVocabulary::CVTerm> terms = cv.getTerms();
	
		// load mappings from mapping file
		String mapping_file = getStringOption_("mapping_file");
		CVMappings mappings;
		CVMappingFile().load(mapping_file, mappings);
		
		//store HTML version of mapping and CV
		if (getStringOption_("html")!="")
		{
			TextFile file;
			file.push_back("<HTML>");
			file.push_back("  <HEAD>");
			file.push_back("    <TITLE>CV mapping file</TITLE>");
			file.push_back("    <SCRIPT language=javascript type='text/javascript'>");
			file.push_back("      function toggleDiv(layer_ref,force_state) ");
			file.push_back("      {");
			file.push_back("        if (document.getElementById(layer_ref).style.display=='none' || force_state=='true')");
			file.push_back("        {");
			file.push_back("          document.getElementById(layer_ref).style.display = 'block';");
			file.push_back("        }");
			file.push_back("        else if (document.getElementById(layer_ref).style.display=='block' || force_state=='false')");
			file.push_back("        {");
			file.push_back("          document.getElementById(layer_ref).style.display = 'none';");
			file.push_back("        }");
			file.push_back("      }");
			file.push_back("    </SCRIPT>");
			file.push_back("  </HEAD>");
			file.push_back("  <BODY>");
			
			//count the number of terms and add button to expend/collaps all terms
			Int term_count = 0;
			for (vector<CVMappingRule>::const_iterator it = mappings.getMappingRules().begin(); it != mappings.getMappingRules().end(); ++it)
			{
				for (vector<CVMappingTerm>::const_iterator tit = it->getCVTerms().begin(); tit != it->getCVTerms().end(); ++tit)
				{
					++term_count;
				}
			}
			String expand_all = "    <a href=\"javascript:toggleDiv('div0','true')";
			String collapse_all = "    <a href=\"javascript:toggleDiv('div0','false')";
			for (Int i=1; i<term_count; ++i)
			{
				expand_all += String(";toggleDiv('div") + i + "','true')";
				collapse_all += String(";toggleDiv('div") + i + "','false')";
			}			
			file.push_back(expand_all + "\">Expand all</a><BR>");
			file.push_back(collapse_all + "\">Collapse all</a>");
			file.push_back("    <TABLE width=100% border=0>");
			term_count = -1;
			for (vector<CVMappingRule>::const_iterator it = mappings.getMappingRules().begin(); it != mappings.getMappingRules().end(); ++it)
			{
				//create rule line
				file.push_back("      <TR><TD colspan=\"2\"><HR></TD></TR>");
				file.push_back(String("      <TR><TD>Identifier:</TD><TD><B>") + it->getIdentifier() + "</B></TD></TR>");
				file.push_back(String("      <TR><TD>Element:</TD><TD><B>") + it->getElementPath() + "</B></TD></TR>");
				if (it->getRequirementLevel()==CVMappingRule::MUST)
				{
					file.push_back("      <TR><TD>Requirement level:</TD><TD><FONT color=\"red\">MUST</FONT></TD></TR>");
				}
				else if (it->getRequirementLevel()==CVMappingRule::SHOULD)
				{
					file.push_back("      <TR><TD>Requirement level:</TD><TD><FONT color=\"orange\">SHOULD</FONT></TD></TR>");
				}
				else if (it->getRequirementLevel()==CVMappingRule::MAY)
				{
					file.push_back("      <TR><TD>Requirement level:</TD><TD><FONT color=\"green\">MAY</FONT></TD></TR>");
				}
				if (it->getCombinationsLogic()==CVMappingRule::AND)
				{
					file.push_back("      <TR><TD>Combination logic:</TD><TD><FONT color=\"red\">AND</FONT></TD></TR>");
				}
				else if (it->getCombinationsLogic()==CVMappingRule::XOR)
				{
					file.push_back("      <TR><TD>Combination logic:</TD><TD><FONT color=\"orange\">XOR</FONT></TD></TR>");
				}
				else if (it->getCombinationsLogic()==CVMappingRule::OR)
				{
					file.push_back("      <TR><TD>Combination logic:</TD><TD><FONT color=\"green\">OR</FONT></TD></TR>");
				}
				
				//create table with terms
				for (vector<CVMappingTerm>::const_iterator tit = it->getCVTerms().begin(); tit != it->getCVTerms().end(); ++tit)
				{
					//create term line
					String term_line = String("      <TR><TD valign=\"top\">Term:</TD><TD>");
					if (tit->getAllowChildren())
					{
						++term_count;
						term_line += String("<a href=\"javascript:toggleDiv('div") + term_count + "','')\" style=\"text-decoration:none\" >+</a> ";
					}
					else
					{
						term_line += String("&nbsp;&nbsp;");
					}
					//add Term accession, name and description (as popup)
					if (cv.exists(tit->getAccession()))
					{
						const ControlledVocabulary::CVTerm& child_term = cv.getTerm(tit->getAccession());

						String description = child_term.description;
						if (child_term.synonyms.size()!=0)
						{
							description += String(" -- Synonyms: '") + child_term.synonyms.concatenate(", ") + "'";
						}
						term_line += "<span title=\"" + description + "\">";
					}
					term_line += tit->getAccession() + " ! " + tit->getTermName();
					if (cv.exists(tit->getAccession()))
					{
						term_line += "</span>";
						//check if term accession and term name correspond to the CV
						const ControlledVocabulary::CVTerm& main_term = cv.getTerm(tit->getAccession());
						if (main_term.name!=tit->getTermName())
						{
								cerr << "Warning: Accession '" << tit->getAccession() << "' and name '" << tit->getTermName() << "' do not match. Name should be '" << main_term.name << "'." << endl;
						}
					}
					//tags
					StringList tags;
					if (!tit->getUseTerm())
					{
						tags.push_back("children only");
					}
					if (tit->getIsRepeatable())
					{
						tags.push_back("repeatable");
					}
					if(cv.exists(tit->getAccession()))
					{
						const ControlledVocabulary::CVTerm& term = cv.getTerm(tit->getAccession());
						if (term.obsolete)
						{
							tags.push_back("<font color=darkred>obsolete</font>");
						}
						if (term.xref_type!=ControlledVocabulary::CVTerm::NONE)
						{
							tags.push_back("value-type=" + ControlledVocabulary::CVTerm::getXRefTypeName(term.xref_type));
						}
						if (term.units.size()>0)
						{
							StringList units;
							for (set<String>::const_iterator u_it=term.units.begin(); u_it!=term.units.end(); ++u_it)
							{
								units.push_back(*u_it + "!" + cv.getTerm(*u_it).name);
							}
							tags.push_back(String("units=") + units.concatenate(","));
						}
						if (term.xref_binary.size()>0)
						{
							StringList types;
							for (StringList::const_iterator u_it=term.xref_binary.begin(); u_it!=term.xref_binary.end(); ++u_it)
							{
								types.push_back(*u_it + "!" + cv.getTerm(*u_it).name);
							}
							tags.push_back(String("binary-array-types=") + types.concatenate(","));
						}
					}
					if (tags.size()!=0)
					{
						term_line += String("<FONT color=\"grey\"> (") + tags.concatenate(", ") + ")</FONT>";
					}
					file.push_back(term_line);
					
					// check whether we need the whole tree, or just the term itself
					if (tit->getAllowChildren())
					{
						file.push_back(String("        <div id=\"div") + term_count + "\" style=\"display: none\">");
						if (cv.exists(tit->getAccession()))
						{
							writeTermTree_(tit->getAccession(), cv, file, 1);
							//BEGIN - THIS IS NEEDED FOR WRITING PARSERS ONLY
							/*
							set<String> allowed_terms;
							cv.getAllChildTerms(allowed_terms, tit->getAccession());
							for (set<String>::const_iterator atit=allowed_terms.begin(); atit!=allowed_terms.end(); ++atit)
							{
								const ControlledVocabulary::CVTerm& child_term = cv.getTerm(*atit);
								String parser_string = String("os << \"&lt;cvParam cvRef=\\\"MS\\\" accession=\\\"") + child_term.id + "\\\" name=\\\"" + child_term.name + "\\\"";
								for (Size i=0; i<child_term.unparsed.size(); ++i)
								{
									//TODO this does not work anymore. The type is now stored as a member
									if (child_term.unparsed[i].hasSubstring("value-type:xsd\\:int") || child_term.unparsed[i].hasSubstring("value-type:xsd\\:float") || child_term.unparsed[i].hasSubstring("value-type:xsd\\:string"))
									{
										parser_string += " value=\\\"\" &lt;&lt; &lt;&lt; \"\\\"";
									}
								}
								parser_string += "/&gt;\\n\";<BR>";
								file.push_back(parser_string);
							}*/
						}
						else
						{
							file.push_back("          &nbsp;&nbsp;&nbsp;- Missing terms, CV not loaded...");
							cerr << "Warning: no child terms for " << tit->getAccession() << " found!" << endl;
						}
						file.push_back("          </div>");
						file.push_back("        </TD></TD></TR>");
					}
				}
			}
			file.push_back("    </TABLE>");
			file.push_back("  </BODY>");
			file.push_back("</HTML>");
			file.store(getStringOption_("html"));
			return EXECUTION_OK;
		}

		// iterator over all mapping rules and store the mentioned terms
		StringList ignore_namespaces = getStringList_("ignore_cv");
		set<String> ignore_cv_list;
		for (StringList::const_iterator it = ignore_namespaces.begin(); it != ignore_namespaces.end(); ++it)
		{
			ignore_cv_list.insert(*it);
		}
		set<String> used_terms;
		for (vector<CVMappingRule>::const_iterator it = mappings.getMappingRules().begin(); it != mappings.getMappingRules().end(); ++it)
		{
			set<String> allowed_terms;
			// iterate over all allowed terms
			for (vector<CVMappingTerm>::const_iterator tit = it->getCVTerms().begin(); tit != it->getCVTerms().end(); ++tit)
			{
				// check whether the term itself it allowed, or only its children
				if (tit->getUseTerm())
				{
					allowed_terms.insert(tit->getAccession());
				}
							
				// check whether we need the whole tree, or just the term itself
				if (tit->getAllowChildren())
				{
					// check whether we want to ignore this term
					if (!(tit->getAccession().has(':') && ignore_cv_list.find(tit->getAccession().prefix(':')) != ignore_cv_list.end()))
					{
						cv.getAllChildTerms(allowed_terms, tit->getAccession());
					}

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



