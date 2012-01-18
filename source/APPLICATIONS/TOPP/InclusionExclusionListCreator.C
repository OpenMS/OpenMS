// -*- mode: C++; tab-width: 2; -*-
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
// $Authors: Alexandra Zerck, Chris Bielow$
// --------------------------------------------------------------------------

//#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/ANALYSIS/TARGETED/InclusionExclusionList.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_InclusionExclusionListCreator InclusionExclusionListCreator

   @brief A tool for creating inclusion and/or exclusion lists for LC-MS/MS.

	<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ InclusionExclusionListCreator \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> - </td>
		</tr>
		<tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
		</tr>
	</table>
	</CENTER>

   Currently this tool can create tab-delimited inclusion or exclusion lists (m/z, RT start, RT stop).
	 The input can either be peptide identifications from previous runs, a feature map or a FASTA-file with proteins.
	 Inclusion and exclusion charges can be specified for FASTA and IdXML input. If no charges are specified in the case of peptide id input, only
   the charge state of the peptide id is in/excluded, otherwise all given charge states are entered to the list.

   The RT window size can be specified in the RT section of the INI file, either as relative window
   with [rt-rel_rt_window_size*rt,rt+rel_rt_window_size*rt] or absolute window.

   The default is RT in minutes, but seconds can also be used (see INI file).

   <B>The command line parameters of this tool are:</B>
   @verbinclude TOPP_InclusionExclusionListCreator.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPInclusionExclusionListCreator
  : public TOPPBase
{
public:
  TOPPInclusionExclusionListCreator()
    : TOPPBase("InclusionExclusionListCreator","Creates inclusion and/or exclusion lists.")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("include", "<file>", "", "Inclusion list input file in FASTA or featureXML format.",false);
    setValidFormats_("include",StringList::create("featureXML,FASTA"));
    registerInputFile_("exclude","<file>", "", "Exclusion list input file in featureXML, IdXML or FASTA format.",false);
    setValidFormats_("exclude",StringList::create("featureXML,IdXML,FASTA"));
    registerOutputFile_("out", "<file>", "", "Output file (tab delimited).");
    registerInputFile_("rt_model","<file>","","RTModel file used for the rt prediction of peptides in FASTA files.",false);
    //in FASTA or featureXML
    registerIntList_("inclusion_charges","<charge>",IntList(),"List containing the charge states to be considered for the inclusion list compounds, space separated.",false);
    setMinInt_("inclusion_charges", 1);
    registerIntList_("exclusion_charges","<charge>",IntList(),"List containing the charge states to be considered for the exclusion list compounds (for idXML and FASTA input), space separated.",false);
    setMinInt_("exclusion_charges", 1);

    //    setValidFormats_("out", StringList::create("TraML"));

    registerSubsection_("algorithm","Inclusion/Exclusion algorithm section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const
  {
    // there is only one subsection: 'algorithm' (s.a) .. and in it belongs the InclusionExclusionList param
    InclusionExclusionList fdc;
    Param tmp;
    tmp.insert("InclusionExclusionList:",fdc.getParameters());
    return tmp;
  }

  ExitCodes main_(int , const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input/output files
    String include(getStringOption_("include"));
    String exclude(getStringOption_("exclude"));
    String out(getStringOption_("out"));

    if (include == "" && exclude == "")
    {
      writeLog_("Error: No input file given.");
      return MISSING_PARAMETERS;
    }
    // currently we can handle only inclusion OR exclusion, will be possible with the traML output
    if (include != "" && exclude != "")
    {
      writeLog_("Error: Currently only inclusion OR exclusion, both will be possible with the traML output coming soon");
      return ILLEGAL_PARAMETERS;
    }

		IntList incl_charges(getIntList_("inclusion_charges"));
    IntList excl_charges(getIntList_("exclusion_charges"));
    String rt_model_file(getStringOption_("rt_model"));

    //-------------------------------------------------------------
    // loading input: inclusion list part
    //-------------------------------------------------------------

		FileHandler fh;
    TargetedExperiment exp;
    Param const& iel_param = getParam_().copy("algorithm:InclusionExclusionList:",true);
    writeDebug_("Parameters passed to InclusionExclusionList", iel_param, 3);

    InclusionExclusionList list;
    list.setParameters(iel_param);


    std::cout << "\n\n\n\n" << iel_param.getValue("RT:unit") << "\n\n";


    if (include != "")
    {
      FileTypes::Type in_type = fh.getType(include);
      std::vector<IncludeExcludeTarget> incl_targets;
      if(in_type == FileTypes::FEATUREXML)
      {
        // load feature map
        FeatureMap<> map;
        FeatureXMLFile().load(include,map);

        if(!incl_charges.empty())
        {
          writeLog_("Warning: 'inclusion_charges' parameter is not honored for featureXML input.");
					return ILLEGAL_PARAMETERS;
        }

        // convert to targeted experiment
				// for traML output
				//            list.loadTargets(map,incl_targets,exp);
				// for tab-delimited output
				try
				{
					list.writeTargets(map,out);
				}
				catch(Exception::UnableToCreateFile)
        {
          writeLog_("Error: Unable to create output file.");
          return CANNOT_WRITE_OUTPUT_FILE;
        }
      }
      else // FASTA format
      {
        if (!File::exists(rt_model_file))
        {
          writeLog_("Error: RT model file required for FASTA input to predict RT elution time.");
          return MISSING_PARAMETERS;
        }
        if(incl_charges.empty())
        {
          writeLog_("Error: Protein sequences for inclusion given, but no charge states specified.");
          return MISSING_PARAMETERS;
        }
        std::vector<FASTAFile::FASTAEntry> entries;
        // load fasta-file
        FASTAFile().load(include,entries);
        // convert to targeted experiment
				// if traML output
        //list.loadTargets(entries,incl_targets,exp,missed_cleavages);
				// if tab-delimited output
				try
				{
					list.writeTargets(entries,out,incl_charges,rt_model_file);
				}
				catch(Exception::UnableToCreateFile)
        {
          writeLog_("Error: Unable to create output file.");
          return CANNOT_WRITE_OUTPUT_FILE;
				}
      }

			//        exp.setIncludeTargets(incl_targets);
    }
    //-------------------------------------------------------------
    // loading input: exclusion list part
    //-------------------------------------------------------------
    if(exclude != "")
    {
      FileTypes::Type ex_type = fh.getType(exclude);
      //        std::vector<IncludeExcludeTarget> excl_targets;
      if(ex_type == FileTypes::FEATUREXML)
      {
        if(!excl_charges.empty())
        {
          writeLog_("Warning: 'exclusion_charges' parameter is not honored for featureXML input.");
          return ILLEGAL_PARAMETERS;
        }

        // load feature map
        FeatureMap<> map;
        FeatureXMLFile().load(exclude, map);

        // convert to targeted experiment if traML output is selected
				//            list.loadTargets(map,excl_targets,exp);
				// else write tab-delimited file directly
				try
				{
					list.writeTargets(map, out);
				}
				catch(Exception::UnableToCreateFile)
        {
          writeLog_("Error: Unable to create output file.");
          return CANNOT_WRITE_OUTPUT_FILE;
        }
      }
      else if(ex_type == FileTypes::IDXML)
      {
        std::vector<PeptideIdentification> pep_ids;
        std::vector<ProteinIdentification> prot_ids;
        IdXMLFile().load(exclude,prot_ids,pep_ids);
				try
				{
					list.writeTargets(pep_ids,out,excl_charges);
				}
				catch(Exception::UnableToCreateFile)
        {
          writeLog_("Error: Unable to create output file.");
          return CANNOT_WRITE_OUTPUT_FILE;
        }
				catch(Exception::InvalidSize)
				{
					writeLog_("Error: Peptide identification contains several hits. Use IDFilter to filter for significant peptide hits.");
					return ILLEGAL_PARAMETERS;
				}
				catch(Exception::MissingInformation)
				{
					writeLog_("Error: Peptide identification contains no RT information.");
					return ILLEGAL_PARAMETERS;
				}
			}
      else   // FASTA format ...
      {
        if (!File::exists(rt_model_file))
        {
          writeLog_("Error: RT model file required for FASTA input to predict RT elution time.");
          return MISSING_PARAMETERS;
        }
        if(excl_charges.empty())
        {
          writeLog_("Error: Protein sequences for exclusion given, but no charge states specified.");
          return MISSING_PARAMETERS;
        }
        std::vector<FASTAFile::FASTAEntry> entries;
        // load fasta-file
        FASTAFile().load(exclude,entries);
        // convert to targeted experiment for traML output
				//            list.loadTargets(entries,excl_targets,exp,missed_cleavages);
				// else for tab-delimited output
				try
				{
					list.writeTargets(entries,out,excl_charges,rt_model_file);
				}
				catch(Exception::UnableToCreateFile)
        {
          writeLog_("Error: Unable to create output file.");
          return CANNOT_WRITE_OUTPUT_FILE;
        }
      }
//      exp.setExcludeTargets(excl_targets);
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    //TraMLFile().store(out, exp);

    return EXECUTION_OK;
  }

};


int main( int argc, const char** argv )
{
  TOPPInclusionExclusionListCreator tool;
  return tool.main(argc,argv);
}

/// @endcond
