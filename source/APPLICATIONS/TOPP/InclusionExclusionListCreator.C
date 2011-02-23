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
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinder </td>
		</tr>
	</table>
	</CENTER>

   Currently this tool can create tab-delimited inclusion or exclusion lists (m/z, RT start, RT stop).
	 The input can either be peptide identifications from previous runs, a feature map or a FASTA-file with proteins.
	 Inclusion and exclusion charges can be specified for FASTA and IdXML input. If no charges are specified in the case of peptide id input, only
   the charge state of the peptide id is in/excluded, otherwise all given charge states are entered to the list.

   The rt window size can be specified via the rel_rt_window_size parameter,
	 then the window is [rt-rel_rt_window_size*rt,rt+rel_rt_window_size*rt]. The default is rt in minutes, set the rt_in_seconds flag to use seconds.

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
    registerInputFile_("include", "<file>", "", "inclusion list input file in fasta or featureXML format.",false);
    setValidFormats_("include",StringList::create("featureXML,fasta"));
    registerInputFile_("exclude","<file>", "", "exclusion list input file in featureXML, IdXML or fasta format.",false);
    setValidFormats_("exclude",StringList::create("featureXML,IdXML,fasta"));
    registerOutputFile_("out", "<file>", "", "output file (tab delimited).");
    //in fasta or featureXML
    registerIntList_("inclusion_charges","<charge>",IntList(),"List containing the charge states to be considered for the inclusion list compounds, space separated.",false);
    setMinInt_("inclusion_charges", 1);
    registerIntList_("exclusion_charges","<charge>",IntList(),"List containing the charge states to be considered for the exclusion list compounds (for idXML and FASTA input), space separated.",false);
    setMinInt_("exclusion_charges", 1);
    registerIntOption_("missed_cleavages","<int>",0,"Number of missed cleavages used for protein digestion.\n",false);
    registerDoubleOption_("rel_rt_window_size","<double>",.05,"The relative factor for the rt_window_size, e.g. the window is calculated as [rt-rt*rel_rt_window_size,rt+rt*rel_rt_window_size].",false);
		setMinFloat_("rel_rt_window_size", 0.0);
    setMaxFloat_("rel_rt_window_size", 10.0);
    registerInputFile_("rt_model","<file>","","RTModel file used for the rt prediction of peptides in fasta files.",false);
		registerFlag_("rt_in_seconds","Create lists with units as seconds instead of minutes (default is 'minutes')");

    
    registerDoubleOption_("merge:mz_tol","<delta m/z>", 10.0, "Two inclusion/exclusion windows are merged when they (almost) overlap in RT (see 'rt_tol') and are close in m/z by this tolerance. Unit of this is defined in 'mz_tol_unit'.",false);
		setMinFloat_("merge:mz_tol", 0.0);
		registerStringOption_("merge:mz_tol_unit", "<unit>", "ppm", "Unit of 'mz_tol'", false);
    setValidStrings_("merge:mz_tol_unit", StringList::create("ppm,Da"));
    registerDoubleOption_("merge:rt_tol","<RT[s]>", 1.1, "Maximal RT delta (in seconds) which would allow two windows in RT to overlap (which causes merging the windows). Two inclusion/exclusion windows are merged when they (almost) overlap in RT and are close in m/z by this tolerance (see 'mz_tol'). Unit of this param is [seconds].",false);
		setMinFloat_("merge:rt_tol", 0.0);
    registerTOPPSubsection_("merge","Options for merging two or more windows into a single window (some vendor instruments do not allow overlap)");

    //    setValidFormats_("out", StringList::create("TraML"));
  }



  ExitCodes main_(int , const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input/output files
    String include(getStringOption_("include"));
    String exclude(getStringOption_("exclude")), out(getStringOption_("out"));

    if(include == "" && exclude == "")
    {
      writeLog_("Error: No input file given.");
      return MISSING_PARAMETERS;
    }
    // currently we can handle only inclusion OR exclusion, will be possible with the traML output
    if(include != "" && exclude != "")
    {
      writeLog_("Error: Currently only inclusion OR exclusion, both will be possible with the traML output coming soon");
      return ILLEGAL_PARAMETERS;
    }

		IntList incl_charges(getIntList_("inclusion_charges"));
    IntList excl_charges(getIntList_("exclusion_charges"));
    Int missed_cleavages(getIntOption_("missed_cleavages"));
    DoubleReal rel_rt_window_size(getDoubleOption_("rel_rt_window_size"));
    String rt_model_file(getStringOption_("rt_model"));
    bool rt_in_seconds(getFlag_("rt_in_seconds"));

    bool mz_tol_as_ppm (getStringOption_("merge:mz_tol_unit") == "ppm");
    DoubleReal mz_tol (getDoubleOption_("merge:mz_tol"));
    DoubleReal rt_tol (getDoubleOption_("merge:rt_tol"));

    //-------------------------------------------------------------
    // loading input: inclusion list part
    //-------------------------------------------------------------

		FileHandler fh;
    TargetedExperiment exp;
    InclusionExclusionList list(rt_tol, mz_tol, mz_tol_as_ppm);
    if(include != "")
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
					list.writeTargets(map,out,rel_rt_window_size,rt_in_seconds);
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
					list.writeTargets(entries,out,incl_charges,rt_model_file,rel_rt_window_size,rt_in_seconds,missed_cleavages);
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
        FeatureXMLFile().load(exclude,map);

        // convert to targeted experiment if traML output is selected
				//            list.loadTargets(map,excl_targets,exp);
				// else write tab-delimited file directly
				try
				{
					list.writeTargets(map,out,rel_rt_window_size,rt_in_seconds);
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
					list.writeTargets(pep_ids,out,rel_rt_window_size,excl_charges,rt_in_seconds);
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
					list.writeTargets(entries,out,excl_charges,rt_model_file,rel_rt_window_size,rt_in_seconds,missed_cleavages);
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
