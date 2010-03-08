// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Alexandra Zerck$
// --------------------------------------------------------------------------

//#include <OpenMS/FORMAT/TraMLFile.h>
//#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
//#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
//#include <OpenMS/ANALYSIS/PRECURSORSELECTION/InclusionExclusionList.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_InclusionExclusionlistCreator InclusionExclusionlistCreator

   @brief A tool for creating inclusion and/or exclusion lists for LC-MS/MS.

   Currently this tool can create tab-delimited exclusion lists (m/z, RT start, RT stop) given
   peptide identifications from previous runs. If no exclusion_charges are specified, only
   the charge state of the peptide id is excluded, otherwise all given charge states are entered to the list.

   The rt window size can be specified via the rel_rt_window_size parameter, then the window is [rt-rel_rt_window_size*rt,rt+rel_rt_window_size*rt] (the rt in the output file is in minutes).

   
   @TODO Support traML...


   <B>The command line parameters of this tool are:</B>
   @verbinclude TOPP_InclusionExclusionlistCreator.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPInclusionExclusionlistCreator
  : public TOPPBase
{
public:
  TOPPInclusionExclusionlistCreator()
    : TOPPBase("InclusionExclusionlistCreator","Creates inclusion and/or exclusion lists.")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    //registerInputFile_("include", "<file>", "", "inclusion list input file in fasta or featureXML format.",false);
    registerInputFile_("exclude","<file>", "", "exclusion list input file in IdXML format.",false);
    //in fasta or featureXML  
    //registerIntList_("inclusion_charges","<charge>",IntList(),"List containing the charge states to be considered for the inclusion list compounds.",false);
    registerIntList_("exclusion_charges","<charge>",IntList(),"List containing the charge states to be considered for the exclusion list compounds, space separated",false);
    //registerIntOption_("missed_cleavages","<int>",0,"Number of missed cleavages used for protein digestion.\n"false);
    registerDoubleOption_("rel_rt_window_size","<double>",.05,"The relative factor for the rt_window_size, e.g. the window is calculated as [rt-rt*rel_rt_window_size,rt+rt*rel_rt_window_size].",false);
    //   registerInputFile_("rt_model","<file>","","RTModel file used for the rt prediction of peptides in fasta files.",false);
    registerOutputFile_("out", "<file>", "", "output file (tab delimited).");
    //    setValidFormats_("out", StringList::create("TraML"));
  }



  ExitCodes main_(int argc, const char** argv)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input/output files
    //String include(getStringOption_("include"));
    String exclude(getStringOption_("exclude")), out(getStringOption_("out"));
    //    IntList incl_charges(getIntList_("inclusion_charges"));
    IntList excl_charges(getIntList_("exclusion_charges"));
    //Int missed_cleavages(getIntOption_("missed_cleavages"));
    DoubleReal rel_rt_window_size(getDoubleOption_("rel_rt_window_size"));
    //String rt_model_path(getStringOption_("rt_model"));
    
    if(/*include == "" &&*/ exclude == "")
      {
        writeLog_("Error: No input file given.");
        return MISSING_PARAMETERS;
      }
//     // currently we can handle only inclusion OR exclusion, will be possible with the traML output
//     if(include != "" && exclude != "")
//       {
//         writeLog_("Error: Currently only inclusion OR exclusion, will be possible with the traML output coming soon");
//         return ILLEGAL_PARAMETERS;
//       }
    
//     //-------------------------------------------------------------
//     // loading input: inclusion list part
//     //-------------------------------------------------------------

     FileHandler fh;
//     TargetedExperiment exp;
//     InclusionExclusionList list;
//     if(include != "")
//       {
//         FileTypes::Type in_type = fh.getType(include);
//         std::vector<IncludeExcludeTarget> incl_targets;
//         if(in_type == FileTypes::FEATUREXML)
//           {
//             // load feature map
//             FeatureMap<> map;
//             FeatureXMLFile().load(include,map);
            
//             // convert to targeted experiment
//             list.loadTargets(map,incl_targets,exp);
//           }
//         else
//           {
//             if(incl_charges.empty())
//               {
//                 writeLog_("Error: Protein sequences for inclusion given, but no charge states specified.");
//                 return MISSING_PARAMETERS;
//               }
//             std::vector<FASTAFile::FASTAEntry> entries;
//             // load fasta-file
//             FASTAFile().load(include,entries);
//             // convert to targeted experiment
//             list.loadTargets(entries,incl_targets,exp,missed_cleavages);
            
//           }
        
//         exp.setIncludeTargets(incl_targets);
//       }
    //-------------------------------------------------------------
    // loading input: exclusion list part
    //-------------------------------------------------------------
    if(exclude != "")
      {
        FileTypes::Type ex_type = fh.getType(exclude);
        //        std::vector<IncludeExcludeTarget> excl_targets;    
//         if(ex_type == FileTypes::FEATUREXML)
//           {
//             // load feature map
//             FeatureMap<> map;
//             FeatureXMLFile().load(exclude,map);
            
//             // convert to targeted experiment
//             list.loadTargets(map,excl_targets,exp);
//           }
//         else 
          if(ex_type == FileTypes::IDXML)
          {
            std::vector<PeptideIdentification> pep_ids;
            std::vector<ProteinIdentification> prot_ids;
            IdXMLFile().load(exclude,prot_ids,pep_ids);
            std::ofstream outs(out.c_str());
            outs.precision(8);
            if (!outs)
              {
                writeLog_("Error: Unable to create output file.");
                return CANNOT_WRITE_OUTPUT_FILE;
               
              }
            std::vector<PeptideIdentification>::const_iterator pep_id_iter = pep_ids.begin();
            for(;pep_id_iter != pep_ids.end();++pep_id_iter)
              {
                if(pep_id_iter->getHits().size() > 1)
                  {
                    writeLog_("Error: Peptide identification contains several hits. Use IDFilter to filter for significant peptide hits.");
                    return ILLEGAL_PARAMETERS;
                  }
                if(!pep_id_iter->metaValueExists("RT"))
                  {
                    writeLog_("Error: Peptide identification contains no RT information.");
                    return ILLEGAL_PARAMETERS;
                  }
                DoubleReal rt = pep_id_iter->getMetaValue("RT");
                DoubleReal rt_start = (rt - rel_rt_window_size * rt) / 60.; // RT in minutes
                if(rt_start < 0.) rt_start = 0.;
                DoubleReal rt_stop = (rt + rel_rt_window_size * rt) / 60.; // RT in minutes
                std::vector<PeptideHit>::const_iterator pep_hit_iter = pep_id_iter->getHits().begin();
                for(;pep_hit_iter != pep_id_iter->getHits().end();++pep_hit_iter)
                  {
                    Int charge = pep_hit_iter->getCharge();
                    bool charge_found = false;
                    for(Size c = 0; c < excl_charges.size();++c)
                      {
                        DoubleReal mz = pep_hit_iter->getSequence().getMonoWeight(Residue::Full,excl_charges[c])/(DoubleReal)excl_charges[c];
                        outs << mz <<"\t"<<rt_start<<"\t"<<rt_stop<<"\n";
                        if(excl_charges[c] == charge)
                          {
                            charge_found = true;
                          }
                      }
                    if(!charge_found)
                      {
                        DoubleReal mz = pep_hit_iter->getSequence().getMonoWeight(Residue::Full,charge)/(DoubleReal)charge;
                        outs << mz <<"\t"<<rt_start<<"\t"<<rt_stop<<"\n";
                      }
                  }
               
              }
            outs.close();
           }
//         else
//           {
//             if(excl_charges.empty())
//               {
//                 writeLog_("Error: Protein sequences for exclusion given, but no charge states specified.");
//                 return MISSING_PARAMETERS;
//               }
//             std::vector<FASTAFile::FASTAEntry> entries;
//             // load fasta-file
//             FASTAFile().load(include,entries);
//             // convert to targeted experiment
//             list.loadTargets(entries,excl_targets,exp,missed_cleavages);
//           }
//         exp.setExcludeTargets(excl_targets);    
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
  TOPPInclusionExclusionlistCreator tool;
  return tool.main(argc,argv);
}

/// @endcond
