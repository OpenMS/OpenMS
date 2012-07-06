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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_ProteinResolver ProteinResolver


  @brief A peptide-centric algorithm for Protein Inference.
<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$  \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinResolver</td>
    </tr>
  </table>
</CENTER>

  @experimental This tool has not been tested thoroughly and might behave not as expected!

  This application is used to... For further information see
  Meyer-Arendt et al. IsoformResolver: A Peptide-Centric Algorithm for Protein Inference (2011)


  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_ProteinResolver.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

/*
The algorithm tries to assign to each Protein its experimantally validated peptide.
Proteins are grouped into ISD groups(in silco derived) and MSD groups(MS/MS derived)
if they have in siilco derived or MS/MS derived peptides in common. Proteins and peptides span a bipartite graph.
There is an edge between a protein node and a peptide node iff the protein contains the peptide.
ISD groups are connected graphs in the forementionend bipartite graph. MSD groups are subgraphs of ISD groups.
At the moment pointers are used a lot which might be wise to change for reasons of safety.
*/

class TOPPProteinResolver
    : public TOPPBase
{


  public:
    TOPPProteinResolver()
      : TOPPBase("ProteinResolver","protein inference",false)
    {

    }

  protected:
    void registerOptionsAndFlags_()
    {
      registerInputFile_("fasta","<file>","","input database file");
      setValidFormats_("fasta",StringList::create("FASTA"));
      registerInputFile_("in","<file>","","input file - holds experimental data");
      setValidFormats_("in",StringList::create("idXML,consensusXML"));
      //registerInputFileList_("in","<files>", StringList(),"input files holding the experimental data",true,false);
      //setValidFormats_("in",StringList::create("idXML,consensusXML"));
      //registerInputFile_("design","<file>","","Text file containing the experimental design. See documentation for specific format requirements",false,false,StringList());

      registerOutputFile_("protein_groups","<file>","","output file. Contains all protein groups");
      registerOutputFile_("peptide_table","<file>","","output file. Contains one peptide per line and all proteins which contain that peptide");
      registerOutputFile_("protein_table","<file>","","output file. Contains one protein per line");

      registerIntOption_("missed_cleavages","<number>",2,"the number of allowed missed cleavages", false);
      setMinInt_("missed_cleavages", 0);
      registerIntOption_("min_length","<number>",6,"minimum length of peptide", false);
      registerStringOption_("enzyme","<string>","Trypsin","the digestion enzyme", false);
      setValidStrings_("enzyme",StringList::create("Trypsin"));
    }


    ExitCodes main_(int , const char**)
    {
      //-------------------------------------------------------------
      // parsing parameters
      //-------------------------------------------------------------
      String fastafile_name = getStringOption_("fasta");
      String input = getStringOption_("in");
      //StringList input = getStringList_("in");
      FileTypes::Type in_type = FileHandler::getType(input);
      //String design = getStringOption_("design");
      String peptide_table_outfile = getStringOption_("peptide_table");
      String protein_groups_outfile = getStringOption_("protein_groups");
      String protein_table_outfile = getStringOption_("protein_table");

      //-------------------------------------------------------------
      // set up enzymatic digestor
      //-------------------------------------------------------------
      EnzymaticDigestion digestor;
      String enzyme_name = getStringOption_("enzyme");
      digestor.setEnzyme(digestor.getEnzymeByName(enzyme_name));
      UInt min_size = getIntOption_("min_length");
      UInt missed_cleavages = getIntOption_("missed_cleavages");
      digestor.setMissedCleavages(missed_cleavages);

      //-------------------------------------------------------------
      // initialize rest
      //-------------------------------------------------------------
      IdXMLFile IdXML_file;
      ConsensusXMLFile consensusXML_file;
      ConsensusMap consensus;
      vector<ProteinIdentification> protein_identifications;
      vector<PeptideIdentification> peptide_identifications;
      bool id;

      //-------------------------------------------------------------
      // reading data input
      //-------------------------------------------------------------
      if(in_type == FileTypes::IDXML)
      {
        IdXML_file.load(input,protein_identifications, peptide_identifications);
        id = true;
      }
      else //consensusXML
      {
        consensusXML_file.load(input,consensus);
        id = false;
      }

      //-------------------------------------------------------------
      // fasta file
      //-------------------------------------------------------------
      FASTAFile file;
      std::vector<FASTAFile::FASTAEntry> protein_data;
      file.load(fastafile_name, protein_data);

      //-------------------------------------------------------------
      // calculation
      //-------------------------------------------------------------
      vector<ProteinResolver::ProteinEntry> protein_nodes;
      protein_nodes.resize(protein_data.size());
      vector<ProteinResolver::PeptideEntry> peptide_nodes;
      vector<ProteinResolver::ISDGroup> isd_groups;
      vector<ProteinResolver::MSDGroup> msd_groups;
      vector<Size> reindexed_proteins, reindexed_peptides;

      ProteinResolver resolver;
      resolver.resolve(isd_groups, msd_groups, protein_data, reindexed_proteins, reindexed_peptides,
                       protein_nodes, peptide_nodes, peptide_identifications, consensus, digestor, id, min_size);


      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
      //TODO make writeProteinGroups optional -- not in mzTAB included
      resolver.writeProteinGroups(isd_groups, msd_groups, protein_groups_outfile);
      resolver.writeProteinTable(protein_nodes, reindexed_proteins, protein_table_outfile);
      //writemzTab(protein_nodes, peptide_nodes, reindexed_protein,reindexed_peptide,protein_table_outfile, peptide_identifications );

      cout<<"Statistics:"<<endl;
      if(in_type == FileTypes::IDXML)
      {
        resolver.writePeptideTable(peptide_nodes, reindexed_peptides, peptide_identifications, peptide_table_outfile);
        resolver.countTargetDecoy(msd_groups, peptide_identifications);
      }
      else //consensusXML
      {
        resolver.writePeptideTable(peptide_nodes, reindexed_peptides, consensus, peptide_table_outfile);
        resolver.countTargetDecoy(msd_groups, consensus);
      }

      cout<<"number of ISD groups: "<<isd_groups.size()<<endl;
      cout<<"number of MSD groups: "<<msd_groups.size()<<endl;
      Size target_peptides = 0;
      Size decoy_peptides = 0;
      Size target_plus_decoy_peptides = 0;
      Size exp_peps = 0;
      for(vector<ProteinResolver::MSDGroup>::iterator msd = msd_groups.begin(); msd < msd_groups.end(); ++ msd)
      {
        target_peptides += msd->number_of_target;
        decoy_peptides += msd->number_of_decoy;
        target_plus_decoy_peptides += msd->number_of_target_plus_decoy;
        exp_peps += msd->peptides.size();
      }
      Real fdr1 = (float)decoy_peptides/(float)(target_peptides+target_plus_decoy_peptides);
      Real fdr2 = (float)(decoy_peptides+target_plus_decoy_peptides)/(float)target_peptides;
      cout<<"number of target peptides = " << target_peptides<<endl;
      cout<<"number of decoy peptides = " << decoy_peptides<<endl;
      cout<<"number of target+decoy peptides = " << target_plus_decoy_peptides<<endl;
      cout<<"number of peptides in MSD groups = " << exp_peps <<endl;
      cout<<"The estimated FDR for protein list is between "<< fdr1 <<" and " << fdr2 <<endl;

      return EXECUTION_OK;
    }
};

int main( int argc, const char** argv )
{
  TOPPProteinResolver tool;
  return tool.main(argc,argv);
}

/// @endcond
