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

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>
#include <OpenMS/ANALYSIS/QUANTITATION/QuantitativeExperimentalDesign.h>

#include <QDir>

using std::list;
using std::vector;
using std::ofstream;
using std::endl;

using namespace OpenMS;


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
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ ProteinResolver \f$ \longrightarrow \f$</td>
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

  <p><b>Remark:</b>
  If one single input file (through @p in) is provided, the parameters @p in_list and @p in_path are ignored.
  If parameter @p in_list is given, @p in_path is ignored. Parameter @p in_path is only considered, if @p in_list and @p in are empty.
  </p>

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_ProteinResolver.cli

  <B>Input</B>

Since the ProteinResolver offers three different input parameters, there are some possibilites how to use this TOPP tool.
  <dl>
      <dt>One single input file (@p in)</dt>
      <dd>The ProteinResolver simply performs the protein inference based on the above mentioned algortihm of Meyer-Arendt et al. (2011) for that specific file.</dd>

      <dt>Multiple files (@p in_list or @p in_path)</dt>
      <dd>
        <ol>
          <li>If no experimental design file is given, all files are treated as in batch processing.</li>
          <li>If an experimental design file is provided, all files that can be mapped to the same experimental design are treated as
              one single input file (simply by merging them before the computation).</li>
         </ol>
      </dd>
  </dl>


  <B>Output</B>

  <p>Four possible outputs are available:

    <dl>
         <dt>Protein groups</dt>
         <dd>For each MSD group, the ISD group, the protein indices, the peptide indeces, the number of peptides in MSD group, the number
              of proteins in ISD and the number of proteins in ISD are written to the output file</dd>
        <dt>Protein table</dt>
        <dd>The resulting text file contains one protein per line<dd>
        <dt>Peptide table</dt>
        <dd>The output file will contain one peptide per line and all proteins which contain that specific peptide</dd>
        <dt>Statistics:</dt>
        <dd>Number of ISD groups, number of MSD groups, number of target peptides, number of decoy peptides,
            number of target and decoy peptides, number of peptides in MSD groups and estimated FDR for protein list.</dd>
    </dl>

    The results for different input files are appended and written into the same output file. In other words, no matter how many input files you have,
    you will end up with one single output file.
  </p>


  <B>Text file format of the quantitative experimental design:</B>

  <p>
    The text file has to be column-based and must contain only one additional line as header.
    The header must specify two specific columns that represents the file name and an identifier for the experimental setup.
    These two header identifier can be defined as parameter and must be unique (default: "File" and "ExperimentalSetting").
    There are four options how the columns can be separated: tabulator, comma, semi-colon and whitespace.

    <i>Example for text file format:</i>

      <CENTER>
        <table>
          <tr>
            <td ALIGN="center" BGCOLOR="#EBEBEB">Slice</td>
            <td ALIGN="center" BGCOLOR="#EBEBEB">FileName</td>
            <td ALIGN="center" BGCOLOR="#EBEBEB">Experiment</td>
          </tr>
          <tr>
            <td ALIGN="center">1</td>
            <td ALIGN="center">SILAC_2_1</td>
            <td ALIGN="center">S1224</td>
          </tr>
          <tr>
            <td ALIGN="center">4</td>
            <td ALIGN="center">SILAC_3_4</td>
            <td ALIGN="center">D1224</td>
          </tr>
          <tr>
            <td ALIGN="center">2</td>
            <td ALIGN="center">SILAC_10_2</td>
            <td ALIGN="center">S1224</td>
          </tr>
          <tr>
            <td ALIGN="center">7</td>
            <td ALIGN="center">SILAC_8_7</td>
            <td ALIGN="center">S1224</td>
          </tr>
        </table>
      </CENTER>

      In this case the default values of the parameters "experiment" and "file" have to be changed to "Experiment", respectively "FileName".
      The separator should be changed if the file is not tab separated.
      Every other column (here: first column) is just ignored. Not every file mentioned in the design file has to be given as input file;
      and every input file that has no match in the design file is ignored for the computation.<br />
      <br />
      <i>Consider the following scenario:</i><br />
      <br />
      <b>Input files:</b> SILAC_2_1.consensusXML, SILAC_3_4.consensusXML, SILAC_10_2.consensusXML and SILAC_8_7_.consensusXML<br />
      <br />
      <b>First step:</b>  Data from SILAC_2_1.consensusXML and SILAC_10_2.consensusXML is merged, because both files can be mapped to the same setting S1224.
      SILAC_8_7_.consensusXML is ignored, since SILAC_8_7_ is no match to SILAC_8_7.<br />
      <br />
      <b>Second step:</b> ProteinResolver computes results for the merged data, and the data from the file SILAC_3_4.<br />
      <br />
      <b>Third step:</b> ProteinResolver writes the results for experimental setting S1224 and D1224 to the same output file.<br />
      <br />

  </p>

	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_ProteinResolver.html
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
      : TOPPBase("ProteinResolver","protein inference",false),
        resolver_params_(), design_params_()
    {
    }

  protected:

    Param resolver_params_; // parameters for ProteinResolver
    Param design_params_; // parameters for QuantitativeExperimentalDesign

    void registerOptionsAndFlags_()
    {
      registerInputFile_("fasta","<file>","","Input database file",true,false);
      setValidFormats_("fasta",StringList::create("FASTA"));
      registerInputFile_("in","<file>","","Input file holding experimental data",false,false);
      setValidFormats_("in",StringList::create("idXML,consensusXML"));
      registerInputFileList_("in_list","<files>", StringList(),"Input files holding the experimental data. Ignored if parameter 'in' is given",false,false);
      setValidFormats_("in_list",StringList::create("idXML,consensusXML"));
      registerStringOption_("in_path","<file>","","Path to idXMLs or consensusXMLs files. Ignored if either 'in' or 'in_list' is given.",false,false);
      registerInputFile_("design","<file>","","Text file containing the experimental design. See documentation for specific format requirements",false,false);

      registerOutputFile_("protein_groups","<file>","","output file. Contains all protein groups");
      registerOutputFile_("peptide_table","<file>","","output file. Contains one peptide per line and all proteins which contain that peptide");
      registerOutputFile_("protein_table","<file>","","output file. Contains one protein per line");
      registerOutputFile_("additional_info","<file>","","output file for additional info",false,true);

      Param temp = ProteinResolver().getParameters();
      registerFullParam_(temp);

      Param temp_ = QuantitativeExperimentalDesign().getParameters();
      registerFullParam_(temp_);
    }


    void writeProteinGroups_(SVOutStream& out, const vector<ProteinResolver::ResolverResult>& result)
    {
      //ISD group descriptor";
      out << "MSD_group" << "ISD_group" << "Protein_indices" << "Peptide_indices" << "#Peptides_MSD" << "#Proteins_ISD" << "ProteinIDs_ISD" << endl;

      for ( vector<ProteinResolver::ResolverResult>::const_iterator iter = result.begin(); iter != result.end(); ++iter)
      {
        const ProteinResolver::ResolverResult& res = *iter;
        const vector<ProteinResolver::ISDGroup> *isd_groups = res.isds;
        const vector<ProteinResolver::MSDGroup> *msd_groups = res.msds;

        for ( vector<ProteinResolver::ISDGroup>::const_iterator isd = isd_groups->begin(); isd != isd_groups->end(); ++isd )
        {
          for ( list<Size>::const_iterator msd_group = isd->msd_groups.begin(); msd_group != isd->msd_groups.end(); ++msd_group)
          {
            const ProteinResolver::MSDGroup *msd = &msd_groups->at(*msd_group);//[*msd_group];
            //Protein group
            out << msd->index;
            out << isd->index;
            //Protein index
            String protein_indices = "";
            for ( list<ProteinResolver::ProteinEntry*>::const_iterator prot = msd->proteins.begin();prot != msd->proteins.end(); ++prot)
            {
              protein_indices += (*prot)->index;//fasta_entry->identifier;
              if (prot != --(msd->proteins.end() ) ) protein_indices += ";";
            }
            out << protein_indices;
            //pep index
            String peptide_indices = "";
            for ( list<ProteinResolver::PeptideEntry*>::const_iterator peps = msd->peptides.begin();peps != msd->peptides.end(); ++peps)
            {
              if (!(*peps)->experimental) continue;
              peptide_indices +=  (*peps)->index;//identifications[(*peps)->peptide_identification].getHits()[(*pep)->peptide_hit].getSequence().toString();
              if ( peps != --(msd->peptides.end() ) ) peptide_indices += ";";
            }
            out << peptide_indices;
            //Peptides in MSD
            out << msd->peptides.size();
            //#prots in ISD
            out << isd->proteins.size();
            //prots in ISD;
            String prots_ISD = "";
            for ( list<ProteinResolver::ProteinEntry*>::const_iterator prot = isd->proteins.begin();prot != isd->proteins.end(); ++prot)
            {
              prots_ISD += (*prot)->fasta_entry->identifier;
              if(prot != --(isd->proteins.end() ) ) prots_ISD += ";";
            }
            out << prots_ISD;
            out << endl;
          }
        }
      }
    }

    void writePeptideTable_(SVOutStream& out, const vector<ProteinResolver::ResolverResult>& result)
    {
      out << "MSD_group" << "ISD_group" << "Protein_indices" << "Protein_ID" << "Peptide_sequence" << "Var_mods" << "Peptide_MW" << "Score" << "Charge" << "RT" << "MZ" << endl;

      UInt counter = result.size();
      for ( vector<ProteinResolver::ResolverResult>::const_iterator iter = result.begin(); iter != result.end(); ++iter,--counter)
      {
        const ProteinResolver::ResolverResult& res = *iter;
        const vector<Size> *reindexed_peptides= res.reindexed_peptides;
        const vector<ProteinResolver::PeptideEntry> *peptides = res.peptide_entries;

        for(vector<Size>::const_iterator pep = reindexed_peptides->begin(); pep != reindexed_peptides->end(); ++pep)
        {
          //MSD and ISD group
          const ProteinResolver::PeptideEntry* peptide_entry = &peptides->at(*pep);
          out << peptide_entry->msd_group;
          out << peptide_entry->isd_group;
          //Protein index
          String protein_indices = "";
          for(list<ProteinResolver::ProteinEntry*>::const_iterator prot = peptide_entry->proteins.begin(); prot != peptide_entry->proteins.end(); ++prot)
          {
            protein_indices += (*prot)->index;
            if( prot != --(peptide_entry->proteins.end())) protein_indices += ";";
          }
          out << protein_indices;
          //Protein ID
          String protein_ID = "";
          for(list<ProteinResolver::ProteinEntry*>::const_iterator prot = peptide_entry->proteins.begin(); prot != peptide_entry->proteins.end(); ++prot)
          {
            protein_ID += (*prot)->fasta_entry->identifier;
            if( prot != --(peptide_entry->proteins.end())) protein_ID += ";";
          }
          out << protein_ID;
          //peptide sequence
          if (res.input_type == ProteinResolver::ResolverResult::PeptideIdent )
          {
            const vector<PeptideIdentification>& identifications =  *res.peptide_identification;
            const PeptideIdentification& pi = ProteinResolver().getPeptideIdentification(identifications, peptide_entry);
            const PeptideHit& ph = ProteinResolver().getPeptideHit(identifications, peptide_entry);
            const AASequence& seq = ph.getSequence();
            out << seq.toUnmodifiedString();
            //var mods TODO
            out << seq.toString();
            //Pep MW
            out << seq.getMonoWeight();
            //score
            out << ph.getScore();
            //charge
            out << ph.getCharge();
            //RT
            out << String(pi.getMetaValue("RT"));
            //MZ
            out << String(pi.getMetaValue("MZ"));
            out << endl;
          }
          else
          {
            const ConsensusMap& consensus = *res.consensus_map;
            const PeptideIdentification& pi = ProteinResolver().getPeptideIdentification(consensus, peptide_entry);
            const PeptideHit& ph = ProteinResolver().getPeptideHit(consensus, peptide_entry);
            const AASequence& seq = ph.getSequence();
            out << seq.toUnmodifiedString();
            //var mods TODO
            out << seq.toString();
            //Pep MW
            out << seq.getMonoWeight();
            //score
            out << ph.getScore();
            //charge
            out << ph.getCharge();
            //RT
            out << String(pi.getMetaValue("RT"));
            //MZ
            out << String(pi.getMetaValue("MZ"));
            out << endl;
          }
        }
      }
    }

    void writeProteinTable_(SVOutStream& out, const vector<ProteinResolver::ResolverResult>& result)
    {
      out << "MSD_group" << "ISD_group" << "Peptide_indices" << "Protein_index" << "Protein_ID" << "#Peptides_per_Protein" << "Prot_MW" << "Coverage" << endl;

      UInt counter = 0;
      for ( vector<ProteinResolver::ResolverResult>::const_iterator iter = result.begin(); iter != result.end(); ++iter,++counter)
      {
        const ProteinResolver::ResolverResult& res = *iter;
        const vector<Size> *reindexed_proteins = res.reindexed_proteins;
        const vector<ProteinResolver::ProteinEntry> *proteins = res.protein_entries;

        for(vector<Size>::const_iterator prot = reindexed_proteins->begin(); prot != reindexed_proteins->end(); ++prot)
        {
          //MSD and ISD group
          const ProteinResolver::ProteinEntry& protein_entry = proteins->at(*prot);
          out << protein_entry.msd_group;
          out << protein_entry.isd_group;
          //peptide indices
          Size pep_counter = 0;
          String peptide_indices = "";
          for(list<ProteinResolver::PeptideEntry*>::const_iterator pep = protein_entry.peptides.begin(); pep != protein_entry.peptides.end(); ++pep)
          {
            if((*pep)->experimental)
            {
              peptide_indices += (*pep)->index;
              ++pep_counter;
              if(pep_counter < protein_entry.number_of_experimental_peptides)
              {
                peptide_indices += ";";
              }
              else
              {
                break;
              }
            }
          }
          out << peptide_indices;
          //Protein identifier
          out <<  protein_entry.index;
          //TODO 1 a or 2* you know what I mean? tmp += prot->typeToString;
          //Protein ID
          out << protein_entry.fasta_entry->identifier;
          //#Peps in prot
          out << protein_entry.number_of_experimental_peptides;
          //Prot MW
          out << protein_entry.weight;
          //coverage
          out << protein_entry.coverage;
          out << endl;
        }
      }
    }

    void writeStatistics_(SVOutStream& out, const vector<ProteinResolver::ResolverResult>& result)
    {
      for ( vector<ProteinResolver::ResolverResult>::const_iterator iter = result.begin(); iter != result.end(); ++iter)
      {
        const ProteinResolver::ResolverResult& res = *iter;
        const vector<ProteinResolver::ISDGroup> *isd_groups = res.isds;
        const vector<ProteinResolver::MSDGroup> *msd_groups = res.msds;

        out << "Number of ISD groups:" << isd_groups->size() << endl;
        out << "Number of MSD groups:" << msd_groups->size() << endl;
        Size target_peptides = 0;
        Size decoy_peptides = 0;
        Size target_plus_decoy_peptides = 0;
        Size exp_peps = 0;

        for(vector<ProteinResolver::MSDGroup>::const_iterator msd = msd_groups->begin(); msd < msd_groups->end(); ++ msd)
        {
          target_peptides += msd->number_of_target;
          decoy_peptides += msd->number_of_decoy;
          target_plus_decoy_peptides += msd->number_of_target_plus_decoy;
          exp_peps += msd->peptides.size();
        }
        Real fdr1 = (float)decoy_peptides/(float)(target_peptides+target_plus_decoy_peptides);
        Real fdr2 = (float)(decoy_peptides+target_plus_decoy_peptides)/(float)target_peptides;
        out << "Number of target peptides:" << target_peptides << endl;
        out << "Number of decoy peptides:" << decoy_peptides << endl;
        out << "Number of target+decoy peptides:" << target_plus_decoy_peptides << endl;
        out << "Number of peptides in MSD groups:" << exp_peps << endl;
        out << "The estimated FDR for protein list is between"<< fdr1 <<"and" << fdr2 << endl;
      }
    }

    const String getBaseName_(String& file_path)
    {
      const String basename = QFileInfo(file_path.toQString()).baseName().toStdString();
      return basename;
    }



    ExitCodes main_(int , const char**)
    {
      //-------------------------------------------------------------
      // parsing parameters
      //-------------------------------------------------------------
      String fastafile_name = getStringOption_("fasta");
      String input = getStringOption_("in");
      StringList input_list = getStringList_("in_list");
      String input_path = getStringOption_("in_path");

      String design = getStringOption_("design");
      String output_pep_table = getStringOption_("peptide_table");
      String ouput_prot_groups = getStringOption_("protein_groups");
      String output_prot_table = getStringOption_("protein_table");
      String output_stats = getStringOption_("additional_info");

      //-------------------------------------------------------------
      // check input parameters
      //-------------------------------------------------------------
      if ( input.empty() && input_list.empty() && input_path.empty() )
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                          "All input options are empty.");
      }

      //-------------------------------------------------------------
      // check output parameters
      //-------------------------------------------------------------
      if ( output_pep_table.empty() && ouput_prot_groups.empty() && output_prot_table .empty() )
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                          "All output options are empty.");
      }

      //-------------------------------------------------------------
      // read fasta file
      //-------------------------------------------------------------
      FASTAFile file;
      std::vector<FASTAFile::FASTAEntry> protein_data;
      file.load(fastafile_name, protein_data);

      //-------------------------------------------------------------
      // set up protein resolver ( parameters )
      //-------------------------------------------------------------
      ProteinResolver resolver;
      resolver_params_ = resolver.getParameters();
      Logger::LogStream nirvana; // avoid parameter update messages
      resolver_params_.update(getParam_(),false,false,nirvana);
      resolver.setParameters(resolver_params_);
      resolver.setProteinData(protein_data);

      //-------------------------------------------------------------
      // initialize rest
      //-------------------------------------------------------------
      IdXMLFile IdXML_file;
      ConsensusXMLFile consensusXML_file;
      ConsensusMap consensus;
      vector<ProteinIdentification> protein_identifications;
      vector<PeptideIdentification> peptide_identifications;


      //-------------------------------------------------------------
      //-------------------------------------------------------------
      // SINGLE INPUT FILE
      //-------------------------------------------------------------
      //-------------------------------------------------------------

      if ( !input.empty() )
      {
        FileTypes::Type in_type = FileHandler::getType(input);

        if(in_type == FileTypes::IDXML)
        {
          //read data
          IdXML_file.load(input,protein_identifications, peptide_identifications);
          //calculation
          resolver.resolveID(peptide_identifications );
        }
        else //consensusXML
        {
          //read data
          consensusXML_file.load(input,consensus);
          //calculation
          resolver.resolveConsensus(consensus);
        }
      }



      //-------------------------------------------------------------
      //-------------------------------------------------------------
      // MULTIPLE INPUT FILES
      //-------------------------------------------------------------
      //-------------------------------------------------------------

      //-------------------------------------------------------------
      // set up quantitative experimental design
      //-------------------------------------------------------------
      bool experimental_design = !design.empty();
      QuantitativeExperimentalDesign designer;
      TextFile design_file;

      // ensure that no single file is provided
      if ( experimental_design && input.empty() )
      {
        // false -> do not trim lines; -1 -> read all lines
        design_file.load(design, false,-1);

        design_params_ = designer.getParameters();
        design_params_.update(getParam_(),false,false,nirvana);
        designer.setParameters(design_params_);
      }


      //-------------------------------------------------------------
      // multiple files given in a list format
      //-------------------------------------------------------------
      //   - without design: batch processing
      //   - with design: files from same experimental setting are
      //                  considered merged before quantitiation
      //-------------------------------------------------------------

      if ( !input_list.empty()  && input.empty() )
      {
        // designer merges files that belong to the same experimental setting.
        // quantitation is performed
        if( experimental_design)
        {
          designer.applyDesign2Resolver(resolver, design_file, input_list);
        }
        // otherwise batch processing
        else
        {
          for ( StringList::Iterator iter = input_list.begin(); iter != input_list.end(); ++iter)
          {
            FileTypes::Type in_type = FileHandler::getType(input);
            if(in_type == FileTypes::IDXML)
            {
              // load ensures that vector are cleared upon loading
              IdXML_file.load(*iter,protein_identifications, peptide_identifications);
              resolver.resolveID(peptide_identifications );
            }
            else
            {
              // load ensures that consensus map is cleared upon loading
              consensusXML_file.load(*iter, consensus);
              resolver.resolveConsensus(consensus);
            }
          }
        }
      }

      //-------------------------------------------------------------
      // multiple files from given path
      //-------------------------------------------------------------
      //   - without design: batch processing
      //   - with design: files from same experimental setting are
      //                  considered merged before quantitiation
      //-------------------------------------------------------------

      if ( !input_path.empty() && input.empty() && input_list.empty() )
      {
        // read file names in 'in_path'
        QDir dir(input_path.toQString());
        QStringList filters;
        filters << "*.idXML" << "*.IdXML" << "*.consensusXML" << ".ConsensusXML";
        dir.setNameFilters(filters);
        dir.setFilter(QDir::Files | QDir::Readable );
        dir.setSorting(QDir::Name | QDir::IgnoreCase);
        QFileInfoList list = dir.entryInfoList();

        StringList input_list;
        for (int i = 0; i < list.size(); ++i)
        {
          QFileInfo fileInfo = list.at(i);
          input_list.push_back( fileInfo.absoluteFilePath() );
        }

        if( experimental_design)
        {
          designer.applyDesign2Resolver(resolver, design_file, input_list);
        }
        // otherwise batch processing
        else
        {
          for ( StringList::Iterator iter = input_list.begin(); iter != input_list.end(); ++iter)
          {
            FileTypes::Type in_type = FileHandler::getType(input);
            if(in_type == FileTypes::IDXML)
            {
              // load ensures that vector are cleared upon loading
              IdXML_file.load(*iter,protein_identifications, peptide_identifications);
              // result is stored in struct and appended to result-vector of the resolver
              resolver.resolveID(peptide_identifications);
            }
            else
            {
              //load ensures that consensus map is cleared upon loading
              consensusXML_file.load(*iter, consensus);
              // result is stored in struct and appended to result-vector of the resolver
              resolver.resolveConsensus(consensus);
            }
          }
        }
      }


      //-------------------------------------------------------------
      // write output files
      //-------------------------------------------------------------
      if ( !ouput_prot_groups.empty() )
      {
        ofstream outstr(ouput_prot_groups.c_str());
        SVOutStream output(outstr, "\t", "_", String::NONE);
        writeProteinGroups_(output, resolver.getResults());
        outstr.close();
      }
      if ( !output_pep_table.empty() )
      {
        ofstream outstr(output_pep_table.c_str());
        SVOutStream output(outstr, "\t", "_", String::NONE);
        writePeptideTable_(output, resolver.getResults());
        outstr.close();
      }
      if ( !output_prot_table .empty() )
      {
        ofstream outstr(output_prot_table .c_str());
        SVOutStream output(outstr, "\t", "_", String::NONE);
        writeProteinTable_(output, resolver.getResults() );
        outstr.close();
      }
      if ( !output_stats.empty())
      {
        ofstream outstr(output_stats.c_str());
        SVOutStream output(outstr, "\t", "_", String::NONE);
        writeStatistics_(output, resolver.getResults());
        outstr.close();
      }

      return EXECUTION_OK;
    }
};

int main( int argc, const char** argv )
{
  TOPPProteinResolver tool;
  return tool.main(argc,argv);
}

/// @endcond
