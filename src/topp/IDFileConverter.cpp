// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers, Clemens Groepl, Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/PercolatorOutfile.h>
#include <OpenMS/FORMAT/ProtXMLFile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_IDFileConverter IDFileConverter

    @brief Converts identification engine file formats.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDFileConverter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> TPP tools: PeptideProphet, ProteinProphet </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> TPP tools: ProteinProphet\n(for conversion from idXML to pepXML) </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> Sequest protein identification engine </td>
        </tr>
    </table>
</CENTER>

IDFileConverter can be used to convert identification results from external tools/pipelines (like TPP, Sequest, Mascot, OMSSA, X! Tandem)
into other (OpenMS-specific) formats.
For search engine results, it might be advisable to use the respective TOPP Adapters (e.g. OMSSAAdapter) to avoid the extra conversion step.

The most simple format accepted is '.tsv': A tab separated text file, which contains one or more peptide sequences per line.
Each line represents one spectrum, i.e. is stored as a PeptideIdentification with one or more PeptideHits.
Lines starting with "#" are ignored by the parser.

Conversion from the TPP file formats pepXML and protXML to OpenMS' idXML is quite comprehensive, to the extent that the original data can be
represented in the simpler idXML format.

In contrast, support for converting from idXML to pepXML is limited. The purpose here is simply to create pepXML files containing the relevant
information for the use of ProteinProphet.

Support for conversion to/from mzIdentML (.mzid) is still experimental and may lose information.

<B>Details on additional parameters:</B>

@p mz_file: \n
Some search engine output files (like pepXML, mascotXML, Sequest .out files) may not contain retention times, only scan numbers. To be able to look up the actual RT values, the raw file has to be provided using the parameter @p mz_file. (If the identification results should be used later to annotate feature maps or consensus maps, it is critical that they contain RT values. See also @ref TOPP_IDMapper.)

@p mz_name: \n
pepXML files can contain results from multiple experiments. However, the idXML format does not support this. The @p mz_name parameter (or @p mz_file, if given) thus serves to define what parts to extract from the pepXML.

@p scan_regex: \n
For Mascot results exported to XML, the scan numbers (used to look up retention times using @p mz_file) should be given in the "pep_scan_title" XML elements, but the format can vary. If the defaults fail to extract the scan numbers, a Perl-style regular expression can be given through the advanced parameter @p scan_regex, and will be used instead. The regular expression should contain a named group "SCAN" matching the scan number or "RT" matching the actual retention time. For example, if the format of the "pep_scan_title" elements is "scan=123", where 123 is the scan number, the expression "scan=(?<SCAN>\\d+)" can be used to extract the number. (However, the format in this example is actually covered by the defaults.)\n
For Percolator tab-delimited output, information is extracted from the "PSMId" column. By default, extraction of scan numbers and charge states is supported for MS-GF+ Percolator results (retention times and precursor m/z values can then be looked up in the raw data via @p mz_file).
In a user-defined regular expression, the named groups "SCAN" (scan number), "CHARGE" (charge state), "RT" (retention time) and "MZ" (precursor m/z) are supported. The parameter @p count_from_zero defines whether scans are counted from zero or from one (default) in the number extracted via "SCAN". If "CHARGE", "RT" and "MZ" are present, it is not necessary to look up any information in the raw data, so @p mz_file is not needed.

Some information about the supported input types:
  @ref OpenMS::MzIdentMLFile "mzIdentML"
  @ref OpenMS::PepXMLFile "pepXML"
  @ref OpenMS::ProtXMLFile "protXML"
  @ref OpenMS::IdXMLFile "idXML"
  @ref OpenMS::MascotXMLFile "mascotXML"
  @ref OpenMS::OMSSAXMLFile "omssaXML"
  @ref OpenMS::XTandemXMLFile "XTandem.xml"
  @ref OpenMS::SequestOutfile "Sequest .out directory"
  @ref OpenMS::PercolatorOutfile "Percolator tab-delimited output"

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_IDFileConverter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_IDFileConverter.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDFileConverter :
  public TOPPBase
{
public:
  TOPPIDFileConverter() :
    TOPPBase("IDFileConverter", "Converts identification engine file formats.", true)
  {
  }

protected:
  void
  registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<path/file>", "",
                       "Input file or directory containing the data to convert. This may be:\n"
                       "- a single file in a multi-purpose XML format (pepXML, protXML, idXML, mzid),\n"
                       "- a single file in a search engine-specific format (Mascot: mascotXML, OMSSA: omssaXML, X! Tandem: xml, Percolator: psms),\n"
                       "- a single text file (tab separated) with one line for all peptide sequences matching a spectrum (top N hits),\n"
                       "- for Sequest results, a directory containing .out files.\n");
    setValidFormats_("in", ListUtils::create<String>("pepXML,protXML,mascotXML,omssaXML,xml,psms,tsv,idXML,mzid"));

    registerOutputFile_("out", "<file>", "", "Output file", true);
    String formats("idXML,mzid,pepXML,FASTA");
    setValidFormats_("out", ListUtils::create<String>(formats));
    registerStringOption_("out_type", "<type>", "", "Output file type (default: determined from file extension)", false);
    setValidStrings_("out_type", ListUtils::create<String>(formats));

    addEmptyLine_();
    registerInputFile_("mz_file", "<file>", "", "[pepXML, Sequest, Mascot, X! Tandem, Percolator only] Retention times will be looked up in this file", false);
    setValidFormats_("mz_file", ListUtils::create<String>("mzML,mzXML,mzData"));
    addEmptyLine_();
    registerStringOption_("mz_name", "<file>", "", "[pepXML only] Experiment filename/path (extension will be removed) to match in the pepXML file ('base_name' attribute). Only necessary if different from 'mz_file'.", false);
    registerFlag_("use_precursor_data", "[pepXML only] Use precursor RTs (and m/z values) from 'mz_file' for the generated peptide identifications, instead of the RTs of MS2 spectra.", false);
    registerFlag_("peptideprophet_analyzed", "[pepXML output only] Write output in the format of a PeptideProphet analysis result. By default a 'raw' pepXML is produced that contains only search engine results.", false);
    registerStringOption_("score_type", "<choice>", PercolatorOutfile::score_type_names[0], "[Percolator only] Which of the Percolator scores to report as 'the' score for a peptide hit", false);
    setValidStrings_("score_type", vector<String>(PercolatorOutfile::score_type_names, PercolatorOutfile::score_type_names + int(PercolatorOutfile::SIZE_OF_SCORETYPE)));

    registerFlag_("ignore_proteins_per_peptide", "[Sequest only] Workaround to deal with .out files that contain e.g. \"+1\" in references column,\n"
                                                 "but do not list extra references in subsequent lines (try -debug 3 or 4)", true);
    registerStringOption_("scan_regex", "<expression>", "", "[Mascot, Percolator only] Regular expression used to extract the scan number or retention time. See documentation for details.", false, true);
    registerFlag_("count_from_zero", "[Percolator only] Scan numbers extracted by 'scan_regex' start counting at zero (default: start at one).", true);
  }

  ExitCodes
  main_(int, const char**)
  {
    //-------------------------------------------------------------
    // general variables and data
    //-------------------------------------------------------------
    FileHandler fh;
    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    const String in = getStringOption_("in");

    ProgressLogger logger;
    logger.setLogType(ProgressLogger::CMD);
    logger.startProgress(0, 1, "Loading...");

    if (File::isDirectory(in))
    {
      const String in_directory = File::absolutePath(in).ensureLastChar('/');
      const String mz_file = getStringOption_("mz_file");
      const bool ignore_proteins_per_peptide = getFlag_("ignore_proteins_per_peptide");

      UInt i = 0;
      FileTypes::Type type;
      MSExperiment<Peak1D> msexperiment;
      // Note: we had issues with leading zeroes, so let us represent scan numbers as Int (next line used to be map<String, float> num_and_rt;)  However, now String::toInt() might throw.
      map<Int, float> num_and_rt;
      vector<String> NativeID;

      // The mz-File (if given)
      if (!mz_file.empty())
      {
        type = fh.getTypeByFileName(mz_file);
        fh.loadExperiment(mz_file, msexperiment, type);

        for (MSExperiment<Peak1D>::Iterator spectra_it = msexperiment.begin(); spectra_it != msexperiment.end(); ++spectra_it)
        {
          String(spectra_it->getNativeID()).split('=', NativeID);
          try
          {
            num_and_rt[NativeID[1].toInt()] = spectra_it->getRT();
            // cout << "num_and_rt: " << NativeID[1] << " = " << NativeID[1].toInt() << " : " << num_and_rt[NativeID[1].toInt()] << endl; // CG debuggging 2009-07-01
          }
          catch (Exception::ConversionError& e)
          {
            writeLog_(String("Error: Cannot read scan number as integer. '") + e.getMessage());
          }
        }
      }

      // Get list of the actual Sequest .out-Files
      StringList in_files;
      if (!File::fileList(in_directory, String("*.out"), in_files))
      {
        writeLog_(String("Error: No .out files found in '") + in_directory + "'. Aborting!");
      }

      // Now get to work ...
      for (vector<String>::const_iterator in_files_it = in_files.begin(); in_files_it != in_files.end(); ++in_files_it)
      {
        vector<PeptideIdentification> peptide_ids_seq;
        ProteinIdentification protein_id_seq;
        vector<double> pvalues_seq;
        vector<String> in_file_vec;

        SequestOutfile sequest_outfile;

        writeDebug_(String("Reading file ") + *in_files_it, 3);

        try
        {
          sequest_outfile.load((String) (in_directory + *in_files_it), peptide_ids_seq, protein_id_seq, 1.0, pvalues_seq, "Sequest", ignore_proteins_per_peptide);

          in_files_it->split('.', in_file_vec);

          for (Size j = 0; j < peptide_ids_seq.size(); ++j)
          {

            // We have to explicitly set the identifiers, because the normal set ones are composed of search engine name and date, which is the same for a bunch of sequest out-files.
            peptide_ids_seq[j].setIdentifier(*in_files_it + "_" + i);

            Int scan_number = 0;
            if (!mz_file.empty())
            {
              try
              {
                scan_number = in_file_vec[2].toInt();
                peptide_ids_seq[j].setRT(num_and_rt[scan_number]);
              }
              catch (Exception::ConversionError& e)
              {
                writeLog_(String("Error: Cannot read scan number as integer. '") + e.getMessage());
              }
              catch (exception& e)
              {
                writeLog_(String("Error: Cannot read scan number as integer. '") + e.what());
              }
              //double real_mz = ( peptide_ids_seq[j].getMZ() - hydrogen_mass )/ (double)peptide_ids_seq[j].getHits()[0].getCharge(); // ???? semantics of mz
              const double real_mz = peptide_ids_seq[j].getMZ() / (double) peptide_ids_seq[j].getHits()[0].getCharge();
              peptide_ids_seq[j].setMZ(real_mz);
            }

            writeDebug_(String("scan: ") + String(scan_number) + String("  RT: ") + String(peptide_ids_seq[j].getRT()) + "  MZ: " + String(peptide_ids_seq[j].getMZ()) + "  Ident: " + peptide_ids_seq[j].getIdentifier(), 4);

            peptide_identifications.push_back(peptide_ids_seq[j]);
          }

          protein_id_seq.setIdentifier(*in_files_it + "_" + i);
          protein_identifications.push_back(protein_id_seq);
          ++i;
        }
        catch (Exception::ParseError& pe)
        {
          writeLog_(pe.getMessage() + String("(file: ") + *in_files_it + ")");
          throw;
        }
        catch (...)
        {
          writeLog_(String("Error reading file: ") + *in_files_it);
          throw;
        }
      }

      writeDebug_("All files processed.", 3);
    } // ! directory
    else
    {
      FileTypes::Type in_type = fh.getType(in);

      if (in_type == FileTypes::PEPXML)
      {
        String exp_name = getStringOption_("mz_file");
        String orig_name =  getStringOption_("mz_name");
        bool use_precursor_data = getFlag_("use_precursor_data");

        if (exp_name.empty())
        {
          PepXMLFile().load(in, protein_identifications,
                            peptide_identifications, orig_name);
        }
        else
        {
          MSExperiment<> exp;
          fh.loadExperiment(exp_name, exp);
          if (!orig_name.empty())
          {
            exp_name = orig_name;
          }
          PepXMLFile().load(in, protein_identifications,
                            peptide_identifications, exp_name, exp,
                            use_precursor_data);
        }
      }
      else if (in_type == FileTypes::IDXML)
      {
        IdXMLFile().load(in, protein_identifications, peptide_identifications);
      }
      else if (in_type == FileTypes::MZIDENTML)
      {
        LOG_WARN << "Converting from mzid: you might experience loss of information depending on the capabilities of the target format." << endl;
        MzIdentMLFile().load(in, protein_identifications, peptide_identifications);
      }
      else if (in_type == FileTypes::PROTXML)
      {
        protein_identifications.resize(1);
        peptide_identifications.resize(1);
        ProtXMLFile().load(in, protein_identifications[0],
                           peptide_identifications[0]);
      }
      else if (in_type == FileTypes::OMSSAXML)
      {
        protein_identifications.resize(1);
        OMSSAXMLFile().load(in, protein_identifications[0],
                            peptide_identifications, true);
      }
      else if (in_type == FileTypes::MASCOTXML)
      {
        String scan_regex = getStringOption_("scan_regex");
        String exp_name = getStringOption_("mz_file");
        MascotXMLFile::RTMapping rt_mapping;
        if (!exp_name.empty())
        {
          PeakMap exp;
          // load only MS2 spectra:
          fh.getOptions().addMSLevel(2);
          fh.loadExperiment(exp_name, exp, FileTypes::MZML, log_type_);
          MascotXMLFile::generateRTMapping(exp.begin(), exp.end(), rt_mapping);
        }
        protein_identifications.resize(1);
        MascotXMLFile().load(in, protein_identifications[0],
                             peptide_identifications, rt_mapping, scan_regex);
      }
      else if (in_type == FileTypes::XML) // X! Tandem
      {
        ProteinIdentification protein_id;
        XTandemXMLFile().load(in, protein_id, peptide_identifications);
        protein_id.setSearchEngineVersion("");
        protein_id.setSearchEngine("XTandem");
        protein_identifications.push_back(protein_id);
        String exp_name = getStringOption_("mz_file");
        if (!exp_name.empty())
        {
          PeakMap exp;
          fh.getOptions().addMSLevel(2);
          fh.loadExperiment(exp_name, exp, FileTypes::MZML, log_type_);
          for (vector<PeptideIdentification>::iterator it = peptide_identifications.begin(); it != peptide_identifications.end(); ++it)
          {
            UInt id = (Int)it->getMetaValue("spectrum_id");
            --id; // native IDs were written 1-based
            if (id < exp.size())
            {
              it->setRT(exp[id].getRT());
              double pre_mz(0.0);
              if (!exp[id].getPrecursors().empty()) pre_mz = exp[id].getPrecursors()[0].getMZ();
              it->setMZ(pre_mz);
              it->removeMetaValue("spectrum_id");
            }
            else
            {
              LOG_ERROR << "XTandem xml: Error: id '" << id << "' not found in peak map!" << endl;
            }
          }
        }
      }
      else if (in_type == FileTypes::PSMS) // Percolator
      {
        String mz_file = getStringOption_("mz_file");
        MSExperiment<> experiment;
        MSExperiment<>* experiment_p = 0;
        if (!mz_file.empty())
        {
          fh.loadExperiment(mz_file, experiment);
          experiment_p = &experiment;
        }
        String score_type = getStringOption_("score_type");
        enum PercolatorOutfile::ScoreType perc_score =
          PercolatorOutfile::getScoreType(score_type);
        String scan_regex = getStringOption_("scan_regex");
        bool count_from_zero = getFlag_("count_from_zero");
        protein_identifications.resize(1);
        PercolatorOutfile().load(in, protein_identifications[0],
                                 peptide_identifications, perc_score, 
                                 scan_regex, count_from_zero, experiment_p);
      }
      else if (in_type == FileTypes::TSV)
      {
        ProteinIdentification protein_id;
        protein_id.setSearchEngineVersion("");
        protein_id.setSearchEngine("XTandem");
        protein_identifications.push_back(protein_id);

        TextFile tf;
        tf.load(in, true, -1, true);
        for (TextFile::Iterator it = tf.begin(); it != tf.end(); ++it)
        {
          it->trim();
          // skip empty and comment lines
          if (it->empty() || it->hasPrefix("#")) continue;

          PeptideIdentification pepid;
          StringList peps;
          it->split('\t', peps, false);
          std::vector<PeptideHit> hits;
          for (StringList::const_iterator sit=peps.begin(); sit != peps.end(); ++sit)
          {
            PeptideHit hit;
            hit.setSequence(AASequence::fromString(*sit));
            hits.push_back(hit);
          }
          pepid.setHits(hits);
          peptide_identifications.push_back(pepid);
        }
      }
      else
      {
        writeLog_("Error: Unknown input file type given. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
    }
    logger.endProgress();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    const String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));
    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(out);
    }
    if (out_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    logger.startProgress(0, 1, "Storing...");
    if (out_type == FileTypes::PEPXML)
    {
      bool peptideprophet_analyzed = getFlag_("peptideprophet_analyzed");
      String mz_file = getStringOption_("mz_file");
      String mz_name = getStringOption_("mz_name");
      PepXMLFile().store(out, protein_identifications, peptide_identifications, mz_file, mz_name, peptideprophet_analyzed);
    }
    else if (out_type == FileTypes::IDXML)
    {
      IdXMLFile().store(out, protein_identifications, peptide_identifications);
    }
    else if (out_type == FileTypes::MZIDENTML)
    {
      MzIdentMLFile().store(out, protein_identifications, peptide_identifications);
    }
    else if (out_type == FileTypes::FASTA)
    {
      Size count = 0;
      ofstream fasta(out.c_str(), ios::out);
      for (Size i = 0; i < peptide_identifications.size(); ++i)
      {
        for (Size l = 0; l < peptide_identifications[i].getHits().size(); ++l)
        {
          const PeptideHit& hit = peptide_identifications[i].getHits()[l];
          String seq = hit.getSequence().toUnmodifiedString();
          std::set<String> prot = hit.extractProteinAccessions();
          fasta << ">" << seq
                << " " << ++count
                << " " << hit.getSequence().toString() 
                << " " << ListUtils::concatenate(StringList(prot.begin(), prot.end()), ";")
                << "\n";
          // FASTA files should have at most 60 characters of sequence info per line
          for (Size j = 0; j < seq.size(); j += 60)
          {
            Size k = min(j + 60, seq.size());
            fasta << seq.substr(j, k - j) << "\n";
          }
        }
      }
      fasta.close();
    }
    else
    {
      writeLog_("Unsupported output file type given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    logger.endProgress();


    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPIDFileConverter tool;
  return tool.main(argc, argv);
}

///@endcond
