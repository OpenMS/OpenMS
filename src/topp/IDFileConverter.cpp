// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers, Clemens Groepl, Chris Bielow, Mathias Walzer,
// Hendrik Weisser
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
#include <OpenMS/FORMAT/TextFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>

#include <boost/math/special_functions/fpclassify.hpp> // for "isnan"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_IDFileConverter IDFileConverter

    @brief Converts peptide/protein identification engine file formats.

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

@p mz_file:@n
Some search engine output files (like pepXML, mascotXML, Sequest .out files) may not contain retention times, only scan numbers or spectrum IDs. To be able to look up the actual RT values, the raw file has to be provided using the parameter @p mz_file. (If the identification results should be used later to annotate feature maps or consensus maps, it is critical that they contain RT values. See also @ref TOPP_IDMapper.)

@p mz_name:@n
pepXML files can contain results from multiple experiments. However, the idXML format does not support this. The @p mz_name parameter (or @p mz_file, if given) thus serves to define what parts to extract from the pepXML.

@p scan_regex:@n
This advanced parameter defines a spectrum reference format via a Perl-style regular expression. The reference format connects search hits to the MS2 spectra that were searched, and may be needed to look up e.g. retention times in the raw data (@p mz_file). See the documentation of class @ref OpenMS::SpectrumLookup "SpectrumLookup" for details on how to specify spectrum reference formats. Note that it is not necessary to look up any information in the raw data if that information can be extracted directly from the spectrum reference, in which case @p mz_file is not needed.@n
For Mascot results exported to (Mascot) XML, scan numbers that can be used to look up retention times (via @p mz_file) should be given in the "pep_scan_title" XML elements, but the format can vary. Some default formats are defined in the Mascot XML reader, but if those fail to extract the scan numbers, @p scan_regex can be used to overwrite the defaults.@n
For pepXML, supplying @p scan_regex may be necessary for files exported from Mascot, but only if the default reference formats (same as for Mascot XML) do not match. The spectrum references to which @p scan_regex is applied are read from the "spectrum" attribute of the "spectrum_query" elements.@n
For Percolator tab-delimited output, information is extracted from the "PSMId" column. By default, extraction of scan numbers and charge states is supported for MS-GF+ Percolator results (retention times and precursor m/z values can then be looked up in the raw data via @p mz_file).@n

Some information about the supported input types:
@li @ref OpenMS::MzIdentMLFile "mzIdentML"
@li @ref OpenMS::IdXMLFile "idXML"
@li @ref OpenMS::PepXMLFile "pepXML"
@li @ref OpenMS::ProtXMLFile "protXML"
@li @ref OpenMS::MascotXMLFile "Mascot XML"
@li @ref OpenMS::OMSSAXMLFile "OMSSA XML"
@li @ref OpenMS::XTandemXMLFile "X! Tandem XML"
@li @ref OpenMS::SequestOutfile "Sequest .out directory"
@li @ref OpenMS::PercolatorOutfile "Percolator tab-delimited output"


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

private:
  bool add_ionmatches_(vector<PeptideIdentification>& peptide_identifications, String filename, double tolerance)
  {
      TheoreticalSpectrumGenerator tg;
      Param tgp(tg.getDefaults());
      tgp.setValue("add_metainfo", "true");
      tgp.setValue("add_losses", "true");
      tgp.setValue("add_precursor_peaks", "true");
      tgp.setValue("add_abundant_immonium_ions", "true");
      tgp.setValue("add_first_prefix_ion", "true");
      tgp.setValue("add_y_ions", "true");
      tgp.setValue("add_b_ions", "true");
      tgp.setValue("add_a_ions", "true");
      tgp.setValue("add_x_ions", "true");
      tg.setParameters(tgp);    

      SpectrumAlignment sa;
      Param sap = sa.getDefaults();
      sap.setValue("tolerance", tolerance, "...");
      sa.setParameters(sap);
      SpectrumAnnotator annot;
      bool ret = true;
      PeakMap expmap;
      SpectrumLookup lookup;
      FileHandler().loadExperiment(filename, expmap);
      lookup.readSpectra(expmap.getSpectra());

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < (SignedSize)peptide_identifications.size(); ++i)
      {
        try
        {
          String ref = peptide_identifications[i].getMetaValue("spectrum_reference");
          Size index = lookup.findByNativeID(ref);
          annot.addIonMatchStatistics(peptide_identifications[i], expmap[index], tg, sa);
        }
        catch (Exception::ElementNotFound&)
        {
#ifdef _OPENMP
#pragma omp critical (IDFileConverter_ERROR)
#endif
          {
            LOG_ERROR << "Error: Failed to look up spectrum - none with corresponding native ID found." << endl;
            ret = false;
          }
        }
      }
    return ret;
  }

protected:
  void registerOptionsAndFlags_() override
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
    registerInputFile_("mz_file", "<file>", "", "[pepXML, Sequest, Mascot, X! Tandem, mzid, Percolator only] Retention times and native spectrum ids (spectrum_references) will be looked up in this file", false);
    setValidFormats_("mz_file", ListUtils::create<String>("mzML,mzXML,mzData"));
    addEmptyLine_();
    registerStringOption_("mz_name", "<file>", "", "[pepXML only] Experiment filename/path (extension will be removed) to match in the pepXML file ('base_name' attribute). Only necessary if different from 'mz_file'.", false);
    registerFlag_("peptideprophet_analyzed", "[pepXML output only] Write output in the format of a PeptideProphet analysis result. By default a 'raw' pepXML is produced that contains only search engine results.", false);
    registerStringOption_("score_type", "<choice>", PercolatorOutfile::score_type_names[0], "[Percolator only] Which of the Percolator scores to report as 'the' score for a peptide hit", false);
    setValidStrings_("score_type", vector<String>(PercolatorOutfile::score_type_names, PercolatorOutfile::score_type_names + int(PercolatorOutfile::SIZE_OF_SCORETYPE)));

    registerFlag_("ignore_proteins_per_peptide", "[Sequest only] Workaround to deal with .out files that contain e.g. \"+1\" in references column,\n"
                                                 "but do not list extra references in subsequent lines (try -debug 3 or 4)", true);
    registerStringOption_("scan_regex", "<expression>", "", "[Mascot, pepXML, Percolator only] Regular expression used to extract the scan number or retention time. See documentation for details.", false, true);
    registerFlag_("no_spectra_data_override", "[+mz_file only] Setting this flag will avoid overriding 'spectra_data' in ProteinIdentifications if mz_file is given and 'spectrum_reference's are added/updated. Use only if you are sure it is absolutely the same mz_file as used for identification.", true);
    registerDoubleOption_("add_ionmatch_annotation", "<tolerance>", 0,"[+mz_file only] Will annotate the contained identifications with their matches in the given mz_file. Will take quite some while. Match tolerance is .4", false, true);
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // general variables and data
    //-------------------------------------------------------------
    FileHandler fh;
    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;
    SpectrumMetaDataLookup lookup;

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    const String in = getStringOption_("in");
    const String mz_file = getStringOption_("mz_file");

    ProgressLogger logger;
    logger.setLogType(ProgressLogger::CMD);
    logger.startProgress(0, 1, "Loading...");

    if (File::isDirectory(in))
    {
      const String in_directory = File::absolutePath(in).ensureLastChar('/');
      const bool ignore_proteins_per_peptide = getFlag_("ignore_proteins_per_peptide");

      UInt i = 0;
      FileTypes::Type type;
      PeakMap msexperiment;
      // Note: we had issues with leading zeroes, so let us represent scan numbers as Int (next line used to be map<String, float> num_and_rt;)  However, now String::toInt() might throw.
      map<Int, float> num_and_rt;
      vector<String> NativeID;

      // The mz-File (if given)
      if (!mz_file.empty())
      {
        type = fh.getTypeByFileName(mz_file);
        fh.loadExperiment(mz_file, msexperiment, type, log_type_, false, false);

        for (PeakMap::Iterator spectra_it = msexperiment.begin(); spectra_it != msexperiment.end(); ++spectra_it)
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
        String mz_name =  getStringOption_("mz_name");
        if (mz_file.empty())
        {
          PepXMLFile().load(in, protein_identifications,
                            peptide_identifications, mz_name);
        }
        else
        {
          PeakMap exp;
          fh.loadExperiment(mz_file, exp, FileTypes::UNKNOWN, log_type_, false,
                            false);
          if (mz_name.empty()) mz_name = mz_file;
          String scan_regex = getStringOption_("scan_regex");
          // we may have to parse Mascot spectrum references in pepXML, too:
          MascotXMLFile::initializeLookup(lookup, exp, scan_regex);
          PepXMLFile().load(in, protein_identifications,
                            peptide_identifications, mz_name, lookup);
        }
      }

      else if (in_type == FileTypes::IDXML)
      {
        IdXMLFile().load(in, protein_identifications, peptide_identifications);
        // get spectrum_references from the mz data, if necessary:
        if (!mz_file.empty())
        {
          SpectrumMetaDataLookup::addMissingSpectrumReferences(
            peptide_identifications, mz_file, false, !getFlag_("no_spectra_data_override"), protein_identifications);

          double add_ions = getDoubleOption_("add_ionmatch_annotation");
          if (add_ions > 0)
          {
            add_ionmatches_(peptide_identifications, mz_file, add_ions);
          }
        }
      }

      else if (in_type == FileTypes::MZIDENTML)
      {
        LOG_WARN << "Converting from mzid: you might experience loss of information depending on the capabilities of the target format." << endl;
        MzIdentMLFile().load(in, protein_identifications,
                             peptide_identifications);

        // get retention times from the mz data, if necessary:
        if (!mz_file.empty())
        {
          SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(
            peptide_identifications, mz_file, false);

          double add_ions = getDoubleOption_("add_ionmatch_annotation");
          if (add_ions > 0)
          {
            add_ionmatches_(peptide_identifications, mz_file, add_ions);
          }
        }
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
        if (!mz_file.empty())
        {
          String scan_regex = getStringOption_("scan_regex");
          PeakMap exp;
          // load only MS2 spectra:
          fh.getOptions().addMSLevel(2);
          fh.loadExperiment(mz_file, exp, FileTypes::MZML, log_type_, false,
                            false);
          MascotXMLFile::initializeLookup(lookup, exp, scan_regex);
        }
        protein_identifications.resize(1);
        MascotXMLFile().load(in, protein_identifications[0], 
                             peptide_identifications, lookup);
      }

      else if (in_type == FileTypes::XML) // X! Tandem
      {
        ProteinIdentification protein_id;
        ModificationDefinitionsSet mod_defs;
        XTandemXMLFile().load(in, protein_id, peptide_identifications,
                              mod_defs);
        protein_id.setSearchEngineVersion("");
        protein_id.setSearchEngine("XTandem");
        protein_identifications.push_back(protein_id);
        if (!mz_file.empty())
        {
          PeakMap exp;
          fh.getOptions().addMSLevel(2);
          fh.loadExperiment(mz_file, exp, FileTypes::MZML, log_type_, false,
                            false);
          for (vector<PeptideIdentification>::iterator it =
                 peptide_identifications.begin(); it !=
                 peptide_identifications.end(); ++it)
          {
            UInt id = (Int)it->getMetaValue("spectrum_id");
            --id; // native IDs were written 1-based
            if (id < exp.size())
            {
              it->setRT(exp[id].getRT());
              double pre_mz(0.0);
              if (!exp[id].getPrecursors().empty())
              {
                pre_mz = exp[id].getPrecursors()[0].getMZ();
              }
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
        String score_type = getStringOption_("score_type");
        enum PercolatorOutfile::ScoreType perc_score =
          PercolatorOutfile::getScoreType(score_type);
        if (!mz_file.empty())
        {
          PeakMap experiment;
          fh.loadExperiment(mz_file, experiment, FileTypes::UNKNOWN, log_type_, false, false);
          lookup.readSpectra(experiment.getSpectra());
        }
        String scan_regex = getStringOption_("scan_regex");
        if (!scan_regex.empty()) lookup.addReferenceFormat(scan_regex);
        protein_identifications.resize(1);
        PercolatorOutfile().load(in, protein_identifications[0],
                                 peptide_identifications, lookup, perc_score);
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
      String mz_name = getStringOption_("mz_name");
      PepXMLFile().store(out, protein_identifications, peptide_identifications,
                         mz_file, mz_name, peptideprophet_analyzed);
    }

    else if (out_type == FileTypes::IDXML)
    {
      IdXMLFile().store(out, protein_identifications, peptide_identifications);
    }

    else if (out_type == FileTypes::MZIDENTML)
    {
      MzIdentMLFile().store(out, protein_identifications,
                            peptide_identifications);
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
          std::set<String> prot = hit.extractProteinAccessionsSet();
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
