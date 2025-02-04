// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers, Clemens Groepl, Chris Bielow, Mathias Walzer,
// Hendrik Weisser
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/PercolatorOutfile.h>
#include <OpenMS/FORMAT/ProtXMLFile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/SYSTEM/File.h>

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
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; IDFileConverter &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
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

IDFileConverter can be used to convert identification results from external tools/pipelines (like TPP, Sequest, Mascot, OMSSA, X! Tandem) into other (OpenMS-specific) formats.
For search engine results, it might be advisable to use the respective TOPP Adapters (e.g. CometAdapter) to avoid the extra conversion step.

The most simple format accepted is '.tsv': A tab separated text file, which contains one or more peptide sequences per line.
Each line represents one spectrum, i.e. is stored as a PeptideIdentification with one or more PeptideHits.
Lines starting with "#" are ignored by the parser.

Conversion from the TPP file formats pepXML and protXML to OpenMS' idXML is quite comprehensive, to the extent that the original data can be
represented in the simpler idXML format.

In contrast, support for converting from idXML to pepXML is limited. The purpose here is simply to create pepXML files containing the relevant
information for the use of ProteinProphet.
We use the following heuristic: if peptideprophet_analyzed is set, we take the scores from the idXML as is and assume
the PeptideHits contain all necessary information. If peptideprophet is not set, we only provide ProteinProphet-compatible
results with probability-based scores (i.e. Percolator with PEP score or scores from IDPosteriorErrorProbability). All
secondary or non-probability main scores will be written as "search_scores" only.

Support for conversion to/from mzIdentML (.mzid) is still experimental and may lose information.

The xquest.xml format is very specific to Protein-Protein Cross-Linking MS (XL-MS) applications and is only considered useful for compatibility
of OpenPepXL / OpenPepXLLF with the xQuest / xProphet / xTract pipeline. It will only have useful output when converting from idXML or mzid containg XL-MS data.

Also supports generation of .mzML files with theoretical spectra from a .FASTA input.

<B>Details on additional parameters:</B>

@p mz_file: @n
Some search engine output files (like pepXML, mascotXML, Sequest .out files) may not contain retention times, only scan numbers or spectrum IDs. To be able to look up the actual RT values, the raw file has to be provided using the parameter @p mz_file. (If the identification results should be used later to annotate feature maps or consensus maps, it is critical that they contain RT values. See also @ref TOPP_IDMapper.)

@p mz_name: @n
pepXML files can contain results from multiple experiments. However, the idXML format does not support this. The @p mz_name parameter (or @p mz_file, if given) thus serves to define what parts to extract from the pepXML.

@p scan_regex: @n
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
    TOPPBase("IDFileConverter", "Converts identification engine file formats.")
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

#pragma omp parallel for
      for (SignedSize i = 0; i < (SignedSize)peptide_identifications.size(); ++i)
      {
        try
        {
          String ref = peptide_identifications[i].getSpectrumReference();
          Size index = lookup.findByNativeID(ref);
          annot.addIonMatchStatistics(peptide_identifications[i], expmap[index], tg, sa);
        }
        catch (Exception::ElementNotFound&)
        {
#pragma omp critical (IDFileConverter_ERROR)
          {
            OPENMS_LOG_ERROR << "Error: Failed to look up spectrum - none with corresponding native ID found." << endl;
            ret = false;
          }
        }
      }
    return ret;
  }

protected:
  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param p(TheoreticalSpectrumGenerator().getDefaults());
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    p.setValue("enzyme", "Trypsin", "Enzym used to digest the fasta proteins");
    p.setValidStrings("enzyme", ListUtils::create<std::string>(all_enzymes));
    p.setValue("missed_cleavages", 0, "Number of allowed missed cleavages while digesting the fasta proteins");
    p.setValue("min_charge", 1, "Minimum charge");
    p.setValue("max_charge", 1, "Maximum charge");
    p.setValue("precursor_charge", 0, "Manually set precursor charge. (default: 0, meaning max_charge + 1 will be used as precursor charge)");
    return p;
  }
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<path/file>", "",
                       "Input file or directory containing the data to convert. This may be:\n"
                       "- a single file in OpenMS database format (.oms),\n"
                       "- a single file in a multi-purpose XML format (.idXML, .mzid, .pepXML, .protXML),\n"
                       "- a single file in a search engine-specific format (Mascot: .mascotXML, OMSSA: .omssaXML, X! Tandem: .xml, Percolator: .psms, xQuest: .xquest.xml),\n"
                       "- a single file in fasta format (can only be used to generate a theoretical mzML),\n"
                       "- a single text file (tab separated) with one line for all peptide sequences matching a spectrum (top N hits),\n"
                       "- for Sequest results, a directory containing .out files.\n");
    setValidFormats_("in", ListUtils::create<String>("oms,idXML,mzid,fasta,pepXML,protXML,mascotXML,omssaXML,xml,psms,tsv,xquest.xml"));

    registerOutputFile_("out", "<file>", "", "Output file", true);
    String formats("oms,idXML,mzid,pepXML,fasta,xquest.xml,mzML");
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
    registerFlag_("no_spectra_data_override", "[+mz_file only] Avoid overriding 'spectra_data' in protein identifications if 'mz_file' is given and 'spectrum_reference's are added/updated. Use only if you are sure it is absolutely the same 'mz_file' as used for identification.", true);
    registerFlag_("no_spectra_references_override", "[+mz_file only] Avoid overriding 'spectrum_reference' in peptide identifications if 'mz_file' is given and a 'spectrum_reference' is already present.", true);
    registerDoubleOption_("add_ionmatch_annotation", "<tolerance>", 0, "[+mz_file only] Annotate the identifications with ion matches from spectra in 'mz_file' using the given tolerance (in Da). This will take quite some time.", false, true);

    registerFlag_("concatenate_peptides", "[FASTA output only] Will concatenate the top peptide hits to one peptide sequence, rather than write a new peptide for each hit.", true);
    registerIntOption_("number_of_hits", "<integer>", 1, "[FASTA output only] Controls how many peptide hits will be exported. A value of 0 or less exports all hits.", false, true);

    registerSubsection_("fasta_to_mzml", "[FASTA input + MzML output only] Parameters used to adjust simulation of the theoretical spectra.");
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
    IdentificationData id_data;

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    const String in = getStringOption_("in");
    const String mz_file = getStringOption_("mz_file");
    FileTypes::Type in_type = FileTypes::UNKNOWN; // set below if 'in' isn't a directory

    const String out = getStringOption_("out");
    FileTypes::Type out_type = FileHandler::getConsistentOutputfileType(out, getStringOption_("out_type"));
    if (out_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

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
        fh.loadExperiment(mz_file, msexperiment, {type}, log_type_, false, false);

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
            writeLogWarn_(String("Error: Cannot read scan number as integer. '") + e.what());
          }
        }
      }

      // Get list of the actual Sequest .out-Files
      StringList in_files;
      if (!File::fileList(in_directory, "*.out", in_files))
      {
        writeLogError_(String("Error: No .out files found in '") + in_directory + "'. Aborting!");
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
                writeLogError_(String("Error: Cannot read scan number as integer. '") + e.what());
              }
              catch (exception& e)
              {
                writeLogError_(String("Error: Cannot read scan number as integer. '") + e.what());
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
          writeLogError_(pe.what() + String("(file: ") + *in_files_it + ")");
          throw;
        }
        catch (...)
        {
          writeLogError_(String("Error reading file: ") + *in_files_it);
          throw;
        }
      }

      writeDebug_("All files processed.", 3);
    } // ! directory
    else
    {
      in_type = fh.getType(in);
      switch (in_type)
      {
      case FileTypes::PEPXML:
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
          fh.loadExperiment(mz_file, exp, {}, log_type_, false,
                            false);
          if (mz_name.empty()) mz_name = mz_file;
          String scan_regex = getStringOption_("scan_regex");
          // we may have to parse Mascot spectrum references in pepXML, too:
          MascotXMLFile::initializeLookup(lookup, exp, scan_regex);
          PepXMLFile().load(in, protein_identifications,
                            peptide_identifications, mz_name, lookup);
        }
      }
      break;

      case FileTypes::IDXML:
      {
        FileHandler().loadIdentifications(in, protein_identifications, peptide_identifications, {FileTypes::IDXML});
        // get spectrum_references from the mz data, if necessary:
        if (!mz_file.empty())
        {
          SpectrumMetaDataLookup::addMissingSpectrumReferences(
            peptide_identifications,
            mz_file,
            false,
            !getFlag_("no_spectra_data_override"),
            !getFlag_("no_spectra_references_override"),
            protein_identifications);

          double add_ions = getDoubleOption_("add_ionmatch_annotation");
          if (add_ions > 0)
          {
            add_ionmatches_(peptide_identifications, mz_file, add_ions);
          }
        }
      }
      break;

      case FileTypes::MZIDENTML:
      {
        OPENMS_LOG_WARN << "Converting from mzid: you might experience loss of information depending on the capabilities of the target format." << endl;
		FileHandler().loadIdentifications(in, protein_identifications,
                             peptide_identifications, {FileTypes::MZIDENTML});

        // get retention times from the mz data, if necessary:
        if (!mz_file.empty())
        {
		  // Add RTs if missing
		  MSExperiment exp;    
		  MzMLFile mzml_file{};
          mzml_file.getOptions().setMetadataOnly(true);
		  mzml_file.load(mz_file, exp); 
          SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(peptide_identifications, exp);

          double add_ions = getDoubleOption_("add_ionmatch_annotation");
          if (add_ions > 0)
          {
            add_ionmatches_(peptide_identifications, mz_file, add_ions);
          }
        }
      }
      break;

      case FileTypes::PROTXML:
      {
        FileHandler().loadIdentifications(in, protein_identifications,
                           peptide_identifications, {FileTypes::PROTXML});
      }
      break;

      case FileTypes::OMSSAXML:
      {
        FileHandler().loadIdentifications(in, protein_identifications,
                            peptide_identifications);
      }
      break;

      case FileTypes::MASCOTXML:
      {
        if (!mz_file.empty())
        {
          String scan_regex = getStringOption_("scan_regex");
          PeakMap exp;
          // load only MS2 spectra:
          fh.getOptions().addMSLevel(2);
          fh.loadExperiment(mz_file, exp, {}, log_type_, false,
                            false);
          MascotXMLFile::initializeLookup(lookup, exp, scan_regex);
        }
        protein_identifications.resize(1);
        MascotXMLFile().load(in, protein_identifications[0],
                             peptide_identifications, lookup);
      }
      break;

      case FileTypes::XML: // X! Tandem
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
          fh.loadExperiment(mz_file, exp, {}, log_type_, false,
                            false);
          for (PeptideIdentification& pep : peptide_identifications)
          {
            UInt id = (Int)pep.getMetaValue("spectrum_id");
            --id; // native IDs were written 1-based
            if (id < exp.size())
            {
              pep.setRT(exp[id].getRT());
              double pre_mz(0.0);
              if (!exp[id].getPrecursors().empty())
              {
                pre_mz = exp[id].getPrecursors()[0].getMZ();
              }
              pep.setMZ(pre_mz);
              pep.removeMetaValue("spectrum_id");
            }
            else
            {
              OPENMS_LOG_ERROR << "XTandem xml: Error: id '" << id << "' not found in peak map!" << endl;
            }
          }
        }
      }
      break;

      case FileTypes::PSMS: // Percolator
      {
        String score_type = getStringOption_("score_type");
        enum PercolatorOutfile::ScoreType perc_score =
          PercolatorOutfile::getScoreType(score_type);
        if (!mz_file.empty())
        {
          PeakMap experiment;
          fh.loadExperiment(mz_file, experiment, {}, log_type_, false, false);
          lookup.readSpectra(experiment.getSpectra());
        }
        String scan_regex = getStringOption_("scan_regex");
        if (!scan_regex.empty()) lookup.addReferenceFormat(scan_regex);
        protein_identifications.resize(1);
        PercolatorOutfile().load(in, protein_identifications[0],
                                 peptide_identifications, lookup, perc_score);
      }
      break;

      case FileTypes::TSV:
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
      break;

      case FileTypes::XQUESTXML:
      {
        FileHandler().loadIdentifications(in, protein_identifications, peptide_identifications, {FileTypes::XQUESTXML});
      }
      break;

      case FileTypes::FASTA:
      {
        // handle out type
        if (out_type != FileTypes::MZML)
        {
          writeLogError_("Error: Illegal output file type given. Fasta can only be converted to an MzML. Aborting!");
          printUsage_();
          return ILLEGAL_PARAMETERS;
        }

        MSExperiment exp;
        TheoreticalSpectrumGenerator tsg;

        // extract parameters and remove non tsg params
        Param p = getParam_().copy("fasta_to_mzml:", true);
        String enzyme = p.getValue("enzyme").toString();
        Int mc = p.getValue("missed_cleavages");
        Int min_charge = p.getValue("min_charge");
        Int max_charge = p.getValue("max_charge");
        Int pc_charge = p.getValue("precursor_charge");
        p.remove("enzyme");
        p.remove("missed_cleavages");
        p.remove("min_charge");
        p.remove("max_charge");
        p.remove("precursor_charge");

        if (min_charge > max_charge)
        {
          writeLogError_("Error: 'fasta_to_mzml:min_charge' must be smaller than or equal to 'fasta_to_mzml:max_charge'.");
          printUsage_();
          return ILLEGAL_PARAMETERS;
        }

        OPENMS_PRECONDITION(pc_charge == 0 || pc_charge >= max_charge, "Error: 'fasta_to_mzml:precursor_charge' must be bigger than or equal to 'fasta_to_mzml:max_charge'.\nSet 'precursor_charge' to '0' to automaticly use 'max_charge' + 1.");

        tsg.setParameters(p);
        ProteaseDigestion digestor;
        digestor.setEnzyme(enzyme);
        digestor.setMissedCleavages(mc);

        // loop through fasta input
        FASTAFile::FASTAEntry entry;
        FASTAFile f;
        f.readStart(in);
        UInt count_catches{};
        while (f.readNext(entry))
        {
          // digest sequence of fasta entry
          vector<AASequence> digested_peptides;
          AASequence seq = AASequence::fromString(entry.sequence);
          digestor.digest(seq, digested_peptides);

          // for each peptide calculate the theoretical spectrum
          for (const auto& peptide : digested_peptides)
          {
            PeakSpectrum spec;

            try
            {
              tsg.getSpectrum(spec, peptide, min_charge, max_charge, pc_charge);
            }
            catch (Exception::InvalidSize())
            {
              ++count_catches;
            }

            exp.addSpectrum(move(spec));
          }
        }
        if (count_catches > 0)
        {
          writeLogWarn_("No spectra were calculated for " + String(count_catches) + " peptides because they were to small for generating a C- or X-ion.");
        }
        logger.endProgress();

        logger.startProgress(0, 1, "Storing...");

        FileHandler().storeExperiment(out, exp, {FileTypes::MZML});

        logger.endProgress();

        return EXECUTION_OK;
      }
      break;

      case FileTypes::OMS:
      {
        OMSFile().load(in, id_data);
        if (out_type != FileTypes::OMS)
        {
          IdentificationDataConverter::exportIDs(id_data, protein_identifications, peptide_identifications);
        }
      }
      break;

      default:
        writeLogError_("Error: Unknown input file type given. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
    }
    logger.endProgress();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    logger.startProgress(0, 1, "Storing...");
    switch (out_type)
    {
    case FileTypes::PEPXML:
    {
      bool peptideprophet_analyzed = getFlag_("peptideprophet_analyzed");
      String mz_name = getStringOption_("mz_name");
      PepXMLFile().store(out, protein_identifications, peptide_identifications,
                         mz_file, mz_name, peptideprophet_analyzed);
    }
    break;

    case FileTypes::IDXML:
      FileHandler().storeIdentifications(out, protein_identifications, peptide_identifications, {FileTypes::IDXML});
      break;

    case FileTypes::MZIDENTML:
      FileHandler().storeIdentifications(out, protein_identifications,
                            peptide_identifications, {FileTypes::MZIDENTML});
      break;

    case FileTypes::XQUESTXML:
      FileHandler().storeIdentifications(out, protein_identifications,
                                  peptide_identifications, {FileTypes::XQUESTXML});
      break;

    case FileTypes::FASTA:
    {
      Size count = 0;
      Int max_hits = getIntOption_("number_of_hits");
      if (max_hits < 1)
      {
        max_hits = INT_MAX;
      }

      bool concat = getFlag_("concatenate_peptides");
      //Because by concatenation of peptides [KR]|P sites will probably be created, peptides starting with 'P' are
      //saved separately and later moved to the beginning of the concatenated sequence.
      //This is done to avoid losing information about the preceding peptides if a peptides starts with 'P'.
      String all_p; //peptides beginning with 'P'
      String all_but_p; //all the others

      FASTAFile f;
      f.writeStart(out);
      FASTAFile::FASTAEntry entry;
      for (const PeptideIdentification& pep_id : peptide_identifications)
      {
        Int curr_hit = 1;
        for (const PeptideHit& hit : pep_id.getHits())
        {
          if (curr_hit > max_hits)
          {
            break;
          }
          ++curr_hit;

          String seq = hit.getSequence().toUnmodifiedString();
          if (concat)
          {
            if (seq[0] == 'P')
            {
              all_p += seq;
            }
            else
            {
              all_but_p += seq;
            }
          }
          else
          {
            std::set<String> prot = hit.extractProteinAccessionsSet();
            entry.sequence = seq;
            entry.identifier = seq;
            entry.description = String(count) + " " + hit.getSequence().toString() + " " + ListUtils::concatenate(StringList(prot.begin(), prot.end()), ";");

            f.writeNext(entry);
            ++count;
          }
        }
      }
      if (concat)
      {
        entry.sequence = all_p + all_but_p;
        entry.identifier = protein_identifications[0].getSearchEngine() + "_" + Constants::UserParam::CONCAT_PEPTIDE;
        entry.description = "";

        f.writeNext(entry);
      }
    }
    break;

    case FileTypes::OMS:
    {
      if (in_type != FileTypes::OMS)
      {
        IdentificationDataConverter::importIDs(id_data, protein_identifications, peptide_identifications);
      }
      OMSFile().store(out, id_data);
    }
    break;

    default:
      writeLogError_("Unsupported output file type given. Aborting!");
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
