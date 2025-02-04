// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MascotInfile.h>
#include <OpenMS/FORMAT/PepXMLFileMascot.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <map>
#include <fstream>
#include <string>

#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QProcess>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MascotAdapter MascotAdapter

@brief Identifies peptides in MS/MS spectra via Mascot.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; MascotAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

This wrapper application serves for getting peptide identifications
for MS/MS spectra. It uses a local installation of the Mascot
server to generate the identifications. A second wrapper (MascotAdapterOnline) is
available which is able to perform identifications by
communicating with a Mascot server over the network. So, it is not
necessary to execute MascotAdapterOnline on the same machine
as Mascot.

The minimal version of Mascot supported with this server is 2.1.

This wrapper can be executed in three different
modes:
<ol>
            <li>
            The whole process of ProteinIdentification via Mascot is executed.
            Inputfile is a mzData file containing the MS/MS spectra
            for which the identifications are to be found. The results
            are written as a idXML output file. This mode is selected
            by default.
            </li>

            <li>
            Only the first part of the ProteinIdentification process is performed.
            This means that the MS/MS data is transformed into Mascot
            Generic Format (mgf) which can be used directly with Mascot.
            Being in the cgi directory of the Mascot directory calling a Mascot
            process should look like the following:

            @code ./nph-mascot.exe 1 -commandline -f outputfilename < inputfilename @endcode

            Consult your Mascot reference manual for further details.

            This mode is selected by the <b>-mascot_in</b> option in the command line.
            </li>

            <li>
            Only the second part of the ProteinIdentification process is performed.
            This means that the outputfile of the Mascot server is
            translated into idXML.

            This mode is selected by the <b>-mascot_out</b> option in the command line.
            </li>
</ol>

<br>
If your Mascot server is installed on the same computer as the
TOPP applications the MascotAdapter can be executed in mode 1.
Otherwise the Mascot engine has to be executed manually assisted
by mode 2 and mode 3. The ProteinIdentification steps then look like:

<ul>
    <li>
        execute MascotAdapter in mode 2
        @code ./MascotAdapter -in mzDataFile -out mascotGenericFormatFile -mascot_in @endcode
    </li>
    <li>
        copy mascotGenericFormatFile to your Mascot server
    </li>
    <li>
        call your Mascot server process:
        @code ./nph-mascot.exe 1 -commandline -f mascotOutFile < mascotGenericFormatFile @endcode
    </li>
    <li>
        call the script to export your outfile in mascot xml
        @code ./export_dat.pl do_export=1 export_format=XML file=mascotOutFile _sigthreshold=0
        _showsubset=1 show_same_sets=1 show_unassigned=0 prot_score=0 pep_exp_z=0 pep_score=0
        pep_homol=0 pep_ident=0 pep_seq=1 show_header=1 show_queries=1 pep_rank=0 > mascotXMLFile @endcode
    </li>
    <li>
        copy mascotXMLFile to the server on which the TOPP applications are installed
    </li>
    <li>
        execute MascotAdapter in mode 3
        @code ./MascotAdapter -in mascotXMLFile -out IdXMLFile -mascot_out @endcode
    </li>
</ul>

<p>
For mode 1 you have to specify the directory in which the Mascot
server is installed. This is done by setting the option <b>mascot_dir</b>
in the ini file. Furthermore you have to specify a folder in which
the user has write permissions. This is done by setting the option
<b>temp_data_directory</b> in the ini file.
Two temporary files will be created in this directory during execution
but deleted at the end of execution.
<br>

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MascotAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MascotAdapter.html

You can specify the Mascot parameters <b>precursor_mass_tolerance</b>
(the peptide mass tolerance), <b>peak_mass_tolerance</b> (the MS/MS tolerance),
<b>taxonomy</b> (restriction to a certain subset of the database), <b>modifications</b>,
<b>variable_modifications</b>, <b>charges</b> (the possible charge variants),
<b>db</b> (database where the peptides are searched in), <b>hits</b> (number of hits),
<b>cleavage</b> (the cleavage enzyme), <b>missed_cleavages</b> (number of missed cleavages)
and <b>mass_type</b> (Monoisotopic or Average) via the ini file.

<br>
Known problems with Mascot server execution:
<ul>
    <li>
    getting error message:
    "FATAL_ERROR: M00327
     The ms-monitor daemon/service is not running, please start it."
    </li>

    <li>
    Possible explanations:
    </li>
    <ul>
        <li>
        Your ms-monitor is really not running => consult your Mascot
                                                                                         reference manual for
                                                                                         details about starting
                                                                                         the Mascot server.
        </li>
        <li>
        (Suppose you have Mascot installed in directory mascot.)
        mascot/data/mascot.control is not writable for the current user.
        This has to be changed. Otherwise you will not be able to
        use the Mascot server via the shell and receive the above error
        message.<br>
        => Change write permissions of the file mascot/data/mascot.control
             such that the current user has write permissions to it.
        </li>
    </ul>
</ul>


@todo This adapter is using antiquated internal methods and needs to be updated! E.g. use MascotGenericFile.h instead of MascotInfile.h....
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPMascotAdapter :
  public TOPPBase
{
public:
  TOPPMascotAdapter() :
    TOPPBase("MascotAdapter", "Annotates MS/MS spectra using Mascot.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file in mzData format.\n"
                                           "Note: In mode 'mascot_out' a Mascot results file (.mascotXML) is read");
    setValidFormats_("in", {"mzData", "mascotXML"});
    registerOutputFile_("out", "<file>", "", "output file in idXML format.\n"
                                             "Note: In mode 'mascot_in' Mascot generic format is written.");
    setValidFormats_("out", {"idXML", "mgf"});
    registerStringOption_("out_type", "<type>", "", "output file type (for TOPPAS)", false, false);
    setValidStrings_("out_type", {"idXML", "mgf"});
    registerStringOption_("instrument", "<i>", "Default", "the instrument that was used to measure the spectra", false);
    registerDoubleOption_("precursor_mass_tolerance", "<tol>", 2.0, "the precursor mass tolerance", false);
    registerDoubleOption_("peak_mass_tolerance", "<tol>", 1.0, "the peak mass tolerance", false);
    registerStringOption_("taxonomy", "<tax>", "All entries", "the taxonomy", false);
    setValidStrings_("taxonomy", ListUtils::create<String>("All entries,. . Archaea (Archaeobacteria),. . Eukaryota (eucaryotes),. . . . Alveolata (alveolates),. . . . . . Plasmodium falciparum (malaria parasite),. . . . . . Other Alveolata,. . . . Metazoa (Animals),. . . . . . Caenorhabditis elegans,. . . . . . Drosophila (fruit flies),. . . . . . Chordata (vertebrates and relatives),. . . . . . . . bony vertebrates,. . . . . . . . . . lobe-finned fish and tetrapod clade,. . . . . . . . . . . . Mammalia (mammals),. . . . . . . . . . . . . . Primates,. . . . . . . . . . . . . . . . Homo sapiens (human),. . . . . . . . . . . . . . . . Other primates,. . . . . . . . . . . . . . Rodentia (Rodents),. . . . . . . . . . . . . . . . Mus.,. . . . . . . . . . . . . . . . . . Mus musculus (house mouse),. . . . . . . . . . . . . . . . Rattus,. . . . . . . . . . . . . . . . Other rodentia,. . . . . . . . . . . . . . Other mammalia,. . . . . . . . . . . . Xenopus laevis (African clawed frog),. . . . . . . . . . . . Other lobe-finned fish and tetrapod clade,. . . . . . . . . . Actinopterygii (ray-finned fishes),. . . . . . . . . . . . Takifugu rubripes (Japanese Pufferfish),. . . . . . . . . . . . Danio rerio (zebra fish),. . . . . . . . . . . . Other Actinopterygii,. . . . . . . . Other Chordata,. . . . . . Other Metazoa,. . . . Dictyostelium discoideum,. . . . Fungi,. . . . . . Saccharomyces Cerevisiae (baker's yeast),. . . . . . Schizosaccharomyces pombe (fission yeast),. . . . . . Pneumocystis carinii,. . . . . . Other Fungi,. . . . Viridiplantae (Green Plants),. . . . . . Arabidopsis thaliana (thale cress),. . . . . . Oryza sativa (rice),. . . . . . Other green plants,. . . . Other Eukaryota,. . Bacteria (Eubacteria),. . . . Actinobacteria (class),. . . . . . Mycobacterium tuberculosis complex,. . . . . . Other Actinobacteria (class),. . . . Firmicutes (gram-positive bacteria),. . . . . . Bacillus subtilis,. . . . . . Mycoplasma,. . . . . . Streptococcus Pneumoniae,. . . . . . Streptomyces coelicolor,. . . . . . Other Firmicutes,. . . . Proteobacteria (purple bacteria),. . . . . . Agrobacterium tumefaciens,. . . . . . Campylobacter jejuni,. . . . . . Escherichia coli,. . . . . . Neisseria meningitidis,. . . . . . Salmonella,. . . . . . Other Proteobacteria,. . . . Other Bacteria,. . Viruses,. . . . Hepatitis C virus,. . . . Other viruses,. . Other (includes plasmids and artificial sequences),. . unclassified,. . Species information unavailable"));
    registerStringList_("modifications", "<mods>", StringList(), "the modifications i.e. Carboxymethyl (C)", false);
    registerStringList_("variable_modifications", "<mods>", StringList(), "the variable modifications i.e. Carboxymethyl (C)", false);
    registerStringList_("charges", "[1+ 2+ ...]", ListUtils::create<String>("1+,2+,3+"), "the different charge states", false);
    registerStringOption_("db", "<name>", "MSDB", "the database to search in", false);
    registerStringOption_("hits", "<num>", "AUTO", "the number of hits to report", false);
    registerStringOption_("cleavage", "<enz>", "Trypsin", "The enzyme descriptor to the enzyme used for digestion. (Trypsin is default, None would be best for peptide input or unspecific digestion, for more please refer to your mascot server).", false);
    setValidStrings_("cleavage", ListUtils::create<String>("Trypsin,Arg-C,Asp-N,Asp-N_ambic,Chymotrypsin,CNBr,CNBr+Trypsin,Formic_acid,Lys-C,Lys-C/P,PepsinA,Tryp-CNBr,TrypChymo,Trypsin/P,V8-DE,V8-E,semiTrypsin,LysC+AspN,None"));
    registerIntOption_("missed_cleavages", "<num>", 0, "number of allowed missed cleavages", false);
    setMinInt_("missed_cleavages", 0);
    registerDoubleOption_("sig_threshold", "<num>", 0.05, "significance threshold", false);
    registerDoubleOption_("pep_homol", "<num>", 1, "peptide homology threshold", false);
    registerDoubleOption_("pep_ident", "<num>", 1, "peptide ident threshold", false);
    registerIntOption_("pep_rank", "<num>", 1, "peptide rank", false);
    registerDoubleOption_("prot_score", "<num>", 1, "protein score", false);
    registerDoubleOption_("pep_score", "<num>", 1, "peptide score", false);
    registerIntOption_("pep_exp_z", "<num>", 1, "peptide expected charge", false);
    registerIntOption_("show_unassigned", "<num>", 1, "show_unassigned", false);
    registerDoubleOption_("first_dim_rt", "<num>", 0, "additional information which is added to every peptide identification as metavalue if set > 0", false);
    registerStringOption_("boundary", "<string>", "", "MIME boundary for mascot output format", false);
    registerStringOption_("mass_type", "<type>", "Monoisotopic", "mass type", false);
    setValidStrings_("mass_type", ListUtils::create<String>("Monoisotopic,Average"));
    registerStringOption_("mascot_directory", "<dir>", "", "the directory in which mascot is located", false);
    registerStringOption_("temp_data_directory", "<dir>", "", "a directory in which some temporary files can be stored", false);
  }

  ExitCodes main_(int, const char **) override
  {
    // path to the log file
    String logfile = "mascot.log";
    // log filestream (as long as the real logfile is not determined yet)
    String inputfile_name;
    String outputfile_name;
    String mascot_infile_name = "tmp.mascot_in";
    String mascot_outfile_name = "tmp_mascot_in.out";
    String mascot_output_name = "tmp_mascot.output";
    String mascot_cgi_dir;
    String mascot_data_dir;
    String call;
    String instrument;
    String taxonomy;
    String mascotXML_file_name = "";
    String pepXML_file_name = "";
    MzDataFile mzdata_infile;
    PeakMap experiment;
    MascotXMLFile mascotXML_file;
    PepXMLFileMascot pepXML_file;
    MascotInfile mascot_infile;
    StringList mods;
    StringList variable_mods;
    ProteinIdentification protein_identification;
    vector<PeptideIdentification> identifications;
    double precursor_mass_tolerance(0);
    double peak_mass_tolerance(0);
    double pep_ident(0), sigthreshold(0), pep_homol(0), prot_score(0), pep_score(0);
    int pep_rank(0), pep_exp_z(0), show_unassigned(0);
    String temp_charge;
    string db;
    string hits;
    string cleavage;
    UInt missed_cleavages;
    string mass_type;

    bool mascot_in = false;
    bool mascot_out = false;
    DateTime date_time;
    String date_time_string;
    String boundary = "";
    map<String, vector<AASequence> > modified_peptides;
    double first_dim_rt = 0;

    date_time.now();
    date_time_string = date_time.get();
    date_time_string.substitute(':', '.');        // Windows does not allow ":" in filenames!
    StringList parts;
    date_time_string.split(' ', parts);

    mascot_infile_name = parts[0] + "_" + parts[1] + "_" + mascot_infile_name;
    mascot_outfile_name = parts[0] + "_" + parts[1] + "_" + mascot_outfile_name;
    mascot_output_name = parts[0] + "_" + parts[1] + "_" + mascot_output_name;
    parts.clear();

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    inputfile_name = getStringOption_("in");
    first_dim_rt = getDoubleOption_("first_dim_rt");


    outputfile_name = getStringOption_("out");

    boundary = getStringOption_("boundary");
    if (!boundary.empty())
    {
      writeDebug_(String("Boundary: ") + boundary, 1);
    }

    FileTypes::Type in_type = FileHandler::getType(inputfile_name);
    FileTypes::Type out_type;
    if (!getStringOption_("out_type").empty())
    {
        out_type = FileTypes::nameToType(getStringOption_("out_type"));
    }
    else
    {
        out_type = FileHandler::getType(outputfile_name);
    }
    
    mascot_out = in_type == FileTypes::MASCOTXML;
    mascot_in = out_type == FileTypes::MGF;
    if (mascot_out && mascot_in)
    {
      writeLogError_("When the input file is a mascotXML, only idXML can be written. When the input is mzData, only MGF is written. Please change the output type accordingly.");
      return ILLEGAL_PARAMETERS;
    }
    
    db = getStringOption_("db");
    hits = getStringOption_("hits");
    cleavage = getStringOption_("cleavage");
    missed_cleavages = getIntOption_("missed_cleavages");
    mass_type = getStringOption_("mass_type");

    sigthreshold = getDoubleOption_("sig_threshold");
    pep_homol = getDoubleOption_("pep_homol");
    pep_ident = getDoubleOption_("pep_ident");
    pep_rank = getIntOption_("pep_rank");
    pep_exp_z = getIntOption_("pep_exp_z");
    show_unassigned = getIntOption_("show_unassigned");
    prot_score = getDoubleOption_("prot_score");
    pep_score = getDoubleOption_("pep_score");

    instrument = getStringOption_("instrument");
    precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
    peak_mass_tolerance = getDoubleOption_("peak_mass_tolerance");
    taxonomy = getStringOption_("taxonomy");

    /// fixed modifications
    mods = getStringList_("modifications");

    /// variable modifications
    variable_mods = getStringList_("variable_modifications");

    /// charges
    parts = getStringList_("charges");
    IntList charges;
    for (String& c : parts)
    {
      if (c.hasPrefix("-") || c.hasSuffix("-"))
      {
        charges.push_back(-1 * (c.remove('-').toInt()));
      }
      else
      {
        charges.push_back(c.remove('+').toInt());
      }
    }
    if (charges.empty())
    {
      writeLogError_("No charge states specified for Mascot search. Aborting!");
      return ILLEGAL_PARAMETERS;
    }


    if (mascot_in)
    {
      mascot_infile_name = outputfile_name;
      writeDebug_("Mascot flag: mascot_in (reads in MzData writes Mascot generic format)", 1);
    }
    else if (mascot_out)
    {
      mascotXML_file_name = inputfile_name;

      writeDebug_("Mascot flag: mascot_out (reads in Mascot results file writes idXML file)", 1);
    }
    else
    {
      writeDebug_("No Mascot flag set: reads in MzData writes idXML file", 1);
    }
    if (!mascot_in && !mascot_out)
    {
      // full pipeline:
      mascot_cgi_dir = getStringOption_("mascot_directory");
      if (mascot_cgi_dir.empty())
      {
        writeLogError_("No Mascot directory specified. Aborting!");
        return ILLEGAL_PARAMETERS;
      }
      writeDebug_(String("Mascot directory: ") + mascot_cgi_dir, 1);
      mascot_cgi_dir += "/cgi/";
      mascot_cgi_dir = QDir(mascot_cgi_dir.toQString()).absolutePath();

      mascot_data_dir = getStringOption_("temp_data_directory");

      if (mascot_data_dir.empty())
      {
        writeLogError_("No temp directory specified. Aborting!");
        return ILLEGAL_PARAMETERS;
      }

      writeDebug_(String("Temp directory: ") + mascot_data_dir, 1);
      mascot_data_dir = QDir(mascot_data_dir.toQString()).absolutePath();

      String tmp = mascot_data_dir + "/" + mascot_outfile_name;
      if (!File::writable(tmp))
      {
        writeLogError_(String(" Could not write in temp data directory: ") + tmp + " Aborting!");
        return ILLEGAL_PARAMETERS;
      }
      mascotXML_file_name = mascot_data_dir + "/" + mascot_outfile_name + ".mascotXML";
      pepXML_file_name = mascot_data_dir + "/" + mascot_outfile_name + ".pepXML";
      writeDebug_(String("mascotXML_file_name: ") + mascotXML_file_name, 1);
      writeDebug_(String("pepXML_file_name: ") + pepXML_file_name, 1);
    }

//  contact_person.setName(getStringOption_("contactName", "unknown"));
//  writeDebug_(String("Contact name: ") + contact_person.getName(), 1);
//
//  contact_person.setInstitution(getStringOption_("contactInstitution", "unknown"));
//  writeDebug_(String("Contact institution: ") + contact_person.getInstitution(), 1);
//
//  contact_person.setContactInfo(getStringOption_("contactInfo"));
//  writeDebug_(String("Contact info: ") + contact_person.getContactInfo(), 1);


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    if (!mascot_out)
    {
      mzdata_infile.setLogType(log_type_);
      mzdata_infile.load(inputfile_name, experiment);

      writeDebug_("read " + String(experiment.size()) + " spectra from mzData file", 1);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

      mascot_infile.setInstrument(instrument);
      mascot_infile.setPrecursorMassTolerance(precursor_mass_tolerance);
      mascot_infile.setPeakMassTolerance(peak_mass_tolerance);
      if (!mods.empty())
      {
        mascot_infile.setModifications(mods);
      }
      if (!variable_mods.empty())
      {
        mascot_infile.setVariableModifications(variable_mods);
      }
      mascot_infile.setTaxonomy(taxonomy);
      mascot_infile.setDB(db);
      mascot_infile.setHits(hits);
      mascot_infile.setCleavage(cleavage);
      mascot_infile.setMissedCleavages(missed_cleavages);
      mascot_infile.setMassType(mass_type);
      mascot_infile.setCharges(charges);
      if (!mascot_in)
      {
#ifdef OPENMS_WINDOWSPLATFORM
        /// @todo test this with a real mascot version for windows
        writeLogWarn_(QString("The windows platform version of this tool has not been tested yet! If you encounter problems,") +
                  QString(" please write to the OpenMS mailing list (open-ms-general@lists.sourceforge.net)"));
#endif

        mascot_infile.store(mascot_data_dir + "/" + mascot_infile_name,
                            experiment,
                            "OpenMS search");
        String tmp = logfile;
        tmp = File::absolutePath(tmp);

        writeDebug_("Searching...", 1);
        // calling the Mascot process
        writeDebug_("The Mascot process created the following output:", 1);

        QProcess qp;
        qp.setWorkingDirectory(mascot_cgi_dir.toQString());
        call = " 1 -commandline -f " +
               mascot_data_dir + "/" + mascot_outfile_name + " < " +
               mascot_data_dir + "/" + mascot_infile_name +
#ifdef OPENMS_WINDOWSPLATFORM
               " > " + tmp;
#else
               " >> " + tmp + ";";
#endif
        writeDebug_("CALLING: nph-mascot.exe" + call + "\nCALL Done!    ", 10);
        Int status = qp.execute("nph-mascot.exe", QStringList() << call.toQString());
        if (status != 0)
        {
          writeLogError_("Mascot server problem. Aborting!(Details can be seen in the logfile: \"" + logfile + "\")");
          QFile(String(mascot_data_dir + "/" + mascot_infile_name).toQString()).remove();
          return EXTERNAL_PROGRAM_ERROR;
        }

#ifdef OPENMS_WINDOWSPLATFORM
        call = String("export_dat.pl ") +
#else
        call =  String("export_dat_2.pl ") +
#endif
               " do_export=1 export_format=XML file=" + mascot_data_dir +
               "/" + mascot_outfile_name + " _sigthreshold=" + String(sigthreshold) + " _showsubset=1 show_same_sets=1 show_unassigned=" + String(show_unassigned) +
               " prot_score=" + String(prot_score) + " query_master=1 search_master=1 protein_master=1 peptide_master=1 pep_exp_z=" + String(pep_exp_z) + " pep_score=" + String(pep_score) +
               " pep_homol=" + String(pep_homol) + " query_title=1 pep_ident=" + String(pep_ident) + " pep_seq=1 report=0 " +
               "show_params=1 _showallfromerrortolerant=1 show_header=1 show_queries=1 pep_rank=" + String(pep_rank) + " > " + mascotXML_file_name +

#ifdef OPENMS_WINDOWSPLATFORM
               " && " + " perl export_dat.pl " +
#else
               ";"    + "./export_dat.pl " +
#endif
               " do_export=1 export_format=pepXML file="  + mascot_data_dir +
               "/" + mascot_outfile_name + " _sigthreshold=" + String(sigthreshold) + " _showsubset=1 show_same_sets=1 show_unassigned=" + String(show_unassigned) +
               " prot_score=" + String(prot_score) + " pep_exp_z=" + String(pep_exp_z) + " pep_score=" + String(pep_score) +
               " pep_homol=" + String(pep_homol) + " pep_ident=" + String(pep_ident) + " pep_seq=1 report=0 " +
               "show_params=1 show_header=1 show_queries=1 pep_rank=" + String(pep_rank) + " > " + pepXML_file_name;
        writeDebug_("CALLING: perl " + call + "\nCALL Done!    ", 10);
        status = qp.execute("perl", QStringList() << call.toQString());

        if (status != 0)
        {
          writeLogError_("Mascot server problem. Aborting!(Details can be seen in the logfile: \"" + logfile + "\")");
          QFile(String(mascot_data_dir + "/" + mascot_infile_name).toQString()).remove();
          QFile(mascotXML_file_name.toQString()).remove();
          QFile(pepXML_file_name.toQString()).remove();
          return EXTERNAL_PROGRAM_ERROR;
        }

      }           // from if(!mascot_in)
      else
      {
        if (!boundary.empty())
        {
          mascot_infile.setBoundary(boundary);
        }
        mascot_infile.store(mascot_infile_name,
                            experiment,
                            "OpenMS search");
      }
    }         // from if(!mascot_out)
    if (!mascot_in)
    {
      SpectrumMetaDataLookup lookup;
      if (mascot_out)
      {
        mascotXML_file.load(mascotXML_file_name, protein_identification,
                            identifications, lookup);
      }
      else
      {
        pepXML_file.load(pepXML_file_name, modified_peptides);
        mascotXML_file.load(mascotXML_file_name, protein_identification,
                            identifications, modified_peptides, lookup);
      }

      if (first_dim_rt > 0)
      {
        for (Size i = 0; i < identifications.size(); ++i)
        {
          identifications[i].setMetaValue("first_dim_rt", first_dim_rt);
        }
      }

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
      vector<ProteinIdentification> protein_identifications;
      protein_identifications.push_back(protein_identification);

      // write all (!) parameters as metavalues to the search parameters
      DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_identifications[0].getSearchParameters(), this->getToolPrefix());

      FileHandler().storeIdentifications(outputfile_name,
                        protein_identifications,
                        identifications,
                        {FileTypes::IDXML});

      // Deletion of temporary Mascot files
      if (!mascot_out)
      {
        QFile(String(mascot_data_dir + "/" + mascot_infile_name).toQString()).remove();
        QFile(String(mascot_data_dir + "/" + mascot_outfile_name).toQString()).remove();
        QFile(mascotXML_file_name.toQString()).remove();
        QFile(pepXML_file_name.toQString()).remove();
      }

    }             // from if(!mascot_in)
    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPMascotAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
