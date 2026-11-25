// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHugePages.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/ProteoformDatabase.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/UniPEFFConverter.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>
#include <cctype>
#include <chrono>
#include <fstream>
#include <set>
#include <unordered_map>
#include <limits>
#if defined(__GLIBC__)
  #include <malloc.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
// We do not want this class to show up in the docu:

namespace
{
  /// the format of a protein database file
  enum class DatabaseFormat
  {
    FASTA,
    PEFF,
    UNIPROT_XML
  };

  /// the format of @p file by its content (gzip: a UniProtKB .xml.gz download; '# PEFF'; an XML declaration or <uniprot>), else by its
  /// name (.peff; .xml); FASTA otherwise
  DatabaseFormat databaseFormat_(const std::string& file)
  {
    std::ifstream in(file, std::ios::binary);
    char magic[2] = {0, 0};
    in.read(magic, 2);
    if (in.gcount() == 2 && (unsigned char)magic[0] == 0x1f && (unsigned char)magic[1] == 0x8b) return DatabaseFormat::UNIPROT_XML;
    in.clear();
    in.seekg(0);
    std::string line;
    while (std::getline(in, line))
    {
      StringUtils::trim(line);
      if (! line.empty()) break;
    }
    if (StringUtils::hasPrefix(line, "# PEFF")) return DatabaseFormat::PEFF;
    if (StringUtils::hasPrefix(line, "<?xml") || StringUtils::hasPrefix(line, "<uniprot")) return DatabaseFormat::UNIPROT_XML;
    const std::string name = StringUtils::toLowered(file);
    if (StringUtils::hasSuffix(name, ".peff")) return DatabaseFormat::PEFF;
    if (StringUtils::hasSuffix(name, ".xml") || StringUtils::hasSuffix(name, ".xml.gz")) return DatabaseFormat::UNIPROT_XML;
    return DatabaseFormat::FASTA;
  }
} // namespace

class TOPPFLASHTnT : public TOPPBase
{
public:
  TOPPFLASHTnT():
      TOPPBase("FLASHTnT",
               "FLASHTnT to generate de novo sequence tags from TDP spectrum and match them against proteome DB for proteoform identification.")
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "",
                       "Input mzML: the deconvolved mzML of FLASHDeconv (-out_mzml), or a raw mzML (profile or centroid; recognised by spectra "
                       "without FLASHDeconv's deconvolved masses). A raw mzML is deconvolved here with FLASHDeconv's algorithm (parameters fd:, "
                       "FLASHDeconv's defaults), and its MSn spectra, as read before deconvolution, serve the raw-spectrum scoring as -in_raw "
                       "would (-in_raw is then not needed and not read): the same results as FLASHDeconv and FLASHTnT -in_raw run one after "
                       "the other, without reading the raw file twice or writing the deconvolved mzML. A Thermo .raw file (recognised by its "
                       "extension) is a raw input too: it is read with Thermo's peak picking, so centroided, as ThermoRawFileParser converts it "
                       "by default (needs OpenMS built with WITH_THERMO_RAW and the .NET 8 runtime). "
                       "A Bruker timsTOF .d directory, or a zipped one (.d.zip), is recognised by its name and read with OpenMS's timsTOF reader (DDA-PASEF or DIA-PASEF): the spectra are summed over the ion mobility scans and centroided, the spectra are numbered by retention time (scan=1, 2, ...), and the ion mobility values are not kept (needs OpenMS built with WITH_OPENTIMS, the default). It is a raw input too.", true);
    setValidFormats_("in", ListUtils::create<std::string>("mzML,raw,d"));

    registerInputFile_("db", "<file>", "",
                       "Protein database. Three formats, recognised from the file content (else from the file name):\n"
                       "- FASTA (e.g. .fasta, .fa): one entry per protein. The file must contain the decoys (accessions with the prefix DECOY_).\n"
                       "- PEFF, PSI Extended FASTA (first line '# PEFF', or .peff; e.g. written by UniPEFF -annotation_identifiers): every annotated "
                       "proteoform is a database entry, i.e. each processed form of a protein (full length, mature chains, without the initiator Met, "
                       "signal, transit peptide or propeptide) with up to -peff_max_mods annotated modification units fixed on it (\\ModResPsi, "
                       "\\ModResUnimod; disulfide bonds from \\DisulfideBond pairs). Variable and blind modifications are searched on top. A target "
                       "without a decoy in the file gets one (its residues shuffled within the processed forms, the ends of each form kept, the fixed "
                       "modifications moved with their residues); decoys in the file are used as given. Sequence variants are not used. Every "
                       "entry that the tags, gapped tags or mass ladders support is extended, and the protein keeps the PrSM of a determined mass, "
                       "then the fewest blind modifications, the best score and the most annotated modifications. With -max_blind_mod_count 0 only "
                       "the annotated (and variable) modifications are used: a stricter search (E. coli lysates: about half the proteoforms, almost "
                       "no entrapment hits).\n"
                       "- UniProtKB XML (.xml or .xml.gz, e.g. a UniProt proteome download in XML format): converted on the fly to PEFF as "
                       "'UniPEFF -annotation_identifiers' writes it, then read as PEFF.\n"
                       "A UniProtKB XML file gives the same results as the PEFF converted from it.",
                       true);

    registerInputFile_("target", "<file>", "",
                       "Target proteins: a FASTA file with the proteins known to be in the sample (e.g. the purified protein or a construct), "
                       "without decoys. Needs a FASTA -db with decoys, which measures chance. A target protein whose sequence is in -db is that "
                       "entry; any other is added to the database for this run with a decoy of its own (its residues shuffled, the first two "
                       "kept). The search itself is unchanged: a target protein is identified by the same rule as any other protein. With "
                       "-native, a target protein that the search did not report in a spectrum is characterized there all the same, from its "
                       "sequence tags and sliding ladders in the spectrum (none of either: no PrSM). The PrSM and proteoform files get two columns: "
                       "TargetProtein (0: not a target protein; 1: target protein reported by the search; 2: target protein characterized on "
                       "request only, i.e. not identified by the spectrum) and DecoyMargin (-native: the candidate score minus the score of the "
                       "best decoy candidate of the same source in the spectrum; the margin an identification needs is in the log). Without "
                       "this option the outputs are unchanged.",
                       false);
    setValidFormats_("target", ListUtils::create<std::string>("fasta"));

    registerStringList_("fixed_mod", "<name>", std::vector<std::string>(),
                        "Fixed modifications, by name with the residue (e.g. 'Carbamidomethyl (C)' for alkylated cysteines, 'TMT6plex (K)') or by Unimod "
                        "accession number with the residue (e.g. '4 (C)'; 'UNIMOD:4 (C)' is the same). "
                        "The residue's mass is changed for the whole search: every occurrence in targets and decoys, in sequence tags (a tag reads "
                        "the modified residue as its letter), mass ladders, the extension, the raw-spectrum scoring and the reported masses. Only "
                        "modifications of a residue anywhere, one per residue (not terminal ones). Written as ProForma global modifications, "
                        "e.g. <[+57.021464|UNIMOD:4]@C>. With a PEFF database an annotated modification of such a residue adds to it.",
                        false);

    registerIntOption_("peff_max_mods", "<number>", 1,
                       "PEFF or UniProtKB XML database: at most this many annotated modification units fixed on an entry (0: the processed forms "
                       "only). A unit is a residue with one or more of its annotated modifications on distinct sites of it, or all disulfide bonds of "
                       "a form.", false);
    setMinInt_("peff_max_mods", 0);
    registerStringList_("peff_mod_masses", "<accession=mass>",
                        ListUtils::create<std::string>("MOD:01819=100.016044,MOD:01893=86.000394,MOD:01614=329.05252,MOD:01956=15.994915,"
                                                  "MOD:00657=14.999666,MOD:00017=skip,MOD:00223=skip,MOD:01976=skip,MOD:01780=skip"),
                        "PEFF or UniProtKB XML database: masses of modification accessions that replace those of the modification database (PSI-MOD, Unimod), or "
                        "'skip' to leave the modification out. The default fills the gaps of the bundled PSI-MOD for E. coli (UniProt ptmlist "
                        "masses). A modification without a mass is left out (counted in the log).",
                        false);

    registerInputFile_("in_raw", "<file>", "",
                       "Raw mzML file (profile or centroid), Thermo .raw file or Bruker .d directory FLASHDeconv was run on, read as -in reads it. With it, every PrSM (targets and decoys) is annotated on its raw "
                       "spectrum (fragment isotope envelopes at every charge, as FLASHAnnotator) before the FDR, and ranked by its extension score "
                       "plus 8 times its raw score. The raw annotation columns are written in the PrSM and proteoform files.",
                       false);
    setValidFormats_("in_raw", ListUtils::create<std::string>("mzML,raw,d"));

    registerOutputFile_("out_pro", "<file>", "", "Output proteoform-level tsv file containing proteoform IDs. ", true);
    setValidFormats_("out_pro", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_prsm", "<file>", "", "Output PrSM-level tsv file containing PrSMs.");
    setValidFormats_("out_prsm", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_tag", "<file>", "", "Output tag-level tsv file containing matched tags.", false);
    setValidFormats_("out_tag", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_protein", "<file>", "", "Output protein-level tsv file containing the candidate proteins and "
                        "their scores, as ranked before proteoform characterization.", false);
    setValidFormats_("out_protein", ListUtils::create<std::string>("tsv"));

    auto tnt_param = FLASHTnTAlgorithm().getDefaults();
    tnt_param.removeAll("tag:");
    tnt_param.removeAll("ex:");

    registerFullParam_(tnt_param);

    //registerSubsection_("tnt", "TnT algorithm parameters");
    registerSubsection_("tag", "Tag algorithm parameters");
    registerSubsection_("ex", "Extension algorithm parameters");
    registerSubsection_("fd", "FLASHDeconv parameters, used when -in is a raw mzML (FLASHDeconv's FD: and SD: parameters; its feature tracing and "
                              "isobaric quantification keep their defaults, and decoy masses and the precursor by the fragments are made with -precursor_FDR)");
  }

  Param getSubsectionDefaults_(const std::string& prefix) const override
  {
    if (prefix == "tag")
    {
      auto tnt_param = FLASHTnTAlgorithm().getDefaults();
      return tnt_param.copy("tag:", true);
    }
    else if (prefix == "ex") {
      auto tnt_param = FLASHTnTAlgorithm().getDefaults();
      return tnt_param.copy("ex:", true);
    }
    else if (prefix == "fd") { return fdDefaults_(); }
    else { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", prefix); }
  }

  // FLASHDeconv's parameters as FLASHTnT exposes them. Kept at their defaults: feature tracing (ft:), isobaric quantification (iq:), the
  // RNA averagine and the beta isotope error tolerance of FLASHDeconv's mass FDR; report_FDR and precursor_FDR are set by FLASHTnT (-precursor_FDR).
  static Param fdDefaults_()
  {
    Param p = FLASHDeconvAlgorithm().getDefaults();
    p.removeAll("ft:");
    p.removeAll("iq:");
    p.remove("use_RNA_averagine");
    p.remove("allowed_isotope_error");
    p.remove("report_FDR");
    p.remove("precursor_FDR");
    return p;
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**) override
  {
    OPENMS_LOG_INFO << "Initializing ... " << endl;
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    std::string in_file = getStringOption_("in");
    const std::string in_db = getStringOption_("db");
    const std::string in_target = getStringOption_("target");

    // Fixed modifications change the residue masses of the search, the in-process amino acid set (ResidueDB), before anything reads them:
    // every mass of the search -- tags, ladders, the extension, the raw annotation, the reported masses -- is that residue's.
    std::string fixed_proforma;
    {
      std::set<char> fixed_residues;
      for (const std::string& name : getStringList_("fixed_mod"))
      {
        if (name.empty()) continue;
        const ResidueModification* mod = nullptr;
        try
        {
          mod = FLASHHelperClasses::resolveModification(name); // a name or a Unimod accession, with its residue
        }
        catch (const Exception::BaseException& e)
        {
          OPENMS_LOG_ERROR << "-fixed_mod " << name << ": " << e.what() << endl;
          return ILLEGAL_PARAMETERS;
        }
        const char origin = mod->getOrigin();
        if (mod->getTermSpecificity() != ResidueModification::ANYWHERE || ! std::isupper((unsigned char)origin) || origin == 'X')
        {
          OPENMS_LOG_ERROR << "-fixed_mod " << name << ": only a modification of a residue anywhere can be fixed (not a terminal one)." << endl;
          return ILLEGAL_PARAMETERS;
        }
        if (! fixed_residues.insert(origin).second)
        {
          OPENMS_LOG_ERROR << "-fixed_mod: two fixed modifications of residue " << origin << "." << endl;
          return ILLEGAL_PARAMETERS;
        }
        Residue* residue = const_cast<Residue*>(ResidueDB::getInstance()->getResidue((unsigned char)origin));
        const double before = residue->getMonoWeight(Residue::Internal);
        residue->setMonoWeight(residue->getMonoWeight(Residue::Full) + mod->getDiffMonoMass());
        std::string accession = mod->getUniModAccession();
        StringUtils::substitute(accession, "UniMod:", "UNIMOD:");
        fixed_proforma += "<[" + ProteoformDatabase::proformaTag(accession, mod->getDiffMonoMass()) + "]@" + std::string(1, origin) + ">";
        OPENMS_LOG_INFO << "Fixed modification " << mod->getFullId() << ": residue " << origin << " " << before << " -> "
                        << residue->getMonoWeight(Residue::Internal) << " Da" << endl;
      }
    }
    FLASHTnTFile::setGlobalFixedModifications(fixed_proforma);
    // variable modifications: checked here, so that a name or an accession the database does not know stops the run with a clear message
    for (const auto& name : getParam_().getValue("var_mod").toStringVector())
    {
      if (StringUtils::trimmed(name).empty()) continue;
      try
      {
        bool protein_n_term = false;
        FLASHHelperClasses::resolveVariableModification(name, protein_n_term);
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "-var_mod " << name << ": " << e.what() << endl;
        return ILLEGAL_PARAMETERS;
      }
    }
    std::string in_raw = getStringOption_("in_raw");

    std::string out_tag_file = getStringOption_("out_tag");
    std::string out_prsm_file = getStringOption_("out_prsm");
    std::string out_pro_file = getStringOption_("out_pro");
    std::string out_protein_file = getStringOption_("out_protein");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    MSExperiment map;
    MzMLFile mzml;

    // reading mzMLs with m/z and rt criteria.

    FLASHHelperClasses::loadExperiment(in_file, map, mzml.getOptions(), log_type_);

    // A raw mzML (no spectrum with FLASHDeconv's deconvolved masses) is deconvolved here, as FLASHDeconv does before writing its
    // deconvolved mzML (FLASHDeconvSpectrumFile::toDeconvolvedMap, uncharged masses), and its MSn spectra are kept for the raw-spectrum
    // scoring as they were read: FLASHDeconv removes low peaks from the spectra it is given, which the scoring of a separate -in_raw read
    // never saw. The spectra are keyed by scan number as FLASHAnnotatorAlgorithm::loadSpectra reads them from the file.
    std::map<int, MSSpectrum> raw_msn_spectra;
    const bool raw_input = std::none_of(map.begin(), map.end(), [](const MSSpectrum& spec) { return spec.metaValueExists("DeconvMassInfo"); });
    if (raw_input)
    {
      OPENMS_LOG_INFO << "-in has no deconvolved masses: deconvolving it with FLASHDeconv's algorithm (fd: parameters)" << endl;
      std::string type_accession = "MS:1000768";
      if (! map.getSourceFiles().empty() && ! map.getSourceFiles()[0].getNativeIDTypeAccession().empty())
      {
        type_accession = map.getSourceFiles()[0].getNativeIDTypeAccession();
      }
      auto scanOf = [&type_accession](const MSSpectrum& spec) {
        std::vector<std::string> native_ids;
        StringUtils::split(spec.getNativeID(), ",", native_ids);
        return native_ids.empty() ? -1 : SpectrumLookup::extractScanNumber(native_ids.back(), type_accession);
      };
      // scan numbers from the native IDs; if an MSn spectrum has none, every spectrum's index + 1 instead (unless its native ID has one)
      const bool from_index = std::any_of(map.begin(), map.end(), [&](const MSSpectrum& spec) { return spec.getMSLevel() >= 2 && scanOf(spec) < 0; });
      for (Size i = 0; i < map.size(); i++)
      {
        const auto& spec = map[i];
        if (spec.getMSLevel() < 2) continue;
        int scan = scanOf(spec);
        if (scan < 0)
        {
          if (! from_index) continue;
          scan = (int)i + 1;
        }
        raw_msn_spectra.emplace(scan, spec); // the first spectrum of a scan number, as the file reader keeps it
      }
      FLASHDeconvAlgorithm fd;
      Param fd_param = getParam_().copy("fd:", true);
      // FLASHDeconv's precursor_FDR (the precursor by the fragments and the decoy precursors of the isolation windows) with -precursor_FDR
      const std::string precursor_fdr = getParam_().getValue("precursor_FDR").toString() == "true" ? "true" : "false";
      fd_param.setValue("report_FDR", "false");
      fd_param.setValue("precursor_FDR", precursor_fdr);
      fd.setParameters(fd_param);
      std::vector<DeconvolvedSpectrum> deconvolved_spectra;
      std::vector<FLASHHelperClasses::MassFeature> deconvolved_features;
      fd.run(map, deconvolved_spectra, deconvolved_features);
      MSExperiment deconvolved_map;
      FLASHDeconvSpectrumFile::toDeconvolvedMap(map, deconvolved_spectra, 0, fd.getTolerances(), deconvolved_map);
      map = std::move(deconvolved_map);
      if (! in_raw.empty()) OPENMS_LOG_WARN << "-in_raw is not read: the raw spectra are those of -in" << endl;
    }

    auto tnt_param = getParam_();
    tnt_param.remove("in");
    tnt_param.remove("db");
    tnt_param.remove("target");
    tnt_param.remove("fixed_mod");
    tnt_param.remove("peff_max_mods");
    tnt_param.remove("peff_mod_masses");
    tnt_param.remove("in_raw");
    tnt_param.removeAll("fd:");
    tnt_param.remove("out_tag");
    tnt_param.remove("out_prsm");
    tnt_param.remove("out_pro");
    tnt_param.remove("out_protein");
    tnt_param.remove("log");
    tnt_param.remove("debug");
    tnt_param.remove("threads");
    tnt_param.remove("no_progress");
    tnt_param.remove("force");
    tnt_param.remove("test");

    double pro_fdr = tnt_param.getValue("pro_fdr");
    // the chimeric columns are written only with the residual search, so that other outputs keep their format
    const bool chimeric = (int)tnt_param.getValue("chimeric_rounds") > 0 || tnt_param.getValue("native").toString() == "true"; // -native searches a residual round
    OPENMS_LOG_INFO << "Finding sequence tags from deconvolved spectra ..." << endl;

    std::vector<FASTAFile::FASTAEntry> fasta_entry;
    const auto db_read_start = std::chrono::steady_clock::now();
    // PEFF (or UniProtKB XML, converted to PEFF): its proteoform entries (ProteoformDatabase); FASTA: its entries are the proteins
    ProteoformDatabase database;
    const DatabaseFormat db_format = databaseFormat_(in_db);
    std::string peff_file = in_db;
    if (db_format == DatabaseFormat::UNIPROT_XML)
    {
      peff_file = TempFiles::getTemporaryFile();
      UniPEFFConverter::Options conversion;
      conversion.annotation_identifiers = true; // the disulfide bond pairs
      const Size converted = UniPEFFConverter::convert(in_db, peff_file, conversion);
      OPENMS_LOG_INFO << "UniProtKB XML database " << in_db << ": " << converted << " entries converted to PEFF" << endl;
    }
    const bool peff = db_format != DatabaseFormat::FASTA;
    if (peff)
    {
      database.loadPEFF(peff_file, getIntOption_("peff_max_mods"), getStringList_("peff_mod_masses"), fasta_entry);
      const auto& st = database.getLoadStatistics();
      Size targets = 0, with_mods = 0;
      for (const auto& p : database.proteins) targets += ! p.decoy;
      for (const auto& e : database.entries) with_mods += ! e.fixed.empty();
      OPENMS_LOG_INFO << "PEFF database: " << st.input_targets << " target and " << st.input_decoys << " decoy proteins in the file, "
                      << st.generated_decoys << " decoys generated; " << targets << " target and " << database.proteins.size() - targets
                      << " decoy proteins, " << database.entries.size() << " entries (" << with_mods << " with fixed modifications; "
                      << st.forms << " processed forms of the targets, peff_max_mods " << getIntOption_("peff_max_mods") << "); "
                      << st.modification_sites << " modification sites and " << st.disulfide_bonds << " disulfide bonds of the targets used"
                      << endl;
      if (st.unpaired_half_cystines > 0)
      {
        OPENMS_LOG_WARN << "PEFF database: " << st.unpaired_half_cystines
                        << " half cystines without a disulfide bond pair are left out (write the PEFF with UniPEFF -annotation_identifiers)" << endl;
      }
      if (st.length_mismatches > 0)
      {
        OPENMS_LOG_WARN << "PEFF database: " << st.length_mismatches << " entries whose \\Length differs from their sequence" << endl;
      }
      if (st.empty_sequences > 0) { OPENMS_LOG_WARN << "PEFF database: " << st.empty_sequences << " entries without residues left out" << endl; }
      if (st.unpaired_decoys > 0)
      {
        OPENMS_LOG_WARN << "PEFF database: " << st.unpaired_decoys << " decoys of the file without a target of their name or accession" << endl;
      }
      if (st.unmirrored_decoys > 0)
      {
        OPENMS_LOG_WARN << "PEFF database: " << st.unmirrored_decoys
                        << " decoys of the file whose entries do not mirror their target's (their number or masses); without decoys in the "
                           "file, mirrored ones are generated"
                        << endl;
      }
      auto siteList = [](const std::map<std::string, Size>& counts, Size& sites) {
        std::string list;
        sites = 0;
        for (const auto& [acc, n] : counts)
        {
          sites += n;
          list += (list.empty() ? "" : ", ") + acc + " (" + StringUtils::toStr(n) + ")";
        }
        return list;
      };
      if (! st.unresolved.empty())
      {
        Size sites = 0;
        const std::string list = siteList(st.unresolved, sites);
        OPENMS_LOG_INFO << "PEFF database: " << sites << " modification sites without a mass left out (peff_mod_masses): " << list << endl;
      }
      if (! st.wrong_residue.empty())
      {
        Size sites = 0;
        const std::string list = siteList(st.wrong_residue, sites);
        OPENMS_LOG_WARN << "PEFF database: " << sites << " modification sites on a residue that is not the modification's origin left out: " << list
                        << endl;
      }
    }
    else
    {
      FASTAFile fasta_file;
      fasta_file.load(in_db, fasta_entry);
    }
    // Target proteins (-target): the entries of -db with their sequence, else added with a decoy of their own (shuffled, the first two
    // residues kept, as the protein decoys of DecoyDatabase with which the native decoy margin was set)
    std::set<std::string> target_proteins;
    if (! in_target.empty())
    {
      if (peff)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "-target needs a FASTA -db (not PEFF or UniProtKB XML)");
      }
      std::unordered_map<std::string, Size> entry_of_sequence, entry_of_identifier;
      Size db_decoys = 0;
      for (Size i = 0; i < fasta_entry.size(); i++)
      {
        const bool decoy = StringUtils::hasPrefix(fasta_entry[i].identifier, "DECOY");
        db_decoys += decoy;
        entry_of_identifier.emplace(fasta_entry[i].identifier, i);
        if (! decoy) entry_of_sequence.emplace(fasta_entry[i].sequence, i);
      }
      if (db_decoys == 0)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "-target needs decoys in -db (accessions with the prefix DECOY_): they measure chance");
      }
      std::vector<FASTAFile::FASTAEntry> targets;
      FASTAFile().load(in_target, targets);
      Size in_db_count = 0, added = 0;
      for (auto& t : targets)
      {
        if (t.sequence.empty()) continue;
        if (StringUtils::hasPrefix(t.identifier, "DECOY"))
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "-target holds a decoy entry: " + t.identifier);
        }
        const auto same = entry_of_sequence.find(t.sequence);
        if (same != entry_of_sequence.end())
        {
          target_proteins.insert(fasta_entry[same->second].identifier);
          in_db_count++;
          continue;
        }
        if (entry_of_identifier.count(t.identifier) > 0 || entry_of_identifier.count("DECOY_" + t.identifier) > 0)
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "-target entry has the identifier of a -db entry with another sequence: " + t.identifier);
        }
        FASTAFile::FASTAEntry decoy("DECOY_" + t.identifier, t.description, t.sequence);
        if (decoy.sequence.size() > 3)
        {
          Math::RandomShuffler shuffler(1);
          shuffler.portable_random_shuffle(decoy.sequence.begin() + 2, decoy.sequence.end());
        }
        entry_of_identifier.emplace(t.identifier, fasta_entry.size());
        entry_of_sequence.emplace(t.sequence, fasta_entry.size());
        target_proteins.insert(t.identifier);
        fasta_entry.push_back(t);
        entry_of_identifier.emplace(decoy.identifier, fasta_entry.size());
        fasta_entry.push_back(decoy);
        added++;
      }
      OPENMS_LOG_INFO << "Target proteins: " << target_proteins.size() << " (" << in_db_count << " in the database, " << added
                      << " added with a decoy each)" << endl;
      if (target_proteins.empty()) { OPENMS_LOG_WARN << "-target holds no protein" << endl; }
    }
    const bool with_targets = ! target_proteins.empty();
    OPENMS_LOG_INFO << "Database read: " << fasta_entry.size() << " entries in "
                    << std::chrono::duration<double>(std::chrono::steady_clock::now() - db_read_start).count() << " s" << endl;
    // Run here!
    OPENMS_LOG_INFO << "Processing : " << in_file << endl;

    fstream out_tagger_stream;
    fstream out_prsm_stream;
    fstream out_pro_stream;
    fstream out_protein_stream;

    if (! out_tag_file.empty())
    {
      out_tagger_stream = fstream(out_tag_file, fstream::out);
      FLASHTnTFile::writeTagHeader(out_tagger_stream);
    }

    if (! out_prsm_file.empty())
    {
      out_prsm_stream = fstream(out_prsm_file, fstream::out); // the header after the run: its IonTypes column depends on the ion types
    }

    if (! out_pro_file.empty())
    {
      out_pro_stream = fstream(out_pro_file, fstream::out);
    }

    if (! out_protein_file.empty())
    {
      out_protein_stream = fstream(out_protein_file, fstream::out);
      FLASHTnTFile::writeProteinHeader(out_protein_stream);
    }

    FLASHTnTAlgorithm tnt;
    tnt.setParameters(tnt_param);
    tnt.setTargetProteins(target_proteins);
    if (raw_input)
    {
      tnt.setRawSpectrumLoader([&raw_msn_spectra](const std::set<int>& scans, std::map<int, MSSpectrum>& spectra) {
        for (int scan : scans)
        {
          const auto it = raw_msn_spectra.find(scan);
          if (it != raw_msn_spectra.end() && spectra.find(scan) == spectra.end()) spectra.emplace(scan, it->second);
        }
      });
    }
    else if (! in_raw.empty())
    {
      tnt.setRawSpectrumLoader([&in_raw](const std::set<int>& scans, std::map<int, MSSpectrum>& spectra) {
        FLASHAnnotatorAlgorithm::loadSpectra(in_raw, scans, spectra);
      });
    }
    tnt.run(map, fasta_entry, peff ? &database : nullptr);
    const auto& proteoform_hits = tnt.getProteoforms();

    OPENMS_LOG_INFO << "FLASHTnT run complete. Now writing the results in output files ..." << endl;
    if (! out_tag_file.empty())
    {
      /// TODO flanking_mass_tol should be calculated within FLASHTagger
      FLASHTnTFile::writeTags(tnt, tnt.getMaxTotalModificationMass(), out_tagger_stream);
      out_tagger_stream.close();
    }
    if (! out_prsm_file.empty())
    {
      FLASHTnTFile::writePrSMHeader(out_prsm_stream, chimeric, with_targets, tnt.getIonTypesFromFragmentation());
      FLASHTnTFile::writePrSMs(proteoform_hits, out_prsm_stream, chimeric, with_targets, tnt.getIonTypesFromFragmentation());
      out_prsm_stream.close();
    }
    if (! out_pro_file.empty())
    {
      FLASHTnTFile::writeProHeader(out_pro_stream, chimeric, with_targets, tnt.getIonTypesFromFragmentation());
      FLASHTnTFile::writeProteoforms(proteoform_hits, out_pro_stream, pro_fdr, chimeric, with_targets, tnt.getIonTypesFromFragmentation());
      out_pro_stream.close();
    }
    if (! out_protein_file.empty())
    {
      FLASHTnTFile::writeProteins(tnt.getProteinHits(), out_protein_stream);
      out_protein_stream.close();
    }

    // the results and the input are freed with all threads instead of by the serial destructors
    tnt.releaseResults();
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < (SignedSize)map.size(); ++i)
    {
      map[i] = MSSpectrum();
    }
    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  useTransparentHugePages(argv);
#if defined(__GLIBC__)
  // Keep freed heap memory in the per-thread arenas instead of trimming it back to the kernel after every
  // work item: trim, refault and heap regrowth serialise on the process-wide mmap lock at high thread counts.
  mallopt(M_TRIM_THRESHOLD, std::numeric_limits<int>::max());
  mallopt(M_TOP_PAD, 256 << 20);
#endif
  TOPPFLASHTnT tool;
  return tool.main(argc, argv);
}
