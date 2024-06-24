// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/SpectrumLookup.h>

#include <QtCore/QFileInfo>
#include <QtCore/QRegularExpression>

#include <iomanip>     // setw

#define HIGH_PRECISION 5
#define LOW_PRECISION 3

using namespace std;

namespace OpenMS
{

  MascotGenericFile::MascotGenericFile() :
    ProgressLogger(), DefaultParamHandler("MascotGenericFile"), mod_group_map_()
  {
    defaults_.setValue("database", "MSDB", "Name of the sequence database");
    defaults_.setValue("search_type", "MIS", "Name of the search type for the query", {"advanced"});
    defaults_.setValidStrings("search_type", {"MIS","SQ","PMF"});
    defaults_.setValue("enzyme", "Trypsin", "The enzyme descriptor to the enzyme used for digestion. (Trypsin is default, None would be best for peptide input or unspecific digestion, for more please refer to your mascot server).");
    defaults_.setValue("instrument", "Default", "Instrument definition which specifies the fragmentation rules");
    defaults_.setValue("missed_cleavages", 1, "Number of missed cleavages allowed for the enzyme");
    defaults_.setMinInt("missed_cleavages", 0);
    defaults_.setValue("precursor_mass_tolerance", 3.0, "Tolerance of the precursor peaks");
    defaults_.setMinFloat("precursor_mass_tolerance", 0.0);
    defaults_.setValue("precursor_error_units", "Da", "Units of the precursor mass tolerance");
    defaults_.setValidStrings("precursor_error_units", {"%","ppm","mmu","Da"});
    defaults_.setValue("fragment_mass_tolerance", 0.3, "Tolerance of the peaks in the fragment spectrum");
    defaults_.setMinFloat("fragment_mass_tolerance", 0.0);
    defaults_.setValue("fragment_error_units", "Da", "Units of the fragment peaks tolerance");
    defaults_.setValidStrings("fragment_error_units", {"mmu","Da"});
    defaults_.setValue("charges", "1,2,3", "Charge states to consider, given as a comma separated list of integers (only used for spectra without precursor charge information)");
    defaults_.setValue("taxonomy", "All entries", "Taxonomy specification of the sequences");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("fixed_modifications", std::vector<std::string>(), "List of fixed modifications, according to UniMod definitions.");
    defaults_.setValidStrings("fixed_modifications", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("variable_modifications", std::vector<std::string>(), "Variable modifications given as UniMod definitions.");
    defaults_.setValidStrings("variable_modifications", ListUtils::create<std::string>(all_mods));

    // special modifications, see "updateMembers_" method below:
    defaults_.setValue("special_modifications", "Cation:Na (DE),Deamidated (NQ),Oxidation (HW),Phospho (ST),Sulfo (ST)", "Modifications with specificity groups that are used by Mascot and have to be treated specially", {"advanced"});
    // list from Mascot 2.4; there's also "Phospho (STY)", but that can be
    // represented using "Phospho (ST)" and "Phospho (Y)"

    defaults_.setValue("mass_type", "monoisotopic", "Defines the mass type, either monoisotopic or average");
    defaults_.setValidStrings("mass_type", {"monoisotopic","average"});
    defaults_.setValue("number_of_hits", 0, "Number of hits which should be returned, if 0 AUTO mode is enabled.");
    defaults_.setMinInt("number_of_hits", 0);
    defaults_.setValue("skip_spectrum_charges", "false", "Sometimes precursor charges are given for each spectrum but are wrong, setting this to 'true' does not write any charge information to the spectrum, the general charge information is however kept.");
    defaults_.setValidStrings("skip_spectrum_charges", {"true","false"});
    defaults_.setValue("decoy", "false", "Set to true if mascot should generate the decoy database.");
    defaults_.setValidStrings("decoy", {"true","false"});

    defaults_.setValue("search_title", "OpenMS_search", "Sets the title of the search.", {"advanced"});
    defaults_.setValue("username", "OpenMS", "Sets the username which is mentioned in the results file.", {"advanced"});
    defaults_.setValue("email", "", "Sets the email which is mentioned in the results file. Note: Some server require that a proper email is provided.");

    // the next section should not be shown to TOPP users
    Param p;
    p.setValue("format", "Mascot generic", "Sets the format type of the peak list, this should not be changed unless you write the header only.", {"advanced"});
    p.setValidStrings("format", {"Mascot generic","mzData (.XML)","mzML (.mzML)"}); // Mascot's HTTP interface supports more, but we don't :)
    p.setValue("boundary", "GZWgAaYKjHFeUaLOLEIOMq", "MIME boundary for parameter header (if using HTTP format)", {"advanced"});
    p.setValue("HTTP_format", "false", "Write header with MIME boundaries instead of simple key-value pairs. For HTTP submission only.", {"advanced"});
    p.setValidStrings("HTTP_format", {"true","false"});
    p.setValue("content", "all", "Use parameter header + the peak lists with BEGIN IONS... or only one of them.", {"advanced"});
    p.setValidStrings("content", {"all","peaklist_only","header_only"});
    defaults_.insert("internal:", p);

    defaultsToParam_();
  }

  MascotGenericFile::~MascotGenericFile() = default;

  void MascotGenericFile::updateMembers_()
  {
    // special cases for specificity groups: OpenMS uses e.g. "Deamidated (N)"
    // and "Deamidated (Q)", but Mascot only understands "Deamidated (NQ)"
    String special_mods = param_.getValue("special_modifications").toString();
    vector<String> mod_groups = ListUtils::create<String>(special_mods);
    for (StringList::const_iterator mod_it = mod_groups.begin();
         mod_it != mod_groups.end(); ++mod_it)
    {
      String mod = mod_it->prefix(' ');
      String residues = mod_it->suffix('(').prefix(')');
      for (String::const_iterator res_it = residues.begin();
           res_it != residues.end(); ++res_it)
      {
        mod_group_map_[mod + " (" + String(*res_it) + ")"] = *mod_it;
      }
    }
  }

  void MascotGenericFile::store(const String& filename, const PeakMap& experiment, bool compact)
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::MGF))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::MGF) + "'");
    }

    if (!File::writable(filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    ofstream os(filename.c_str());
    store(os, filename, experiment, compact);
    os.close();
  }

  void MascotGenericFile::store(ostream& os, const String& filename, const PeakMap& experiment, bool compact)
  {
    // stream formatting may get changed, so back up:
    const ios_base::fmtflags old_flags = os.flags();
    const streamsize old_precision = os.precision();

    store_compact_ = compact;
    if (param_.getValue("internal:content") != "peaklist_only")
      writeHeader_(os);
    if (param_.getValue("internal:content") != "header_only")
      writeMSExperiment_(os, filename, experiment);

    // reset formatting:
    os.flags(old_flags);
    os.precision(old_precision);
  }

  void MascotGenericFile::writeParameterHeader_(const String& name, ostream& os)
  {
    if (param_.getValue("internal:HTTP_format") == "true")
    {
      os << "--" << param_.getValue("internal:boundary") << "\n" << "Content-Disposition: form-data; name=\"" << name << "\"" << "\n\n";
    }
    else
    {
      os << name << "=";
    }
  }

  void MascotGenericFile::writeModifications_(const vector<String>& mods,
                                              ostream& os, bool variable_mods)
  {
    String tag = variable_mods ? "IT_MODS" : "MODS";
    // @TODO: remove handling of special cases when a general solution for
    // specificity groups in UniMod is implemented (ticket #387)
    set<String> filtered_mods;
    for (StringList::const_iterator it = mods.begin(); it != mods.end(); ++it)
    {
      map<String, String>::iterator pos = mod_group_map_.find(*it);
      if (pos == mod_group_map_.end())
      {
        filtered_mods.insert(*it);
      }
      else
      {
        filtered_mods.insert(pos->second);
      }
    }
    for (set<String>::const_iterator it = filtered_mods.begin();
         it != filtered_mods.end(); ++it)
    {
      writeParameterHeader_(tag, os);
      os << *it << "\n";
    }
  }

  void MascotGenericFile::writeHeader_(ostream& os)
  {
    // search title
    if (param_.getValue("search_title") != "")
    {
      writeParameterHeader_("COM", os);
      os << param_.getValue("search_title") << "\n";
    }

    // user name
    writeParameterHeader_("USERNAME", os);
    os << param_.getValue("username") << "\n";

    // email
    if (!param_.getValue("email").toString().empty())
    {
      writeParameterHeader_("USEREMAIL", os);
      os << param_.getValue("email") << "\n";
    }

    // format
    writeParameterHeader_("FORMAT", os); // make sure this stays within the first 5 lines of the file, since we use it to recognize our own MGF files in case their file suffix is not MGF
    os << param_.getValue("internal:format") << "\n";

    // precursor mass tolerance unit : Da
    writeParameterHeader_("TOLU", os);
    os << param_.getValue("precursor_error_units") << "\n";

    // ion mass tolerance unit : Da
    writeParameterHeader_("ITOLU", os);
    os << param_.getValue("fragment_error_units") << "\n";

    // format version
    writeParameterHeader_("FORMVER", os);
    os << "1.01" << "\n";

    // db name
    writeParameterHeader_("DB", os);
    os << param_.getValue("database") << "\n";

    // decoys
    if (param_.getValue("decoy").toBool() == true)
    {
      writeParameterHeader_("DECOY", os);
      os << 1 << "\n";
    }

    // search type
    writeParameterHeader_("SEARCH", os);
    os << param_.getValue("search_type") << "\n";

    // number of peptide candidates in the list
    writeParameterHeader_("REPORT", os);
    UInt num_hits((UInt)param_.getValue("number_of_hits"));
    if (num_hits != 0)
    {
      os << param_.getValue("number_of_hits") << "\n";
    }
    else
    {
      os << "AUTO" << "\n";
    }

    // cleavage enzyme
    writeParameterHeader_("CLE", os);
    os << param_.getValue("enzyme") << "\n";

    // average/monoisotopic
    writeParameterHeader_("MASS", os);
    os << param_.getValue("mass_type") << "\n";

    // fixed modifications
    StringList fixed_mods = ListUtils::toStringList<std::string>(param_.getValue("fixed_modifications"));
    writeModifications_(fixed_mods, os);

    // variable modifications
    StringList var_mods = ListUtils::toStringList<std::string>(param_.getValue("variable_modifications"));
    writeModifications_(var_mods, os, true);

    // instrument
    writeParameterHeader_("INSTRUMENT", os);
    os << param_.getValue("instrument") << "\n";

    // missed cleavages
    writeParameterHeader_("PFA", os);
    os << param_.getValue("missed_cleavages") << "\n";

    // precursor mass tolerance
    writeParameterHeader_("TOL", os);
    os << param_.getValue("precursor_mass_tolerance") << "\n";

    // ion mass tolerance_
    writeParameterHeader_("ITOL", os);
    os << param_.getValue("fragment_mass_tolerance") << "\n";

    // taxonomy
    writeParameterHeader_("TAXONOMY", os);
    os << param_.getValue("taxonomy") << "\n";

    // charge
    writeParameterHeader_("CHARGE", os);
    os << param_.getValue("charges") << "\n";
  }

  void MascotGenericFile::writeSpectrum(ostream& os, const PeakSpectrum& spec, const String& filename, const String& native_id_type_accession)
  {
    Precursor precursor;
    if (!spec.getPrecursors().empty())
    {
      precursor = spec.getPrecursors()[0];
    }
    if (spec.getPrecursors().size() > 1)
    {
      cerr << "Warning: The spectrum written to Mascot file has more than one precursor. The first precursor is used!\n";
    }
    if (spec.size() >= 10000)
    {
      String msg = "Spectrum to be written as MGF has " + String(spec.size()) +
        " peaks; the upper limit is 10,000. Only centroided data is allowed - this is most likely profile data.";
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       msg);
    }
    double mz(precursor.getMZ()), rt(spec.getRT());

    if (mz == 0)
    {
      // retention time
      cout << "No precursor m/z information for spectrum with rt " << rt
           << " present, skipping spectrum!\n";
    }
    else
    {
      os << "\n";
      os << "BEGIN IONS\n";
      if (!store_compact_)
      {
        // if a TITLE is available, it was (most likely) parsed from an MGF
        // or generated to be written out in an MGF
        if (spec.metaValueExists("TITLE"))
        {
          os << "TITLE=" << spec.getMetaValue("TITLE") << "\n";
        }
        else
        {
          os << "TITLE=" << precisionWrapper(mz) << "_" << precisionWrapper(rt)
             << "_" << spec.getNativeID() << "_" << filename << "\n";
        }
        os << "PEPMASS=" << precisionWrapper(mz) <<  "\n";
        os << "RTINSECONDS=" << precisionWrapper(rt) << "\n";
        if (native_id_type_accession == "UNKNOWN")
        {
          os << "SCANS=" << spec.getNativeID().substr(spec.getNativeID().find_last_of("=")+1) << "\n";
        }
        else
        {
          os << "SCANS=" << SpectrumLookup::extractScanNumber(spec.getNativeID(), native_id_type_accession) << "\n";
        }
      }
      else
      {
        // if a TITLE is available, it was (most likely) parsed from an MGF
        // or generated to be written out in an MGF
        if (spec.metaValueExists("TITLE"))
        {
          os << "TITLE=" << spec.getMetaValue("TITLE") << "\n";
        }
        else
        {
          os << "TITLE=" << fixed << setprecision(HIGH_PRECISION) << mz << "_"
             << setprecision(LOW_PRECISION) << rt << "_"
             << spec.getNativeID() << "_" << filename << "\n";
        }
        os << "PEPMASS=" << setprecision(HIGH_PRECISION) << mz << "\n";
        os << "RTINSECONDS=" << setprecision(LOW_PRECISION) << rt << "\n";
        if (native_id_type_accession == "UNKNOWN")
        {
          os << "SCANS=" << spec.getNativeID().substr(spec.getNativeID().find_last_of("=")+1) << "\n";
        }
        else
        {
          os << "SCANS=" << SpectrumLookup::extractScanNumber(spec.getNativeID(), native_id_type_accession) << "\n";
        }
      }

      int charge(precursor.getCharge());

      if (charge != 0)
      {
        bool skip_spectrum_charges(param_.getValue("skip_spectrum_charges").toBool());
        if (!skip_spectrum_charges)
        {
          String cs = charge < 0 ? "-" : "+";
          os << "CHARGE=" << charge << cs << "\n";
        }
      }

      if (!store_compact_)
      {
        for (PeakSpectrum::const_iterator it = spec.begin(); it != spec.end(); ++it)
        {
          os << precisionWrapper(it->getMZ()) << " "
             << precisionWrapper(it->getIntensity()) << "\n";
        }
      }
      else
      {
        for (PeakSpectrum::const_iterator it = spec.begin(); it != spec.end(); ++it)
        {
          PeakSpectrum::PeakType::IntensityType intensity = it->getIntensity();
          if (intensity == 0.0)
          {
            continue; // skip zero-intensity peaks
          }
          os << fixed << setprecision(HIGH_PRECISION) << it->getMZ() << " "
             << setprecision(LOW_PRECISION) << intensity << "\n";
        }
      }
      os << "END IONS\n";
    }
  }

  std::pair<String, String> MascotGenericFile::getHTTPPeakListEnclosure(const String& filename) const
  {
    std::pair<String, String> r;
    r.first = String("--" + (std::string)param_.getValue("internal:boundary") + "\n" + R"(Content-Disposition: form-data; name="FILE"; filename=")" + filename + "\"\n\n");
    r.second = String("\n\n--" + (std::string)param_.getValue("internal:boundary") + "--\n");
    return r;
  }

  void MascotGenericFile::writeMSExperiment_(ostream& os, const String& filename, const PeakMap& experiment)
  {

    std::pair<String, String> enc = getHTTPPeakListEnclosure(filename);
    if (param_.getValue("internal:HTTP_format").toBool())
    {
      os << enc.first;
    }

    QFileInfo fileinfo(filename.c_str());
    QString filtered_filename = fileinfo.completeBaseName();
    filtered_filename.remove(QRegularExpression("[^a-zA-Z0-9]"));


    String native_id_type_accession;
    const vector<SourceFile>& sourcefiles = experiment.getSourceFiles();
    if (sourcefiles.empty())
    {
      OPENMS_LOG_WARN << "MascotGenericFile: no native ID accession." << endl;
      native_id_type_accession = "UNKNOWN";
    }
    else
    {
      native_id_type_accession = experiment.getSourceFiles()[0].getNativeIDTypeAccession();
      if (native_id_type_accession.empty())
      {
        OPENMS_LOG_WARN << "MascotGenericFile: empty native ID accession." << endl;
        native_id_type_accession = "UNKNOWN";
      }
    }


    this->startProgress(0, experiment.size(), "storing mascot generic file");
    for (Size i = 0; i < experiment.size(); i++)
    {
      this->setProgress(i);
      if (experiment[i].getMSLevel() == 2)
      {
        writeSpectrum(os, experiment[i], filtered_filename, native_id_type_accession);
      }
      else if (experiment[i].getMSLevel() == 0)
      {
        OPENMS_LOG_WARN << "MascotGenericFile: MSLevel is set to 0, ignoring this spectrum!" << "\n";
      }
    }
    // close file
    if (param_.getValue("internal:HTTP_format").toBool())
    {
      os << enc.second;
    }
    this->endProgress();
  }

} // namespace OpenMS
