// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <QFileInfo>
#include <QtCore/QRegExp>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#define HIGH_PRECISION 8
#define LOW_PRECISION 6

using namespace std;

namespace OpenMS
{

  MascotGenericFile::MascotGenericFile() :
    ProgressLogger(), DefaultParamHandler("MascotGenericFile"), mod_group_map_()
  {
    defaults_.setValue("database", "MSDB", "Name of the sequence database");
    defaults_.setValue("search_type", "MIS", "Name of the search type for the query", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("search_type", ListUtils::create<String>("MIS,SQ,PMF"));
    defaults_.setValue("enzyme", "Trypsin", "The enzyme descriptor to the enzyme used for digestion. (Trypsin is default, None would be best for peptide input or unspecific digestion, for more please refer to your mascot server).");
    defaults_.setValue("instrument", "Default", "Instrument definition which specifies the fragmentation rules");
    defaults_.setValue("missed_cleavages", 1, "Number of missed cleavages allowed for the enzyme");
    defaults_.setMinInt("missed_cleavages", 0);
    defaults_.setValue("precursor_mass_tolerance", 3.0, "Tolerance of the precursor peaks");
    defaults_.setMinFloat("precursor_mass_tolerance", 0.0);
    defaults_.setValue("precursor_error_units", "Da", "Units of the precursor mass tolerance");
    defaults_.setValidStrings("precursor_error_units", ListUtils::create<String>("%,ppm,mmu,Da"));
    defaults_.setValue("fragment_mass_tolerance", 0.3, "Tolerance of the peaks in the fragment spectrum");
    defaults_.setMinFloat("fragment_mass_tolerance", 0.0);
    defaults_.setValue("fragment_error_units", "Da", "Units of the fragment peaks tolerance");
    defaults_.setValidStrings("fragment_error_units", ListUtils::create<String>("mmu,Da"));
    defaults_.setValue("charges", "1,2,3", "Allowed charge states, given as a comma separated list of integers");
    defaults_.setValue("taxonomy", "All entries", "Taxonomy specification of the sequences");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    defaults_.setValue("fixed_modifications", ListUtils::create<String>(""), "List of fixed modifications, according to UniMod definitions.");
    defaults_.setValidStrings("fixed_modifications", all_mods);
    defaults_.setValue("variable_modifications", ListUtils::create<String>(""), "Variable modifications given as UniMod definitions.");
    defaults_.setValidStrings("variable_modifications", all_mods);

    // special modifications, see "updateMembers_" method below:
    defaults_.setValue("special_modifications", "Cation:Na (DE),Deamidated (NQ),Oxidation (HW),Phospho (ST),Sulfo (ST)", "Modifications with specificity groups that are used by Mascot and have to be treated specially", ListUtils::create<String>("advanced"));
    // list from Mascot 2.4; there's also "Phospho (STY)", but that can be
    // represented using "Phospho (ST)" and "Phospho (Y)"

    defaults_.setValue("mass_type", "monoisotopic", "Defines the mass type, either monoisotopic or average");
    defaults_.setValidStrings("mass_type", ListUtils::create<String>("monoisotopic,average"));
    defaults_.setValue("number_of_hits", 0, "Number of hits which should be returned, if 0 AUTO mode is enabled.");
    defaults_.setMinInt("number_of_hits", 0);
    defaults_.setValue("skip_spectrum_charges", "false", "Sometimes precursor charges are given for each spectrum but are wrong, setting this to 'true' does not write any charge information to the spectrum, the general charge information is however kept.");
    defaults_.setValidStrings("skip_spectrum_charges", ListUtils::create<String>("true,false"));

    defaults_.setValue("search_title", "OpenMS_search", "Sets the title of the search.", ListUtils::create<String>("advanced"));
    defaults_.setValue("username", "OpenMS", "Sets the username which is mentioned in the results file.", ListUtils::create<String>("advanced"));
    defaults_.setValue("email", "", "Sets the email which is mentioned in the results file. Note: Some server require that a proper email is provided.");

    // the next section should not be shown to TOPP users
    Param p;
    p.setValue("format", "Mascot generic", "Sets the format type of the peak list, this should not be changed unless you write the header only.", ListUtils::create<String>("advanced"));
    p.setValidStrings("format", ListUtils::create<String>("Mascot generic,mzData (.XML),mzML (.mzML)")); // Mascot's HTTP interface supports more, but we don't :)
    p.setValue("boundary", "GZWgAaYKjHFeUaLOLEIOMq", "MIME boundary for parameter header (if using HTTP format)", ListUtils::create<String>("advanced"));
    p.setValue("HTTP_format", "false", "Write header with MIME boundaries instead of simple key-value pairs. For HTTP submission only.", ListUtils::create<String>("advanced"));
    p.setValidStrings("HTTP_format", ListUtils::create<String>("true,false"));
    p.setValue("content", "all", "Use parameter header + the peak lists with BEGIN IONS... or only one of them.", ListUtils::create<String>("advanced"));
    p.setValidStrings("content", ListUtils::create<String>("all,peaklist_only,header_only"));
    defaults_.insert("internal:", p);

    defaultsToParam_();
  }

  MascotGenericFile::~MascotGenericFile()
  {
  }

  void MascotGenericFile::updateMembers_()
  {
    // special cases for specificity groups: OpenMS uses e.g. "Deamidated (N)"
    // and "Deamidated (Q)", but Mascot only understands "Deamidated (NQ)"
    String special_mods = param_.getValue("special_modifications");
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
    if (!File::writable(filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    ofstream os(filename.c_str());
    store(os, filename, experiment, compact);
    os.close();
  }

  void MascotGenericFile::store(ostream& os, const String& filename, const PeakMap& experiment, bool compact)
  {
    const streamsize precision = os.precision(); // may get changed, so back-up
    
    store_compact_ = compact;
    if (param_.getValue("internal:content") != "peaklist_only")
      writeHeader_(os);
    if (param_.getValue("internal:content") != "header_only")
      writeMSExperiment_(os, filename, experiment);

    os.precision(precision); // reset precision
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
    writeParameterHeader_("FORMAT", os);    // make sure this stays within the first 5 lines of the file, since we use it to recognize our own MGF files in case their file suffix is not MGF
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

    //db name
    writeParameterHeader_("DB", os);
    os << param_.getValue("database") << "\n";

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

    //cleavage enzyme
    writeParameterHeader_("CLE", os);
    os << param_.getValue("enzyme") << "\n";

    //average/monoisotopic
    writeParameterHeader_("MASS", os);
    os << param_.getValue("mass_type") << "\n";

    //fixed modifications
    vector<String> fixed_mods = param_.getValue("fixed_modifications");
    writeModifications_(fixed_mods, os);

    //variable modifications
    vector<String> var_mods = param_.getValue("variable_modifications");
    writeModifications_(var_mods, os, true);

    //instrument
    writeParameterHeader_("INSTRUMENT", os);
    os << param_.getValue("instrument") << "\n";

    //missed cleavages
    writeParameterHeader_("PFA", os);
    os << param_.getValue("missed_cleavages") << "\n";

    //precursor mass tolerance
    writeParameterHeader_("TOL", os);
    os << param_.getValue("precursor_mass_tolerance") << "\n";

    //ion mass tolerance_
    writeParameterHeader_("ITOL", os);
    os << param_.getValue("fragment_mass_tolerance") << "\n";

    //taxonomy
    writeParameterHeader_("TAXONOMY", os);
    os << param_.getValue("taxonomy") << "\n";

    //charge
    writeParameterHeader_("CHARGE", os);
    os << param_.getValue("charges") << "\n";
  }

  void MascotGenericFile::writeSpectrum_(ostream& os, const PeakSpectrum& spec, const String& filename)
  {
    Precursor precursor;
    if (spec.getPrecursors().size() > 0)
    {
      precursor = spec.getPrecursors()[0];
    }
    if (spec.getPrecursors().size() > 1)
    {
      cerr << "Warning: The spectrum written to Mascot file has more than one precursor. The first precursor is used!\n";
    }
    if (spec.size() >= 10000)
    {
      throw Exception::InvalidValue(
        __FILE__, __LINE__, __PRETTY_FUNCTION__, "Spectrum to be written as "
        "MGF has more than 10.000 peaks, which is the maximum upper limit. "
        "Only centroided data is allowed. This is most likely raw data.", 
        String(spec.size()));
    }
    double mz(precursor.getMZ()), rt(spec.getRT());

    if (mz == 0)
    {
      //retention time
      cout << "No precursor m/z information for spectrum with rt " << rt 
           << " present, skipping spectrum!\n";
    }
    else
    {
      os << "\n";
      os << "BEGIN IONS\n";
      if (!store_compact_)
      {
        os << "TITLE=" << precisionWrapper(mz) << "_" << precisionWrapper(rt) 
           << "_" << spec.getNativeID() << "_" << filename << "\n";
        os << "PEPMASS=" << precisionWrapper(mz) <<  "\n";
        os << "RTINSECONDS=" << precisionWrapper(rt) << "\n";
      }
      else
      {
        os << "TITLE=" << setprecision(HIGH_PRECISION) << mz << "_" 
           << setprecision(LOW_PRECISION) << rt << "_" 
           << spec.getNativeID() << "_" << filename << "\n";
        os << "PEPMASS=" << setprecision(HIGH_PRECISION) << mz << "\n";
        os << "RTINSECONDS=" << setprecision(LOW_PRECISION) << rt << "\n";
      }

      int charge(precursor.getCharge());

      if (charge != 0)
      {
        bool skip_spectrum_charges(param_.getValue("skip_spectrum_charges").toBool());
        if (!skip_spectrum_charges)
        {
          os << "CHARGE=" << charge << "\n";
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
          if (intensity == 0.0) continue; // skip zero-intensity peaks
          os << setprecision(HIGH_PRECISION) << it->getMZ() << " " 
             << setprecision(LOW_PRECISION) << intensity << "\n";
        }
      }
      os << "END IONS\n";
    }
  }

  std::pair<String, String> MascotGenericFile::getHTTPPeakListEnclosure(const String& filename) const
  {
    std::pair<String, String> r;
    r.first = String("--" + String(param_.getValue("internal:boundary")) + "\n" + "Content-Disposition: form-data; name=\"FILE\"; filename=\"" + filename + "\"\n\n");
    r.second = String("\n\n--" + String(param_.getValue("internal:boundary")) + "--\n");
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
    filtered_filename.remove(QRegExp("[^a-zA-Z0-9]"));
    this->startProgress(0, experiment.size(), "storing mascot generic file");
    for (Size i = 0; i < experiment.size(); i++)
    {
      this->setProgress(i);
      if (experiment[i].getMSLevel() == 2)
      {
        writeSpectrum_(os, experiment[i], filtered_filename);
      }
      else if (experiment[i].getMSLevel() == 0)
      {
        LOG_WARN << "MascotGenericFile: MSLevel is set to 0, ignoring this spectrum!" << "\n";
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
