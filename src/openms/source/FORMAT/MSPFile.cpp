// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MSPFile.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <regex>
#include <map>

using namespace std;

namespace OpenMS
{
  MSPFile::MSPFile() :
    DefaultParamHandler("MSPFile")
  {
    defaults_.setValue("parse_headers", "false", "Flag whether header information should be parsed an stored for each spectrum");
    vector<std::string> parse_strings{"true","false"};
    defaults_.setValidStrings("parse_headers", parse_strings);
    defaults_.setValue("parse_peakinfo", "true", "Flag whether the peak annotation information should be parsed and stored for each peak");
    defaults_.setValidStrings("parse_peakinfo", parse_strings);
    defaults_.setValue("parse_firstpeakinfo_only", "true", "Flag whether only the first (default for 1:1 correspondence in SpecLibSearcher) or all peak annotation information should be parsed and stored for each peak.");
    defaults_.setValidStrings("parse_firstpeakinfo_only", parse_strings);
    defaults_.setValue("instrument", "", "If instrument given, only spectra of these type of instrument (Inst= in header) are parsed");
    defaults_.setValidStrings("instrument", {"","it","qtof","toftof"});

    defaultsToParam_();
  }

  MSPFile::MSPFile(const MSPFile & rhs) = default;

  MSPFile & MSPFile::operator=(const MSPFile & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
    }
    return *this;
  }

  MSPFile::~MSPFile() = default;

  void MSPFile::load(const String & filename, vector<PeptideIdentification> & ids, PeakMap & exp)
  {
    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    if (!File::readable(filename))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    // groups everything inside the shortest pair of parentheses
    const std::regex rex(R"(\((.*?)\))");
    // matches 2+ whitespaces or tabs or returns "   ", "\t", "\r"
    // Note: this is a hack because one of the encountered formats has single whitespaces in quotes.
    // TODO choose a format during construction of the class. If we actually knew how to call and define them.
    const std::regex ws_rex(R"(\s{2,}|\t|\r)");

    exp.reset();

    //set DocumentIdentifier
    exp.setLoadedFileType(filename);
    exp.setLoadedFilePath(filename);

    String line;
    ifstream is(filename.c_str());

    std::map<String, String> modname_to_unimod;
    modname_to_unimod["Pyro-glu"] = "Gln->pyro-Glu";
    modname_to_unimod["CAM"] = "Carbamidomethyl";
    modname_to_unimod["AB_old_ICATd8"] = "ICAT-D:2H(8)";
    modname_to_unimod["AB_old_ICATd0"] = "ICAT-D";

    bool parse_headers(param_.getValue("parse_headers").toBool());
    bool parse_peakinfo(param_.getValue("parse_peakinfo").toBool());
    bool parse_firstpeakinfo_only(param_.getValue("parse_firstpeakinfo_only").toBool());
    std::string instrument((std::string)param_.getValue("instrument"));
    bool inst_type_correct(true);
    [[maybe_unused]] bool spectrast_format(false); // TODO: implement usage
    Size spectrum_number = 0;

    PeakSpectrum spec;

    // line number counter
    Size line_number = 0;

    while (getline(is, line))
    {
      ++line_number;

      if (line.hasPrefix("Name:"))
      {
        vector<String> split, split2;
        line.split(' ', split);
        split[1].split('/', split2);
        String peptide = split2[0];
        // in newer NIST versions, the charge is followed by the modification(s) e.g. "_1(0,A,Acetyl)"
        vector<String> split3;
        split2[1].split('_', split3);
        Int charge = split3[0].toInt();

        // remove modifications inside the peptide string, since it is also defined in 'Mods=' comment
        peptide = std::regex_replace(peptide, rex, "");
        PeptideIdentification id;
        id.insertHit(PeptideHit(0, 0, charge, AASequence::fromString(peptide)));
        ids.push_back(id);
        inst_type_correct = true;
      }
      else if (line.hasPrefix("MW:"))
      {
        vector<String> split;
        line.split(' ', split);
        if (split.size() == 2)
        {
          UInt charge = ids.back().getHits().begin()->getCharge();
          spec.getPrecursors().resize(1);
          spec.getPrecursors()[0].setMZ((split[1].toDouble() + (double)charge * Constants::PROTON_MASS_U) / (double)charge);
        }
      }
      else if (line.hasPrefix("Comment:"))
      {
        // slow, but we need the modifications from the header and the instrument type
        vector<String> split;
        line.split(' ', split);
        for (vector<String>::const_iterator it = split.begin(); it != split.end(); ++it)
        {
          if (!inst_type_correct)
          {
            break;
          }
          if (!instrument.empty() && it->hasPrefix("Inst="))
          {
            String inst_type = it->suffix('=');
            if (instrument != inst_type)
            {
              inst_type_correct = false;
              ids.erase(--ids.end());
            }
            break;
          }

          if (it->hasPrefix("Mods=") && *it != "Mods=0")
          {
            String mods = it->suffix('=');
            // e.g. Mods=2/7,K,Carbamyl/22,K,Carbamyl
            vector<String> mod_split;
            mods.split('/', mod_split);
            if (mod_split.size() <= 1) // e.g. Mods=2(0,A,Acetyl)(11,M,Oxidation)
            {
              mod_split.clear();
              mod_split.emplace_back(mods.prefix('('));
              Size sz = mod_split[0].toInt();
              std::smatch sm;
              std::string::const_iterator cit = mods.cbegin();
              // go through all pairs of parentheses
              while (std::regex_search(cit, mods.cend(), sm, rex) && mod_split.size()-1 <= sz)
              {
                if (sm.size() == 2) // 2 = match
                {
                  mod_split.emplace_back(sm[1].str());
                }
                // set cit to after match
                cit = sm[0].second;
              }
            }
            AASequence peptide = ids.back().getHits().begin()->getSequence();
            for (Size i = 1; i <= (UInt)mod_split[0].toInt(); ++i)
            {
              vector<String> single_mod;
              mod_split[i].split(',', single_mod);

              String mod_name = single_mod[2];
              if (modname_to_unimod.find(mod_name) != modname_to_unimod.end())
              {
                mod_name = modname_to_unimod[mod_name];
              }
              String origin  = single_mod[1];
              Size position = single_mod[0].toInt();

              //cerr << "MSP modification: " << origin << " " << mod_name << " " << position << "\n";

              if (position > 0 && position < peptide.size() - 1)
              {
                peptide.setModification(position, mod_name);
              }
              else if (position == 0)
              {
                // we must decide whether this can be a terminal mod
                try
                {
                  peptide.setNTerminalModification(mod_name);
                }
                catch (Exception::ElementNotFound& /*e*/)
                {
                  peptide.setModification(position, mod_name);
                }
              }
              else if (position == peptide.size() - 1)
              {
                // we must decide whether this can be a terminal mod
                try
                {
                  peptide.setCTerminalModification(mod_name);
                }
                catch (Exception::ElementNotFound& /*e*/)
                {
                  peptide.setModification(position, mod_name);
                }
              }
              else
              {
                cerr << "MSPFile: Error: ignoring modification: '" << line << "' in line " << line_number << "\n";
              }
            }
            vector<PeptideHit> hits(ids.back().getHits());
            hits.begin()->setSequence(peptide);
            ids.back().setHits(hits);
          }
        }

        if (parse_headers && inst_type_correct)
        {
          parseHeader_(line, spec);
        }
      }
      else if (line.hasPrefix("Num peaks:") || line.hasPrefix("NumPeaks:"))
      {
        if (line.hasPrefix("NumPeaks:"))
        {
          spectrast_format = true;
        }

        if (!inst_type_correct)
        {
          while (getline(is, line) && ++line_number && !line.empty() && isdigit(line[0]))
          {
          }
        }
        else
        {
          PeptideHit& hitToAnnotate = ids.back().getHits()[0];
          vector<PeptideHit::PeakAnnotation> annots;
          while (getline(is, line) && ++line_number && !line.empty() && isdigit(line[0]))
          {
            std::sregex_token_iterator iter(line.begin(),
                                            line.end(),
                                            ws_rex,
                                            -1);
            std::sregex_token_iterator end;
            if (iter == end)
            {
              throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          line, R"(not <mz><tab/spaces><intensity><tab/spaces>"<annotation>"<tab/spaces>"<comment>" in line )" + String(line_number));
            }
            Peak1D peak;
            float mz = String(iter->str()).toFloat();
            peak.setMZ(mz);
            ++iter;
            if (iter == end)
            {
              throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          line, R"(not <mz><tab/spaces><intensity><tab/spaces>"<annotation>"<tab/spaces>"<comment>" in line )" + String(line_number));
            }
            float ity = String(iter->str()).toFloat();
            peak.setIntensity(ity);
            ++iter;
            if (parse_peakinfo && iter != end)
            {
              //e.g. "b32-H2O^3/0.11,y19-H2O^2/0.26"
              String annot = iter->str();
              annot = annot.unquote();
              if (annot.hasPrefix("?"))  //"? 2/2 0.6" or "?i 2/2 0.6" whatever i means, it will be lost here
              {
                annots.push_back(PeptideHit::PeakAnnotation{"?", 0, mz, ity});
              }
              else
              {
                if (annot.has(' ')) annot = annot.prefix(' '); // in case of different format "b8/-0.07,y9-46/-0.01 2/2 32.4" we only need the first part
                StringList splitstr;
                annot.split(',',splitstr);
                for (auto& str : splitstr)
                {
                  String splitstrprefix = str.prefix('/');
                  int charge = 1;
                  StringList splitstr2;
                  splitstrprefix.split('^', splitstr2);
                  if (splitstr2.size() > 1)
                  {
                    charge = splitstr2[1].toInt();
                  }
                  annots.push_back(PeptideHit::PeakAnnotation{splitstr2[0], charge, mz, ity});
                  if (parse_firstpeakinfo_only) break;
                }
              }
            }
            else if (parse_peakinfo)
            {
              throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          line, "Requested reading peak info but no annotation found for line " + String(line_number));
            }
            spec.push_back(peak);
          }
          hitToAnnotate.setPeakAnnotations(annots);
          spec.setNativeID(String("index=") + spectrum_number);
          exp.addSpectrum(spec);
          // clear spectrum
          spec.clear(true);
        }
        spectrum_number++;
      }
    }
  }

  void MSPFile::parseHeader_(const String & header, PeakSpectrum & spec)
  {
    // first header from std_protein of NIST spectra DB
    // Spec=Consensus Pep=Tryptic Fullname=R.AAANFFSASCVPCADQSSFPK.L/2 Mods=0 Parent=1074.480 Inst=it Mz_diff=0.500 Mz_exact=1074.4805 Mz_av=1075.204 Protein="TRFE_BOVIN" Organism="Protein Standard" Se=2^X23:ex=3.1e-008/1.934e-005,td=5.14e+007/2.552e+019,sd=0/0,hs=45.8/5.661,bs=6.3e-021,b2=1.2e-015,bd=5.87e+020^O22:ex=3.24e-005/0.0001075,td=304500/5.909e+297,pr=3.87e-007/1.42e-006,bs=1.65e-301,b2=1.25e-008,bd=1.3e+299 Sample=1/bovine-serotransferrin_cam,23,26 Nreps=23/34 Missing=0.3308/0.0425 Parent_med=1074.88/0.23 Max2med_orig=22.1/9.5 Dotfull=0.618/0.029 Dot_cons=0.728/0.040 Unassign_all=0.161 Unassigned=0.000 Dotbest=0.70 Naa=21 DUScorr=2.3/3.8/0.61 Dottheory=0.86 Pfin=4.3e+010 Probcorr=1 Tfratio=8e+008 Specqual=0.0
    vector<String> split;
    header.split(' ', split);

    for (vector<String>::const_iterator it = split.begin(); it != split.end(); ++it)
    {
      vector<String> split2;
      String tmp = *it;
      tmp.trim();
      tmp.split('=', split2);
      if (split2.size() == 2)
      {
        spec.setMetaValue(split2[0], split2[1]);
      }
    }
  }

  //TODO adapt store to write new? format
  void MSPFile::store(const String & filename, const PeakMap & exp) const
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::MSP))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::MSP) + "'");
    }

    if (!File::writable(filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    ofstream out(filename.c_str());

    for (const MSSpectrum& it : exp)
    {
      if (!it.getPeptideIdentifications().empty() && !it.getPeptideIdentifications().begin()->getHits().empty())
      {
        PeptideHit hit = *it.getPeptideIdentifications().begin()->getHits().begin();
        String peptide;
        for (const Residue& pit : hit.getSequence())
        {
          if (pit.isModified() && pit.getOneLetterCode() == "M" &&
              fabs(pit.getModification()->getDiffFormula().getMonoWeight() - 16.0) < 0.01)
          {
            peptide += "M(O)"; // TODO why are we writing specifically only oxidations?
          }
          else
          {
            peptide += pit.getOneLetterCode();
          }
        }
        out << "Name: " << peptide << "/" << hit.getCharge() << "\n";
        out << "MW: " << hit.getSequence().getMonoWeight() << "\n";
        out << "Comment:";

        // modifications
        // e.g. 2/9,C,Carbamidomethyl/12,C,Carbamidomethyl
        Size num_mods(0);
        vector<String> modifications;
        if (hit.getSequence().hasNTerminalModification())
        {
          String mod = hit.getSequence().getNTerminalModificationName();
          ++num_mods;
          String modification = "0," + hit.getSequence().begin()->getOneLetterCode() + "," + mod;
          modifications.push_back(modification);
        }

        // @improvement improve writing support (Andreas)
        UInt pos(0);
        for (AASequence::ConstIterator pit = hit.getSequence().begin(); pit != hit.getSequence().end(); ++pit, ++pos)
        {
          if (!pit->isModified())
          {
            continue;
          }

          String mod = pit->getModificationName();
          String res = pit->getOneLetterCode();
          ++num_mods;
          String modification = String(pos) + "," + res + "," + mod;
          modifications.push_back(modification);
        }

        String mods;
        mods.concatenate(modifications.begin(), modifications.end(), "/");
        if (!mods.empty())
        {
          out << " Mods=" << String(num_mods)  << "/" << mods;
        }
        else
        {
          out << " Mods=0";
        }
        out << " Inst=it\n";         // @improvement write instrument type, protein...and other information
        out << "Num peaks: " << it.size() << "\n";

        // normalize to 10,000
        PeakSpectrum rich_spec = it;
        double max_int(0);
        for (const Peak1D& sit : rich_spec)
        {
          if (sit.getIntensity() > max_int)
          {
            max_int = sit.getIntensity();
          }
        }

        if (max_int != 0)
        {
          for (Peak1D& sit : rich_spec)
          {
            sit.setIntensity(sit.getIntensity() / max_int * 10000.0);
          }
        }
        else
        {
          cerr << "MSPFile: spectrum contains only zero intensities!" << endl;
        }

        int ion_name = -1;
        for (Size k = 0; k < rich_spec.getStringDataArrays().size(); k++)
        {
          if (rich_spec.getStringDataArrays()[k].getName() == Constants::UserParam::IonNames)
          {
            ion_name = (int)k;
            break;
          }
        }

        Size k = 0;
        for (const Peak1D& sit : rich_spec)
        {
          out << sit.getPosition()[0] << "\t" << sit.getIntensity() << "\t";
          if (ion_name >= 0)
          {
            out << "\"" << rich_spec.getStringDataArrays()[ion_name][k] << "\"";
            k++;
          }
          else
          {
            out << "\"?\"";
          }
          out << "\n";
        }
        out << "\n";
      }
    }
  }

} // namespace OpenMS
