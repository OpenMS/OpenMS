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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  template<class T>   // primary template
  bool extractName(T& value, const std::string& header_name,
                   const std::vector<std::string>& tmp_line,
                   const std::map<std::string, int>& header_dict)
  {
    auto tmp = header_dict.find( header_name );
    if (tmp != header_dict.end())
    {
      value = tmp_line[ tmp->second ];
      // perform cleanup
      value = value.remove('"');
      value = value.remove('\'');
      value = value.remove(',');
      return true;
    }
    return false;
  }

  template<>   // specialization for int
  bool extractName<int>(int& value, const std::string& header_name,
                        const std::vector<std::string>& tmp_line,
                        const std::map<std::string, int>& header_dict)
  {
    auto tmp = header_dict.find( header_name );
    if (tmp != header_dict.end())
    {
      value = String(tmp_line[ tmp->second ]).toInt();
      return true;
    }
    return false;
  }

  template<>   // specialization for double
  bool extractName<double>(double& value, const std::string& header_name,
                        const std::vector<std::string>& tmp_line,
                        const std::map<std::string, int>& header_dict)
  {
    auto tmp = header_dict.find( header_name );
    if (tmp != header_dict.end())
    {
      value = String(tmp_line[ tmp->second ]).toDouble();
      return true;
    }
    return false;
  }

  template<>   // specialization for bool
  bool extractName<bool>(bool& value, const std::string& header_name,
                        const std::vector<std::string>& tmp_line,
                        const std::map<std::string, int>& header_dict)
  {
    auto tmp = header_dict.find( header_name );
    if (tmp != header_dict.end())
    {
      auto str_value = tmp_line[ tmp->second ];
      if (str_value == "1" || str_value == "TRUE") value = true;
      else if (str_value == "0" || str_value == "FALSE") value = false;
      else return false;

      // all went well, we set the value and can return
      return true;
    }
    return false;
  }

  TransitionTSVFile::TransitionTSVFile() :
    DefaultParamHandler("TransitionTSVFile")
  {
    defaults_.setValue("retentionTimeInterpretation", "iRT", "How to interpret the provided retention time (the retention time column can either be interpreted to be in iRT, minutes or seconds)", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("retentionTimeInterpretation", ListUtils::create<String>("iRT,seconds,minutes"));
    defaults_.setValue("override_group_label_check", "false", "Override an internal check that assures that all members of the same PeptideGroupLabel have the same PeptideSequence (this ensures that only different isotopic forms of the same peptide can be grouped together in the same label group). Only turn this off if you know what you are doing.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("override_group_label_check", ListUtils::create<String>("true,false"));
    defaults_.setValue("force_invalid_mods", "false", "Force reading even if invalid modifications are encountered (OpenMS may not recognize the modification)", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("force_invalid_mods", ListUtils::create<String>("true,false"));

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  TransitionTSVFile::~TransitionTSVFile()
  {
  }

  void TransitionTSVFile::updateMembers_()
  {
    retentionTimeInterpretation_ = param_.getValue("retentionTimeInterpretation");
    override_group_label_check_ = param_.getValue("override_group_label_check").toBool();
    force_invalid_mods_ = param_.getValue("force_invalid_mods").toBool();
  }

  const char* TransitionTSVFile::strarray_[] =
  {
    "PrecursorMz",
    "ProductMz",
    "PrecursorCharge",
    "ProductCharge",
    "LibraryIntensity",
    "NormalizedRetentionTime",
    "PeptideSequence",
    "ModifiedPeptideSequence",
    "PeptideGroupLabel",
    "LabelType",
    "CompoundName",
    "SumFormula",
    "SMILES",
    "ProteinId",
    "UniprotId",
    "FragmentType",
    "FragmentSeriesNumber",
    "Annotation",
    "CollisionEnergy",
    "TransitionGroupId",
    "TransitionId",
    "Decoy",
    "DetectingTransition",
    "IdentifyingTransition",
    "QuantifyingTransition"
  };

  const std::vector<std::string> TransitionTSVFile::header_names_(strarray_, strarray_ + 25);

  void TransitionTSVFile::getTSVHeader_(const std::string& line, char& delimiter,
                                          std::vector<std::string> header, std::map<std::string, int>& header_dict)
  {
    std::string tmp;

    int nr_delimiters = 3;
    Size min_header_size = 8;
    const char possibleDelimiters[3] = {',', ';', '\t'};

    for (int i = 0; i < nr_delimiters; i++)
    {
      std::stringstream lineStream(line);
      delimiter = possibleDelimiters[i];
      while (std::getline(lineStream, tmp, delimiter))
      {
        header.push_back(tmp);
      }
      if (header.size() >= min_header_size)
      {
        break; // found the delimiter, got the correct header
      }
      header.clear();
    }

    for (Size i = 0; i < header.size(); i++)
    {
      header_dict[header[i]] = i;
    }
    char txt_delimiter = delimiter;
    if (txt_delimiter == '\t')
    {
      txt_delimiter = 't';
    }

    // could not determine the delimiter correctly
    if (header.size() < min_header_size)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Determined your csv/tsv file to have delimiter '" + (String)txt_delimiter + 
          "', but the parsed header has only " + (String)header.size() + " fields instead of the minimal " + 
          (String)min_header_size + ". Please check your input file.");
    }
  }

  void TransitionTSVFile::readUnstructuredTSVInput_(const char* filename, FileTypes::Type filetype, std::vector<TSVTransition>& transition_list)
  {
    std::ifstream data(filename);
    std::string   line;
    std::string   tmp;

    // read header
    std::vector<std::string>   tmp_line;
    std::vector<std::string>   header;
    std::map<std::string, int> header_dict;
    char delimiter = ',';

    // SpectraST MRM Files do not have a header
    if (FileTypes::typeToName(filetype) == "mrm")
    {
      delimiter = '\t';

      header_dict["SpectraSTBestSample"] = 0;
      header_dict["SpectraSTmaxNumUsed/totNumUsed"] = 1;
      header_dict["SpectraSTpI"] = 2;
      header_dict["PrecursorMz"] = 3;
      header_dict["SpectraSTRetentionTime"] = 4;
      header_dict["ProductMz"] = 5;
      header_dict["LibraryIntensity"] = 6;
      header_dict["SpectraSTAnnotation"] = 7;
      header_dict["FragmentCharge"] = 8;
      header_dict["SpectraSTFullPeptideName"] = 9;
      header_dict["SpectraSTUnknown"] = 10;
      header_dict["SpectraSTNumberOfProteinsMappedTo"] = 11;
      header_dict["ProteinName"] = 12;
    }
    // Read header for TSV input
    else
    {
      std::getline(data, line);

      getTSVHeader_(line, delimiter, header, header_dict);
    }

    bool spectrast_legacy = false; // we will check below if SpectraST was run in legacy (<5.0) mode or if the RT normalization was forgotten.
    int cnt = 0;
    while (std::getline(data, line))
    {
      line.push_back(delimiter); // avoid losing last column if it is empty
      std::stringstream lineStream(line);

      while (std::getline(lineStream, tmp, delimiter))
      {
        tmp_line.push_back(tmp);
      }
      cnt++;

#ifdef TRANSITIONTSVREADER_TESTING
      for (Size i = 0; i < tmp_line.size(); i++)
      {
        std::cout << "line " << i << " " << tmp_line[i] << std::endl;
      }

      for (std::map<std::string, int>::iterator iter = header_dict.begin(); iter != header_dict.end(); ++iter)
      {
        std::cout << "header " << iter->first << " " << iter->second << std::endl;
      }
#endif

      if (tmp_line.size() != header_dict.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Error reading the file on line " + String(cnt) + ": length of the header and length of the line" +
                                         " do not match: " + String(tmp_line.size()) + " != " + String(header_dict.size()));
      }

      TSVTransition mytransition;
      bool skip_transition = false; // skip unannotated transitions in SpectraST MRM files

      //// Required columns (they are guaranteed to be present, see getTSVHeader_)
      // PrecursorMz
      mytransition.precursor = String(tmp_line[header_dict["PrecursorMz"]]).toDouble();

      // ProductMz
      if (!extractName<double>(mytransition.product, "ProductMz", tmp_line, header_dict) &&
          !extractName<double>(mytransition.product, "FragmentMz", tmp_line, header_dict)) // Spectronaut
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Expected a header named ProductMz or FragmentMz but found none");
      }

      // LibraryIntensity
      if (!extractName<double>(mytransition.library_intensity, "LibraryIntensity", tmp_line, header_dict) &&
          !extractName<double>(mytransition.library_intensity, "RelativeFragmentIntensity", tmp_line, header_dict)) // Spectronaut
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Expected a header named LibraryIntensity or RelativeFragmentIntensity but found none");
      }

      //// Additional columns for both proteomics and metabolomics
      // NormalizedRetentionTime
      if (!extractName<double>(mytransition.rt_calibrated, "RetentionTimeCalculatorScore", tmp_line, header_dict) && // Skyline
          !extractName<double>(mytransition.rt_calibrated, "iRT", tmp_line, header_dict) && // Spectronaut
          !extractName<double>(mytransition.rt_calibrated, "NormalizedRetentionTime", tmp_line, header_dict) &&
          !extractName<double>(mytransition.rt_calibrated, "RetentionTime", tmp_line, header_dict) &&
          !extractName<double>(mytransition.rt_calibrated, "Tr_recalibrated", tmp_line, header_dict))
      {
        if (header_dict.find("SpectraSTRetentionTime") != header_dict.end())
        {
          spectrastRTExtract(tmp_line[header_dict["SpectraSTRetentionTime"]], mytransition.rt_calibrated, spectrast_legacy);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "Expected a header named RetentionTime, NormalizedRetentionTime, iRT, RetentionTimeCalculatorScore, Tr_recalibrated or SpectraSTRetentionTime but found none");
        }
      }

      // PrecursorCharge
      !extractName(mytransition.precursor_charge, "PrecursorCharge", tmp_line, header_dict) &&
      !extractName(mytransition.precursor_charge, "Charge", tmp_line, header_dict); // charge is assumed to be the charge of the precursor

      !extractName(mytransition.fragment_type, "FragmentType", tmp_line, header_dict) &&
      !extractName(mytransition.fragment_type, "FragmentIonType", tmp_line, header_dict); // Skyline

      !extractName(mytransition.fragment_charge, "FragmentCharge", tmp_line, header_dict) &&
      !extractName(mytransition.fragment_charge, "ProductCharge", tmp_line, header_dict);

      !extractName<int>(mytransition.fragment_nr, "FragmentSeriesNumber", tmp_line, header_dict) &&
      !extractName<int>(mytransition.fragment_nr, "FragmentNumber", tmp_line, header_dict) &&
      !extractName<int>(mytransition.fragment_nr, "FragmentIonOrdinal", tmp_line, header_dict);

      extractName<double>(mytransition.fragment_mzdelta, "FragmentMzDelta", tmp_line, header_dict);
      extractName<int>(mytransition.fragment_modification, "FragmentModification", tmp_line, header_dict);

      //// Proteomics
      !extractName(mytransition.ProteinName, "ProteinName", tmp_line, header_dict) &&
      !extractName(mytransition.ProteinName, "ProteinId", tmp_line, header_dict); // Spectronaut

      extractName(mytransition.peptide_group_label, "PeptideGroupLabel", tmp_line, header_dict);

      extractName(mytransition.label_type, "LabelType", tmp_line, header_dict);

      !extractName(mytransition.PeptideSequence, "PeptideSequence", tmp_line, header_dict) &&
      !extractName(mytransition.PeptideSequence, "Sequence", tmp_line, header_dict) && // Skyline
      !extractName(mytransition.PeptideSequence, "StrippedSequence", tmp_line, header_dict); // Spectronaut

      !extractName(mytransition.FullPeptideName, "FullUniModPeptideName", tmp_line, header_dict) &&
      !extractName(mytransition.FullPeptideName, "FullPeptideName", tmp_line, header_dict) &&
      !extractName(mytransition.FullPeptideName, "ModifiedSequence", tmp_line, header_dict) && // Spectronaut
      !extractName(mytransition.FullPeptideName, "ModifiedPeptideSequence", tmp_line, header_dict);

      //// IPF
      !extractName<bool>(mytransition.detecting_transition, "detecting_transition", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.detecting_transition, "DetectingTransition", tmp_line, header_dict);
      !extractName<bool>(mytransition.identifying_transition, "identifying_transition", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.identifying_transition, "IdentifyingTransition", tmp_line, header_dict);
      !extractName<bool>(mytransition.quantifying_transition, "quantifying_transition", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.quantifying_transition, "QuantifyingTransition", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.quantifying_transition, "Quantitative", tmp_line, header_dict); // Skyline

      //// Targeted Metabolomics
      !extractName(mytransition.CompoundName, "CompoundName", tmp_line, header_dict) &&
      !extractName(mytransition.CompoundName, "CompoundId", tmp_line, header_dict);
      extractName(mytransition.SumFormula, "SumFormula", tmp_line, header_dict);
      extractName(mytransition.SMILES, "SMILES", tmp_line, header_dict);

      //// Meta
      extractName(mytransition.Annotation, "Annotation", tmp_line, header_dict);
      // UniprotId
      !extractName(mytransition.uniprot_id, "UniprotId", tmp_line, header_dict) &&
      !extractName(mytransition.uniprot_id, "UniprotID", tmp_line, header_dict);
      if (mytransition.uniprot_id == "NA") mytransition.uniprot_id = "";

      !extractName<double>(mytransition.CE, "CE", tmp_line, header_dict) &&
      !extractName<double>(mytransition.CE, "CollisionEnergy", tmp_line, header_dict);

      // Decoy
      !extractName<bool>(mytransition.decoy, "decoy", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.decoy, "Decoy", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.decoy, "IsDecoy", tmp_line, header_dict);

      if (header_dict.find("SpectraSTAnnotation") != header_dict.end())
      {
        skip_transition = spectrastAnnotationExtract(tmp_line[header_dict["SpectraSTAnnotation"]], mytransition);
      }

      //// Generate Group IDs
      // SpectraST
      if (FileTypes::typeToName(filetype) == "mrm")
      {
        std::vector<String> substrings;
        String(tmp_line[header_dict["SpectraSTFullPeptideName"]]).split("/", substrings);
        AASequence peptide = AASequence::fromString(substrings[0]);

        mytransition.FullPeptideName = peptide.toString();
        mytransition.PeptideSequence = peptide.toUnmodifiedString();
        mytransition.precursor_charge = substrings[1];

        mytransition.transition_name = String(cnt);

        mytransition.group_id = mytransition.FullPeptideName + String("_") + String(mytransition.precursor_charge);
      }
      // Generate transition_group_id and transition_name if not defined
      else
      {
        // Use TransitionId if available, else generate from attributes
        if (!extractName(mytransition.transition_name, "transition_name", tmp_line, header_dict) &&
            !extractName(mytransition.transition_name, "TransitionName", tmp_line, header_dict) &&
            !extractName(mytransition.transition_name, "TransitionId", tmp_line, header_dict))
        {
          mytransition.transition_name = String(cnt);
        }

        // Use TransitionGroupId if available, else generate from attributes
        if (!extractName(mytransition.group_id, "transition_group_id", tmp_line, header_dict) &&
            !extractName(mytransition.group_id, "TransitionGroupId", tmp_line, header_dict) &&
            !extractName(mytransition.group_id, "TransitionGroupName", tmp_line, header_dict))
        {
          mytransition.group_id = AASequence::fromString(mytransition.FullPeptideName).toString() + String("_") + String(mytransition.precursor_charge);
        }
      }

      cleanupTransitions_(mytransition);

      if (!skip_transition)
      {
        transition_list.push_back(mytransition);
      }

#ifdef TRANSITIONTSVREADER_TESTING
      std::cout << mytransition.precursor << std::endl;
      std::cout << mytransition.product << std::endl;
      std::cout << mytransition.rt_calibrated << std::endl;
      std::cout << mytransition.transition_name << std::endl;
      std::cout << mytransition.CE << std::endl;
      std::cout << mytransition.library_intensity << std::endl;
      std::cout << mytransition.group_id << std::endl;
      std::cout << mytransition.decoy << std::endl;
      std::cout << mytransition.PeptideSequence << std::endl;
      std::cout << mytransition.ProteinName << std::endl;
      std::cout << mytransition.Annotation << std::endl;
      std::cout << mytransition.FullPeptideName << std::endl;
      std::cout << mytransition.precursor_charge << std::endl;
      std::cout << mytransition.peptide_group_label << std::endl;
      std::cout << mytransition.fragment_charge << std::endl;
      std::cout << mytransition.fragment_nr << std::endl;
      std::cout << mytransition.fragment_mzdelta << std::endl;
      std::cout << mytransition.fragment_modification << std::endl;
      std::cout << mytransition.fragment_type << std::endl;
      std::cout << mytransition.uniprot_id << std::endl;
#endif

      tmp_line.clear();
    }

    if (spectrast_legacy && retentionTimeInterpretation_ == "iRT")
    {
      std::cout << "Warning: SpectraST was not run in RT normalization mode but the converted list was interpreted to have iRT units. Check whether you need to adapt the parameter -algorithm:retentionTimeInterpretation. You can ignore this warning if you used a legacy SpectraST 4.0 file." << std::endl;

    }
  }

  void TransitionTSVFile::spectrastRTExtract(const String str_inp, double & value, bool & spectrast_legacy)
  {
    // If SpectraST was run in RT normalization mode, the retention time is annotated as following: "3887.50(57.30)"
    // 3887.50 refers to the non-normalized RT of the individual or consensus run, and 57.30 refers to the normalized
    // iRT.
    size_t start_position = str_inp.find("(");
    if (start_position != std::string::npos)
    {
      ++start_position;
      size_t end_position = str_inp.find(")");
      if (end_position != std::string::npos)
      {
        value = String(str_inp.substr(start_position, end_position - start_position)).toDouble();
      }
    }
    else
    {
      // SpectraST was run without RT Normalization mode
      spectrast_legacy = true;
      value = str_inp.toDouble();
    }
  }

  bool TransitionTSVFile::spectrastAnnotationExtract(const String str_inp, TSVTransition & mytransition)
  {
    // Parses SpectraST fragment ion annotations
    // Example: y13^2/0.000,b16-18^2/-0.013,y7-45/0.000
    // Important: m2:8 are not yet supported! See SpectraSTPeakList::annotateInternalFragments for further information
    mytransition.Annotation = str_inp;

    std::vector<String> all_fragment_annotations;
    str_inp.split(",", all_fragment_annotations);

    if (all_fragment_annotations[0].find("[") == std::string::npos && // non-unique peak annotation
        all_fragment_annotations[0].find("]") == std::string::npos && // non-unique peak annotation
        all_fragment_annotations[0].find("I") == std::string::npos && // immonium ion
        all_fragment_annotations[0].find("p") == std::string::npos && // precursor ion
        all_fragment_annotations[0].find("i") == std::string::npos && // isotope ion
        all_fragment_annotations[0].find("m") == std::string::npos &&
        all_fragment_annotations[0].find("?") == std::string::npos
        )
    {
      std::vector<String> best_fragment_annotation_with_deviation;
      all_fragment_annotations[0].split("/", best_fragment_annotation_with_deviation);
      String best_fragment_annotation = best_fragment_annotation_with_deviation[0];

      if (best_fragment_annotation.find("^") != std::string::npos)
      {
        std::vector<String> best_fragment_annotation_charge;
        best_fragment_annotation.split("^", best_fragment_annotation_charge);
        mytransition.fragment_charge = String(best_fragment_annotation_charge[1]);
        best_fragment_annotation = best_fragment_annotation_charge[0];
      }
      else
      {
        mytransition.fragment_charge = 1; // assume 1 (most frequent charge state)
      }

      if (best_fragment_annotation.find("-") != std::string::npos)
      {
        std::vector<String> best_fragment_annotation_modification;
        best_fragment_annotation.split("-", best_fragment_annotation_modification);
        mytransition.fragment_type = best_fragment_annotation_modification[0].substr(0, 1);
        mytransition.fragment_nr = String(best_fragment_annotation_modification[0].substr(1)).toInt();
        mytransition.fragment_modification = -1 * String(best_fragment_annotation_modification[1]).toInt();

      }
      else if (best_fragment_annotation.find("+") != std::string::npos)
      {
        std::vector<String> best_fragment_annotation_modification;
        best_fragment_annotation.split("+", best_fragment_annotation_modification);
        mytransition.fragment_type = best_fragment_annotation_modification[0].substr(0, 1);
        mytransition.fragment_nr = String(best_fragment_annotation_modification[0].substr(1)).toInt();
        mytransition.fragment_modification = String(best_fragment_annotation_modification[1]).toInt();
      }
      else
      {
        mytransition.fragment_type = best_fragment_annotation.substr(0, 1);
        mytransition.fragment_nr = String(best_fragment_annotation.substr(1)).toInt();
        mytransition.fragment_modification = 0;
      }

      mytransition.fragment_mzdelta = String(best_fragment_annotation_with_deviation[1]).toDouble();
    }
    else
    {
      // The fragment ion could not be annotated and will likely not be used for detection transitions;
      // we thus skip it and reduce the size of the output TraML.
      return true;
    }

    return false;
  }

  void TransitionTSVFile::cleanupTransitions_(TSVTransition& mytransition)
  {
    // deal with FullPeptideNames like PEPTIDE/2
    std::vector<String> substrings;
    mytransition.FullPeptideName.split("/", substrings);
    if (substrings.size() == 2)
    {
      mytransition.FullPeptideName = substrings[0];
      mytransition.precursor_charge = substrings[1];
    }
  }

  void TransitionTSVFile::TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenMS::TargetedExperiment& exp)
  {
    // For the CV terms, see
    // http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo

    typedef std::vector<OpenMS::TargetedExperiment::Compound> CompoundVectorType;

    CompoundVectorType compounds;
    PeptideVectorType peptides;
    ProteinVectorType proteins;

    std::map<String, int> peptide_map;
    std::map<String, int> compound_map;
    std::map<String, int> protein_map;

    resolveMixedSequenceGroups_(transition_list);

    Size progress = 0;
    startProgress(0, transition_list.size(), "conversion to internal data representation");
    for (std::vector<TSVTransition>::iterator tr_it = transition_list.begin(); tr_it != transition_list.end(); ++tr_it)
    {

      ReactionMonitoringTransition rm_trans;
      createTransition_(tr_it, rm_trans);
      exp.addTransition(rm_trans);

      // check whether we need a new peptide
      if (peptide_map.find(tr_it->group_id) == peptide_map.end() && 
          compound_map.find(tr_it->group_id) == compound_map.end() )
      {
        // should we make a peptide or a compound ?
        if (tr_it->isPeptide())
        {
          OpenMS::TargetedExperiment::Peptide peptide;
          createPeptide_(tr_it, peptide);
          peptides.push_back(peptide);
          peptide_map[peptide.id] = 0;
        }
        else 
        {
          OpenMS::TargetedExperiment::Compound compound;
          createCompound_(tr_it, compound);
          compounds.push_back(compound);
          compound_map[compound.id] = 0;
        }
      }

      // check whether we need a new protein
      if (tr_it->isPeptide() && protein_map.find(tr_it->ProteinName) == protein_map.end())
      {
        OpenMS::TargetedExperiment::Protein protein;
        createProtein_(tr_it, protein);
        proteins.push_back(protein);
        protein_map[tr_it->ProteinName] = 0;
      }

      setProgress(progress++);
    }
    endProgress();

    exp.setCompounds(compounds);
    exp.setPeptides(peptides);
    exp.setProteins(proteins);
  }

  void TransitionTSVFile::TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenSwath::LightTargetedExperiment& exp)
  {
    std::map<String, int> compound_map;
    std::map<String, int> protein_map;

    resolveMixedSequenceGroups_(transition_list);

    Size progress = 0;
    startProgress(0, transition_list.size(), "conversion to internal data representation");
    for (std::vector<TSVTransition>::iterator tr_it = transition_list.begin(); tr_it != transition_list.end(); ++tr_it)
    {
      OpenSwath::LightTransition transition;
      transition.transition_name  = tr_it->transition_name;
      transition.peptide_ref  = tr_it->group_id;
      transition.library_intensity  = tr_it->library_intensity;
      transition.precursor_mz  = tr_it->precursor;
      transition.product_mz  = tr_it->product;
      transition.fragment_charge = 0; // use zero for charge that is not set
      if (!tr_it->fragment_charge.empty() && tr_it->fragment_charge != "NA")
      {
        transition.fragment_charge = tr_it->fragment_charge.toInt();
      }

      transition.decoy = tr_it->decoy;
      transition.detecting_transition = tr_it->detecting_transition;
      transition.identifying_transition = tr_it->identifying_transition;
      transition.quantifying_transition = tr_it->quantifying_transition;

      exp.transitions.push_back(transition);

      // check whether we need a new compound
      if (compound_map.find(tr_it->group_id) == compound_map.end())
      {
        OpenSwath::LightCompound compound;
        if (tr_it->isPeptide())
        {
          OpenMS::TargetedExperiment::Peptide tramlpeptide;
          createPeptide_(tr_it, tramlpeptide);
          OpenSwathDataAccessHelper::convertTargetedCompound(tramlpeptide, compound);
        }
        else
        {
          OpenMS::TargetedExperiment::Compound tramlcompound;
          createCompound_(tr_it, tramlcompound);
          OpenSwathDataAccessHelper::convertTargetedCompound(tramlcompound, compound);
        }
        exp.compounds.push_back(compound);
        compound_map[compound.id] = 0;
      }

      // check whether we need a new protein
      if (tr_it->isPeptide() && protein_map.find(tr_it->ProteinName) == protein_map.end())
      {
        OpenSwath::LightProtein protein;
        protein.id = tr_it->ProteinName;
        protein.sequence = "";
        exp.proteins.push_back(protein);
        protein_map[tr_it->ProteinName] = 0;
      }
      setProgress(progress++);
    }
    endProgress();
  }

  void TransitionTSVFile::resolveMixedSequenceGroups_(std::vector<TransitionTSVFile::TSVTransition>& transition_list)
  {
    // Create temporary map by group label
    std::map< String, std::vector<TSVTransition*> > label_transition_map;
    for (std::vector<TSVTransition>::iterator tr_it = transition_list.begin(); tr_it != transition_list.end(); ++tr_it)
    {
      if (!tr_it->peptide_group_label.empty() )
      {
        label_transition_map[tr_it->peptide_group_label].push_back(& (*tr_it));
      }
    }

    // Iterate through all the group labels and perform sanity check whether the peptide sequence is the same for all of them
    for (std::map< String, std::vector<TSVTransition*> >::iterator pep_it = label_transition_map.begin(); pep_it != label_transition_map.end(); ++pep_it)
    {
      String curr_sequence;
      if (!pep_it->second.empty() )
      { 
        curr_sequence = (*pep_it->second.begin())->PeptideSequence;
      }

      for (std::vector<TSVTransition*>::iterator tr_it = pep_it->second.begin(); tr_it != pep_it->second.end(); ++tr_it)
      {
        // Sanity check: different peptide sequence in the same peptide label group means that something is probably wrong ...
        if (!curr_sequence.empty() && (*tr_it)->PeptideSequence != curr_sequence)
        {

          if (override_group_label_check_)
          {
            // We wont fix it but give out a warning
            LOG_WARN << "Warning: Found multiple peptide sequences for peptide label group " << pep_it->first << 
              //" found multiple peptide sequences: " << curr_sequence << " and " << (*tr_it)->PeptideSequence << 
              ". Since 'override_group_label_check' is on, nothing will be changed." << std::endl;
          }
          else
          {
            // Lets fix it and inform the user
            LOG_WARN << "Warning: Found multiple peptide sequences for peptide label group " << pep_it->first << 
              //" found multiple peptide sequences: " << curr_sequence << " and " << (*tr_it)->PeptideSequence << 
              ". This is most likely an error and to fix this, a new peptide label group will be inferred - " << 
              "to override this decision, please use the override_group_label_check parameter." << std::endl;
            (*tr_it)->peptide_group_label = (*tr_it)->group_id;
          }

        }
      }
    }

  }

  void TransitionTSVFile::createTransition_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::ReactionMonitoringTransition& rm_trans)
  {
    // the following attributes will be stored as meta values (userParam):
    // - annotation (as by SpectraST)
    // the following attributes will be stored as CV values (CV):
    // - collision energy
    // - library intensity (product ion intensity)
    // - decoy / target transition (binary MS:1002007 or MS:1002008)
    // the following attributes will be stored as attributes:
    // - id (native id)
    // the following attributes will be stored in sub-tags:
    // - Precursor:
    //   * target precursor mass isolation window [Q1] (CV Param)
    // - Product:
    //   * charge state (CV Param)
    //   * target product mass isolation window [Q3] (CV Param)
    //   - Interpretation (only best)
    //     * Fragment number (number in series) (CV Param)
    //     * Fragment type (which series) (CV Param)

    rm_trans.setNativeID(tr_it->transition_name);
    rm_trans.setPrecursorMZ(tr_it->precursor);
    rm_trans.setProductMZ(tr_it->product);
    if (tr_it->isPeptide())
    {
      rm_trans.setPeptideRef(tr_it->group_id);
    }
    else
    {
      rm_trans.setCompoundRef(tr_it->group_id);
    }

    rm_trans.setLibraryIntensity(tr_it->library_intensity);
    if (!tr_it->fragment_charge.empty() && tr_it->fragment_charge != "NA")
    {
      OpenMS::ReactionMonitoringTransition::Product p = rm_trans.getProduct();
      p.setChargeState(tr_it->fragment_charge.toInt());
      rm_trans.setProduct(p);
    }

    // add interpretation
    OpenMS::ReactionMonitoringTransition::Product p = rm_trans.getProduct();
    TargetedExperiment::Interpretation interpretation;

    // check if we have any information about the interpretation
    bool interpretation_set = false;
    if (tr_it->fragment_nr != -1 ||
        tr_it->fragment_mzdelta != -1 ||
        tr_it->fragment_modification < 0 ||
        tr_it->fragment_type != "" )
    {
      interpretation_set = true;
    }

    if (tr_it->fragment_nr != -1)
    {
      interpretation.rank = 1; // we only store the best interpretation
    }

    if (tr_it->fragment_nr != -1)
    {
      interpretation.ordinal = tr_it->fragment_nr;
    }

    if (tr_it->fragment_mzdelta != -1)
    {
      CVTerm frag_mzdelta;
      frag_mzdelta.setCVIdentifierRef("MS");
      frag_mzdelta.setAccession("MS:1000904");
      frag_mzdelta.setName("product ion m/z delta");
      frag_mzdelta.setValue(tr_it->fragment_mzdelta);
      interpretation.addCVTerm(frag_mzdelta);
    }

    if (tr_it->fragment_modification < 0)
    {
      CVTerm frag_loss;
      frag_loss.setCVIdentifierRef("MS");
      frag_loss.setAccession("MS:1001524");
      frag_loss.setName("fragment neutral loss");
      frag_loss.setValue(tr_it->fragment_modification);
      interpretation.addCVTerm(frag_loss);
    }

    // figure out which fragment it is
    if (tr_it->fragment_type == "v")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001237");
      ion.setName("frag: v ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "w")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001238");
      ion.setName("frag: w ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "x")
    {
      interpretation.iontype = TargetedExperiment::IonType::XIon;
    }
    else if (tr_it->fragment_type == "y")
    {
      interpretation.iontype = TargetedExperiment::IonType::YIon;
    }
    else if (tr_it->fragment_type == "z")
    {
      interpretation.iontype = TargetedExperiment::IonType::ZIon;
    }
    else if (tr_it->fragment_type == "a")
    {
      interpretation.iontype = TargetedExperiment::IonType::AIon;
    }
    else if (tr_it->fragment_type == "b")
    {
      interpretation.iontype = TargetedExperiment::IonType::BIon;
    }
    else if (tr_it->fragment_type == "c")
    {
      interpretation.iontype = TargetedExperiment::IonType::CIon;
    }
    else if (tr_it->fragment_type == "d")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001236");
      ion.setName("frag: d ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "unknown")
    {
      // unknown means that we should write CV Term "1001240"
      interpretation.iontype = TargetedExperiment::IonType::NonIdentified;
    }
    else if (tr_it->fragment_type == "")
    {
      // empty means that we have no information whatsoever
      interpretation.iontype = TargetedExperiment::IonType::Unannotated;
    }
    else
    {
      interpretation.iontype = TargetedExperiment::IonType::NonIdentified;
    }

    // dont add empty interpretations
    if (interpretation_set) 
    {
      p.addInterpretation(interpretation);
    }
    rm_trans.setProduct(p);

    // add collision energy
    if (tr_it->CE > 0.0)
    {
      CVTerm CE;
      CE.setCVIdentifierRef("MS");
      CE.setAccession("MS:1000045"); // collision energy
      CE.setName("collision energy");
      CE.setValue(tr_it->CE);
      rm_trans.addCVTerm(CE);
    }

    if (!tr_it->decoy)
    {
      rm_trans.setDecoyTransitionType(ReactionMonitoringTransition::TARGET);
    }
    else
    {
      rm_trans.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);
    }

    if (!tr_it->Annotation.empty())
    {
      rm_trans.setMetaValue("annotation", tr_it->Annotation);
    }
    if (tr_it->detecting_transition) {rm_trans.setDetectingTransition(true);}
    else if (!tr_it->detecting_transition) {rm_trans.setDetectingTransition(false);}

    if (tr_it->identifying_transition) {rm_trans.setIdentifyingTransition(true);}
    else if (!tr_it->identifying_transition) {rm_trans.setIdentifyingTransition(false);}

    if (tr_it->quantifying_transition) {rm_trans.setQuantifyingTransition(true);}
    else if (!tr_it->quantifying_transition) {rm_trans.setQuantifyingTransition(false);}
  }

  void TransitionTSVFile::createProtein_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Protein& protein)
  {
    // the following attributes will be stored as CV values (CV):
    // - uniprot accession number (if available)
    // the following attributes will be stored as attributes:
    // - id
    protein.id = tr_it->ProteinName;

    if (!tr_it->uniprot_id.empty())
    {
      // accession numbers
      CVTerm acc;
      OpenMS::DataValue dtype(tr_it->uniprot_id);
      acc.setCVIdentifierRef("MS");
      acc.setAccession("MS:1000885"); // Accession number for a specific protein in a database.
      acc.setName("protein accession");
      acc.setValue(dtype);
      protein.addCVTerm(acc);
    }
  }

  void TransitionTSVFile::interpretRetentionTime_(std::vector<TargetedExperiment::RetentionTime>& retention_times, const OpenMS::DataValue rt_value)
  {
    TargetedExperiment::RetentionTime retention_time;
    retention_time.setRT(rt_value);
    if (retentionTimeInterpretation_ == "iRT")
    {
      retention_time.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::IRT;
      // no unit, since it is iRT (normalized RT)
    }
    else if (retentionTimeInterpretation_ == "seconds" || retentionTimeInterpretation_ == "minutes")
    {
      retention_time.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::LOCAL;
      if (retentionTimeInterpretation_ == "seconds")
      {
        retention_time.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
      }
      else if (retentionTimeInterpretation_ == "minutes")
      {
        retention_time.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::MINUTE;
      }
    }
    retention_times.push_back(retention_time);
  }

  void TransitionTSVFile::createPeptide_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Peptide& peptide)
  {

    // the following attributes will be stored as meta values (userParam):
    //  - full_peptide_name (full unimod peptide name)
    // the following attributes will be stored as CV values (CV):
    // - retention time
    // - charge state
    // - group label
    // the following attributes will be stored as attributes:
    // - id
    // - sequence

    peptide.id = tr_it->group_id;
    peptide.sequence = tr_it->PeptideSequence;

    // per peptide user params
    peptide.setMetaValue("full_peptide_name", tr_it->FullPeptideName);
    if (!tr_it->label_type.empty())
    {
      peptide.setMetaValue("LabelType", tr_it->label_type);
    }

    // per peptide CV terms
    peptide.setPeptideGroupLabel(tr_it->peptide_group_label);
    if (!tr_it->precursor_charge.empty() && tr_it->precursor_charge != "NA")
    {
      peptide.setChargeState(tr_it->precursor_charge.toInt());
    }

    // add retention time for the peptide
    std::vector<TargetedExperiment::RetentionTime> retention_times;
    OpenMS::DataValue rt_value(tr_it->rt_calibrated);
    interpretRetentionTime_(retention_times, rt_value);
    peptide.rts = retention_times;

    // Try to parse full UniMod string including modifications. If we fail, we
    // can force reading and only parse the "naked" sequence.
    std::vector<TargetedExperiment::Peptide::Modification> mods;
    AASequence aa_sequence;
    try
    {
      aa_sequence = AASequence::fromString(tr_it->FullPeptideName);
    } catch (Exception::InvalidValue & e)
    {
      if (force_invalid_mods_)
      {
        std::cout << "Warning while reading file: " << e.what() << std::endl;
        aa_sequence = AASequence::fromString(tr_it->PeptideSequence);
      }
      else 
      {
        std::cerr << "Error while reading file (use force_invalid_mods to override): " << e.what() << std::endl;
        throw;
      }
    }

    std::vector<String> tmp_proteins;
    tmp_proteins.push_back(tr_it->ProteinName);
    peptide.protein_refs = tmp_proteins;

    // check if the naked peptide sequence is equal to the unmodified AASequence
    if (peptide.sequence != aa_sequence.toUnmodifiedString())
    {
      if (force_invalid_mods_)
      {
        // something is wrong, return and do not try and add any modifications
        return;
      }
      LOG_WARN << "Warning: The peptide sequence " << peptide.sequence << " and the full peptide name " << aa_sequence << 
        " are not equal. Please check your input." << std::endl;
      LOG_WARN << "(use force_invalid_mods to override)" << std::endl;
    }

    // Unfortunately, we cannot store an AASequence here but have to work with
    // the TraML modification object.
    // In TraML, the modification the AA starts with residue 1 but the
    // OpenMS objects start with zero -> we start counting with zero here
    // and the TraML handler will add 1 when storing the file.
    if (std::string::npos == tr_it->FullPeptideName.find("["))
    {
      if (aa_sequence.hasNTerminalModification())
      {
        const ResidueModification& rmod = *(aa_sequence.getNTerminalModification());
        addModification_(mods, -1, rmod);
      }
      if (aa_sequence.hasCTerminalModification())
      {
        const ResidueModification& rmod = *(aa_sequence.getCTerminalModification());
        addModification_(mods, aa_sequence.size(), rmod);
      }
      for (Size i = 0; i != aa_sequence.size(); i++)
      {
        if (aa_sequence[i].isModified())
        {
          const ResidueModification& rmod = *(aa_sequence.getResidue(i).getModification());
          addModification_(mods, i, rmod);
        }
      }
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Error, could not parse modifications on " + tr_it->FullPeptideName +
                                       ". Please use unimod / freetext identifiers like PEPT(Phosphorylation)IDE(UniMod:27)A.");
    }

    peptide.mods = mods;

    OPENMS_POSTCONDITION(aa_sequence.toUnmodifiedString() == peptide.sequence,
                         (String("Internal error: the sequences of the naked and modified peptide sequence are unequal(")
                          + aa_sequence.toUnmodifiedString() + " != " + peptide.sequence).c_str())
  }

  void TransitionTSVFile::createCompound_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Compound& compound)
  {

    // the following attributes will be stored as meta values (userParam):
    //  - CompoundName (name of the compound)
    // the following attributes will be stored as CV values (CV):
    // - label type
    // the following attributes will be stored as attributes:
    // - retention time
    // - charge state
    // - SMILES
    // - id

    compound.id = tr_it->group_id;

    compound.molecular_formula = tr_it->SumFormula;
    compound.smiles_string = tr_it->SMILES;
    compound.setMetaValue("CompoundName", tr_it->CompoundName);

    // does this apply to compounds as well?
    if (!tr_it->label_type.empty())
    {
      compound.setMetaValue("LabelType", tr_it->label_type);
    }

    if (!tr_it->precursor_charge.empty() && tr_it->precursor_charge != "NA")
    {
      compound.setChargeState(tr_it->precursor_charge.toInt());
    }

    // add retention time for the compound
    std::vector<TargetedExperiment::RetentionTime> retention_times;
    OpenMS::DataValue rt_value(tr_it->rt_calibrated);
    interpretRetentionTime_(retention_times, rt_value);
    compound.rts = retention_times;
  }

  void TransitionTSVFile::addModification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
                                             int location, const ResidueModification& rmod)
  {
    TargetedExperiment::Peptide::Modification mod;
    mod.location = location;
    mod.mono_mass_delta = rmod.getDiffMonoMass();
    mod.avg_mass_delta = rmod.getDiffAverageMass();
    mod.unimod_id = rmod.getUniModRecordId();
    mods.push_back(mod);
  }

  TransitionTSVFile::TSVTransition TransitionTSVFile::convertTransition_(const ReactionMonitoringTransition* it, OpenMS::TargetedExperiment& targeted_exp)
  {
    TSVTransition mytransition;
    mytransition.precursor = it->getPrecursorMZ();
    mytransition.product = it->getProductMZ();
    mytransition.rt_calibrated = -1;
    mytransition.fragment_type = "";
    mytransition.fragment_nr = -1;
    mytransition.fragment_charge = "NA";

    if (!it->getPeptideRef().empty())
    {
      const OpenMS::TargetedExperiment::Peptide& pep = targeted_exp.getPeptideByRef(it->getPeptideRef());
      mytransition.group_id = it->getPeptideRef();

#ifdef TRANSITIONTSVREADER_TESTING
      LOG_DEBUG << "Peptide rts empty " <<
      pep.rts.empty()  << " or no cv term " << pep.getRetentionTime() << std::endl;
#endif

      if (pep.hasRetentionTime())
      {
        mytransition.rt_calibrated = pep.getRetentionTime();
      }

      mytransition.PeptideSequence = pep.sequence;
      mytransition.ProteinName = "NA";
      mytransition.uniprot_id = "NA";
      if (!pep.protein_refs.empty())
      {
        const OpenMS::TargetedExperiment::Protein& prot = targeted_exp.getProteinByRef(pep.protein_refs[0]);
        mytransition.ProteinName = prot.id;
        if (prot.hasCVTerm("MS:1000885"))
        {
          mytransition.uniprot_id = prot.getCVTerms()["MS:1000885"][0].getValue().toString();
        }
      }

      mytransition.FullPeptideName = TargetedExperimentHelper::getAASequence(pep).toUniModString();

      mytransition.precursor_charge = "NA";
      if (pep.hasCharge())
      {
        mytransition.precursor_charge = String(pep.getChargeState());
      }
      mytransition.peptide_group_label = "NA";
      if (pep.getPeptideGroupLabel() != "")
      {
        mytransition.peptide_group_label = pep.getPeptideGroupLabel();
      }
      if (pep.metaValueExists("LabelType"))
      {
        mytransition.label_type = pep.getMetaValue("LabelType").toString();
      }
    }
    else if (!it->getCompoundRef().empty())
    {
      const OpenMS::TargetedExperiment::Compound& compound = targeted_exp.getCompoundByRef(it->getCompoundRef());
      mytransition.group_id = it->getCompoundRef();

      if (compound.hasRetentionTime())
      {
        mytransition.rt_calibrated = compound.getRetentionTime();
      }

      mytransition.precursor_charge = "NA";
      if (compound.hasCharge())
      {
        mytransition.precursor_charge = String(compound.getChargeState());
      }

      // get metabolomics specific terms
      mytransition.SumFormula = compound.molecular_formula;
      mytransition.SMILES = compound.smiles_string;
      if (compound.metaValueExists("CompoundName"))
      {
        mytransition.CompoundName = compound.getMetaValue("CompoundName");
      }
    }
    else
    {
      // Error?
    }

    if (it->isProductChargeStateSet())
    {
      mytransition.fragment_charge = String(it->getProductChargeState());
    }

    const ReactionMonitoringTransition::Product & product = it->getProduct();
    for (std::vector<TargetedExperiment::Interpretation>::const_iterator
        int_it = product.getInterpretationList().begin(); int_it !=
        product.getInterpretationList().end(); ++int_it)
    {
      // only report first / best interpretation
      if (int_it->rank == 1 || product.getInterpretationList().size() == 1)
      {
        if (int_it->ordinal != 0) mytransition.fragment_nr = int_it->ordinal;

        switch (int_it->iontype)
        {
          case Residue::AIon:
            mytransition.fragment_type = "a";
            break;
          case Residue::BIon:
            mytransition.fragment_type = "b";
            break;
          case Residue::CIon:
            mytransition.fragment_type = "c";
            break;
          case Residue::XIon:
            mytransition.fragment_type = "x";
            break;
          case Residue::YIon:
            mytransition.fragment_type = "y";
            break;
          case Residue::ZIon:
            mytransition.fragment_type = "z";
            break;
          case Residue::Precursor:
            mytransition.fragment_type = "prec";
            break;
          case Residue::BIonMinusH20:
            mytransition.fragment_type = "b-H20";
            break;
          case Residue::YIonMinusH20:
            mytransition.fragment_type = "y-H20";
            break;
          case Residue::BIonMinusNH3:
            mytransition.fragment_type = "b-NH3";
            break;
          case Residue::YIonMinusNH3:
            mytransition.fragment_type = "y-NH3";
            break;
          case Residue::NonIdentified:
            mytransition.fragment_type = "unknown";
            break;
          case Residue::Unannotated:
            // means no annotation and no input cvParam - to write out a cvParam, use Residue::NonIdentified
            mytransition.fragment_type = "";
            break;
          // invalid values
          case Residue::Full: break;
          case Residue::Internal: break;
          case Residue::NTerminal: break;
          case Residue::CTerminal: break;
          case Residue::SizeOfResidueType:
            break;
        }
      }
    }

    mytransition.transition_name = it->getNativeID();
    mytransition.CE = -1;
    if (it->hasCVTerm("MS:1000045"))
    {
      mytransition.CE = it->getCVTerms()["MS:1000045"][0].getValue().toString().toDouble();
    }
    mytransition.library_intensity = -1;
    if (it->getLibraryIntensity() > -100)
    {
      mytransition.library_intensity = it->getLibraryIntensity();
    }
    mytransition.decoy = false;
    if (it->getDecoyTransitionType() == ReactionMonitoringTransition::TARGET)
    {
      mytransition.decoy = false;
    }
    else if (it->getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
    {
      mytransition.decoy = true;
    }
    mytransition.Annotation = "NA";
    if (it->metaValueExists("annotation"))
    {
      mytransition.Annotation = it->getMetaValue("annotation").toString();
    }
    if (it->metaValueExists("Peptidoforms"))
    {
      it->getMetaValue("Peptidoforms").toString().split('|', mytransition.peptidoforms);
    }
    mytransition.detecting_transition = it->isDetectingTransition();
    mytransition.identifying_transition = it->isIdentifyingTransition();
    mytransition.quantifying_transition = it->isQuantifyingTransition();

    return(mytransition);
  }

  void TransitionTSVFile::writeTSVOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> mytransitions;

    Size progress = 0;
    startProgress(0, targeted_exp.getTransitions().size(), "writing OpenSWATH Transition List TSV file");
    for (Size i = 0; i < targeted_exp.getTransitions().size(); i++)
    {
      TransitionTSVFile::TSVTransition mytransition = convertTransition_(&targeted_exp.getTransitions()[i],targeted_exp);
      mytransitions.push_back(mytransition);
      setProgress(progress++);
    }
    endProgress();

    // start writing
    std::ofstream os(filename);
    os.precision(writtenDigits(double()));
    for (Size i = 0; i < header_names_.size(); i++)
    {
      os << header_names_[i];
      if (i != header_names_.size() - 1)
      {
        os << "\t";
      }
    }
    os << std::endl;

    for (std::vector<TSVTransition>::iterator it = mytransitions.begin(); it != mytransitions.end(); ++it)
    {

      String line;
      line +=
        (String)it->precursor                + "\t"
        + (String)it->product                  + "\t"
        + (String)it->precursor_charge         + "\t"
        + (String)it->fragment_charge          + "\t"
        + (String)it->library_intensity        + "\t"
        + (String)it->rt_calibrated            + "\t"
        + (String)it->PeptideSequence          + "\t"
        + (String)it->FullPeptideName          + "\t"
        + (String)it->peptide_group_label      + "\t"
        + (String)it->label_type               + "\t"
        + (String)it->CompoundName             + "\t"
        + (String)it->SumFormula               + "\t"
        + (String)it->SMILES                   + "\t"
        + (String)it->ProteinName              + "\t"
        + (String)it->uniprot_id               + "\t"
        + (String)it->fragment_type            + "\t"
        + (String)it->fragment_nr              + "\t"
        + (String)it->Annotation               + "\t"
        + (String)it->CE                       + "\t"
        + (String)it->group_id                 + "\t"
        + (String)it->transition_name          + "\t"
        + (String)it->decoy                    + "\t"
        + (String)it->detecting_transition     + "\t"
        + (String)it->identifying_transition   + "\t"
        + (String)it->quantifying_transition;

      os << line << std::endl;

    }
    os.close();
  }

  // public methods
  void TransitionTSVFile::convertTargetedExperimentToTSV(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    if (targeted_exp.containsInvalidReferences())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Your input file contains invalid references, cannot process file.");
    }
    writeTSVOutput_(filename, targeted_exp);
  }

  void TransitionTSVFile::convertTSVToTargetedExperiment(const char* filename, FileTypes::Type filetype, OpenMS::TargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> transition_list;
    readUnstructuredTSVInput_(filename, filetype, transition_list);
    TSVToTargetedExperiment_(transition_list, targeted_exp);
  }

  void TransitionTSVFile::convertTSVToTargetedExperiment(const char* filename, FileTypes::Type filetype, OpenSwath::LightTargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> transition_list;
    readUnstructuredTSVInput_(filename, filetype, transition_list);
    TSVToTargetedExperiment_(transition_list, targeted_exp);
  }

  void TransitionTSVFile::validateTargetedExperiment(OpenMS::TargetedExperiment& targeted_exp)
  {
    if (targeted_exp.containsInvalidReferences())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Invalid input, contains duplicate or invalid references");
    }
  }

}

