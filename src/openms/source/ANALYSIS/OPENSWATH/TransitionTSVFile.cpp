// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <utility>

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
      value = value.remove('\''); //SPWFIXME is this gonna be an issue with modomics IDs that have single quotes?
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
    if (tmp != header_dict.end() && !String(tmp_line[ tmp->second ]).empty())
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
    auto tmp = header_dict.find(header_name);
    if (tmp != header_dict.end() && !String(tmp_line[ tmp->second ]).empty())
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
    if (tmp != header_dict.end() && !String(tmp_line[ tmp->second ]).empty())
    {
      OpenMS::String str_value = tmp_line[ tmp->second ];
      if (str_value == "1" || str_value.toUpper() == "TRUE") value = true;
      else if (str_value == "0" || str_value.toUpper() == "FALSE") value = false;
      else return false;

      // all went well, we set the value and can return
      return true;
    }
    return false;
  }

  TransitionTSVFile::TransitionTSVFile() :
    DefaultParamHandler("TransitionTSVFile")
  {
    defaults_.setValue("retentionTimeInterpretation", "iRT", "How to interpret the provided retention time (the retention time column can either be interpreted to be in iRT, minutes or seconds)", {"advanced"});
    defaults_.setValidStrings("retentionTimeInterpretation", {"iRT","seconds","minutes"});
    defaults_.setValue("override_group_label_check", "false", "Override an internal check that assures that all members of the same PeptideGroupLabel have the same PeptideSequence (this ensures that only different isotopic forms of the same peptide can be grouped together in the same label group). Only turn this off if you know what you are doing.", {"advanced"});
    defaults_.setValidStrings("override_group_label_check", {"true","false"});
    defaults_.setValue("force_invalid_mods", "false", "Force reading even if invalid modifications are encountered (OpenMS may not recognize the modification)", {"advanced"});
    defaults_.setValidStrings("force_invalid_mods", {"true","false"});

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  TransitionTSVFile::~TransitionTSVFile() = default;

  void TransitionTSVFile::updateMembers_()
  {
    retentionTimeInterpretation_ = param_.getValue("retentionTimeInterpretation").toString();
    override_group_label_check_ = param_.getValue("override_group_label_check").toBool();
    force_invalid_mods_ = param_.getValue("force_invalid_mods").toBool();
  }

  const std::vector<std::string> TransitionTSVFile::header_names_ = 
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
    "Adducts",
    "ProteinId",
    "UniprotId",
    "GeneName",
    "FragmentType",
    "FragmentSeriesNumber",
    "Annotation",
    "CollisionEnergy",
    "PrecursorIonMobility",
    "TransitionGroupId",
    "TransitionId",
    "Decoy",
    "DetectingTransition",
    "IdentifyingTransition",
    "QuantifyingTransition",
    "Peptidoforms",
    "OligoName",
    "NuctideSequence",
    "FullNuctideName",
//    "Nucleoforms"
  };


  void TransitionTSVFile::getTSVHeader_(const std::string& line, char& delimiter, std::map<std::string, int>& header_dict) const
  {
    std::string tmp;
    std::vector<std::string> header;
    Size min_header_size = 4;
    const std::array<char, 3> possible_delimiters = {'\t', ';', ','};

    for (auto possible_delimiter : possible_delimiters)
    {
      header.clear();
      std::stringstream lineStream(line);
      delimiter = possible_delimiter;
      while (std::getline(lineStream, tmp, delimiter))
      {
        String tmp2(tmp);
        tmp2 = tmp2.remove('"');
        tmp2 = tmp2.remove('\'');
        tmp2 = tmp2.remove(',');
        header.push_back(tmp2);
      }
      if (header.size() >= min_header_size)
      {
        break; // found the delimiter, got the correct header
      }
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
    std::map<std::string, int> header_dict;
    char delimiter = ',';

    // SpectraST MRM Files do not have a header
    if (filetype == FileTypes::MRM)
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
      TextFile::getLine(data, line);
      getTSVHeader_(line, delimiter, header_dict);
    }

    bool spectrast_legacy = false; // we will check below if SpectraST was run in legacy (<5.0) mode or if the RT normalization was forgotten.
    int cnt = 0;
    while (TextFile::getLine(data, line)) // make sure line endings are handled correctly
    {
      line.push_back(delimiter); // avoid losing last column if it is empty
      std::stringstream lineStream(line);

      while (std::getline(lineStream, tmp, delimiter)) // default getline is fine here, we only want to split the line
      {
        tmp_line.push_back(tmp);
      }
      cnt++;

      #ifdef TRANSITIONTSVREADER_TESTING
      for (Size i = 0; i < tmp_line.size(); i++)
      {
        std::cout << "line " << i << " " << tmp_line[i] << std::endl;
      }

      for (const auto& iter : header_dict)
      {
        std::cout << "header " << iter.first << " " << iter.second << std::endl;
      }
      #endif

      if (tmp_line.size() != header_dict.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Error reading the file on line " + String(cnt) + ": length of the header and length of the line" +
                                         " do not match: " + String(tmp_line.size()) + " != " + String(header_dict.size()));
      }


      // Check for ambiguous entity type: more than one of ProteinId, OligoName, CompoundName not empty
      // isPeptide treats both empty and "NA" as not set.
      // TODO: This breaks a number of tests right now.
      /*
      int non_empty_count = 0;
      String protein_id, oligo_id, compound_name;
      extractName(protein_id, "FullPeptideName", tmp_line, header_dict);
      extractName(oligo_id, "FullNuctideName", tmp_line, header_dict);
      extractName(compound_name, "CompoundName", tmp_line, header_dict);

      if (!protein_id.empty() || !(protein_id == "NA")) non_empty_count++;
      if (!oligo_id.empty() || !(oligo_id == "NA")) non_empty_count++;
      if (!compound_name.empty() || !(compound_name == "NA")) non_empty_count++;
      if (non_empty_count > 1)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Ambiguous entity type: more than one of ProteinId, OligoName, and CompoundName is not empty in the same row: " + String(cnt));
      }
      */

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
          !extractName<double>(mytransition.library_intensity, "RelativeIntensity", tmp_line, header_dict) && // Spectronaut
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
      void(!extractName(mytransition.precursor_charge, "PrecursorCharge", tmp_line, header_dict) &&
      !extractName(mytransition.precursor_charge, "Charge", tmp_line, header_dict)); // charge is assumed to be the charge of the precursor

      void(!extractName(mytransition.fragment_type, "FragmentType", tmp_line, header_dict) &&
      !extractName(mytransition.fragment_type, "FragmentIonType", tmp_line, header_dict)); // Skyline

      void(!extractName(mytransition.fragment_charge, "FragmentCharge", tmp_line, header_dict) &&
      !extractName(mytransition.fragment_charge, "ProductCharge", tmp_line, header_dict));

      void(!extractName<int>(mytransition.fragment_nr, "FragmentSeriesNumber", tmp_line, header_dict) &&
      !extractName<int>(mytransition.fragment_nr, "FragmentNumber", tmp_line, header_dict) &&
      !extractName<int>(mytransition.fragment_nr, "FragmentIonOrdinal", tmp_line, header_dict));

      void(extractName<double>(mytransition.drift_time, "PrecursorIonMobility", tmp_line, header_dict));
      void(extractName<double>(mytransition.fragment_mzdelta, "FragmentMzDelta", tmp_line, header_dict));
      void(extractName<int>(mytransition.fragment_modification, "FragmentModification", tmp_line, header_dict));

      //// Proteomics
      extractName(mytransition.GeneName, "GeneName", tmp_line, header_dict);

      String proteins;
      void(!extractName(proteins, "ProteinName", tmp_line, header_dict) &&
      !extractName(proteins, "ProteinId", tmp_line, header_dict)); // Spectronaut
      if (proteins != "NA" && !proteins.empty())
      {
        proteins.split(';', mytransition.ProteinName);
      }

      void(extractName(mytransition.peptide_group_label, "PeptideGroupLabel", tmp_line, header_dict));

      void(extractName(mytransition.label_type, "LabelType", tmp_line, header_dict));

      void(!extractName(mytransition.PeptideSequence, "PeptideSequence", tmp_line, header_dict) &&
      !extractName(mytransition.PeptideSequence, "Sequence", tmp_line, header_dict) && // Skyline
      !extractName(mytransition.PeptideSequence, "StrippedSequence", tmp_line, header_dict)); // Spectronaut

      void(!extractName(mytransition.FullPeptideName, "FullUniModPeptideName", tmp_line, header_dict) &&
      !extractName(mytransition.FullPeptideName, "FullPeptideName", tmp_line, header_dict) &&
      !extractName(mytransition.FullPeptideName, "ModifiedSequence", tmp_line, header_dict) && // Spectronaut
      !extractName(mytransition.FullPeptideName, "ModifiedPeptideSequence", tmp_line, header_dict));

      //// Nucleotides

      void(!extractName(mytransition.NuctideSequence, "NuctideSequence", tmp_line, header_dict));
      void(!extractName(mytransition.FullNuctideName, "FullNuctideName", tmp_line, header_dict));

      String oligos;
      void(!extractName(oligos, "OligoName", tmp_line, header_dict));
      if (oligos != "NA" && !oligos.empty())
      {
        oligos.split(';', mytransition.OligoName);
      }

      void(!extractName(proteins, "ProteinName", tmp_line, header_dict) &&
      !extractName(proteins, "ProteinId", tmp_line, header_dict)); // Spectronaut
      if (proteins != "NA" && !proteins.empty())
      {
        proteins.split(';', mytransition.ProteinName);
      }


      //// IPF
      String peptidoforms;
      void(!extractName<bool>(mytransition.detecting_transition, "detecting_transition", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.detecting_transition, "DetectingTransition", tmp_line, header_dict));

      void(!extractName<bool>(mytransition.identifying_transition, "identifying_transition", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.identifying_transition, "IdentifyingTransition", tmp_line, header_dict));

      void(!extractName<bool>(mytransition.quantifying_transition, "quantifying_transition", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.quantifying_transition, "QuantifyingTransition", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.quantifying_transition, "Quantitative", tmp_line, header_dict)); // Skyline

      void(extractName(peptidoforms, "Peptidoforms", tmp_line, header_dict));
      peptidoforms.split('|', mytransition.peptidoforms);

      //// Targeted Metabolomics
      void(extractName(mytransition.CompoundName, "CompoundName", tmp_line, header_dict));
      void(extractName(mytransition.SumFormula, "SumFormula", tmp_line, header_dict));
      void(extractName(mytransition.SMILES, "SMILES", tmp_line, header_dict));
      void(extractName(mytransition.Adducts, "Adducts", tmp_line, header_dict));

      //// Meta
      void(extractName(mytransition.Annotation, "Annotation", tmp_line, header_dict));
      
      // UniprotId
      String uniprot_ids;
      void(!extractName(uniprot_ids, "UniprotId", tmp_line, header_dict) &&
      !extractName(uniprot_ids, "UniprotID", tmp_line, header_dict));
      if (uniprot_ids != "NA" && !uniprot_ids.empty())
      {
        uniprot_ids.split(';', mytransition.uniprot_id);
      }

      void(!extractName<double>(mytransition.CE, "CE", tmp_line, header_dict) &&
      !extractName<double>(mytransition.CE, "CollisionEnergy", tmp_line, header_dict));

      // Decoy
      void(!extractName<bool>(mytransition.decoy, "decoy", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.decoy, "Decoy", tmp_line, header_dict) &&
      !extractName<bool>(mytransition.decoy, "IsDecoy", tmp_line, header_dict));

      if (header_dict.find("SpectraSTAnnotation") != header_dict.end())
      {
        skip_transition = spectrastAnnotationExtract(tmp_line[header_dict["SpectraSTAnnotation"]], mytransition);
      }

      //// Generate Group IDs
      // SpectraST
      if (filetype == FileTypes::MRM)
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

  void TransitionTSVFile::spectrastRTExtract(const String& str_inp, double & value, bool & spectrast_legacy)
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

  bool TransitionTSVFile::spectrastAnnotationExtract(const String& str_inp, TSVTransition & mytransition)
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
    NuctideVectorType nuctides;
    OligoVectorType oligos;

    std::map<String, int> compound_map;
    std::map<String, int> peptide_map;
    std::map<String, int> protein_map;
    std::map<String, int> nuctide_map;
    std::map<String, int> oligo_map;

    resolveMixedSequenceGroups_(transition_list);

    Size progress = 0;
    startProgress(0, transition_list.size(), "conversion to internal data representation");
    for (auto tr_it = transition_list.begin(); tr_it != transition_list.end(); ++tr_it)
    {
      ReactionMonitoringTransition rm_trans;
      createTransition_(tr_it, rm_trans);
      exp.addTransition(rm_trans);

      // check whether we need a new peptide
      if (peptide_map.find(tr_it->group_id) == peptide_map.end() &&
          compound_map.find(tr_it->group_id) == compound_map.end() &&
          nuctide_map.find(tr_it->group_id) == nuctide_map.end())
      {
        // should we make a peptide a nuctide or a compound ?
        OpenSwath::TransType type = tr_it->getType();
        if (type == OpenSwath::TransType::PEPTIDE)
        {
          OpenMS::TargetedExperiment::Peptide peptide;
          createPeptide_(tr_it, peptide);
          peptides.push_back(peptide);
          peptide_map[peptide.id] = 0;
        }
        else if (type == OpenSwath::TransType::NUCTIDE)
        {
          OpenMS::TargetedExperiment::Nuctide nuctide;
          createNuctide_(tr_it, nuctide);
          nuctides.push_back(nuctide);
          nuctide_map[nuctide.id] = 0;
        }
        else
        {
          OpenMS::TargetedExperiment::Compound compound;
          createCompound_(tr_it, compound);
          compounds.push_back(compound);
          compound_map[compound.id] = 0;
        }
      }

      // check whether we need new proteins
      if (tr_it->getType() == OpenSwath::TransType::PEPTIDE)
      {
        for (size_t i = 0; i < tr_it->ProteinName.size(); ++i)
        {
          if (protein_map.find(tr_it->ProteinName[i]) == protein_map.end())
          {
            OpenMS::TargetedExperiment::Protein protein;
            String protein_name = tr_it->ProteinName[i];
            String uniprot_id = "";
            if (tr_it->uniprot_id.size() == tr_it->ProteinName.size())
            {
              uniprot_id = tr_it->uniprot_id[i];
            }
            createProtein_(protein_name, uniprot_id, protein);
            proteins.push_back(protein);
            protein_map[tr_it->ProteinName[i]] = 0;
          }
        }
      }
      else if ((tr_it->getType() == OpenSwath::TransType::NUCTIDE))
      {
        for (size_t i = 0; i < tr_it->OligoName.size(); ++i)
        {
          if (oligo_map.find(tr_it->OligoName[i]) == oligo_map.end())
          {
            OpenMS::TargetedExperiment::Oligo oligo;
            String oligo_name = tr_it->OligoName[i];
            createOligo_(oligo_name, oligo);
            oligos.push_back(oligo);
            oligo_map[tr_it->OligoName[i]] = 0;
          }
        }
      }

    setProgress(progress++);
  }
    endProgress();

    exp.setCompounds(compounds);
    exp.setPeptides(peptides);
    exp.setProteins(proteins);
    exp.setNuctides(nuctides);
    exp.setOligos(oligos);

    OPENMS_POSTCONDITION(exp.getTransitions().size() == transition_list.size(), "Input and output list need to have equal size.")
  }

  void TransitionTSVFile::TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenSwath::LightTargetedExperiment& exp)
  {
    std::map<String, int> compound_map;
    std::map<String, int> protein_map;
    std::map<String, int> oligo_map;

    resolveMixedSequenceGroups_(transition_list);

    Size progress = 0;
    startProgress(0, transition_list.size(), "conversion to internal data representation");
    for (auto tr_it = transition_list.cbegin(); tr_it != transition_list.cend(); ++tr_it)
    {
      OpenSwath::LightTransition transition;
      transition.transition_name  = tr_it->transition_name;
      transition.transition_ref  = tr_it->group_id;
      transition.library_intensity  = tr_it->library_intensity;
      transition.precursor_mz  = tr_it->precursor;
      transition.product_mz  = tr_it->product;
      transition.precursor_im = tr_it->drift_time;
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
        if (tr_it->getType()==OpenSwath::TransType::PEPTIDE)
        {
          OpenMS::TargetedExperiment::Peptide tramlpeptide;
          createPeptide_(tr_it, tramlpeptide);
          OpenSwathDataAccessHelper::convertTargetedPeptide(tramlpeptide, compound);
        }
        else if (tr_it->getType()==OpenSwath::TransType::NUCTIDE)
        {
          OpenMS::TargetedExperiment::Nuctide tramlnuctide;
          createNuctide_(tr_it, tramlnuctide);
          OpenSwathDataAccessHelper::convertTargetedNuctide(tramlnuctide, compound);
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

      // check whether we need new proteins
      for (Size i = 0; i < tr_it->ProteinName.size(); ++i)
      {
        if (tr_it->getType()==OpenSwath::TransType::PEPTIDE && protein_map.find(tr_it->ProteinName[i]) == protein_map.end())
        {
          OpenSwath::LightProtein protein;
          protein.id = tr_it->ProteinName[i];
          protein.sequence = "";
          exp.proteins.push_back(protein);
          protein_map[tr_it->ProteinName[i]] = 0;
        }
        else if ((tr_it->getType()==OpenSwath::TransType::NUCTIDE) && oligo_map.find(tr_it->OligoName[i]) == oligo_map.end())
        {
          OpenSwath::LightOligo oligo;
          oligo.id = tr_it->OligoName[i];
          oligo.sequence = "";
          exp.oligos.push_back(oligo);
          oligo_map[tr_it->OligoName[i]] = 0;
        }
      }

      setProgress(progress++);
    }
    endProgress();

    OPENMS_POSTCONDITION(exp.transitions.size() == transition_list.size(), "Input and output list need to have equal size.")
  }

  void TransitionTSVFile::resolveMixedSequenceGroups_(std::vector<TransitionTSVFile::TSVTransition>& transition_list) const
  {
    // Create temporary map by group label
    std::map<String, std::vector<TSVTransition*> > label_transition_map;
    for (auto & tr_it : transition_list)
    {
      if (!tr_it.peptide_group_label.empty())
      {
        label_transition_map[tr_it.peptide_group_label].push_back(&tr_it);
      }
    }

    // Iterate through all the group labels and perform sanity check whether
    // the peptide sequence is the same for all of them
    for (auto & pep_it : label_transition_map)
    {
      String curr_sequence;
      if (!pep_it.second.empty())
      { 
        curr_sequence = (*pep_it.second.begin())->PeptideSequence;
      }

      for (auto & tr_it : pep_it.second)
      {
        // Sanity check: different peptide sequence in the same peptide label
        // group means that something is probably wrong ...
        if (!curr_sequence.empty() && tr_it->PeptideSequence != curr_sequence)
        {
          if (override_group_label_check_)
          {
            // We wont fix it but give out a warning
            OPENMS_LOG_WARN << "Warning: Found multiple peptide sequences for peptide label group " << pep_it.first << 
              ". Since 'override_group_label_check' is on, nothing will be changed." << std::endl;
          }
          else
          {
            // Lets fix it and inform the user
            OPENMS_LOG_WARN << "Warning: Found multiple peptide sequences for peptide label group " << pep_it.first << 
              ". This is most likely an error and to fix this, a new peptide label group will be inferred - " << 
              "to override this decision, please use the override_group_label_check parameter." << std::endl;
            tr_it->peptide_group_label = tr_it->group_id;
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

    rm_trans.setTransGroupRef(tr_it->group_id,
                         tr_it->getType());

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
        !tr_it->fragment_type.empty() )
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
    else if (tr_it->fragment_type.empty())
    {
      // empty means that we have no information whatsoever
      interpretation.iontype = TargetedExperiment::IonType::Unannotated;
    }
    else
    {
      interpretation.iontype = TargetedExperiment::IonType::NonIdentified;
    }

    // don't add empty interpretations
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

    rm_trans.setDetectingTransition(tr_it->detecting_transition);
    rm_trans.setIdentifyingTransition(tr_it->identifying_transition);
    rm_trans.setQuantifyingTransition(tr_it->quantifying_transition);

    if (!tr_it->peptidoforms.empty())
    {
      rm_trans.setMetaValue("Peptidoforms", ListUtils::concatenate(tr_it->peptidoforms, "|"));
    }
  }

  void TransitionTSVFile::createProtein_(String protein_name, const String& uniprot_id, OpenMS::TargetedExperiment::Protein& protein)
  {
    // the following attributes will be stored as CV values (CV):
    // - uniprot accession number (if available)
    // the following attributes will be stored as attributes:
    // - id
    protein.id = std::move(protein_name);

    if (!uniprot_id.empty())
    {
      // accession numbers
      CVTerm acc;
      OpenMS::DataValue dtype(uniprot_id);
      acc.setCVIdentifierRef("MS");
      acc.setAccession("MS:1000885"); // Accession number for a specific protein in a database.
      acc.setName("protein accession");
      acc.setValue(dtype);
      protein.addCVTerm(acc);
    }
  }

  void TransitionTSVFile::createOligo_(String oligo_name, OpenMS::TargetedExperiment::Oligo& oligo)
  {
    // the following attributes will be stored as attributes:
    // - id
    oligo.id = std::move(oligo_name);
  }

  void TransitionTSVFile::interpretRetentionTime_(std::vector<TargetedExperiment::RetentionTime>& retention_times, const OpenMS::DataValue& rt_value)
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

  void TransitionTSVFile::createPeptide_(std::vector<TSVTransition>::const_iterator tr_it, OpenMS::TargetedExperiment::Peptide& peptide)
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
    if (!tr_it->GeneName.empty())
    {
      peptide.setMetaValue("GeneName", tr_it->GeneName);
    }
    if (!tr_it->SumFormula.empty())
    {
      peptide.setMetaValue("SumFormula", tr_it->SumFormula);
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

    // add ion mobility drift time
    if (tr_it->drift_time >= 0.0)
    {
      peptide.setDriftTime(tr_it->drift_time);
    }

    // Try to parse full UniMod string including modifications. If the string
    // is not parseable (e.g. contains invalid modifications), we can force
    // reading and only parse the "naked" sequence.
    // Note: If the user did not provide a modified sequence string, we will
    // fall back to the "naked" sequence by default.
    std::vector<TargetedExperiment::Peptide::Modification> mods;
    AASequence aa_sequence;
    String sequence = tr_it->FullPeptideName;
    if (sequence.empty()) sequence = tr_it->PeptideSequence;
    try
    {
      aa_sequence = AASequence::fromString(sequence);
    } catch (Exception::InvalidValue & e)
    {
      if (force_invalid_mods_)
      {
        // fallback: parse the "naked" peptide sequence which should always work
        OPENMS_LOG_DEBUG << "Invalid sequence when parsing '" << tr_it->FullPeptideName << "'" << std::endl;
        aa_sequence = AASequence::fromString(tr_it->PeptideSequence);
      }
      else
      {
        OPENMS_LOG_DEBUG << "Invalid sequence when parsing '" << tr_it->FullPeptideName << "'" << std::endl;
        std::cerr << "Error while reading file (use 'force_invalid_mods' parameter to override): " << e.what() << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Invalid input, cannot parse: " + tr_it->FullPeptideName);
      }
    }

    peptide.protein_refs = tr_it->ProteinName;

    // check if the naked peptide sequence is equal to the unmodified AASequence
    if (peptide.sequence != aa_sequence.toUnmodifiedString())
    {
      if (force_invalid_mods_)
      {
        // something is wrong, return and do not try and add any modifications
        return;
      }
      OPENMS_LOG_WARN << "Warning: The peptide sequence " << peptide.sequence << " and the full peptide name " << aa_sequence << 
        " are not equal. Please check your input." << std::endl;
      OPENMS_LOG_WARN << "(use force_invalid_mods to override)" << std::endl;
    }

    // Unfortunately, we cannot store an AASequence here but have to work with
    // the TraML modification object.
    // In TraML, the modification the AA starts with residue 1 but the
    // OpenMS objects start with zero -> we start counting with zero here
    // and the TraML handler will add 1 when storing the file.
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

    peptide.mods = mods;

    OPENMS_POSTCONDITION(aa_sequence.toUnmodifiedString() == peptide.sequence,
                         (String("Internal error: the sequences of the naked and modified peptide sequence are unequal(")
                          + aa_sequence.toUnmodifiedString() + " != " + peptide.sequence).c_str())
  }

    void TransitionTSVFile::createNuctide_(std::vector<TSVTransition>::const_iterator tr_it, OpenMS::TargetedExperiment::Nuctide& nuctide)
  {
    // the following attributes will be stored as meta values (userParam):
    //  - full_nuctide_name
    // the following attributes will be stored as CV values (CV):
    // - retention time
    // - charge state
    // - group label
    // the following attributes will be stored as attributes:
    // - id
    // - sequence

    nuctide.id = tr_it->group_id;
    nuctide.sequence = tr_it->NuctideSequence;

    // per nuctide user params
    nuctide.setMetaValue("full_nuctide_name", tr_it->FullNuctideName);
    if (!tr_it->label_type.empty())
    {
      nuctide.setMetaValue("LabelType", tr_it->label_type);
    }
    if (!tr_it->SumFormula.empty())
    {
      nuctide.setMetaValue("SumFormula", tr_it->SumFormula);
    }

    // per peptide CV terms
    nuctide.setNuctideGroupLabel(tr_it->nuctide_group_label);
    if (!tr_it->precursor_charge.empty() && tr_it->precursor_charge != "NA")
    {
      nuctide.setChargeState(tr_it->precursor_charge.toInt());
    }

    // add retention time for the nuctide
    std::vector<TargetedExperiment::RetentionTime> retention_times;
    OpenMS::DataValue rt_value(tr_it->rt_calibrated);
    interpretRetentionTime_(retention_times, rt_value);
    nuctide.rts = retention_times;

    // add ion mobility drift time
    if (tr_it->drift_time >= 0.0)
    {
      nuctide.setDriftTime(tr_it->drift_time);
    }

    NASequence na_sequence;
    String sequence = tr_it->FullPeptideName;
    if (sequence.empty()) sequence = tr_it->PeptideSequence;
    try
    {
      na_sequence = NASequence::fromString(sequence);
    } catch (Exception::InvalidValue & e)
    {
      if (force_invalid_mods_)
      {
        // fallback: parse the "naked" peptide sequence which should always work
        OPENMS_LOG_DEBUG << "Invalid sequence when parsing '" << tr_it->FullNuctideName << "'" << std::endl;
        na_sequence = NASequence::fromString(tr_it->NuctideSequence);
      }
      else
      {
        OPENMS_LOG_DEBUG << "Invalid sequence when parsing '" << tr_it->FullNuctideName << "'" << std::endl;
        std::cerr << "Error while reading file (use 'force_invalid_mods' parameter to override): " << e.what() << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Invalid input, cannot parse: " + tr_it->FullNuctideName);
      }
    }

    nuctide.oligo_refs = tr_it->OligoName;
  }

  void TransitionTSVFile::createCompound_(std::vector<TSVTransition>::const_iterator tr_it, OpenMS::TargetedExperiment::Compound& compound)
  {
    // the following attributes will be stored as meta values (userParam):
    //  - CompoundName (name of the compound)
    //  - Adducts (adduct associated to the compound)
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
    if (!tr_it->Adducts.empty()) compound.setMetaValue("Adducts", tr_it->Adducts);

    // does this apply to compounds as well?
    if (!tr_it->label_type.empty())
    {
      compound.setMetaValue("LabelType", tr_it->label_type);
    }

    // add ion mobility drift time
    if (tr_it->drift_time >= 0.0)
    {
      compound.setDriftTime(tr_it->drift_time);
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
                                           int location,
                                           const ResidueModification& rmod)
  {
    TargetedExperiment::Peptide::Modification mod;
    mod.location = location;
    mod.mono_mass_delta = rmod.getDiffMonoMass();
    mod.avg_mass_delta = rmod.getDiffAverageMass();
    mod.unimod_id = rmod.getUniModRecordId(); // NOTE: will be -1 if not found in UniMod (e.g. user-defined modifications)
    mods.push_back(mod);
  }

  TransitionTSVFile::TSVTransition TransitionTSVFile::convertTransition_(const ReactionMonitoringTransition* it, OpenMS::TargetedExperiment& targeted_exp)
  {
    // Constants for better maintainability and readability
    static const String NA_VALUE = "NA";
    static const double INVALID_TIME = -1.0;
    static const double INVALID_CE = -1.0;
    static const int INVALID_FRAGMENT_NR = -1;
    static const double MIN_LIBRARY_INTENSITY = -100.0;
    static const String CV_COLLISION_ENERGY = "MS:1000045";
    static const String CV_PROTEIN_ACCESSION = "MS:1000885";
    
    // Fragment type mapping for O(1) lookup performance instead of switch statement
    static const std::unordered_map<TargetedExperiment::IonType, String> FRAGMENT_TYPE_MAP = {
      {TargetedExperiment::IonType::AIon, "a"}, {TargetedExperiment::IonType::BIon, "b"}, {TargetedExperiment::IonType::CIon, "c"},
      {TargetedExperiment::IonType::XIon, "x"}, {TargetedExperiment::IonType::YIon, "y"}, {TargetedExperiment::IonType::ZIon, "z"},
      {TargetedExperiment::IonType::Zp1Ion, "z."}, {TargetedExperiment::IonType::Zp2Ion, "z'"}, {TargetedExperiment::IonType::Precursor, "prec"},
      {TargetedExperiment::IonType::BIonMinusH20, "b-H20"}, {TargetedExperiment::IonType::YIonMinusH20, "y-H20"},
      {TargetedExperiment::IonType::BIonMinusNH3, "b-NH3"}, {TargetedExperiment::IonType::YIonMinusNH3, "y-NH3"},
      {TargetedExperiment::IonType::NonIdentified, "unknown"}, {TargetedExperiment::IonType::Unannotated, ""}
    };

    // Input validation
    if (!it)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Invalid transition pointer provided");
    }

    // Initialize transition with safe default values
    TSVTransition mytransition;
    mytransition.precursor = it->getPrecursorMZ();
    mytransition.product = it->getProductMZ();
    mytransition.rt_calibrated = INVALID_TIME;
    mytransition.fragment_type = "";
    mytransition.fragment_nr = INVALID_FRAGMENT_NR;
    mytransition.fragment_charge = NA_VALUE;

    if (!it->getTransGroupRef().empty())
    {
      const String& trans_ref = it->getTransGroupRef();
      mytransition.group_id = trans_ref;

      try
      {
        if (it->getTransGroupType() == OpenSwath::TransType::PEPTIDE)
        {
          const OpenMS::TargetedExperiment::Peptide& pep = targeted_exp.getPeptideByRef(trans_ref);

          #ifdef TRANSITIONTSVREADER_TESTING
          OPENMS_LOG_DEBUG << "Peptide rts empty " << pep.rts.empty()
                           << " or no cv term " << pep.getRetentionTime() << std::endl;
          #endif

          // Set retention time
          if (pep.hasRetentionTime())
          {
            mytransition.rt_calibrated = pep.getRetentionTime();
          }
          
          // Set sequences and names
          mytransition.PeptideSequence = pep.sequence;
          mytransition.FullPeptideName = TargetedExperimentHelper::getAASequence(pep).toUniModString();
          mytransition.GeneName = NA_VALUE;
          
          // Process protein references efficiently
          if (!pep.protein_refs.empty())
          {
            mytransition.ProteinName.reserve(pep.protein_refs.size());
            mytransition.uniprot_id.reserve(pep.protein_refs.size());
            
            for (const auto& prot_ref : pep.protein_refs)
            {
              const OpenMS::TargetedExperiment::Protein& prot = targeted_exp.getProteinByRef(prot_ref);
              mytransition.ProteinName.push_back(prot.id);
              if (prot.hasCVTerm(CV_PROTEIN_ACCESSION))
              {
                mytransition.uniprot_id.push_back(prot.getCVTerms().at(CV_PROTEIN_ACCESSION)[0].getValue().toString());
              }
            }
          }
          
          // Set common properties using ternary operators for cleaner code
          mytransition.drift_time = (pep.getDriftTime() >= 0.0) ? pep.getDriftTime() : INVALID_TIME;
          mytransition.precursor_charge = pep.hasCharge() ? String(pep.getChargeState()) : NA_VALUE;
          mytransition.peptide_group_label = !pep.getPeptideGroupLabel().empty() ? pep.getPeptideGroupLabel() : NA_VALUE;
          
          // Set meta values
          if (pep.metaValueExists("LabelType"))
          {
            mytransition.label_type = pep.getMetaValue("LabelType").toString();
          }
          if (pep.metaValueExists("GeneName"))
          {
            mytransition.GeneName = pep.getMetaValue("GeneName").toString();
          }
        }
        else if (it->getTransGroupType() == OpenSwath::TransType::NUCTIDE)
        {
          const OpenMS::TargetedExperiment::Nuctide& nuc = targeted_exp.getNuctideByRef(trans_ref);

          #ifdef TRANSITIONTSVREADER_TESTING
          OPENMS_LOG_DEBUG << "Nucleotide rts empty " << nuc.rts.empty()
                           << " or no cv term " << nuc.getRetentionTime() << std::endl;
          #endif

          // Set retention time
          if (nuc.hasRetentionTime())
          {
            mytransition.rt_calibrated = nuc.getRetentionTime();
          }
          
          // Set sequences
          mytransition.NuctideSequence = nuc.sequence;
          mytransition.FullNuctideName = TargetedExperimentHelper::getNASequence(nuc).toString();
          
          // Process oligo references efficiently
          if (!nuc.oligo_refs.empty())
          {
            mytransition.OligoName.reserve(nuc.oligo_refs.size());
            mytransition.uniprot_id.reserve(nuc.oligo_refs.size());
            
            for (const auto& oligo_ref : nuc.oligo_refs)
            {
              const OpenMS::TargetedExperiment::Oligo& oligo = targeted_exp.getOligoByRef(oligo_ref);
              mytransition.OligoName.push_back(oligo.id);
              if (oligo.hasCVTerm(CV_PROTEIN_ACCESSION))
              {
                mytransition.uniprot_id.push_back(oligo.getCVTerms().at(CV_PROTEIN_ACCESSION)[0].getValue().toString());
              }
            }
          }
          
          // Set common properties - eliminates the duplication from original code
          mytransition.drift_time = (nuc.getDriftTime() >= 0.0) ? nuc.getDriftTime() : INVALID_TIME;
          mytransition.precursor_charge = nuc.hasCharge() ? String(nuc.getChargeState()) : NA_VALUE;
          mytransition.nuctide_group_label = !nuc.getNuctideGroupLabel().empty() ? nuc.getNuctideGroupLabel() : NA_VALUE;
          
          if (nuc.metaValueExists("LabelType"))
          {
            mytransition.label_type = nuc.getMetaValue("LabelType").toString();
          }
        }
        else if (it->getTransGroupType() == OpenSwath::TransType::COMPOUND)
        {
          const OpenMS::TargetedExperiment::Compound& compound = targeted_exp.getCompoundByRef(trans_ref);

          if (compound.hasRetentionTime())
          {
            mytransition.rt_calibrated = compound.getRetentionTime();
          }

          // Set common properties
          mytransition.drift_time = (compound.getDriftTime() >= 0.0) ? compound.getDriftTime() : INVALID_TIME;
          mytransition.precursor_charge = compound.hasCharge() ? String(compound.getChargeState()) : NA_VALUE;

          // Set metabolomics specific terms
          mytransition.SumFormula = compound.molecular_formula;
          mytransition.SMILES = compound.smiles_string;
          
          if (compound.metaValueExists("CompoundName"))
          {
            mytransition.CompoundName = compound.getMetaValue("CompoundName");
          }
          if (compound.metaValueExists("Adducts"))
          {
            mytransition.Adducts = compound.getMetaValue("Adducts");
          }
        }
        else
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "ReactionMonitoringTransition has unsupported transition type.");
        }
      }
      catch (const Exception::InvalidValue& e)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Referenced entity '" + trans_ref + "' not found in targeted experiment: " + e.what());
      }
    }

    // Set fragment charge if available
    if (it->isProductChargeStateSet())
    {
      mytransition.fragment_charge = String(it->getProductChargeState());
    }

    // Process fragment interpretation with efficient lookup table
    const auto& product = it->getProduct();
    for (const auto& interpretation : product.getInterpretationList())
    {
      // Only report first/best interpretation
      if (interpretation.rank == 1 || product.getInterpretationList().size() == 1)
      {
        if (interpretation.ordinal != 0)
        {
          mytransition.fragment_nr = interpretation.ordinal;
        }

        // Use O(1) lookup instead of large switch statement
        auto fragment_type_it = FRAGMENT_TYPE_MAP.find(interpretation.iontype);
        if (fragment_type_it != FRAGMENT_TYPE_MAP.end())
        {
          mytransition.fragment_type = fragment_type_it->second;
        }
        break; // Only process first valid interpretation
      }
    }

    // Set remaining properties efficiently
    mytransition.transition_name = it->getNativeID();
    
    mytransition.CE = INVALID_CE;
    if (it->hasCVTerm(CV_COLLISION_ENERGY))
    {
      mytransition.CE = it->getCVTerms().at(CV_COLLISION_ENERGY)[0].getValue().toString().toDouble();
    }
    
    mytransition.library_intensity = (it->getLibraryIntensity() > MIN_LIBRARY_INTENSITY) ?
      it->getLibraryIntensity() : INVALID_TIME;
    
    // Simplified decoy logic
    mytransition.decoy = (it->getDecoyTransitionType() == ReactionMonitoringTransition::DECOY);
    
    mytransition.Annotation = it->metaValueExists("annotation") ?
      it->getMetaValue("annotation").toString() : NA_VALUE;
      
    if (it->metaValueExists("Peptidoforms"))
    {
      String(it->getMetaValue("Peptidoforms")).split('|', mytransition.peptidoforms);
    }
    
    mytransition.detecting_transition = it->isDetectingTransition();
    mytransition.identifying_transition = it->isIdentifyingTransition();
    mytransition.quantifying_transition = it->isQuantifyingTransition();

    return mytransition;
  }

  void TransitionTSVFile::writeTSVOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> mytransitions;

    Size progress = 0;
    startProgress(0, targeted_exp.getTransitions().size(), "writing OpenSWATH Transition List TSV file");
    for (const auto& tr : targeted_exp.getTransitions())
    {
      mytransitions.push_back(convertTransition_(&tr, targeted_exp));
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

    for (const auto& it : mytransitions)
    {
      String line;
      line +=
        (String)it.precursor                + "\t"
        + (String)it.product                  + "\t"
        + (String)it.precursor_charge         + "\t"
        + (String)it.fragment_charge          + "\t"
        + (String)it.library_intensity        + "\t"
        + (String)it.rt_calibrated            + "\t"
        + (String)it.PeptideSequence          + "\t"
        + (String)it.FullPeptideName          + "\t"
        + (String)it.peptide_group_label      + "\t"
        + (String)it.label_type               + "\t"
        + (String)it.CompoundName             + "\t"
        + (String)it.SumFormula               + "\t"
        + (String)it.SMILES                   + "\t"
        + (String)it.Adducts                  + "\t"
        + ListUtils::concatenate(it.ProteinName, ";")              + "\t"
        + ListUtils::concatenate(it.uniprot_id, ";")               + "\t"
        + (String)it.GeneName                 + "\t"
        + (String)it.fragment_type            + "\t"
        + (String)it.fragment_nr              + "\t"
        + (String)it.Annotation               + "\t"
        + (String)it.CE                       + "\t"
        + (String)it.drift_time               + "\t"
        + (String)it.group_id                 + "\t"
        + (String)it.transition_name          + "\t"
        + (String)it.decoy                    + "\t"
        + (String)it.detecting_transition     + "\t"
        + (String)it.identifying_transition   + "\t"
        + (String)it.quantifying_transition   + "\t"
        + ListUtils::concatenate(it.peptidoforms, "|") + "\t"
        + ListUtils::concatenate(it.OligoName, "|") + "\t"
        + (String)it.NuctideSequence + "\t"
        + (String)it.FullNuctideName
        ;

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

  void TransitionTSVFile::validateTargetedExperiment(const OpenMS::TargetedExperiment& targeted_exp)
  {
    if (targeted_exp.containsInvalidReferences())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Invalid input, contains duplicate or invalid references");
    }
  }

}

