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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  TransitionTSVReader::TransitionTSVReader() :
    DefaultParamHandler("TransitionTSVReader")
  {
    defaults_.setValue("retentionTimeInterpretation", "iRT", "How to interpret the provided retention time (the retention time column can either be interpreted to be in iRT, minutes or seconds)", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("retentionTimeInterpretation", ListUtils::create<String>("iRT,seconds,minutes"));
    defaults_.setValue("override_group_label_check", "false", "Override an internal check that assures that all members of the same PeptideGroupLabel have the same PeptideSequence (this ensures that only different isotopic forms of the same peptide can be grouped together in the same label group). Only turn this off if you know what you are doing.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("override_group_label_check", ListUtils::create<String>("true,false"));

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  TransitionTSVReader::~TransitionTSVReader()
  {
  }

  void TransitionTSVReader::updateMembers_()
  {
    retentionTimeInterpretation_ = param_.getValue("retentionTimeInterpretation");
    override_group_label_check_ = param_.getValue("override_group_label_check").toBool();
  }

  const char* TransitionTSVReader::strarray_[] =
  {
    "PrecursorMz",
    "ProductMz",
    "Tr_recalibrated",
    "transition_name",
    "CE",
    "LibraryIntensity",
    "transition_group_id",
    "decoy",
    "PeptideSequence",
    "ProteinName",
    "Annotation",
    "FullUniModPeptideName",
    "MissedCleavages",
    "Replicates",
    "NrModifications",
    "PrecursorCharge",
    "ProductCharge",
    "PeptideGroupLabel",
    "LabelType",
    "UniprotID",
    "detecting_transition",
    "identifying_transition",
    "quantifying_transition"
  };

  const std::vector<std::string> TransitionTSVReader::header_names_(strarray_, strarray_ + 23);

  void TransitionTSVReader::getTSVHeader_(const std::string& line, char& delimiter,
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
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Determined your csv/tsv file to have delimiter '" + (String)txt_delimiter + "', but the parsed header has only " + (String)header.size() + " fields instead of the minimal " + (String)min_header_size + ". Please check your input file.");
    }

    int requiredFields[8] = { 0, 1, 3, 5, 6, 7, 8, 9};
    /*
     * required fields:
     *

      "PrecursorMz",
      "ProductMz",
      "transition_name",
      "LibraryIntensity",
      "transition_group_id",
      "decoy",
      "PeptideSequence",
      "ProteinName"
    */
    for (int i = 0; i < 8; i++)
    {
      if (header_dict.find(header_names_[requiredFields[i]]) == header_dict.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "I determined that your your csv/tsv file has the delimiter '" + (String)txt_delimiter +
                                         "'.\nBut the parsed header does not have the required field \"" + (String)header_names_[requiredFields[i]] + "\". Please check your input file.");
      }
    }

  }

  void TransitionTSVReader::readUnstructuredTSVInput_(const char* filename, FileTypes::Type filetype, std::vector<TSVTransition>& transition_list)
  {
    std::ifstream data(filename);
    std::string   line;
    std::string   tmp;

    // read header
    std::vector<std::string>   tmp_line;
    std::vector<std::string>   header;
    std::map<std::string, int> header_dict;
    char delimiter = ',';

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
      header_dict["ProductCharge"] = 8;
      header_dict["SpectraSTFullPeptideName"] = 9;
      header_dict["SpectraSTUnknown"] = 10;
      header_dict["SpectraSTNumberOfProteinsMappedTo"] = 11;
      header_dict["ProteinName"] = 12;
    }
    else
    {
      std::getline(data, line);

      getTSVHeader_(line, delimiter, header, header_dict);
    }

    bool spectrast_legacy = 0; // we will check below if SpectraST was run in legacy (<5.0) mode or if the RT normalization was forgotten.
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
      for (int i = 0; i < tmp_line.size(); i++)
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
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                         "Error reading the file on line " + String(cnt) + ": length of the header and length of the line" +
                                         " do not match: " + String(tmp_line.size()) + " != " + String(header_dict.size()));
      }

      TSVTransition mytransition;

      // Required columns (they are guaranteed to be present, see getTSVHeader_)
      mytransition.precursor                    =                      String(tmp_line[header_dict["PrecursorMz"]]).toDouble();
      mytransition.product                      =                      String(tmp_line[header_dict["ProductMz"]]).toDouble();
      mytransition.library_intensity            =                      String(tmp_line[header_dict["LibraryIntensity"]]).toDouble();
      mytransition.ProteinName                  =                             tmp_line[header_dict["ProteinName"]];
      mytransition.CE                           =  -1.0;
      mytransition.decoy                        =   0;
      mytransition.product_charge              =  -1;
      mytransition.fragment_nr                  =  -1;
      mytransition.fragment_mzdelta             =  -1;
      mytransition.fragment_modification        =   0;
      mytransition.detecting_transition         =   true;
      mytransition.identifying_transition       =   false;
      mytransition.quantifying_transition       =   true;

      if (FileTypes::typeToName(filetype) == "mrm")
      {
        std::vector<String> substrings;
        String(tmp_line[header_dict["SpectraSTFullPeptideName"]]).split("/", substrings);
        AASequence peptide = AASequence::fromString(substrings[0]);

        mytransition.FullPeptideName = peptide.toString();
        mytransition.PeptideSequence = peptide.toUnmodifiedString();
        mytransition.precursor_charge = substrings[1].toInt();

        mytransition.transition_name = String(cnt) + ("_") + String(tmp_line[header_dict["ProteinName"]]) + String("_") + mytransition.FullPeptideName + String("_") + String(tmp_line[header_dict["PrecursorMz"]]) + "_" + String(tmp_line[header_dict["ProductMz"]]);
        mytransition.group_id =  String(tmp_line[header_dict["ProteinName"]]) + String("_") + mytransition.FullPeptideName + String("_") + String(mytransition.precursor_charge);
      }
      else
      {
        mytransition.transition_name              =                             tmp_line[header_dict["transition_name"]];
        mytransition.group_id                     =                             tmp_line[header_dict["transition_group_id"]];
        mytransition.PeptideSequence              =                             tmp_line[header_dict["PeptideSequence"]];
        mytransition.precursor_charge             =  -1;
      }

      if (header_dict.find("RetentionTime") != header_dict.end())
      {
        mytransition.rt_calibrated = String(tmp_line[header_dict["RetentionTime"]]).toDouble();
      }
      else if (header_dict.find("Tr_recalibrated") != header_dict.end())
      {
        mytransition.rt_calibrated = String(tmp_line[header_dict["Tr_recalibrated"]]).toDouble();
      }
      else if (header_dict.find("SpectraSTRetentionTime") != header_dict.end())
      {
        // If SpectraST was run in RT normalization mode, the retention time is annotated as following: "3887.50(57.30)"
        // 3887.50 refers to the non-normalized RT of the individual or consensus run, and 57.30 refers to the normalized
        // iRT.
        size_t start_position = tmp_line[header_dict["SpectraSTRetentionTime"]].find("(");
        if (start_position != std::string::npos)
        {
          ++start_position;
          size_t end_position = tmp_line[header_dict["SpectraSTRetentionTime"]].find(")");
          if (end_position != std::string::npos)
          {
            mytransition.rt_calibrated = String(tmp_line[header_dict["SpectraSTRetentionTime"]].substr(start_position, end_position - start_position)).toDouble();
          }
        }
        else
        {
          // SpectraST was run without RT Normalization mode
          spectrast_legacy = 1;
          mytransition.rt_calibrated = String(tmp_line[header_dict["SpectraSTRetentionTime"]]).toDouble();
        }
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                         "Expected a header named RetentionTime, Tr_recalibrated or SpectraSTRetentionTime but found none");
      }

      if (header_dict.find("Annotation") != header_dict.end())
      {
        mytransition.Annotation = tmp_line[header_dict["Annotation"]];
      }
      if (header_dict.find("CE") != header_dict.end())
      {
        mytransition.CE = String(tmp_line[header_dict["CE"]]).toDouble();
      }
      else if (header_dict.find("CollisionEnergy") != header_dict.end())
      {
        mytransition.CE = String(tmp_line[header_dict["CollisionEnergy"]]).toDouble();
      }

      if (header_dict.find("decoy") != header_dict.end())
      {
        mytransition.decoy                        =                      String(tmp_line[header_dict["decoy"]]).toInt();
      }
      if (header_dict.find("detecting_transition") != header_dict.end())
      {
        if  (String(tmp_line[header_dict["detecting_transition"]]) == "1") { mytransition.detecting_transition = true; }
        else if (String(tmp_line[header_dict["detecting_transition"]]) == "0") { mytransition.detecting_transition = false; }
      }
      if (header_dict.find("identifying_transition") != header_dict.end())
      {
        if  (String(tmp_line[header_dict["identifying_transition"]]) == "1") { mytransition.identifying_transition = true; }
        else if (String(tmp_line[header_dict["identifying_transition"]]) == "0") { mytransition.identifying_transition = false; }
      }
      if (header_dict.find("quantifying_transition") != header_dict.end())
      {
        if  (String(tmp_line[header_dict["quantifying_transition"]]) == "1") { mytransition.quantifying_transition = true; }
        else if (String(tmp_line[header_dict["quantifying_transition"]]) == "0") { mytransition.quantifying_transition = false; }
      }
      if (header_dict.find("FullUniModPeptideName") != header_dict.end())
      {
        mytransition.FullPeptideName              =                             tmp_line[header_dict["FullUniModPeptideName"]];
      }
      else if (header_dict.find("FullPeptideName") != header_dict.end())
      {
        // previously, only FullPeptideName was used and not FullUniModPeptideName
        mytransition.FullPeptideName              =                             tmp_line[header_dict["FullPeptideName"]];
      }
      if (header_dict.find("PrecursorCharge") != header_dict.end())
      {
        mytransition.precursor_charge             =                      String(tmp_line[header_dict["PrecursorCharge"]]).toInt();
      }
      else if (header_dict.find("Charge") != header_dict.end())
      {
        // charge is assumed to be the charge of the precursor
        mytransition.precursor_charge             =                      String(tmp_line[header_dict["Charge"]]).toInt();
      }

      if (header_dict.find("PeptideGroupLabel") != header_dict.end()) 
      {
        mytransition.peptide_group_label = tmp_line[header_dict["PeptideGroupLabel"]];
      }
      if (header_dict.find("LabelType") != header_dict.end())
      {
        mytransition.label_type                   =                             tmp_line[header_dict["LabelType"]];
      }
      if (header_dict.find("UniprotID") != header_dict.end())
      {
        if (tmp_line[header_dict["UniprotID"]] != "NA")
        {
          mytransition.uniprot_id                   =                             tmp_line[header_dict["UniprotID"]];
        }
      }
      if (header_dict.find("FragmentType") != header_dict.end())
      {
        mytransition.fragment_type                =                             tmp_line[header_dict["FragmentType"]];
      }
      if (header_dict.find("ProductCharge") != header_dict.end())
      {
        mytransition.product_charge              =                      String(tmp_line[header_dict["ProductCharge"]]).toInt();
      }
      if (header_dict.find("FragmentSeriesNumber") != header_dict.end())
      {
        mytransition.fragment_nr                  =                      String(tmp_line[header_dict["FragmentSeriesNumber"]]).toInt();
      }
      if (header_dict.find("FragmentMzDelta") != header_dict.end())
      {
        mytransition.fragment_mzdelta = String(tmp_line[header_dict["FragmentMzDelta"]]).toInt();
      }
      if (header_dict.find("FragmentModification") != header_dict.end())
      {
        mytransition.fragment_modification = String(tmp_line[header_dict["FragmentModification"]]).toInt();
      }
      if (header_dict.find("SpectraSTAnnotation") != header_dict.end())
      {
        // Parses SpectraST fragment ion annotations
        // Example: y13^2/0.000,b16-18^2/-0.013,y7-45/0.000
        // Important: m2:8 are not yet supported! See SpectraSTPeakList::annotateInternalFragments for further information
        mytransition.Annotation = tmp_line[header_dict["SpectraSTAnnotation"]];

        std::vector<String> all_fragment_annotations;
        String(tmp_line[header_dict["SpectraSTAnnotation"]]).split(",", all_fragment_annotations);

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
            mytransition.product_charge = String(best_fragment_annotation_charge[1]).toInt();
            best_fragment_annotation = best_fragment_annotation_charge[0];
          }
          else
          {
            mytransition.product_charge = 1;
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
      }

      cleanupTransitions_(mytransition);

      transition_list.push_back(mytransition);

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
      std::cout << mytransition.product_charge << std::endl;
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

  void TransitionTSVReader::cleanupTransitions_(TSVTransition& mytransition)
  {
    mytransition.transition_name  = mytransition.transition_name.remove('"');
    mytransition.transition_name  = mytransition.transition_name.remove('\'');

    mytransition.PeptideSequence  = mytransition.PeptideSequence.remove('"');
    mytransition.PeptideSequence  = mytransition.PeptideSequence.remove('\'');

    mytransition.ProteinName  = mytransition.ProteinName.remove('"');
    mytransition.ProteinName  = mytransition.ProteinName.remove('\'');
    mytransition.ProteinName  = mytransition.ProteinName.remove(',');

    mytransition.Annotation = mytransition.Annotation.remove('"');
    mytransition.Annotation = mytransition.Annotation.remove('\'');

    mytransition.FullPeptideName = mytransition.FullPeptideName.remove('"');
    mytransition.FullPeptideName = mytransition.FullPeptideName.remove('\'');

    mytransition.group_id  = mytransition.group_id.remove('"');
    mytransition.group_id  = mytransition.group_id.remove('\'');

    mytransition.peptide_group_label  = mytransition.peptide_group_label.remove('"');
    mytransition.peptide_group_label  = mytransition.peptide_group_label.remove('\'');

    mytransition.label_type  = mytransition.label_type.remove('"');
    mytransition.label_type  = mytransition.label_type.remove('\'');

    mytransition.fragment_type  = mytransition.fragment_type.remove('"');
    mytransition.fragment_type  = mytransition.fragment_type.remove('\'');

    mytransition.uniprot_id  = mytransition.uniprot_id.remove('"');
    mytransition.uniprot_id  = mytransition.uniprot_id.remove('\'');

    // deal with FullPeptideNames like PEPTIDE/2
    std::vector<String> substrings;
    mytransition.FullPeptideName.split("/", substrings);
    if (substrings.size() == 2)
    {
      mytransition.FullPeptideName = substrings[0];
      mytransition.precursor_charge = substrings[1].toInt();
    }
  }

  void TransitionTSVReader::TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenMS::TargetedExperiment& exp)
  {
    // For the CV terms, see
    // http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo

    PeptideVectorType peptides;
    ProteinVectorType proteins;

    std::map<String, int> peptide_map;
    std::map<String, int> protein_map;

    resolveMixedSequenceGroups_(transition_list);

    Size progress = 0;
    startProgress(0, transition_list.size(), "converting to TraML format");
    for (std::vector<TSVTransition>::iterator tr_it = transition_list.begin(); tr_it != transition_list.end(); ++tr_it)
    {

      ReactionMonitoringTransition rm_trans;
      createTransition_(tr_it, rm_trans);
      exp.addTransition(rm_trans);

      // check whether we need a new peptide
      if (peptide_map.find(tr_it->group_id) == peptide_map.end())
      {
        OpenMS::TargetedExperiment::Peptide peptide;
        createPeptide_(tr_it, peptide);
        peptides.push_back(peptide);
        peptide_map[peptide.id] = 0;
      }

      // check whether we need a new protein
      if (protein_map.find(tr_it->ProteinName) == protein_map.end())
      {
        OpenMS::TargetedExperiment::Protein protein;
        createProtein_(tr_it, protein);


        proteins.push_back(protein);
        protein_map[tr_it->ProteinName] = 0;
      }

      setProgress(progress++);
    }
    endProgress();

    exp.setPeptides(peptides);
    exp.setProteins(proteins);
  }

  void TransitionTSVReader::TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenSwath::LightTargetedExperiment& exp)
  {
    std::map<String, int> peptide_map;
    std::map<String, int> protein_map;

    resolveMixedSequenceGroups_(transition_list);

    OpenMS::TargetedExperiment::Peptide tramlpeptide;

    Size progress = 0;
    startProgress(0, transition_list.size(), "converting to Transition List Format");
    for (std::vector<TSVTransition>::iterator tr_it = transition_list.begin(); tr_it != transition_list.end(); ++tr_it)
    {
      OpenSwath::LightTransition transition;
      transition.transition_name  = tr_it->transition_name;
      transition.peptide_ref  = tr_it->group_id;
      transition.library_intensity  = tr_it->library_intensity;
      transition.precursor_mz  = tr_it->precursor;
      transition.product_mz  = tr_it->product;
      transition.product_charge  = tr_it->product_charge;
      if (tr_it->decoy == 0)
      {
        transition.decoy = false;
      }
      else
      {
        transition.decoy = true;
      }

      transition.detecting_transition = tr_it->detecting_transition;
      transition.identifying_transition = tr_it->identifying_transition;
      transition.quantifying_transition = tr_it->quantifying_transition;

      exp.transitions.push_back(transition);

      // check whether we need a new peptide
      if (peptide_map.find(tr_it->group_id) == peptide_map.end())
      {
        OpenSwath::LightPeptide peptide;
        createPeptide_(tr_it, tramlpeptide);
        OpenSwathDataAccessHelper::convertTargetedPeptide(tramlpeptide, peptide);
        exp.peptides.push_back(peptide);
        peptide_map[peptide.id] = 0;
      }

      // check whether we need a new protein
      if (protein_map.find(tr_it->ProteinName) == protein_map.end())
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

  void TransitionTSVReader::resolveMixedSequenceGroups_(std::vector<TransitionTSVReader::TSVTransition>& transition_list)
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

  void TransitionTSVReader::createTransition_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::ReactionMonitoringTransition& rm_trans)
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
    rm_trans.setPeptideRef(tr_it->group_id);
    rm_trans.setLibraryIntensity(tr_it->library_intensity);
    if (tr_it->product_charge != -1)
    {
      OpenMS::ReactionMonitoringTransition::Product p = rm_trans.getProduct();
      p.setChargeState(tr_it->product_charge);
      rm_trans.setProduct(p);
    }

    // add interpretation
    OpenMS::ReactionMonitoringTransition::Product p = rm_trans.getProduct();
    CVTermList interpretation;

    if (tr_it->fragment_nr != -1)
    {
      CVTerm rank;
      rank.setCVIdentifierRef("MS");
      rank.setAccession("MS:1000926");
      rank.setName("product interpretation rank");
      rank.setValue(1); // we only store the best interpretation
      interpretation.addCVTerm(rank);
    }

    if (tr_it->fragment_nr != -1)
    {
      CVTerm frag_nr;
      frag_nr.setCVIdentifierRef("MS");
      frag_nr.setAccession("MS:1000903");
      frag_nr.setName("product ion series ordinal");
      frag_nr.setValue(tr_it->fragment_nr);
      interpretation.addCVTerm(frag_nr);
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
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001228");
      ion.setName("frag: x ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "y")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001220");
      ion.setName("frag: y ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "z")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001230");
      ion.setName("frag: z ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "a")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001229");
      ion.setName("frag: a ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "b")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001224");
      ion.setName("frag: b ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "c")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001231");
      ion.setName("frag: c ion");
      interpretation.addCVTerm(ion);
    }
    else if (tr_it->fragment_type == "d")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001236");
      ion.setName("frag: d ion");
      interpretation.addCVTerm(ion);
    }
    else
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001240");
      ion.setName("non-identified ion");
      interpretation.addCVTerm(ion);
    }

    p.addInterpretation(interpretation);
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

    if (tr_it->decoy == 0)
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
    if (tr_it->detecting_transition) {rm_trans.setMetaValue("detecting_transition", "true");}
    else if (!tr_it->detecting_transition) {rm_trans.setMetaValue("detecting_transition", "false");}

    if (tr_it->identifying_transition) {rm_trans.setMetaValue("identifying_transition", "true");}
    else if (!tr_it->identifying_transition) {rm_trans.setMetaValue("identifying_transition", "false");}

    if (tr_it->quantifying_transition) {rm_trans.setMetaValue("quantifying_transition", "true");}
    else if (!tr_it->quantifying_transition) {rm_trans.setMetaValue("quantifying_transition", "false");}
  }

  void TransitionTSVReader::createProtein_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Protein& protein)
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

  void TransitionTSVReader::createPeptide_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Peptide& peptide)
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
    if (tr_it->precursor_charge != -1)
    {
      peptide.setChargeState(tr_it->precursor_charge);
    }

    // add retention time for the peptide
    std::vector<TargetedExperiment::RetentionTime> retention_times;
    OpenMS::DataValue rt_value(tr_it->rt_calibrated);

    if (retentionTimeInterpretation_ == "iRT")
    {
      TargetedExperiment::RetentionTime retention_time;

      {
        CVTerm rt;
        rt.setCVIdentifierRef("MS");
        rt.setAccession("MS:1000896"); // normalized RT
        rt.setName("normalized retention time");
        rt.setValue(rt_value);
        retention_time.addCVTerm(rt);
      }

      {
        CVTerm rt;
        rt.setCVIdentifierRef("MS");
        rt.setAccession("MS:1002005"); // iRT
        rt.setName("iRT retention time normalization standard");
        retention_time.addCVTerm(rt);
      }

      retention_times.push_back(retention_time);
    }
    else if (retentionTimeInterpretation_ == "seconds" || retentionTimeInterpretation_ == "minutes")
    {
      TargetedExperiment::RetentionTime retention_time;

      {
        CVTerm rt;

        CVTerm::Unit u;
        if (retentionTimeInterpretation_ == "seconds")
        {
          u.accession = "UO:0000010";
          u.name = "second";
          u.cv_ref = "UO";
        }
        else if (retentionTimeInterpretation_ == "minutes")
        {
          u.accession = "UO:0000031";
          u.name = "minute";
          u.cv_ref = "UO";
        }

        rt.setCVIdentifierRef("MS");
        rt.setAccession("MS:1000895"); // local RT
        rt.setName("local retention time");
        rt.setValue(rt_value);
        rt.setUnit(u);
        retention_time.addCVTerm(rt);
      }

      retention_times.push_back(retention_time);
    }
    peptide.rts = retention_times;

    // try to parse it and get modifications out
    // TODO: at this point we could check whether the modification is actually valid
    // aas.setModification(it->location, "UniMod:" + mo->getAccession().substr(7));
    std::vector<TargetedExperiment::Peptide::Modification> mods;

    AASequence aa_sequence = AASequence::fromString(tr_it->FullPeptideName);

    ModificationsDB* mod_db = ModificationsDB::getInstance();

    // in TraML, the modification the AA starts with residue 1 but the
    // OpenMS objects start with zero -> we start counting with zero here
    // and the TraML handler will add 1 when storing the file.
    if (std::string::npos == tr_it->FullPeptideName.find("["))
    {
      if (!aa_sequence.getNTerminalModification().empty())
      {
        ResidueModification rmod = mod_db->getTerminalModification(aa_sequence.getNTerminalModification(), ResidueModification::N_TERM);
        addModification_(mods, -1, rmod, aa_sequence.getNTerminalModification());
      }
      if (!aa_sequence.getCTerminalModification().empty())
      {
        ResidueModification rmod = mod_db->getTerminalModification(aa_sequence.getCTerminalModification(), ResidueModification::C_TERM);
        addModification_(mods, aa_sequence.size(), rmod, aa_sequence.getCTerminalModification());
      }
      for (Size i = 0; i != aa_sequence.size(); i++)
      {
        if (aa_sequence[i].isModified())
        {
          // search the residue in the modification database (if the sequence is valid, we should find it)
          TargetedExperiment::Peptide::Modification mod;
          ResidueModification rmod = mod_db->getModification(aa_sequence.getResidue(i).getOneLetterCode(),
                                                             aa_sequence.getResidue(i).getModification(), ResidueModification::ANYWHERE);
          addModification_(mods, i, rmod, aa_sequence.getResidue(i).getModification());
        }
      }
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "Error, could not parse modifications on " + tr_it->FullPeptideName +
                                       ". Please use unimod / freetext identifiers like PEPT(Phosphorylation)IDE(UniMod:27)A.");
    }

    peptide.mods = mods;

    std::vector<String> tmp_proteins;
    tmp_proteins.push_back(tr_it->ProteinName);
    peptide.protein_refs = tmp_proteins;

    OPENMS_POSTCONDITION(aa_sequence.toUnmodifiedString() == peptide.sequence,
                         (String("Internal error: the sequences of the naked and modified peptide sequence are unequal(")
                          + aa_sequence.toUnmodifiedString() + " != " + peptide.sequence).c_str())
  }

  void TransitionTSVReader::addModification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
                                             int location, ResidueModification& rmod, const String& name)
  {
    TargetedExperiment::Peptide::Modification mod;
    String unimod_str = rmod.getUniModAccession();
    mod.location = location;
    mod.mono_mass_delta = rmod.getDiffMonoMass();
    mod.avg_mass_delta = rmod.getDiffAverageMass();
    // CV term with the full unimod accession number and name
    CVTerm unimod_name;
    unimod_name.setCVIdentifierRef("UNIMOD");
    unimod_name.setAccession(unimod_str.toUpper());
    unimod_name.setName(name);
    mod.addCVTerm(unimod_name);
    mods.push_back(mod);
  }

  void TransitionTSVReader::writeTSVOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> mytransitions;
    //for (const std::vector<ReactionMonitoringTransition>::iterator it = targeted_exp.getTransitions().begin(); it != targeted_exp.getTransitions().end(); it++)

    Size progress = 0;
    startProgress(0, targeted_exp.getTransitions().size(), "converting to OpenSWATH transition TSV format");
    for (Size i = 0; i < targeted_exp.getTransitions().size(); i++)
    {
      // get the current transition and try to find the corresponding chromatogram
      const ReactionMonitoringTransition* it = &targeted_exp.getTransitions()[i];

      TSVTransition mytransition;

      const OpenMS::TargetedExperiment::Peptide& pep = targeted_exp.getPeptideByRef(it->getPeptideRef());

      mytransition.precursor = it->getPrecursorMZ();
      mytransition.product = it->getProductMZ();
      mytransition.rt_calibrated = -1;

#ifdef TRANSITIONTSVREADER_TESTING
      std::cout << "Peptide rts empty " <<
        pep.rts.empty()  << " or no cv term " << pep.rts[0].hasCVTerm("MS:1000896") << std::endl;
#endif

      if (!pep.rts.empty() && pep.rts[0].hasCVTerm("MS:1000896"))
      {
        mytransition.rt_calibrated = pep.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
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
      mytransition.group_id = it->getPeptideRef();
      mytransition.decoy = 0;
      if (it->getDecoyTransitionType() == ReactionMonitoringTransition::TARGET)
      {
        mytransition.decoy = 0;
      }
      else if (it->getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
      {
        mytransition.decoy = 1;
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
      mytransition.Annotation = "NA";
      if (it->metaValueExists("annotation"))
      {
        mytransition.Annotation = it->getMetaValue("annotation").toString();
      }
      if (it->metaValueExists("detecting_transition"))
      {
        if (it->getMetaValue("detecting_transition").toBool())
        {
          mytransition.detecting_transition = true;
        }
        else if (!it->getMetaValue("detecting_transition").toBool())
        {
          mytransition.detecting_transition = false;
        }
      }
      else
      {
        mytransition.detecting_transition = true;
      }
      if (it->metaValueExists("identifying_transition"))
      {
        if (it->getMetaValue("identifying_transition").toBool())
        {
          mytransition.identifying_transition = true;
        }
        else if (!it->getMetaValue("identifying_transition").toBool())
        {
          mytransition.identifying_transition = false;
        }
      }
      else
      {
        mytransition.identifying_transition = false;
      }
      if (it->metaValueExists("quantifying_transition"))
      {
        if (it->getMetaValue("quantifying_transition").toBool())
        {
          mytransition.quantifying_transition = true;
        }
        else if (!it->getMetaValue("quantifying_transition").toBool())
        {
          mytransition.quantifying_transition = false;
        }
      }
      else
      {
        mytransition.quantifying_transition = true;
      }
      mytransition.FullPeptideName = "";
      {
        // Instead of relying on the full_peptide_name, rather look at the actual modifications!
        OpenSwath::LightPeptide lightpep;
        OpenSwathDataAccessHelper::convertTargetedPeptide(pep, lightpep);
        for (int loc = -1; loc <= (int)lightpep.sequence.size(); loc++)
        {
          if (loc > -1 && loc < (int)lightpep.sequence.size())
          {
            mytransition.FullPeptideName += lightpep.sequence[loc];
          }
          // C-terminal and N-terminal modifications may be at positions -1 or lightpep.sequence
          for (Size modloc = 0; modloc < lightpep.modifications.size(); modloc++)
          {
            if (lightpep.modifications[modloc].location == loc)
            {
              mytransition.FullPeptideName += "(" + lightpep.modifications[modloc].unimod_id + ")";
            }
          }
        }
      }
      mytransition.precursor_charge = -1;
      if (pep.getChargeState() > 0)
      {
        mytransition.precursor_charge = pep.getChargeState();
      }
      mytransition.product_charge = -1;
      if (it->getProductChargeState() > 0)
      {
        mytransition.product_charge = it->getProductChargeState();
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
        + (String)it->rt_calibrated            + "\t"
        + (String)it->transition_name          + "\t"
        + (String)it->CE                       + "\t"
        + (String)it->library_intensity        + "\t"
        + (String)it->group_id                 + "\t"
        + (String)it->decoy                    + "\t"
        + (String)it->PeptideSequence          + "\t"
        + (String)it->ProteinName              + "\t"
        + (String)it->Annotation               + "\t"
        + (String)it->FullPeptideName          + "\t"
        + (String)0                            + "\t"
        + (String)0                            + "\t"
        + (String)0                            + "\t"
        + (String)it->precursor_charge         + "\t"
        + (String)it->product_charge          + "\t"
        + (String)it->peptide_group_label      + "\t"
        + (String)it->label_type               + "\t"
        + (String)it->uniprot_id               + "\t"
        + (String)it->detecting_transition     + "\t"
        + (String)it->identifying_transition   + "\t"
        + (String)it->quantifying_transition;

      os << line << std::endl;

    }
    os.close();
  }

  // public methods
  void TransitionTSVReader::convertTargetedExperimentToTSV(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    writeTSVOutput_(filename, targeted_exp);
  }

  void TransitionTSVReader::convertTSVToTargetedExperiment(const char* filename, FileTypes::Type filetype, OpenMS::TargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> transition_list;
    readUnstructuredTSVInput_(filename, filetype, transition_list);
    TSVToTargetedExperiment_(transition_list, targeted_exp);

  }

  void TransitionTSVReader::convertTSVToTargetedExperiment(const char* filename, FileTypes::Type filetype, OpenSwath::LightTargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> transition_list;
    readUnstructuredTSVInput_(filename, filetype, transition_list);
    TSVToTargetedExperiment_(transition_list, targeted_exp);
  }

  void TransitionTSVReader::validateTargetedExperiment(OpenMS::TargetedExperiment& targeted_exp)
  {
    // check that all proteins ids are unique
    std::map<String, int> unique_protein_map;
    for (ProteinVectorType::const_iterator prot_it = targeted_exp.getProteins().begin(); prot_it != targeted_exp.getProteins().end(); ++prot_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_protein_map.find(prot_it->id) != unique_protein_map.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Found duplicate protein id (must be unique): " + String(prot_it->id));
      }
      unique_protein_map[prot_it->id] = 0;
    }

    // check that all peptide ids are unique
    std::map<String, int> unique_peptide_map;
    for (PeptideVectorType::const_iterator pep_it = targeted_exp.getPeptides().begin(); pep_it != targeted_exp.getPeptides().end(); ++pep_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_peptide_map.find(pep_it->id) != unique_peptide_map.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Found duplicate peptide id (must be unique): " + String(pep_it->id));
      }
      unique_peptide_map[pep_it->id] = 0;
    }

    // check that all transition ids are unique
    std::map<String, int> unique_transition_map;
    for (TransitionVectorType::const_iterator tr_it = targeted_exp.getTransitions().begin(); tr_it != targeted_exp.getTransitions().end(); ++tr_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_transition_map.find(tr_it->getNativeID()) != unique_transition_map.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Found duplicate transition id (must be unique): " + String(tr_it->getNativeID()));
      }
      unique_transition_map[tr_it->getNativeID()] = 0;
    }
  }

}
