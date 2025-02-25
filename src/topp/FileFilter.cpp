// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Lars Nilse, Chris Bielow, Hendrik Brauer $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>

#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/COMPARISON/ZhangSimilarityScore.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <memory>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FileFilter FileFilter

@brief Extracts portions of the data from an mzML, featureXML or consensusXML file.
<center>
<table>
    <tr>
        <th ALIGN = "center"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=2> &rarr; FileFilter &rarr;</td>
        <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool yielding output @n in mzML, featureXML @n or consensusXML format</td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool that profits on reduced input </td>
    </tr>

</table>
</center>
With this tool it is possible to extract m/z, retention time and intensity ranges from an input file
and to write all data that lies within the given ranges to an output file.

Depending on the input file type, additional specific operations are possible:
- mzML
    - extract spectra of a certain MS level
    - filter by signal-to-noise estimation
    - filter by scan mode of the spectra
    - filter by scan polarity of the spectra
- remove MS2 scans whose precursor matches identifications (from an idXML file in 'id:blacklist')
- featureXML
    - filter by feature charge
    - filter by feature size (number of subordinate features)
    - filter by overall feature quality
- consensusXML
    - filter by size (number of elements in consensus features)
    - filter by consensus feature charge
    - filter by map (extracts specified maps and re-evaluates consensus centroid)@n e.g. FileFilter -map 2 3 5 -in file1.consensusXML -out file2.consensusXML@n If a single map is specified, the feature itself can be extracted.@n e.g. FileFilter -map 5 -in file1.consensusXML -out file2.featureXML
- featureXML / consensusXML:
- remove items with a certain meta value annotation. Allowing for >, < and = comparisons. List types are compared by length, not content. Integer, Double and String are compared using their build-in operators.
    - filter sequences, e.g. "LYSNLVER" or the modification "(Phospho)"@n e.g. FileFilter -id:sequences_whitelist Phospho -in file1.consensusXML -out file2.consensusXML
    - filter accessions, e.g. "sp|P02662|CASA1_BOVIN"
    - remove features with annotations
    - remove features without annotations
    - remove unassigned peptide identifications
    - filter id with best score of features with multiple peptide identifications@n e.g. FileFilter -id:remove_unannotated_features -id:remove_unassigned_ids -id:keep_best_score_id -in file1.featureXML -out file2.featureXML
    - remove features with id clashes (different sequences mapped to one feature)

The priority of the id-flags is (decreasing order): remove_annotated_features / remove_unannotated_features -> remove_clashes -> keep_best_score_id -> sequences_whitelist / accessions_whitelist

MS2 and higher spectra can be filtered according to precursor m/z (see 'peak_options:pc_mz_range'). This flag can be combined with 'rt' range to filter precursors by RT and m/z.
If you want to extract an MS1 region with untouched MS2 spectra included, you will need to split the dataset by MS level, then use the 'mz' option for MS1 data and 'peak_options:pc_mz_range' for MS2 data. Afterwards merge the two files again. RT can be filtered at any step.

@note For filtering peptide/protein identification data, see the @ref TOPP_IDFilter tool.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FileFilter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FileFilter.html

For the parameters of the S/N algorithm section see the class documentation there: @n
@ref OpenMS::SignalToNoiseEstimatorMedian "peak_options:sn"@n

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileFilter :
  public TOPPBase
{
public:

  TOPPFileFilter() :
    TOPPBase("FileFilter", "Extracts or manipulates portions of data from peak, feature or consensus-feature files.")
  {
  }

private:

  static bool sequenceIsWhiteListed_(const AASequence& peptide_hit_sequence, 
                                     const StringList& whitelist, 
                                     const String& sequence_comparison_method) 
  {
    const String& sequence_str = peptide_hit_sequence.toString();
    const String& sequence_unmodified_str = peptide_hit_sequence.toUnmodifiedString();
    if (sequence_comparison_method == "substring") 
    {
      for (const String & s : whitelist)
      {
        if (sequence_str.hasSubstring(s) || sequence_unmodified_str.hasSubstring(s))
        {
          return true;
        }
      }
    } 
    else if (sequence_comparison_method == "exact")
    {
      for (const String & s : whitelist)
      {
       if (sequence_str == s || sequence_unmodified_str == s)
       {
         return true;
       }
      } 
    }
    else 
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid sequence comparison method given: '" + sequence_comparison_method + "'");
    }
    return false;
  }

  static void replacePrecursorCharge(MSExperiment& e, int charge_in, int charge_out)
  {
    for (auto& s : e.getSpectra())
    {
      for (auto& p : s.getPrecursors())
      {
        if (p.getCharge() == charge_in) { p.setCharge(charge_out); }
      }
    }
  }

  static bool checkPeptideIdentification_(BaseFeature& feature,
                                          const bool remove_annotated_features,
                                          const bool remove_unannotated_features,
                                          const StringList& sequences,
                                          const String& sequence_comparison_method,
                                          const StringList& accessions,
                                          const bool keep_best_score_id,
                                          const bool remove_clashes)
  {
    //flag: remove_annotated_features and non-empty peptideIdentifications
    if (remove_annotated_features && !feature.getPeptideIdentifications().empty())
    {
      return false;
    }
    //flag: remove_unannotated_features and no peptideIdentifications
    if (remove_unannotated_features && feature.getPeptideIdentifications().empty())
    {
      return false;
    }
    //flag: remove_clashes
    if (remove_clashes && !feature.getPeptideIdentifications().empty())
    {
      String temp = feature.getPeptideIdentifications().begin()->getHits().begin()->getSequence().toString();
      //loop over all peptideIdentifications
      for (const PeptideIdentification& pep : feature.getPeptideIdentifications())
      {
        //loop over all peptideHits
        for (const PeptideHit& pep_hit : pep.getHits())
        {
          if (pep_hit.getSequence().toString() != temp)
          {
            return false;
          }
        }
      }
    }
    //flag: keep_best_score_id
    if (keep_best_score_id && !feature.getPeptideIdentifications().empty())
    {
      PeptideIdentification temp = feature.getPeptideIdentifications().front();
      //loop over all peptideIdentifications
      for (const PeptideIdentification& pep : feature.getPeptideIdentifications())
      {
        //loop over all peptideHits
        for (const PeptideHit& pep_hit : pep.getHits())
        {
          if ((pep.isHigherScoreBetter() && pep_hit.getScore() > temp.getHits().front().getScore()) ||
              (!pep.isHigherScoreBetter() && pep_hit.getScore() < temp.getHits().front().getScore()))
          {
            temp = pep;
          }
        }
      }
      feature.setPeptideIdentifications(vector<PeptideIdentification>(1, temp));
      // not filtering sequences or accessions
      if (sequences.empty() && accessions.empty())
      {
        return true;
      }
    }
    //flag: sequences or accessions
    if (!sequences.empty() || !accessions.empty())
    {
      bool sequen = false;
      bool access = false;
      //loop over all peptideIdentifications
      for (const PeptideIdentification& pep_id : feature.getPeptideIdentifications())
      {
        //loop over all peptideHits
        for (const PeptideHit& pep_hit : pep_id.getHits())
        {
          if (sequenceIsWhiteListed_(pep_hit.getSequence(), sequences, sequence_comparison_method)) 
          {
            sequen = true;
          }
          
          //loop over all accessions of the peptideHits
          set<String> protein_accessions = pep_hit.extractProteinAccessionsSet();
          for (set<String>::const_iterator p_acc_it = protein_accessions.begin(); p_acc_it != protein_accessions.end(); ++p_acc_it)
          {
            //loop over all accessions entries of the StringList
            for (StringList::const_iterator acc_it = accessions.begin(); acc_it != accessions.end(); ++acc_it)
            {
              if (p_acc_it->hasSubstring(*acc_it))
              {
                access = true;
              }
            }
          }
        }
      }
      if (!sequences.empty() && !accessions.empty())
      {
        return sequen && access;
      }
      if (!sequences.empty())
      {
        return sequen;
      }
      else
      {
        return access;
      }
    }
    return true;
  }

protected:

  typedef PeakMap MapType;

  void registerOptionsAndFlags_() override
  {
    std::vector<String> formats = ListUtils::create<String>("mzML,featureXML,consensusXML");

    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", formats);

    registerStringOption_("in_type", "<type>", "", "Input file type -- default: determined from file extension or content", false);
    setValidStrings_("in_type", formats);

    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", formats);

    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content", false);
    setValidStrings_("out_type", formats);

    registerStringOption_("rt", "[min]:[max]", ":", "Retention time range to extract", false);
    registerStringOption_("mz", "[min]:[max]", ":", "m/z range to extract (applies to ALL ms levels!)", false);
    registerStringOption_("int", "[min]:[max]", ":", "Intensity range to extract", false);

    registerFlag_("sort", "Sorts the output according to RT and m/z.");

    registerTOPPSubsection_("peak_options", "Peak data options");
    registerDoubleOption_("peak_options:sn", "<s/n ratio>", 0, "Write peaks with S/N > 'sn' values only", false);
    registerIntList_("peak_options:rm_pc_charge", "i j ...", IntList(), "Remove MS(2) spectra with these precursor charges. All spectra without precursor are kept!", false);
    registerStringOption_("peak_options:pc_mz_range", "[min]:[max]", ":", "MSn (n>=2) precursor filtering according to their m/z value. Do not use this flag in conjunction with 'mz', unless you want to actually remove peaks in spectra (see 'mz'). RT filtering is covered by 'rt' and compatible with this flag.", false);
    registerDoubleList_("peak_options:pc_mz_list", "mz_1 mz_2 ...", DoubleList(), "List of m/z values. If a precursor window covers ANY of these values, the corresponding MS/MS spectrum will be kept.", false);
    registerIntList_("peak_options:level", "i j ...", ListUtils::create<Int>("1,2,3"), "MS levels to extract", false);
    registerFlag_("peak_options:sort_peaks", "Sorts the peaks according to m/z");
    registerFlag_("peak_options:no_chromatograms", "No conversion to space-saving real chromatograms, e.g. from SRM scans");
    registerFlag_("peak_options:remove_chromatograms", "Removes chromatograms stored in a file");
    registerFlag_("peak_options:remove_empty", "Removes spectra and chromatograms without peaks.");
    registerStringOption_("peak_options:mz_precision", "32 or 64", 64, "Store base64 encoded m/z data using 32 or 64 bit precision", false);
    setValidStrings_("peak_options:mz_precision", ListUtils::create<String>("32,64"));
    registerStringOption_("peak_options:int_precision", "32 or 64", 32, "Store base64 encoded intensity data using 32 or 64 bit precision", false);
    setValidStrings_("peak_options:int_precision", ListUtils::create<String>("32,64"));
    registerStringOption_("peak_options:indexed_file", "true or false", "true", "Whether to add an index to the file when writing", false);
    setValidStrings_("peak_options:indexed_file", ListUtils::create<String>("true,false"));

    registerStringOption_("peak_options:zlib_compression", "true or false", "false", "Whether to store data with zlib compression (lossless compression)", false);
    setValidStrings_("peak_options:zlib_compression", ListUtils::create<String>("true,false"));

    registerTOPPSubsection_("peak_options:numpress", "Numpress compression for peak data");
    registerStringOption_("peak_options:numpress:masstime", "<compression_scheme>", "none", "Apply MS Numpress compression algorithms in m/z or rt dimension (recommended: linear)", false);
    setValidStrings_("peak_options:numpress:masstime", MSNumpressCoder::NamesOfNumpressCompression, (int)MSNumpressCoder::SIZE_OF_NUMPRESSCOMPRESSION);
    registerDoubleOption_("peak_options:numpress:lossy_mass_accuracy", "<error>", -1.0, "Desired (absolute) m/z accuracy for lossy compression (e.g. use 0.0001 for a mass accuracy of 0.2 ppm at 500 m/z, default uses -1.0 for maximal accuracy).", false, true);
    registerStringOption_("peak_options:numpress:intensity", "<compression_scheme>", "none", "Apply MS Numpress compression algorithms in intensity dimension (recommended: slof or pic)", false);
    setValidStrings_("peak_options:numpress:intensity", MSNumpressCoder::NamesOfNumpressCompression, (int)MSNumpressCoder::SIZE_OF_NUMPRESSCOMPRESSION);
    registerStringOption_("peak_options:numpress:float_da", "<compression_scheme>", "none", "Apply MS Numpress compression algorithms for the float data arrays (recommended: slof or pic)", false);
    setValidStrings_("peak_options:numpress:float_da", MSNumpressCoder::NamesOfNumpressCompression, (int)MSNumpressCoder::SIZE_OF_NUMPRESSCOMPRESSION);

    registerTOPPSubsection_("spectra", "Remove spectra or select spectra (removing all others) with certain properties");
    registerFlag_("spectra:remove_zoom", "Remove zoom (enhanced resolution) scans");

    registerStringOption_("spectra:remove_mode", "<mode>", "", "Remove scans by scan mode", false);
    setValidStrings_("spectra:remove_mode", InstrumentSettings::NamesOfScanMode, (int)InstrumentSettings::SIZE_OF_SCANMODE);

    addEmptyLine_();
    registerStringOption_("spectra:remove_activation", "<activation>", "", "Remove MSn scans where any of its precursors features a certain activation method", false);
    setValidStrings_("spectra:remove_activation", Precursor::NamesOfActivationMethod, (int)Precursor::SIZE_OF_ACTIVATIONMETHOD);

    registerStringOption_("spectra:remove_collision_energy", "[min]:[max]", ":", "Remove MSn scans with a collision energy in the given interval", false);
    registerStringOption_("spectra:remove_isolation_window_width", "[min]:[max]", ":", "Remove MSn scans whose isolation window width is in the given interval", false);

    addEmptyLine_();
    registerFlag_("spectra:select_zoom", "Select zoom (enhanced resolution) scans");
    registerStringOption_("spectra:select_mode", "<mode>", "", "Selects scans by scan mode\n", false);
    setValidStrings_("spectra:select_mode", InstrumentSettings::NamesOfScanMode, (int)InstrumentSettings::SIZE_OF_SCANMODE);
    registerStringOption_("spectra:select_activation", "<activation>", "", "Retain MSn scans where any of its precursors features a certain activation method", false);
    setValidStrings_("spectra:select_activation", Precursor::NamesOfActivationMethod, (int)Precursor::SIZE_OF_ACTIVATIONMETHOD);
    registerStringOption_("spectra:select_collision_energy", "[min]:[max]", ":", "Select MSn scans with a collision energy in the given interval", false);
    registerStringOption_("spectra:select_isolation_window_width", "[min]:[max]", ":", "Select MSn scans whose isolation window width is in the given interval", false);

    addEmptyLine_();
    registerStringOption_("spectra:select_polarity", "<polarity>", "", "Retain MSn scans with a certain scan polarity", false);
    setValidStrings_("spectra:select_polarity", IonSource::NamesOfPolarity, (int)IonSource::SIZE_OF_POLARITY);

    registerTOPPSubsection_("spectra:blackorwhitelist", "Black or white listing of of MS2 spectra by spectral similarity");
    registerInputFile_("spectra:blackorwhitelist:file", "<file>", "",   "Input file containing MS2 spectra that should be retained or removed from the mzML file!\n"
                                                                        "Matching tolerances are taken from 'spectra:blackorwhitelist:similarity_threshold|rt|mz' options.\n", false);
    setValidFormats_("spectra:blackorwhitelist:file", ListUtils::create<String>("mzML"));
    registerDoubleOption_("spectra:blackorwhitelist:similarity_threshold", "<similarity>", -1, "Similarity threshold when matching MS2 spectra. (-1 = disabled).", false);
    registerDoubleOption_("spectra:blackorwhitelist:rt", "tolerance", 0.01, "Retention tolerance [s] when matching precursor positions. (-1 = disabled)", false);
    registerDoubleOption_("spectra:blackorwhitelist:mz", "tolerance", 0.01, "m/z tolerance [Th] when matching precursor positions. (-1 = disabled)", false);
    registerStringOption_("spectra:blackorwhitelist:use_ppm_tolerance", "", "false", "If ppm tolerance should be used. Otherwise Da are used.", false, false);
    registerStringOption_("spectra:blackorwhitelist:blacklist", "", "true", "True: remove matched MS2. False: retain matched MS2 spectra. Other levels are kept", false, false);
    setValidStrings_("spectra:blackorwhitelist:blacklist", ListUtils::create<String>("false,true"));
    setMinFloat_("spectra:blackorwhitelist:similarity_threshold", -1.0);
    setMaxFloat_("spectra:blackorwhitelist:similarity_threshold", 1.0);

    registerStringOption_("spectra:replace_pc_charge", "in_charge:out_charge", ":", "Replaces in_charge with out_charge in all precursors.", false, false);
    addEmptyLine_();
    registerTOPPSubsection_("feature", "Feature data options");
    registerStringOption_("feature:q", "[min]:[max]", ":", "Overall quality range to extract [0:1]", false);

    addEmptyLine_();
    registerTOPPSubsection_("consensus", "Consensus feature data options");
    registerIntList_("consensus:map", "i j ...", ListUtils::create<Int>(""), "Non-empty list of maps to be extracted from a consensus (indices are 0-based).", false);
    registerFlag_("consensus:map_and", "Consensus features are kept only if they contain exactly one feature from each map (as given above in 'map')");

    // black and white listing
    registerTOPPSubsection_("consensus:blackorwhitelist", "Black or white listing of of MS2 spectra by consensus features");
    registerStringOption_("consensus:blackorwhitelist:blacklist", "", "true", "True: remove matched MS2. False: retain matched MS2 spectra. Other levels are kept", false, false);
    setValidStrings_("consensus:blackorwhitelist:blacklist", ListUtils::create<String>("false,true"));

    registerInputFile_("consensus:blackorwhitelist:file", "<file>", "", "Input file containing consensus features whose corresponding MS2 spectra should be removed from the mzML file!\n"
                                                                        "Matching tolerances are taken from 'consensus:blackorwhitelist:rt' and 'consensus:blackorwhitelist:mz' options.\n"
                                                                        "If consensus:blackorwhitelist:maps is specified, only these will be used.\n", false);
    setValidFormats_("consensus:blackorwhitelist:file", ListUtils::create<String>("consensusXML"));
    registerIntList_("consensus:blackorwhitelist:maps", "i j ...", ListUtils::create<Int>(""), "Maps used for black/white list filtering", false);

    registerDoubleOption_("consensus:blackorwhitelist:rt", "tolerance", 60.0, "Retention tolerance [s] for precursor to consensus feature position", false);
    registerDoubleOption_("consensus:blackorwhitelist:mz", "tolerance", 0.01, "m/z tolerance [Th] for precursor to consensus feature position", false);
    registerStringOption_("consensus:blackorwhitelist:use_ppm_tolerance", "", "false", "If ppm tolerance should be used. Otherwise Da are used.", false, false);

    setValidStrings_("consensus:blackorwhitelist:use_ppm_tolerance", ListUtils::create<String>("false,true"));

    setMinFloat_("consensus:blackorwhitelist:rt", 0);
    setMinFloat_("consensus:blackorwhitelist:mz", 0);

    addEmptyLine_();
    registerTOPPSubsection_("f_and_c", "Feature & Consensus data options");
    registerStringOption_("f_and_c:charge", "[min]:[max]", ":", "Charge range to extract", false);
    registerStringOption_("f_and_c:size", "[min]:[max]", ":", "Size range to extract", false);
    registerStringList_("f_and_c:remove_meta", "<name> 'lt|eq|gt' <value>", StringList(), "Expects a 3-tuple (=3 entries in the list), i.e. <name> 'lt|eq|gt' <value>; the first is the name of meta value, followed by the comparison operator (equal, less or greater) and the value to compare to. All comparisons are done after converting the given value to the corresponding data value type of the meta value (for lists, this simply compares length, not content!)!", false);
    registerFlag_("f_and_c:remove_hull", "Remove hull from features.", false);

    addEmptyLine_();
    // XXX: Change description
    registerTOPPSubsection_("id", "ID options. The Priority of the id-flags is: remove_annotated_features / remove_unannotated_features -> remove_clashes -> keep_best_score_id -> sequences_whitelist  / accessions_whitelist");
    registerFlag_("id:remove_clashes", "Remove features with id clashes (different sequences mapped to one feature)", true);
    registerFlag_("id:keep_best_score_id", "in case of multiple peptide identifications, keep only the id with best score");
    registerStringList_("id:sequences_whitelist", "<sequence>", StringList(), "Keep only features containing whitelisted substrings, e.g. features containing LYSNLVER or the modification (Oxidation). To control comparison method used for whitelisting, see 'id:sequence_comparison_method'.", false);
    registerStringOption_("id:sequence_comparison_method", "substring|exact", "substring", "Comparison method used to determine if a feature is whitelisted.", false, true);
    registerStringList_("id:accessions_whitelist", "<accessions>", StringList(), "keep only features with white listed accessions, e.g. sp|P02662|CASA1_BOVIN", false);
    // XXX: Proper description of this parameter.
    setValidStrings_("id:sequence_comparison_method", ListUtils::create<String>("substring,exact"));
    registerFlag_("id:remove_annotated_features", "Remove features with annotations");
    registerFlag_("id:remove_unannotated_features", "Remove features without annotations");
    registerFlag_("id:remove_unassigned_ids", "Remove unassigned peptide identifications");
    registerInputFile_("id:blacklist", "<file>", "", "Input file containing MS2 identifications whose corresponding MS2 spectra should be removed from the mzML file!\n"
                                                     "Matching tolerances are taken from 'id:rt' and 'id:mz' options.\n"
                                                     "This tool will require all IDs to be matched to an MS2 spectrum, and quit with error otherwise. Use 'id:blacklist_imperfect' to allow for mismatches.", false);
    setValidFormats_("id:blacklist", ListUtils::create<String>("idXML"));
    registerDoubleOption_("id:rt", "tolerance", 0.1, "Retention tolerance [s] for precursor to id position", false);
    registerDoubleOption_("id:mz", "tolerance", 0.001, "m/z tolerance [Th] for precursor to id position", false);
    setMinFloat_("id:rt", 0);
    setMinFloat_("id:mz", 0);
    registerFlag_("id:blacklist_imperfect", "Allow for mismatching precursor positions (see 'id:blacklist')");


    addEmptyLine_();
    registerSubsection_("algorithm", "S/N algorithm section");

  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    SignalToNoiseEstimatorMedian<MapType::SpectrumType> sn;
    Param tmp;
    tmp.insert("SignalToNoise:", sn.getParameters());
    return tmp;
  }

  bool checkMetaOk(const MetaInfoInterface& mi, const StringList& meta_info)
  {
    if (!mi.metaValueExists(meta_info[0]))
    {
      return true; // not having the meta value means passing the test
    }
    DataValue v_data = mi.getMetaValue(meta_info[0]);
    DataValue v_user;
    if (v_data.valueType() == DataValue::STRING_VALUE)
    {
      v_user = String(meta_info[2]);
    }
    else if (v_data.valueType() == DataValue::INT_VALUE)
    {
      v_user = String(meta_info[2]).toInt();
    }
    else if (v_data.valueType() == DataValue::DOUBLE_VALUE)
    {
      v_user = String(meta_info[2]).toDouble();
    }
    else if (v_data.valueType() == DataValue::STRING_LIST)
    {
      v_user = (StringList)ListUtils::create<String>(meta_info[2]);
    }
    else if (v_data.valueType() == DataValue::INT_LIST) 
    {
      v_user = ListUtils::create<Int>(meta_info[2]);
    }
    else if (v_data.valueType() == DataValue::DOUBLE_LIST) 
    {
      v_user = ListUtils::create<double>(meta_info[2]);
    }
    else if (v_data.valueType() == DataValue::EMPTY_VALUE)
    {
      v_user = DataValue::EMPTY;
    }
    if (meta_info[1] == "lt")
    {
      return !(v_data < v_user);
    }
    else if (meta_info[1] == "eq")
    {
      return !(v_data == v_user);
    }
    else if (meta_info[1] == "gt")
    {
      return !(v_data > v_user);
    }
    else
    {
      writeLogError_("Internal Error. Meta value filtering got invalid comparison operator ('" + meta_info[1] + "'), which should have been caught before! Aborting!");
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Illegal meta value filtering operator!");
    }
  }

  ExitCodes main_(int, const char**) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input file name and type
    String in = getStringOption_("in");
    FileHandler fh;

    FileTypes::Type in_type = fh.getType(in);
    //only use flag in_type, if the in_type cannot be determined by file
    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = FileTypes::nameToType(getStringOption_("in_type"));
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    //output file name and type
    String out = getStringOption_("out");

    FileTypes::Type out_type = fh.getTypeByFileName(out);

    //only use flag out_type, if the out_type cannot be determined by file
    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = FileTypes::nameToType(getStringOption_("out_type"));
      writeDebug_(String("Output file type: ") + FileTypes::typeToName(out_type), 2);
    }
    //use in_type as out_type, if out_type cannot be determined by file or out_type flag
    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = in_type;
      writeDebug_(String("Output file type: ") + FileTypes::typeToName(out_type), 2);
    }

    bool no_chromatograms(getFlag_("peak_options:no_chromatograms"));

    //ranges
    double mz_l, mz_u, rt_l, rt_u, it_l, it_u, charge_l, charge_u, size_l, size_u, q_l, q_u, pc_left, pc_right, select_collision_l, remove_collision_l, select_collision_u, remove_collision_u, select_isolation_width_l, remove_isolation_width_l, select_isolation_width_u, remove_isolation_width_u, replace_pc_charge_in, replace_pc_charge_out;

    //initialize ranges
    mz_l = rt_l = it_l = charge_l = size_l = q_l = pc_left = select_collision_l = remove_collision_l = select_isolation_width_l = remove_isolation_width_l = replace_pc_charge_in = -1 * numeric_limits<double>::max();
    mz_u = rt_u = it_u = charge_u = size_u = q_u = pc_right = select_collision_u = remove_collision_u = select_isolation_width_u = remove_isolation_width_u = replace_pc_charge_out = numeric_limits<double>::max();

    String rt = getStringOption_("rt");
    String mz = getStringOption_("mz");
    String pc_mz_range = getStringOption_("peak_options:pc_mz_range");
    String it = getStringOption_("int");
    IntList levels = getIntList_("peak_options:level");
    IntList maps = getIntList_("consensus:map");
    double sn = getDoubleOption_("peak_options:sn");
    String charge = getStringOption_("f_and_c:charge");
    String size = getStringOption_("f_and_c:size");
    String q = getStringOption_("feature:q");
    String remove_collision_energy = getStringOption_("spectra:remove_collision_energy");
    String select_collision_energy = getStringOption_("spectra:select_collision_energy");
    String remove_isolation_width = getStringOption_("spectra:remove_isolation_window_width");
    String select_isolation_width = getStringOption_("spectra:select_isolation_window_width");
    String replace_pc_charge = getStringOption_("spectra:replace_pc_charge");

    int mz32 = getStringOption_("peak_options:mz_precision").toInt();
    int int32 = getStringOption_("peak_options:int_precision").toInt();
    bool indexed_file = getStringOption_("peak_options:indexed_file") == "true";
    bool zlib_compression = getStringOption_("peak_options:zlib_compression") == "true";

    //-----------------------------------
    // MS Numpress options
    //-----------------------------------
    MSNumpressCoder::NumpressConfig npconfig_mz, npconfig_int, npconfig_fda;
    npconfig_mz.estimate_fixed_point = true; // critical
    npconfig_int.estimate_fixed_point = true; // critical
    npconfig_fda.estimate_fixed_point = true; // critical
    // npconfig_mz.numpressErrorTolerance = getDoubleOption_("peak_options:numpress:masstime_error");
    // npconfig_int.numpressErrorTolerance = getDoubleOption_("peak_options:numpress:intensity_error");
    // npconfig_fda.numpressErrorTolerance = getDoubleOption_("peak_options:numpress:intensity_error");
    npconfig_mz.setCompression(getStringOption_("peak_options:numpress:masstime"));
    npconfig_int.setCompression(getStringOption_("peak_options:numpress:intensity"));
    npconfig_fda.setCompression(getStringOption_("peak_options:numpress:float_da"));
    double mass_acc = getDoubleOption_("peak_options:numpress:lossy_mass_accuracy");
    npconfig_mz.linear_fp_mass_acc = mass_acc; // set the desired mass accuracy

    //-----------------------------------
    // ID-filtering parameters
    //-----------------------------------
    bool remove_annotated_features = getFlag_("id:remove_annotated_features");
    bool remove_unannotated_features = getFlag_("id:remove_unannotated_features");
    bool remove_unassigned_ids = getFlag_("id:remove_unassigned_ids");
    StringList sequences = getStringList_("id:sequences_whitelist");
    String sequence_comparison_method = getStringOption_("id:sequence_comparison_method");
    StringList accessions = getStringList_("id:accessions_whitelist");
    bool keep_best_score_id = getFlag_("id:keep_best_score_id");
    bool remove_clashes = getFlag_("id:remove_clashes");
    
    bool remove_hulls = getFlag_("f_and_c:remove_hull");

    // convert bounds to numbers
    try
    {
      //rt
      parseRange_(rt, rt_l, rt_u);
      //mz
      parseRange_(mz, mz_l, mz_u);
      //mz precursor
      parseRange_(pc_mz_range, pc_left, pc_right);
      //int
      parseRange_(it, it_l, it_u);
      //charge (features only)
      parseRange_(charge, charge_l, charge_u);
      //size (features and consensus features only)
      parseRange_(size, size_l, size_u);
      //overall quality (features only)
      parseRange_(q, q_l, q_u);
      //remove collision energy
      parseRange_(remove_collision_energy, remove_collision_l, remove_collision_u);
      //select collision energy
      parseRange_(select_collision_energy, select_collision_l, select_collision_u);
      //remove isolation window width
      parseRange_(remove_isolation_width, remove_isolation_width_l, remove_isolation_width_u);
      //select isolation window width
      parseRange_(select_isolation_width, select_isolation_width_l, select_isolation_width_u);
      //parse precursor charge from in to out
      parseRange_(replace_pc_charge, replace_pc_charge_in, replace_pc_charge_out);
    }
    catch (Exception::ConversionError& ce)
    {
      writeLogError_(String("Error: Invalid boundary given: ") + ce.what() + ". Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // sort by RT and m/z
    bool sort = getFlag_("sort");
    writeDebug_("Sorting output data: " + String(sort), 3);

    // handle remove_meta
    StringList meta_info = getStringList_("f_and_c:remove_meta");
    bool remove_meta_enabled = (!meta_info.empty());
    if (remove_meta_enabled && meta_info.size() != 3)
    {
      writeLogError_("Error: Param 'f_and_c:remove_meta' has invalid number of arguments. Expected 3, got " + String(meta_info.size()) + ". Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    if (remove_meta_enabled && !(meta_info[1] == "lt" || meta_info[1] == "eq" || meta_info[1] == "gt"))
    {
      writeLogError_("Error: Param 'f_and_c:remove_meta' has invalid second argument. Expected one of 'lt', 'eq' or 'gt'. Got '" + meta_info[1] + "'. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (in_type == FileTypes::MZML)
    {
      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      FileHandler f;
      f.getOptions().setRTRange(DRange<1>(rt_l, rt_u));
      f.getOptions().setMZRange(DRange<1>(mz_l, mz_u));
      f.getOptions().setIntensityRange(DRange<1>(it_l, it_u));
      f.getOptions().setMSLevels(levels);

      // set precision options
      if (mz32 == 32)
      { 
        f.getOptions().setMz32Bit(true); 
      } 
      else if (mz32 == 64) 
      { 
        f.getOptions().setMz32Bit(false);
      }
      if (int32 == 32)
      {
        f.getOptions().setIntensity32Bit(true); 
      } 
      else if (int32 == 64)
      { 
        f.getOptions().setIntensity32Bit(false);
      }

      // set writing index (e.g. indexedmzML)
      f.getOptions().setWriteIndex(indexed_file);
      f.getOptions().setCompression(zlib_compression);
      // numpress compression
      f.getOptions().setNumpressConfigurationMassTime(npconfig_mz);
      f.getOptions().setNumpressConfigurationIntensity(npconfig_int);
      f.getOptions().setNumpressConfigurationFloatDataArray(npconfig_fda);

      MapType exp;
      f.loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

      // remove spectra with meta values:
      if (remove_meta_enabled)
      {
        MapType exp_tmp;
        for (MapType::ConstIterator it = exp.begin(); it != exp.end(); ++it)
        {
          if (checkMetaOk(*it, meta_info)) exp_tmp.addSpectrum(*it);
        }
        exp.clear(false);
        exp.getSpectra().insert(exp.begin(), exp_tmp.begin(), exp_tmp.end());
      }


      if (!no_chromatograms)
      {
        // convert the spectra chromatograms to real chromatograms
        ChromatogramTools chrom_tools;
        chrom_tools.convertSpectraToChromatograms(exp, true);
      }

      bool remove_chromatograms(getFlag_("peak_options:remove_chromatograms"));
      if (remove_chromatograms)
      {
        exp.setChromatograms(vector<MSChromatogram >());
      }

      bool remove_empty = getFlag_("peak_options:remove_empty");
      if (remove_empty)
      {
        auto& spectra = exp.getSpectra();
        spectra.erase(
          remove_if(spectra.begin(), spectra.end(), [](const MSSpectrum & s){ return s.empty();} )
          ,spectra.end());
        auto& chroms = exp.getChromatograms();
        chroms.erase(
          remove_if(chroms.begin(), chroms.end(), [](const MSChromatogram & c){ return c.empty();} )
          ,chroms.end());
      }

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

      // remove forbidden precursor charges
      IntList rm_pc_charge = getIntList_("peak_options:rm_pc_charge");
      if (!rm_pc_charge.empty())
      {
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), HasPrecursorCharge<MapType::SpectrumType>(rm_pc_charge, false)), exp.end());
      }

      // remove precursors out of certain m/z range for all spectra with a precursor (MS2 and above)
      if (!pc_mz_range.empty())
      {
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), InPrecursorMZRange<MapType::SpectrumType>(pc_left, pc_right, true)), exp.end());
      }

      // keep MS/MS spectra whose precursors cover at least of the given m/z values
      std::vector<double> vec_mz = getDoubleList_("peak_options:pc_mz_list");
      if (!vec_mz.empty())
      {
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), IsInIsolationWindow<MapType::SpectrumType>(vec_mz, true)), exp.end());
      }


      // remove by scan mode (might be a lot of spectra)
      String remove_mode = getStringOption_("spectra:remove_mode");
      if (!remove_mode.empty())
      {
        writeDebug_("Removing mode: " + remove_mode, 3);
        for (Size i = 0; i < InstrumentSettings::SIZE_OF_SCANMODE; ++i)
        {
          if (InstrumentSettings::NamesOfScanMode[i] == remove_mode)
          {
            exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), HasScanMode<MapType::SpectrumType>((InstrumentSettings::ScanMode)i)), exp.end());
          }
        }
      }

      //select by scan mode (might be a lot of spectra)
      String select_mode = getStringOption_("spectra:select_mode");
      if (!select_mode.empty())
      {
        writeDebug_("Selecting mode: " + select_mode, 3);
        for (Size i = 0; i < InstrumentSettings::SIZE_OF_SCANMODE; ++i)
        {
          if (InstrumentSettings::NamesOfScanMode[i] == select_mode)
          {
            exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), HasScanMode<MapType::SpectrumType>((InstrumentSettings::ScanMode)i, true)), exp.end());
          }
        }
      }

      //remove by activation mode (might be a lot of spectra)
      String remove_activation = getStringOption_("spectra:remove_activation");
      if (!remove_activation.empty())
      {
        writeDebug_("Removing scans with activation mode: " + remove_activation, 3);
        for (Size i = 0; i < Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
        {
          if (Precursor::NamesOfActivationMethod[i] == remove_activation)
          {
            exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), HasActivationMethod<MapType::SpectrumType>(ListUtils::create<String>(remove_activation))), exp.end());
          }
        }
      }

      //select by activation mode
      String select_activation = getStringOption_("spectra:select_activation");
      if (!select_activation.empty())
      {
        writeDebug_("Selecting scans with activation mode: " + select_activation, 3);
        for (Size i = 0; i < Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
        {
          if (Precursor::NamesOfActivationMethod[i] == select_activation)
          {
            exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), HasActivationMethod<MapType::SpectrumType>(ListUtils::create<String>(select_activation), true)), exp.end());
          }
        }
      }

      //select by scan polarity
      String select_polarity = getStringOption_("spectra:select_polarity");
      if (!select_polarity.empty())
      {
        writeDebug_("Selecting polarity: " + select_polarity, 3);
        for (Size i = 0; i < IonSource::SIZE_OF_POLARITY; ++i)
        {
          if (IonSource::NamesOfPolarity[i] == select_polarity)
          {
            exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), HasScanPolarity<MapType::SpectrumType>((IonSource::Polarity)i, true)), exp.end());
          }
        }
      }

      //remove zoom scans (might be a lot of spectra)
      if (getFlag_("spectra:remove_zoom"))
      {
        writeDebug_("Removing zoom scans", 3);
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), IsZoomSpectrum<MapType::SpectrumType>()), exp.end());
      }

      if (getFlag_("spectra:select_zoom"))
      {
        writeDebug_("Selecting zoom scans", 3);
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), IsZoomSpectrum<MapType::SpectrumType>(true)), exp.end());
      }

      //remove based on collision energy
      if (remove_collision_l != -1 * numeric_limits<double>::max() || remove_collision_u != numeric_limits<double>::max())
      {
        writeDebug_(String("Removing collision energy scans in the range: ") + remove_collision_l + ":" + remove_collision_u, 3);
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), IsInCollisionEnergyRange<PeakMap::SpectrumType>(remove_collision_l, remove_collision_u)), exp.end());
      }
      if (select_collision_l != -1 * numeric_limits<double>::max() || select_collision_u != numeric_limits<double>::max())
      {
        writeDebug_(String("Selecting collision energy scans in the range: ") + select_collision_l + ":" + select_collision_u, 3);
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), IsInCollisionEnergyRange<PeakMap::SpectrumType>(select_collision_l, select_collision_u, true)), exp.end());
      }

      //remove based on isolation window size
      if (remove_isolation_width_l != -1 * numeric_limits<double>::max() || remove_isolation_width_u != numeric_limits<double>::max())
      {
        writeDebug_(String("Removing isolation windows with width in the range: ") + remove_isolation_width_l + ":" + remove_isolation_width_u, 3);
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), IsInIsolationWindowSizeRange<PeakMap::SpectrumType>(remove_isolation_width_l, remove_isolation_width_u)), exp.end());
      }
      if (select_isolation_width_l != -1 * numeric_limits<double>::max() || select_isolation_width_u != numeric_limits<double>::max())
      {
        writeDebug_(String("Selecting isolation windows with width in the range: ") + select_isolation_width_l + ":" + select_isolation_width_u, 3);
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), IsInIsolationWindowSizeRange<PeakMap::SpectrumType>(select_isolation_width_l, select_isolation_width_u, true)), exp.end());
      }

      // reannoate precursor charge if both range values are set
      if (replace_pc_charge_in != -1 * numeric_limits<double>::max() && replace_pc_charge_out != numeric_limits<double>::max())
      {
        replacePrecursorCharge(exp, (int)replace_pc_charge_in, (int)replace_pc_charge_out);
      }        

      //remove empty scans
      exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), IsEmptySpectrum<MapType::SpectrumType>()), exp.end());

      //sort
      if (sort)
      {
        exp.sortSpectra(true);
        if (getFlag_("peak_options:sort_peaks"))
        {
          OPENMS_LOG_INFO << "Info: Using 'peak_options:sort_peaks' in combination with 'sort' is redundant, since 'sort' implies 'peak_options:sort_peaks'." << std::endl;
        }
      }
      else if (getFlag_("peak_options:sort_peaks"))
      {
        for (Size i = 0; i < exp.size(); ++i)
        {
          exp[i].sortByPosition();
        }
      }

      // calculate S/N values and delete data points below S/N threshold
      if (sn > 0)
      {
        SignalToNoiseEstimatorMedian<MapType::SpectrumType> snm;
        Param const& dc_param = getParam_().copy("algorithm:SignalToNoise:", true);
        snm.setParameters(dc_param);
        for (auto& spec : exp)
        {
          snm.init(spec);
          for (Size i = 0; i != spec.size(); ++i)
          {
            if (snm.getSignalToNoise(i) < sn) spec[i].setIntensity(0);
          }
          spec.erase(remove_if(spec.begin(), spec.end(), InIntensityRange<MapType::PeakType>(1, numeric_limits<MapType::PeakType::IntensityType>::max(), true)), spec.end());
        }
      }

      //
      String id_blacklist = getStringOption_("id:blacklist");
      if (!id_blacklist.empty())
      {
        OPENMS_LOG_INFO << "Filtering out MS2 spectra from raw file using blacklist ..." << std::endl;
        bool blacklist_imperfect = getFlag_("id:blacklist_imperfect");

        int ret = filterByBlackList(exp, id_blacklist, blacklist_imperfect, getDoubleOption_("id:rt"), getDoubleOption_("id:mz"));
        if (ret != EXECUTION_OK)
        {
          return (ExitCodes)ret;
        }
      }

      // check if filtering by consensus feature is enabled
      String consensus_blackorwhitelist = getStringOption_("consensus:blackorwhitelist:file");
      if (!consensus_blackorwhitelist.empty())
      {
        OPENMS_LOG_INFO << "Filtering MS2 spectra from raw file using consensus features ..." << std::endl;
        IntList il = getIntList_("consensus:blackorwhitelist:maps");
        set<UInt64> maps(il.begin(), il.end());
        double rt_tol = getDoubleOption_("consensus:blackorwhitelist:rt");
        double mz_tol = getDoubleOption_("consensus:blackorwhitelist:mz");
        bool is_ppm = getStringOption_("consensus:blackorwhitelist:use_ppm_tolerance") == "false" ? false : true;
        bool is_blacklist = getStringOption_("consensus:blackorwhitelist:blacklist") == "true" ? true : false;
        
        ConsensusMap consensus_map;
        FileHandler cxml_file;
        cxml_file.loadConsensusFeatures(consensus_blackorwhitelist, consensus_map);
        consensus_map.sortByMZ();

        int ret = filterByBlackOrWhiteList(is_blacklist, exp, consensus_map, rt_tol, mz_tol, is_ppm, maps);
        if (ret != EXECUTION_OK)
        { 
          return (ExitCodes)ret; 
        }
      }

      // filter spectra if they occur in spectra:blackorwhitelist:file 
      // (determined by comparing rt/mz/similarity)
      String lib_file_name = getStringOption_("spectra:blackorwhitelist:file");
      if (!lib_file_name.empty())
      {
        OPENMS_LOG_INFO << "Filtering MS2 spectra based on precursor rt, mz, and spectral similarity ..." << std::endl;
        double tol_rt = getDoubleOption_("spectra:blackorwhitelist:rt");
        double tol_mz = getDoubleOption_("spectra:blackorwhitelist:mz");
        double tol_sim = getDoubleOption_("spectra:blackorwhitelist:similarity_threshold");
        bool is_ppm = getStringOption_("spectra:blackorwhitelist:use_ppm_tolerance") == "true" ? true : false; 
        bool is_blacklist = getStringOption_("spectra:blackorwhitelist:blacklist") == "true" ? true : false;

        PeakMap lib_file;
        FileHandler().loadExperiment(lib_file_name, lib_file, {FileTypes::MZML});

        int ret = filterByBlackOrWhiteList(is_blacklist, exp, lib_file, tol_rt, tol_mz, tol_sim, is_ppm);
        if (ret != EXECUTION_OK)
        { 
          return (ExitCodes)ret;
        }
      }



      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FILTERING));
      f.storeExperiment(out, exp,{FileTypes::MZML}, log_type_);
    }
    else if (in_type == FileTypes::FEATUREXML || in_type == FileTypes::CONSENSUSXML)
    {
      bool meta_ok = true; // assume true by default (as meta might not be checked below)

      if (in_type == FileTypes::FEATUREXML)
      {
        //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------

        FeatureMap feature_map;
        FileHandler f;
        //f.setLogType(log_type_);
        // this does not work yet implicitly - not supported by FeatureXMLFile
        f.getFeatOptions().setRTRange(DRange<1>(rt_l, rt_u));
        f.getFeatOptions().setMZRange(DRange<1>(mz_l, mz_u));
        f.getFeatOptions().setIntensityRange(DRange<1>(it_l, it_u));
        f.loadFeatures(in, feature_map);


        //-------------------------------------------------------------
        // calculations
        //-------------------------------------------------------------

        //copy all properties
        FeatureMap map_sm = feature_map;
        //.. but delete feature information
        map_sm.clear(false);

        // only keep charge ch_l:ch_u   (WARNING: feature files without charge information have charge=0, see Ctor of KERNEL/Feature.h)
        for (Feature& fm : feature_map)
        {
          bool const rt_ok = f.getFeatOptions().getRTRange().encloses(DPosition<1>(fm.getRT()));
          bool const mz_ok = f.getFeatOptions().getMZRange().encloses(DPosition<1>(fm.getMZ()));
          bool const int_ok = f.getFeatOptions().getIntensityRange().encloses(DPosition<1>(fm.getIntensity()));
          bool const charge_ok = ((charge_l <= fm.getCharge()) && (fm.getCharge() <= charge_u));
          bool const size_ok = ((size_l <= fm.getSubordinates().size()) && (fm.getSubordinates().size() <= size_u));
          bool const q_ok = ((q_l <= fm.getOverallQuality()) && (fm.getOverallQuality() <= q_u));
          if (remove_hulls) fm.getConvexHulls().clear();

          if (rt_ok && mz_ok && int_ok && charge_ok && size_ok && q_ok)
          {
            if (remove_meta_enabled)
            {
              meta_ok = checkMetaOk(fm, meta_info);
            }
            bool const annotation_ok = checkPeptideIdentification_(fm, remove_annotated_features, remove_unannotated_features, sequences, sequence_comparison_method, accessions, keep_best_score_id, remove_clashes);
            if (annotation_ok && meta_ok) map_sm.push_back(fm);
          }
        }
        //delete unassigned PeptideIdentifications
        if (remove_unassigned_ids)
        {
          map_sm.getUnassignedPeptideIdentifications().clear();
        }
        //update minimum and maximum position/intensity
        map_sm.updateRanges();

        // sort if desired
        if (sort)
        {
          map_sm.sortByPosition();
        }

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------

        //annotate output with data processing info
        addDataProcessing_(map_sm, getProcessingInfo_(DataProcessing::FILTERING));

        f.storeFeatures(out, map_sm, {FileTypes::FEATUREXML});
      }
      else if (in_type == FileTypes::CONSENSUSXML)
      {
        //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------

        ConsensusMap consensus_map;
        FileHandler f;
        //f.setLogType(log_type_);
        f.getOptions().setRTRange(DRange<1>(rt_l, rt_u));
        f.getOptions().setMZRange(DRange<1>(mz_l, mz_u));
        f.getOptions().setIntensityRange(DRange<1>(it_l, it_u));
        f.loadConsensusFeatures(in, consensus_map);

        //-------------------------------------------------------------
        // calculations
        //-------------------------------------------------------------

        // copy all properties
        ConsensusMap consensus_map_filtered = consensus_map;
        //.. but delete feature information
        consensus_map_filtered.resize(0);

        for (ConsensusFeature& cm : consensus_map)
        {
          const bool charge_ok = ((charge_l <= cm.getCharge()) && (cm.getCharge() <= charge_u));
          const bool size_ok = ((cm.size() >= size_l) && (cm.size() <= size_u));

          if (charge_ok && size_ok)
          {
            // this is expensive, so evaluate after everything else passes the test
            if (remove_meta_enabled)
            {
              meta_ok = checkMetaOk(cm, meta_info);
            }
            const bool annotation_ok = checkPeptideIdentification_(cm, remove_annotated_features, remove_unannotated_features, sequences, sequence_comparison_method, accessions, keep_best_score_id, remove_clashes);
            if (annotation_ok && meta_ok) consensus_map_filtered.push_back(cm);
          }
        }
        //delete unassigned PeptideIdentifications
        if (remove_unassigned_ids)
        {
          consensus_map_filtered.getUnassignedPeptideIdentifications().clear();
        }
        //update minimum and maximum position/intensity
        consensus_map_filtered.updateRanges();

        // sort if desired
        if (sort)
        {
          consensus_map_filtered.sortByPosition();
        }

        if (out_type == FileTypes::FEATUREXML)
        {
          if (maps.size() == 1) // When extracting a feature map from a consensus map, only one map ID should be specified. Hence 'maps' should contain only one integer.
          {
            FeatureMap feature_map_filtered;
            FileHandler ff;

            for (ConsensusMap::Iterator cm_it = consensus_map_filtered.begin(); cm_it != consensus_map_filtered.end(); ++cm_it)
            {

              for (ConsensusFeature::HandleSetType::const_iterator fh_iter = cm_it->getFeatures().begin(); fh_iter != cm_it->getFeatures().end(); ++fh_iter)
              {
                if ((int)fh_iter->getMapIndex() == maps[0])
                {
                  Feature feature;
                  feature.setRT(fh_iter->getRT());
                  feature.setMZ(fh_iter->getMZ());
                  feature.setIntensity(fh_iter->getIntensity());
                  feature.setCharge(fh_iter->getCharge());
                  feature_map_filtered.push_back(feature);
                }
              }
            }

            //-------------------------------------------------------------
            // writing output
            //-------------------------------------------------------------

            //annotate output with data processing info
            addDataProcessing_(feature_map_filtered, getProcessingInfo_(DataProcessing::FILTERING));

            feature_map_filtered.applyMemberFunction(&UniqueIdInterface::setUniqueId);

            ff.storeFeatures(out, feature_map_filtered, {FileTypes::FEATUREXML});
          }
          else
          {
            writeLogError_("Error: When extracting a feature map from a consensus map, only one map ID should be specified. The 'map' parameter contains more than one. Aborting!");
            printUsage_();
            return ILLEGAL_PARAMETERS;
          }
        }
        else if (out_type == FileTypes::CONSENSUSXML)
        {
          // generate new consensuses with features that appear in the 'maps' list
          ConsensusMap cm_new; // new consensus map

          for (IntList::iterator map_it = maps.begin(); map_it != maps.end(); ++map_it)
          {
            cm_new.getColumnHeaders()[*map_it].filename = consensus_map_filtered.getColumnHeaders()[*map_it].filename;
            cm_new.getColumnHeaders()[*map_it].size = consensus_map_filtered.getColumnHeaders()[*map_it].size;
            cm_new.getColumnHeaders()[*map_it].unique_id = consensus_map_filtered.getColumnHeaders()[*map_it].unique_id;
          }

          cm_new.setProteinIdentifications(consensus_map_filtered.getProteinIdentifications());

          const bool and_connective = getFlag_("consensus:map_and");
          for (ConsensusFeature& cm : consensus_map_filtered) // iterate over consensuses in the original consensus map
          {
            ConsensusFeature consensus_feature_new(cm); // new consensus feature
            consensus_feature_new.clear();

            ConsensusFeature::HandleSetType::const_iterator fh_it = cm.getFeatures().begin();
            ConsensusFeature::HandleSetType::const_iterator fh_it_end = cm.getFeatures().end();
            for (; fh_it != fh_it_end; ++fh_it) // iterate over features in consensus
            {
              if (ListUtils::contains(maps, fh_it->getMapIndex()))
              {
                consensus_feature_new.insert(*fh_it);
              }
            }

            if ((!consensus_feature_new.empty() && !and_connective) || (consensus_feature_new.size() == maps.size() && and_connective)) // add the consensus to the consensus map only if it is non-empty
            {
              consensus_feature_new.computeConsensus(); // evaluate position of the consensus
              cm_new.push_back(consensus_feature_new);
            }
          }

          // assign unique ids
          cm_new.applyMemberFunction(&UniqueIdInterface::setUniqueId);

          //-------------------------------------------------------------
          // writing output
          //-------------------------------------------------------------

          if (maps.empty())
          {
            //annotate output with data processing info
            addDataProcessing_(consensus_map_filtered, getProcessingInfo_(DataProcessing::FILTERING));

            f.storeConsensusFeatures(out, consensus_map_filtered, {FileTypes::CONSENSUSXML});
          }
          else
          {
            //annotate output with data processing info
            addDataProcessing_(cm_new, getProcessingInfo_(DataProcessing::FILTERING));

            f.storeConsensusFeatures(out, cm_new, {FileTypes::CONSENSUSXML});
          }
        }
      }
      else
      {
        writeLogError_("Error: Unknown output file type given. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
    }
    else
    {
      writeLogError_("Error: Unknown input file type given. Aborting!");
      printUsage_();
      return INCOMPATIBLE_INPUT_DATA;
    }

    return EXECUTION_OK;
  }

  ExitCodes filterByBlackList(MapType& exp, const String& id_blacklist, bool blacklist_imperfect, double rt_tol, double mz_tol)
  {
    vector<ProteinIdentification> protein_ids;
    vector<PeptideIdentification> peptide_ids;
    FileHandler().loadIdentifications(id_blacklist, protein_ids, peptide_ids);

    // translate idXML entries into something more handy
    typedef std::vector<Peak2D> IdType;
    IdType ids; // use Peak2D since it has sorting operators already
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      if (!(peptide_ids[i].hasRT() && peptide_ids[i].hasMZ()))
      {
        OPENMS_LOG_ERROR << "Identifications given in 'id:blacklist' are missing RT and/or MZ coordinates. Cannot do blacklisting without. Quitting." << std::endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
      Peak2D p;
      p.setRT(peptide_ids[i].getRT());
      p.setMZ(peptide_ids[i].getMZ());
      ids.push_back(p);
    }

    std::sort(ids.begin(), ids.end(), Peak2D::RTLess());

    set<Size> blacklist_idx;
    set<Size> ids_covered;
    for (Size i = 0; i != exp.size(); ++i)
    {
      if (exp[i].getMSLevel() == 2)
      {
        if (!exp[i].getPrecursors().empty())
        {
          double pc_rt = exp[i].getRT();
          double pc_mz = exp[i].getPrecursors()[0].getMZ();

          IdType::iterator p_low = std::lower_bound(ids.begin(), ids.end(), pc_rt - rt_tol, Peak2D::RTLess());
          IdType::iterator p_high = std::lower_bound(ids.begin(), ids.end(), pc_rt + rt_tol, Peak2D::RTLess());

          // if precursor is out of the whole range, then p_low==p_high == (begin()||end())
          // , thus the following loop will not run
          for (IdType::iterator id_it = p_low; id_it != p_high; ++id_it) // RT already checked.. now check m/z
          {
            if (pc_mz - mz_tol < id_it->getMZ() && id_it->getMZ() < pc_mz + mz_tol)
            {
              blacklist_idx.insert(i);
              ids_covered.insert(std::distance(ids.begin(), id_it));
              // no break, since we might cover more IDs here
            }
          }
        }
      }
    }

    OPENMS_LOG_INFO << "Removing " << blacklist_idx.size() << " MS2 spectra." << endl;
    if (ids_covered.size() != ids.size())
    {
      if (!blacklist_imperfect)
      {
        OPENMS_LOG_ERROR << "Covered only " << ids_covered.size() << "/" << ids.size() << " IDs. Check if your input files (raw + ids) match and if your tolerances ('rt' and 'mz') are set properly.\n"
                  << "If you are sure unmatched ids are ok, set the 'id:blacklist_imperfect' flag!" << std::endl;
        return UNEXPECTED_RESULT;
      }
      else
      {
        OPENMS_LOG_WARN << "Covered only " << ids_covered.size() << "/" << ids.size() << " IDs. Check if your input files (raw + ids) match and if your tolerances ('rt' and 'mz') are set properly.\n"
                 << "Remove the 'id:blacklist_imperfect' flag of you want this to be an error!" << std::endl;
      }
    }


    PeakMap exp2 = exp;
    exp2.clear(false);

    for (Size i = 0; i != exp.size(); ++i)
    {
      if (blacklist_idx.find(i) ==
          blacklist_idx.end())
      {
        exp2.addSpectrum(exp[i]);
      }
    }

    exp = exp2;
    return EXECUTION_OK;
  }

  ExitCodes filterByBlackOrWhiteList(bool is_blacklist, MapType& exp, const ConsensusMap& consensus_map, double rt_tol, double mz_tol, bool unit_ppm, std::set<UInt64> map_ids)
  {
    std::vector<Peak2D> feature_pos;
    // if map_id are specified, only use these for blacklisting
    for (const ConsensusFeature& cm : consensus_map)
    {
      for (const FeatureHandle& fh : cm)
      {
        UInt64 map_index = fh.getMapIndex();
        if (map_ids.empty() || map_ids.find(map_index) != map_ids.end())
        {
          Peak2D p;
          p.setMZ(fh.getMZ());
          p.setRT(fh.getRT());
          feature_pos.push_back(p);
        }
      }
    }

    // sort by rt to use binary search
    std::sort(feature_pos.begin(), feature_pos.end(), Peak2D::RTLess());
    set<Size> list_idx;
    for (Size i = 0; i != exp.size(); ++i)
    {
      if (exp[i].getMSLevel() == 2)
      {
        if (!exp[i].getPrecursors().empty())
        {
          double pc_mz = exp[i].getPrecursors()[0].getMZ();
          double pc_rt = exp[i].getRT(); // use rt of MS2

          std::vector<Peak2D>::iterator p_low = std::lower_bound(feature_pos.begin(), feature_pos.end(), pc_rt - rt_tol, Peak2D::RTLess());
          std::vector<Peak2D>::iterator p_high = std::lower_bound(feature_pos.begin(), feature_pos.end(), pc_rt + rt_tol, Peak2D::RTLess());

          double mz_tol_da = unit_ppm ? pc_mz * 1e-6 * mz_tol : mz_tol;

          // if precursor is out of the whole range, then p_low==p_high == (begin()||end())
          // , thus the following loop will not run
          for (std::vector<Peak2D>::iterator f_it = p_low; f_it != p_high; ++f_it) // RT already checked.. now check m/z
          {
            if (pc_mz - mz_tol_da < f_it->getMZ() && f_it->getMZ() < pc_mz + mz_tol_da)
            {
              list_idx.insert(i);
              // no break, since we might cover more features here
            }
          }
        }
      }
    }

    // create new experiment
    PeakMap exp2;
    exp2.getExperimentalSettings() = (ExperimentalSettings)exp.getExperimentalSettings(); // copy meta data

    for (Size i = 0; i != exp.size(); ++i)
    {
      // don't need to sort list as it is increasing
      if (is_blacklist)
      {
        // blacklist: add all spectra not contained in list
        if (list_idx.find(i) == list_idx.end())
        {
          exp2.addSpectrum(exp[i]);
        }
      }
      else   // whitelist: add all non MS2 spectra, and MS2 only if in list
      {
        if (exp[i].getMSLevel() != 2 || list_idx.find(i) != list_idx.end())
        {
          exp2.addSpectrum(exp[i]);
        }
      }
    }

    exp = exp2;
    return EXECUTION_OK;
  }

  ExitCodes filterByBlackOrWhiteList(bool is_blacklist, PeakMap& exp, const PeakMap& lib_file, double rt_tol, double mz_tol, double sim_tol,  bool unit_ppm)
  {
    const bool enable_mz_check = (mz_tol >= 0);
    const bool enable_rt_check = (rt_tol >= 0);
    const bool enable_sim_check = (sim_tol > -1);

    auto comp_function= std::unique_ptr<PeakSpectrumCompareFunctor>(new (ZhangSimilarityScore));

    set<Size> list_idx;

    for (auto const & lib_spectrum : lib_file)
    {
      if (!lib_spectrum.getPrecursors().empty())
      {
        // extract precursor positions from query file
        double lib_mz = lib_spectrum.getPrecursors()[0].getMZ();
        double lib_rt = lib_spectrum.getRT();

        // look-up matching spectra in input file (TODO: use KD-tree)
        int exp_index = -1;
        for (auto const & exp_spectrum : exp)
        {
          // keep track of current spectrum index
          ++exp_index;

          // TODO: extend to other MS levels and multiple precursors
          if (exp_spectrum.getMSLevel() != 2 || exp_spectrum.getPrecursors().empty())
          { 
            continue;
          }

          // skip if m/z's don't match
          const double pc_mz = exp_spectrum.getPrecursors()[0].getMZ();
          const double mz_tol_da = unit_ppm ? pc_mz * 1e-6 * mz_tol : mz_tol;
          if (enable_mz_check && fabs(pc_mz - lib_mz) > mz_tol_da)
          {
            continue; 
          }

          // skip if rt's don't match
          const double pc_rt = exp_spectrum.getRT();
          if (enable_rt_check && fabs(pc_rt - lib_rt) > rt_tol)
          {
            continue;
          }

          // skip if not similar enough
          if (enable_sim_check && (*comp_function)(exp_spectrum, lib_spectrum) < sim_tol)
          {
            continue;
          }

          writeDebug_("Similarity score: " + String((*comp_function)(exp_spectrum, lib_spectrum)), 10);

          // we have matching spectra
          list_idx.insert(exp_index); 
        }
      }
    }

    // create new experiment
    PeakMap exp2 = exp; // copy meta data
    exp2.clear(false); // clear spectra

    for (Size i = 0; i != exp.size(); ++i)
    {
      // don't need to sort list as it is increasing
      if (is_blacklist)
      {
        // blacklist: add all spectra not contained in list
        if (list_idx.find(i) == list_idx.end())
        {
          exp2.addSpectrum(exp[i]);
        }
      }
      else   // whitelist: add all non-MS2 spectra + matched MS2 spectra
      {
        if (exp[i].getMSLevel() != 2 || list_idx.find(i) != list_idx.end())
        {
          exp2.addSpectrum(exp[i]);
        }
      }
    }

    exp = exp2;
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFileFilter tool;
  return tool.main(argc, argv);
}

/// @endcond
