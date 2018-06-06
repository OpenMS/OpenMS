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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>

#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>


#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_ProteomicsLFQ
 **/

class UTILProteomicsLFQ :
  public TOPPBase
{
public:
  UTILProteomicsLFQ() :
    TOPPBase("ProteomicsLFQ", "A standard proteomics LFQ pipeline.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file list>", StringList(), "input files");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFileList_("ids", "<file list>", StringList(), "unfiltered identifications");
    setValidFormats_("ids", ListUtils::create<String>("idXML,mzId"));
    registerInputFile_("design", "<file>", "", "design file", false);
    setValidFormats_("design", ListUtils::create<String>("tsv"));

    registerOutputFile_("out", "<file>", "", "output mzTab file");
    setValidFormats_("out", ListUtils::create<String>("tsv")); // TODO: add file extension for mzTab

    /// TODO: think about export of quality control files (qcML?)

    Param pp_defaults = PeakPickerHiRes().getDefaults();
    Param ff_defaults = FeatureFinderIdentificationAlgorithm().getDefaults();
    Param ma_defaults = MapAlignmentAlgorithmIdentification().getDefaults();
    Param fl_defaults = FeatureGroupingAlgorithmKD().getDefaults();
    Param pep_defaults = Math::PosteriorErrorProbabilityModel().getParameters();
    //Param pi_defaults = ProteinInferenceAlgorithmXX().getDefaults();
    Param pq_defaults = PeptideAndProteinQuant().getDefaults();

    Param combined;
    combined.insert("Centroiding:", pp_defaults);
    combined.insert("Peptide Quantification:", ff_defaults);
    combined.insert("Alignment:", ma_defaults);
    combined.insert("Linking:", fl_defaults);
    // combined.insert("Protein Inference:", pi_defaults);
    combined.insert("Protein Quantification:", pq_defaults);
    combined.insert("Posterior Error Probability:", pep_defaults);

    registerFullParam_(combined);
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parameter handling
    //-------------------------------------------------------------

    // TODO: check handling of single MS file and n-MS files of a single run

    // Read tool parameters
    StringList in = getStringList_("in");
    String out = getStringOption_("out");
    StringList in_ids = getStringList_("ids");
    String design_file = getStringOption_("design");

    if (in.size() != in_ids.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, 
        OPENMS_PRETTY_FUNCTION, "Number of id and spectra files don't match.");
    }
  
    // map between mzML file and corresponding id file
    map<String, String> mzfile2idfile;
    for (Size i = 0; i != in.size(); ++i)
    {
      const String& in_abs_path = File::absolutePath(in[i]);
      const String& id_abs_path = File::absolutePath(in_ids[i]);
      mzfile2idfile[in_abs_path] = id_abs_path;
      LOG_DEBUG << "Spectra: " << in[i] << "\t Ids: " << in_ids[i] << endl;
      if (!File::exists(in[i]))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, 
          OPENMS_PRETTY_FUNCTION, "Spectra file '" + in[i] + "' does not exist.");
      }

      if (!File::exists(in_ids[i]))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, 
          OPENMS_PRETTY_FUNCTION, "Id file '" + in[i] + "' does not exist.");
      }
    }
 
    ExperimentalDesign design;
    if (!design_file.empty())
    { 
      design = ExperimentalDesignFile::load(design_file, false);
    }
    else
    {  // default to unfractionated design
      ExperimentalDesign::MSFileSection msfs;
      Size count{1};
      for (String & s : in)
      {
        ExperimentalDesign::MSFileSectionEntry e;
        e.fraction = 1;
        e.fraction_group = count;
        e.label = 1;
        e.path = s;
        e.sample = count;
      }
    }
    std::map<unsigned int, std::vector<String> > frac2ms = design.getFractionToMSFilesMapping();

    Param pp_param = getParam_().copy("Centroiding:", true);
    writeDebug_("Parameters passed to PeakPickerHiRes algorithm", pp_param, 3);
    PeakPickerHiRes pp;
    pp.setLogType(log_type_);
    pp.setParameters(pp_param);

    Param ff_param = getParam_().copy("Peptide Quantification:", true);
    writeDebug_("Parameters passed to FeatureFinderIdentification algorithm", ff_param, 3);

    Param ma_param = getParam_().copy("Alignment:", true);
    writeDebug_("Parameters passed to MapAlignmentAlgorithmIdentification algorithm", ma_param, 3);

    Param fl_param = getParam_().copy("Linking:", true);
    writeDebug_("Parameters passed to FeatureGroupingAlgorithmKD algorithm", fl_param, 3);

    Param pep_param = getParam_().copy("Posterior Error Probability:", true);
    writeDebug_("Parameters passed to PEP algorithm", pep_param, 3);
    Math::PosteriorErrorProbabilityModel PEP_model;
    // PEP_model.setLogType(log_type_); TODO: add to PEP
    PEP_model.setParameters(pep_param);

    // TODO: inference parameter

    Param pq_param = getParam_().copy("Protein Quantification:", true);
    writeDebug_("Parameters passed to PeptideAndProteinQuant algorithm", pq_param, 3);
    PeptideAndProteinQuant quantifier;
    //quantifier.setLogType(log_type_);
    quantifier.setParameters(pq_param);

    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param com_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to both sub-algorithms (mtd and epd)", com_param, 3);

    //-------------------------------------------------------------
    // Loading input
    //-------------------------------------------------------------
    ConsensusMap consensus;
    for (auto const ms_files : frac2ms) // for each fraction->ms file(s)
    {
      vector<FeatureMap> feature_maps;
      ConsensusMap consensus_fraction;      

      // debug output
      LOG_DEBUG << "Processing fraction number: " << ms_files.first << endl; 
      LOG_DEBUG << "Files: " << endl; 
      
      for (String const & mz_file : ms_files.second) // for each MS file
      {
        LOG_DEBUG << mz_file << endl;          
      }

      //TODO: check if we want to parallelize that
      for (String const & mz_file : ms_files.second) // for each MS file
      {
        // TODO: check if s is part of in 
      
        // load raw file
        MzMLFile mzML_file;
        mzML_file.setLogType(log_type_);

        // setup for fwhm estimation
        String in_option = getStringOption_("in");
        MSExperiment ms_peakmap;
        std::vector<Int> ms_level(1, 1);
        mzML_file.getOptions().setMSLevels(ms_level);
        mzML_file.load(in_option, ms_peakmap);

        MassTraceDetection mt_ext;
        mtd_param.insert("", com_param);
        mtd_param.remove("chrom_fwhm");
        mt_ext.setParameters(mtd_param);
        std::vector<MassTrace> m_traces;
        mt_ext.run(ms_peakmap, m_traces, 1000);

        double median_fwhm;
        std::vector<double> fwhm_1000;

        for (auto &m : m_traces)
        {
          if (m.getSize() == 0) continue;

          m.updateMeanMZ();
          m.updateWeightedMZsd();

          // FWHM
          double fwhm = m.estimateFWHM(true); // smoothed or not?
          fwhm_1000.push_back(fwhm);
        }

        median_fwhm = std::accumulate(fwhm_1000.begin(), fwhm_1000.end(), fwhm_1000.size());

        std::cout << "Calculated median fwhm: " << median_fwhm << std::endl;

        PeakMap ms_centroided;
        {  // create scope for raw data so it is properly freed (Note: clear() is not sufficient)
          PeakMap ms_raw;
          mzML_file.load(mz_file, ms_raw);

          if (ms_raw.empty())
          {
            LOG_WARN << "The given file does not contain any spectra.";
            return INCOMPATIBLE_INPUT_DATA;
          }

          // remove MS2 peak data and check if spectra are sorted
          for (Size i = 0; i < ms_raw.size(); ++i)
          {
            if (ms_raw[i].getMSLevel() == 2)
            {
              ms_raw[i].clear(false);
            }
            if (!ms_raw[i].isSorted())
            {
              ms_raw[i].sortByPosition();
              writeLog_("Info: Sorted peaks by m/z.");
            }
          }

          //-------------------------------------------------------------
          // Centroiding of MS1
          //-------------------------------------------------------------
          pp.pickExperiment(ms_raw, ms_centroided, true);
        }

        // writing picked mzML files for data submission
        // annotate output with data processing info
        // TODO: how to store picked files? by specifying a folder? or by output files that match in number to input files
        // TODO: overwrite primaryMSRun with picked mzML name (for submission)
        // mzML_file.store(OUTPUTFILENAME, ms_centroided);
        // TODO: free all MS2 spectra (to release memory!)

        //-------------------------------------------------------------
        // Posterior Error Probability calculation
        //-------------------------------------------------------------
        
        // TODO: load id file for this MS file

        vector<ProteinIdentification> protein_ids;
        vector<PeptideIdentification> peptide_ids;
        IdXMLFile().load(mzfile2idfile[mz_file], protein_ids, peptide_ids);

        map<String, vector<vector<double> > > all_scores = Math::PosteriorErrorProbabilityModel::extractAndTransformScores(
          protein_ids, 
          peptide_ids, 
          false, 
          true,   
          true,  
          0.05);
                
        for (auto & score : all_scores)  // for all search engine scores (should only be 1)
        {
          vector<String> engine_info;
          score.first.split(',', engine_info);
          String engine = engine_info[0];
          Int charge = (engine_info.size() == 2) ? engine_info[1].toInt() : -1;

          // fit to score vector
          bool return_value = PEP_model.fit(score.second[0]);
          if (!return_value) 
          {
            writeLog_("Unable to fit data. Algorithm did not run through for the following search engine: " + engine);
          }

          bool unable_to_fit_data(true), data_might_not_be_well_fit(true);
          Math::PosteriorErrorProbabilityModel::updateScores(
            PEP_model,
            engine,
            charge,
            true, // prob_correct 
            false, // split_charge
            protein_ids,
            peptide_ids,
            unable_to_fit_data,
            data_might_not_be_well_fit);

          if (unable_to_fit_data)
          {
            writeLog_(String("Unable to fit data for search engine: ") + engine);
          }
          else if (data_might_not_be_well_fit) 
          {
            writeLog_(String("Data might not be well fitted for search engine: ") + engine);
          }
        }        

        //-------------------------------------------------------------
        // Feature detection
        //-------------------------------------------------------------   
        vector<ProteinIdentification> ext_protein_ids;
        vector<PeptideIdentification> ext_peptide_ids;

        // create empty feature map and annotate MS file
        FeatureMap fm;
        StringList sl;
        sl.push_back(mz_file);
        fm.setPrimaryMSRunPath(sl);
        feature_maps.push_back(fm);

        FeatureFinderIdentificationAlgorithm ff;
        ff.getMSData().swap(ms_centroided);     
        ff.getProgressLogger().setLogType(log_type_);
        ff.setParameters(ff_param);
        ff.run(peptide_ids, protein_ids, ext_peptide_ids, ext_protein_ids, feature_maps.back());

        // TODO: think about external ids ;)
        // TODO: free parts of feature map not needed for further processing (e.g., subfeatures...)
      }

      //-------------------------------------------------------------
      // Align all features of this fraction
      //-------------------------------------------------------------
      vector<TransformationDescription> transformations;
   
      //TODO: check if we need to set reference
      Size reference_index(0);
      MapAlignmentAlgorithmIdentification aligner;
      aligner.setLogType(log_type_);
      aligner.setParameters(ma_param);
      aligner.align(feature_maps, transformations, reference_index);


      // find model parameters:
      Param model_params = TOPPMapAlignerBase::getModelDefaults("b_spline");
      String model_type = model_params.getValue("type");
      model_params = model_params.copy(model_type + ":", true);
      for (TransformationDescription & t : transformations)
      {
        LOG_DEBUG << "Using " << t.getDataPoints().size() << " points in fit." << endl; 
        if (t.getDataPoints().size() > 10)
        {
          t.fitModel(model_type, model_params);
        }
      }

      // Apply transformations
      for (Size i = 0; i < feature_maps.size(); ++i)
      {
        MapAlignmentTransformer::transformRetentionTimes(feature_maps[i],
          transformations[i]);
      }

      //-------------------------------------------------------------
      // Link all features of this fraction
      //-------------------------------------------------------------
      LOG_DEBUG << "Linking: " << feature_maps.size() << " features" << endl;
      FeatureGroupingAlgorithmKD linker;
      linker.setLogType(log_type_);
      linker.setParameters(fl_param);
      linker.group(feature_maps, consensus_fraction);
      ConsensusXMLFile().store("debug_fraction_" + String(ms_files.first) +  ".consensusXML", consensus_fraction);
      LOG_DEBUG << "to produce a consensus map with " << consensus_fraction.getColumnHeaders().size() << " columns." << endl;

      //-------------------------------------------------------------
      // ID conflict resolution
      //-------------------------------------------------------------
      IDConflictResolverAlgorithm::resolve(consensus_fraction);

      //-------------------------------------------------------------
      // ConsensusMap normalization
      //-------------------------------------------------------------
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(consensus_fraction, 
        ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, 
        "", 
        "");

      // append consensus map calculated for this fraction number
      consensus.appendColumns(consensus_fraction);

      // end of scope of fraction related data
    }

    // TODO: FileMerger merge ids (here? or already earlier? filtered?)
    // TODO: check if it makes sense to integrate SVT imputation algorithm (branch)

    //-------------------------------------------------------------
    // Protein inference
    //-------------------------------------------------------------
    // TODO: ProteinInference on merged ids (how merged?)
    // TODO: Output coverage on protein and group level
    ProteinIdentification infered_protein_groups;
    vector<PeptideIdentification> infered_peptides;

    // TODO: maybe check if some consensus ID algorithms are applicable

    //-------------------------------------------------------------
    // Protein quantification and export to mzTab
    //-------------------------------------------------------------
    // TODO: ProteinQuantifier on (merged?) consensusXML (with 1% FDR?) + inference ids (unfiltered?)? 
    if (infered_protein_groups.getIndistinguishableProteins().empty())
    {
      throw Exception::MissingInformation(
       __FILE__, 
       __LINE__, 
       OPENMS_PRETTY_FUNCTION, 
       "No information on indistinguishable protein groups found.");
    }

    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides(infered_peptides); 
    quantifier.quantifyProteins(infered_protein_groups);

    //-------------------------------------------------------------
    // Export of MzTab file as final output
    //-------------------------------------------------------------

    // Annotate quants to protein(groups) for easier export in mzTab
    auto const & protein_quants = quantifier.getProteinResults();
    PeptideAndProteinQuant::annotateQuantificationsToProteins(protein_quants, infered_protein_groups);
    vector<ProteinIdentification>& proteins = consensus.getProteinIdentifications();
    proteins.insert(proteins.begin(), infered_protein_groups); // insert inference information as first protein identification

    // Fill MzTab with meta data and quants annotated in identification data structure
    const bool report_unmapped(true);
    const bool report_unidentified_features(false);
    MzTab m = MzTab::exportConsensusMapToMzTab(consensus, String("null"), report_unidentified_features, report_unmapped);
    MzTabFile().store(out, m);    

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  UTILProteomicsLFQ tool;
  return tool.main(argc, argv);
}

/// @endcond

