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
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>
#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>


#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

#include <OpenMS/KERNEL/ConversionHelper.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>

#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_ProteomicsLFQ

  TODO: add test for

  ./bin/ProteomicsLFQ \
-in ../share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML ../share/OpenMS/examples/FRACTIONS/BSA1_F2.mzML ../share/OpenMS/examples/FRACTIONS/BSA2_F1.mzML ../share/OpenMS/examples/FRACTIONS/BSA2_F2.mzML ../share/OpenMS/examples/FRACTIONS/BSA3_F1.mzML ../share/OpenMS/examples/FRACTIONS/BSA3_F2.mzML \
-ids ../share/OpenMS/examples/FRACTIONS/BSA1_F1.idXML ../share/OpenMS/examples/FRACTIONS/BSA1_F2.idXML ../share/OpenMS/examples/FRACTIONS/BSA2_F1.idXML ../share/OpenMS/examples/FRACTIONS/BSA2_F2.idXML ../share/OpenMS/examples/FRACTIONS/BSA3_F1.idXML ../share/OpenMS/examples/FRACTIONS/BSA3_F2.idXML \
-design ../share/OpenMS/examples/FRACTIONS/BSA_design.tsv \
-Alignment:max_rt_shift 0 \
-fasta ../share/OpenMS/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta \
-out BSA.mzTab -threads 4 -debug 667 


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
    registerInputFileList_("in", "<file list>", StringList(), "Input files");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFileList_("ids", "<file list>", StringList(), 
      "Identifications filtered at PSM level (e.g., q-vaue < 0.01)."
      "And annotated with PEP as main score.\n"
      "We suggest using:\n"
      "1. PercolatorAdapter tool (score_type = 'q-value', -post-processing-tdc)\n"
      "2. FalseDiscoveryRate (FDR:PSM = 0.01)\n"
      "3. IDScoreSwitcher (-old_score q-value -new_score MS:1001493 -new_score_orientation lower_better -new_score_type)\n"
      "To obtain well calibrated PEPs and an inital reduction of PSMs\n"
      "ID files must be provided in same order as spectra files.");
    setValidFormats_("ids", ListUtils::create<String>("idXML,mzId"));

    registerInputFile_("design", "<file>", "", "design file", false);
    setValidFormats_("design", ListUtils::create<String>("tsv"));

    registerInputFile_("fasta", "<file>", "", "fasta file", false);
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out", "<file>", "", "output mzTab file");
    setValidFormats_("out", ListUtils::create<String>("mzTab"));

    registerStringOption_("targeted_only", "<option>", "false", "Only ID based quantification.", false, true);
    setValidStrings_("targeted_only", ListUtils::create<String>("true,false"));

    registerStringOption_("mass_recalibration", "<option>", "true", "Mass recalibration.", false, true);
    setValidStrings_("mass_recalibration", ListUtils::create<String>("true,false"));


    /// TODO: think about export of quality control files (qcML?)

    Param pp_defaults = PeakPickerHiRes().getDefaults();
    Param ffm_defaults = FeatureFindingMetabo().getDefaults();
    Param ffi_defaults = FeatureFinderIdentificationAlgorithm().getDefaults();
    Param ma_defaults = MapAlignmentAlgorithmIdentification().getDefaults();
//    Param fl_defaults = FeatureGroupingAlgorithmKD().getDefaults();
    Param fl_defaults = FeatureGroupingAlgorithmQT().getDefaults();
    //Param pi_defaults = ProteinInferenceAlgorithmXX().getDefaults();
    Param pq_defaults = PeptideAndProteinQuant().getDefaults();

    Param combined;
    combined.insert("Centroiding:", pp_defaults);
    combined.insert("Seeds:", ffm_defaults);
    combined.insert("PeptideQuantification:", ffi_defaults);
    combined.insert("Alignment:", ma_defaults);
    combined.insert("Linking:", fl_defaults);
    // combined.insert("Protein Inference:", pi_defaults);
    combined.insert("ProteinQuantification:", pq_defaults);

    registerFullParam_(combined);
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parameter handling
    //-------------------------------------------------------------    

    // Read tool parameters
    StringList in = getStringList_("in");
    String out = getStringOption_("out");
    StringList in_ids = getStringList_("ids");
    String design_file = getStringOption_("design");
    String in_db = getStringOption_("fasta");

    // Validate parameters
    if (in.size() != in_ids.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, 
        OPENMS_PRETTY_FUNCTION, "Number of id and spectra files don't match.");
    }

    //-------------------------------------------------------------
    // Experimental design: read or generate default
    //-------------------------------------------------------------      
    ExperimentalDesign design;
    if (!design_file.empty())
    { // load from file
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
        msfs.push_back(e);
      }      
      design.setMSFileSection(msfs);
    }
    std::map<unsigned int, std::vector<String> > frac2ms = design.getFractionToMSFilesMapping();

    for (auto & f : frac2ms)
    {
      writeDebug_("Fraction " + String(f.first) + ":", 10);
      for (const auto & s : f.second)
      {
        writeDebug_("MS file: " + String(s), 10);
      }
    }

    // Map between mzML file and corresponding id file
    // Here we currently assume that these are provided in the exact same order.
    // In the future, we could warn or reorder them based on the annotated primaryMSRunPath in the ID file.
    map<String, String> mzfile2idfile;
    for (Size i = 0; i != in.size(); ++i)
    {
      const String& in_abs_path = File::absolutePath(in[i]);
      const String& id_abs_path = File::absolutePath(in_ids[i]);
      mzfile2idfile[in_abs_path] = id_abs_path;      
      writeDebug_("Spectra: " + in[i] + "\t Ids: " + in_ids[i],  1);
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

    Param pp_param = getParam_().copy("Centroiding:", true);
    writeDebug_("Parameters passed to PeakPickerHiRes algorithm", pp_param, 3);


    Param ma_param = getParam_().copy("Alignment:", true);
    writeDebug_("Parameters passed to MapAlignmentAlgorithmIdentification algorithm", ma_param, 3);

    Param fl_param = getParam_().copy("Linking:", true);
    writeDebug_("Parameters passed to FeatureGroupingAlgorithmQT algorithm", fl_param, 3);

    Param pep_param = getParam_().copy("Posterior Error Probability:", true);
    writeDebug_("Parameters passed to PEP algorithm", pep_param, 3);

    // TODO: inference parameter

    Param pq_param = getParam_().copy("ProteinQuantification:", true);
    writeDebug_("Parameters passed to PeptideAndProteinQuant algorithm", pq_param, 3);


    Param com_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to both sub-algorithms (mtd and epd)", com_param, 3);

    //-------------------------------------------------------------
    // Loading input
    //-------------------------------------------------------------
    ConsensusMap consensus;
    double median_fwhm(0);

    for (auto const &ms_files : frac2ms) // for each fraction->ms file(s)
    {
      vector<FeatureMap> feature_maps;
      ConsensusMap consensus_fraction;
      const Size fraction = ms_files.first;

      // debug output
      writeDebug_("Processing fraction number: " + String(fraction) + "\nFiles: ",  1);
      for (String const & mz_file : ms_files.second) { writeDebug_(mz_file,  1); }

      // for each MS file
      Size fraction_group{1};
      for (String const & mz_file : ms_files.second)
      {     
        PeakMap ms_centroided;
        { // create scope for raw data so it is properly freed (Note: clear() is not sufficient)
          // load raw file
          MzMLFile mzML_file;
          mzML_file.setLogType(log_type_);

          PeakMap ms_raw;
          mzML_file.load(mz_file, ms_raw);
          ms_raw.clearMetaDataArrays();

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
              ms_raw[i].clear(false);  // delete MS2 peaks
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
          PeakPickerHiRes pp;
          pp.setLogType(log_type_);
          pp.setParameters(pp_param);
          pp.pickExperiment(ms_raw, ms_centroided, true);
           
          //-------------------------------------------------------------
          // HighRes Precursor Mass Correction
          //-------------------------------------------------------------
          std::vector<double> deltaMZs, mzs, rts;
          std::set<Size> corrected_to_highest_intensity_peak = PrecursorCorrection::correctToHighestIntensityMS1Peak(
            ms_centroided, 
            0.01, // check if we can estimate this from data (here it is given in m/z not ppm)
            deltaMZs, 
            mzs, 
            rts
            );      
          writeLog_("Info: Corrected " + String(corrected_to_highest_intensity_peak.size()) + " precursors.");
          if (!deltaMZs.empty())
          {
            vector<double> deltaMZs_ppm, deltaMZs_ppmabs;
            for (Size i = 0; i != deltaMZs.size(); ++i)
            {
              deltaMZs_ppm.push_back(Math::getPPM(mzs[i], mzs[i] + deltaMZs[i]));
              deltaMZs_ppmabs.push_back(Math::getPPMAbs(mzs[i], mzs[i] + deltaMZs[i]));
            }

            double median = Math::median(deltaMZs_ppm.begin(), deltaMZs_ppm.end());
            double MAD =  Math::MAD(deltaMZs_ppm.begin(), deltaMZs_ppm.end(), median);
            double median_abs = Math::median(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end());
            double MAD_abs = Math::MAD(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end(), median_abs);
            writeLog_("Precursor correction:\n  median = " + String(median) + " ppm  MAD = " + String(MAD)
                     +"\n  median (abs.) = " + String(median_abs) + " ppm  MAD = " + String(MAD_abs));
          }
        }

        // writing picked mzML files for data submission
        // annotate output with data processing info
        // TODO: how to store picked files? by specifying a folder? or by output files that match in number to input files
        // TODO: overwrite primaryMSRun with picked mzML name (for submission)
        // mzML_file.store(OUTPUTFILENAME, ms_centroided);
        // TODO: free all MS2 spectra (to release memory!)

        vector<ProteinIdentification> protein_ids;
        vector<PeptideIdentification> peptide_ids;
        const String& mz_file_abs_path = File::absolutePath(mz_file);
        const String& id_file_abs_path = File::absolutePath(mzfile2idfile[mz_file_abs_path]);
        IdXMLFile().load(id_file_abs_path, protein_ids, peptide_ids);

        if (protein_ids.size() != 1)
        {
          LOG_FATAL_ERROR << "Exactly one protein identification runs must be annotated in " << id_file_abs_path << endl;
          return ExitCodes::INCOMPATIBLE_INPUT_DATA;     
        }

        // annotate experimental design
        StringList id_msfile_ref;
        protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);
        if (id_msfile_ref.empty())
        {
          LOG_DEBUG << "MS run path not set in ID file." << endl;
        }
        else
        {
          // TODO: we could add a check (e.g., matching base name) here
          id_msfile_ref.clear();
        }                
        id_msfile_ref.push_back(mz_file);
        protein_ids[0].setPrimaryMSRunPath(id_msfile_ref);
        protein_ids[0].setMetaValue("fraction_group", fraction_group);
        protein_ids[0].setMetaValue("fraction", fraction);

        // update identifers to make them unique 
        // fixes some bugs related to users splitting the original mzML and id files before running the analysis
        // in that case these files might have the same identifier
        const String old_identifier = protein_ids[0].getIdentifier();
        const String new_identifier = old_identifier + "_" + String(fraction_group) + "F" + String(fraction);
        protein_ids[0].setIdentifier(new_identifier);
        for (auto & p : peptide_ids)
        {
          if (p.getIdentifier() == old_identifier)
          {
            p.setIdentifier(new_identifier);
          }
          else
          {
            LOG_WARN << "Peptide ID identifier found not present in the protein ID" << endl;
          }
        }

        // reannotate spectrum references
        SpectrumMetaDataLookup::addMissingSpectrumReferences(
          peptide_ids, 
          mz_file_abs_path,
          true);

        //-------------------------------------------------------------
        // Internal Calibration of spectra peaks and precursor peaks with high-confidence IDs
        // TODO: check if this improves targeted extraction
        //-------------------------------------------------------------
        if (getStringOption_("mass_recalibration") == "true")
        {
          InternalCalibration ic;
          ic.setLogType(log_type_);
          ic.fillCalibrants(peptide_ids, 25.0); // >25 ppm maximum deviation defines an outlier TODO: check if we need to adapt this
          MZTrafoModel::MODELTYPE md = MZTrafoModel::QUADRATIC; // TODO: check if it makes sense to choose the quadratic model
          bool use_RANSAC = (md == MZTrafoModel::LINEAR || md == MZTrafoModel::QUADRATIC);
          Size RANSAC_initial_points = (md == MZTrafoModel::LINEAR) ? 2 : 3;
          Math::RANSACParam p(RANSAC_initial_points, 70, 10, 30, true); // TODO: check defaults (taken from tool)
          MZTrafoModel::setRANSACParams(p);
          // these limits are a little loose, but should prevent grossly wrong models without burdening the user with yet another parameter.
          MZTrafoModel::setCoefficientLimits(25.0, 25.0, 0.5); 

          IntList ms_level = {1};
          double rt_chunk = 300.0; // 5 minutes
          String qc_residual_path, qc_residual_png_path;
          if (debug_level_ >= 1)
          {
            const String & id_basename = File::basename(id_file_abs_path);
            qc_residual_path = id_basename + "qc_residuals.tsv";
            qc_residual_png_path = id_basename + "qc_residuals.png";
          } 

          if (!ic.calibrate(ms_centroided, ms_level, md, rt_chunk, use_RANSAC, 
                        10.0,
                        5.0, 
                        "",                      
                        "",
                        qc_residual_path,
                        qc_residual_png_path,
                        "Rscript"))
          {
            LOG_WARN << "\nCalibration failed. See error message above!" << std::endl;          
          }
        }
        

        //////////////////////////////////////////
        // Chromatographic parameter estimation
        //////////////////////////////////////////
        MassTraceDetection mt_ext;
        Param mtd_param = mt_ext.getParameters();
        writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

        std::vector<MassTrace> m_traces;
        mt_ext.run(ms_centroided, m_traces, 1000);

        std::vector<double> fwhm_1000;
        for (auto &m : m_traces)
        {
          if (m.getSize() == 0) continue;
          m.updateMeanMZ();
          m.updateWeightedMZsd();
          double fwhm = m.estimateFWHM(false);
          fwhm_1000.push_back(fwhm);
        }

        median_fwhm = Math::median(fwhm_1000.begin(), fwhm_1000.end());

        LOG_INFO << "Median FWHM: " << median_fwhm << std::endl;

        //-------------------------------------------------------------
        // Feature detection
        //-------------------------------------------------------------   
        vector<ProteinIdentification> ext_protein_ids;
        vector<PeptideIdentification> ext_peptide_ids;

 
        ///////////////////////////////////////////////
        // Run MTD before FFM

        // create empty feature map and annotate MS file
        FeatureMap seeds;

        StringList sl;
        sl.push_back(mz_file);
        seeds.setPrimaryMSRunPath(sl);

        if (getStringOption_("targeted_only") == "false")
        {
          std::vector<MassTrace> m_traces_full;
          mt_ext.run(ms_centroided, m_traces_full);

          std::vector<MassTrace> splitted_mtraces;
          ElutionPeakDetection epdet;
          epdet.detectPeaks(m_traces_full, splitted_mtraces);

          FeatureFindingMetabo ffm;
          Param ffm_param = getParam_().copy("Seeds:", true);
          ffm_param.setValue("mz_scoring_13C", "true");
          ffm_param.setValue("isotope_filtering_model", "peptides");
          ffm_param.setValue("remove_single_traces", "true");
          ffm_param.setValue("chrom_fwhm", median_fwhm);
          ffm_param.setValue("charge_lower_bound", 2);
          ffm_param.setValue("charge_upper_bound", 5);
          ffm_param.setValue("report_chromatograms", "false");
          ffm.setLogType(log_type_);
          ffm.setParameters(ffm_param);
          std::vector<std::vector< OpenMS::MSChromatogram > > chromatograms;
          writeDebug_("Parameters passed to FeatureFindingMetabo algorithm", ffm_param, 3);

          ffm.run(splitted_mtraces, seeds, chromatograms);
          LOG_INFO << "Using " << seeds.size() << " seeds from untargeted feature extraction." << endl;
        }

        /////////////////////////////////////////////////
        // Run FeatureFinderIdentification

        FeatureMap fm;
        StringList feature_msfile_ref;
        feature_msfile_ref.push_back(mz_file);
        fm.setPrimaryMSRunPath(feature_msfile_ref);
        feature_maps.push_back(fm);

        FeatureFinderIdentificationAlgorithm ffi;
        ffi.getMSData().swap(ms_centroided);
        ffi.getProgressLogger().setLogType(log_type_);

        Param ffi_param = getParam_().copy("PeptideQuantification:", true);
        ffi_param.setValue("detect:peak_width", 5.0 * median_fwhm);
        ffi.setParameters(ffi_param);
        writeDebug_("Parameters passed to FeatureFinderIdentification algorithm", ffi_param, 3);

        ffi.run(peptide_ids, protein_ids, ext_peptide_ids, ext_protein_ids, feature_maps.back(), seeds);
        
        if (debug_level_ > 666)
        {
          FeatureXMLFile().store("debug_fraction_" + String(ms_files.first) + "_" + String(fraction_group) + ".featureXML", feature_maps.back());
        }

        // TODO: think about external ids ;) / maybe for technical replicates?

        // TODO: free parts of feature map not needed for further processing (e.g., subfeatures...)
        ++fraction_group;
      }

      //-------------------------------------------------------------
      // Align all features of this fraction
      //-------------------------------------------------------------
      if (feature_maps.size() > 1) // do we have several maps to align / link?
      {
        vector<TransformationDescription> transformations;
    
        // Determine reference from data, otherwise a change in order of input files
        // leads to slightly different results
        const int reference_index(-1); // set no reference (determine from data)
        MapAlignmentAlgorithmIdentification aligner;
        aligner.setLogType(log_type_);
        aligner.setParameters(ma_param);
        aligner.align(feature_maps, transformations, reference_index);

        // find model parameters:
        Param model_params = TOPPMapAlignerBase::getModelDefaults("b_spline");
        String model_type = model_params.getValue("type");
        model_params = model_params.copy(model_type + ":", true);
        
        vector<TransformationDescription::TransformationStatistics> alignment_stats;
        for (TransformationDescription & t : transformations)
        {
          writeDebug_("Using " + String(t.getDataPoints().size()) + " points in fit.", 1); 
          if (t.getDataPoints().size() > 10)
          {
            t.fitModel(model_type, model_params);
          }
          t.printSummary(LOG_DEBUG);
          alignment_stats.emplace_back(t.getStatistics());
        }

        // determine maximum RT shift after transformation that includes all high confidence IDs 
        using TrafoStat = TransformationDescription::TransformationStatistics;
        double max_alignment_diff = std::max_element(alignment_stats.begin(), alignment_stats.end(),
                [](TrafoStat a, TrafoStat b) 
                { return a.percentiles_after[100] > b.percentiles_after[100]; })->percentiles_after[100];
        // sometimes, very good alignments might lead to bad overall performance. Choose 2 minutes as minimum.
         max_alignment_diff = std::max(max_alignment_diff, 120.0);


        // Apply transformations
        for (Size i = 0; i < feature_maps.size(); ++i)
        {
          MapAlignmentTransformer::transformRetentionTimes(feature_maps[i],
            transformations[i]);
        }
        addDataProcessing_(consensus_fraction,
                       getProcessingInfo_(DataProcessing::ALIGNMENT));
        //-------------------------------------------------------------
        // Link all features of this fraction
        //-------------------------------------------------------------
        writeDebug_("Linking: " + String(feature_maps.size()) + " features.", 1);
        // grouping tolerance = max alignment error + median FWHM
        FeatureGroupingAlgorithmQT linker;
        fl_param.setValue("distance_RT:max_difference", 2.0 * max_alignment_diff + 2.0 * median_fwhm);
        fl_param.setValue("distance_MZ:max_difference", 10.0);
        fl_param.setValue("distance_MZ:unit", "ppm");
        fl_param.setValue("distance_MZ:weight", 5.0);
        fl_param.setValue("distance_intensity:weight", 0.1); 

/*      FeatureGroupingAlgorithmKD linker;
        fl_param.setValue("warp:rt_tol", 2.0 * max_alignment_diff + 2.0 * median_fwhm);
        fl_param.setValue("link:rt_tol", 2.0 * max_alignment_diff + 2.0 * median_fwhm);
        fl_param.setValue("link:mz_tol", 10.0);
        fl_param.setValue("mz_unit", "ppm");
  */      
        linker.setParameters(fl_param);      
        linker.group(feature_maps, consensus_fraction);
        addDataProcessing_(consensus_fraction,
                       getProcessingInfo_(DataProcessing::FEATURE_GROUPING));

      }
      else // only one feature map
      {
        MapConversion::convert(0, feature_maps.back(), consensus_fraction);                           
      }

      ////////////////////////////////////////////////////////////
      // Annotate experimental design in consensus map
      ////////////////////////////////////////////////////////////
      Size j(0);
      // for each MS file (as provided in the experimental design)
      for (String const & mz_file : ms_files.second) 
      {
        const Size fraction = ms_files.first;
        const Size fraction_group = j + 1;
        consensus_fraction.getColumnHeaders()[j].label = "label-free";
        consensus_fraction.getColumnHeaders()[j].filename = mz_file;
        consensus_fraction.getColumnHeaders()[j].unique_id = feature_maps[j].getUniqueId();
        consensus_fraction.getColumnHeaders()[j].setMetaValue("fraction", fraction);
        consensus_fraction.getColumnHeaders()[j].setMetaValue("fraction_group", fraction_group);
        ++j;
      }

      // assign unique ids
      consensus_fraction.applyMemberFunction(&UniqueIdInterface::setUniqueId);

      // sort list of peptide identifications in each consensus feature by map index
      consensus_fraction.sortPeptideIdentificationsByMapIndex();

      if (debug_level_ >= 666)
      {
        ConsensusXMLFile().store("debug_fraction_" + String(ms_files.first) +  ".consensusXML", consensus_fraction);
        writeDebug_("to produce a consensus map with: " + String(consensus_fraction.getColumnHeaders().size()) + " columns.", 1);
      }

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

    consensus.sortByPosition();
    consensus.sortPeptideIdentificationsByMapIndex();

    if (debug_level_ >= 666)
    {
      ConsensusXMLFile().store("debug_after normalization.consensusXML", consensus);
    }



    // TODO: FileMerger merge ids (here? or already earlier? filtered?)
    // TODO: check if it makes sense to integrate SVT imputation algorithm (branch)

    //-------------------------------------------------------------
    // Protein inference
    //-------------------------------------------------------------
    // TODO: ProteinInference on merged ids (how merged?)
    // TODO: Output coverage on protein (and group level?)


    //-------------------------------------------------------------
    // Protein inference
    //-------------------------------------------------------------
    // TODO: Implement inference
    vector<PeptideIdentification> infered_peptides;
    vector<ProteinIdentification> infered_protein_groups(1, ProteinIdentification());

    // first handle unique peptides and indistinguishable groups
    infered_peptides = consensus.getUnassignedPeptideIdentifications();
    for (ConsensusFeature & c : consensus)
    {
      infered_peptides.insert(infered_peptides.end(),
                      c.getPeptideIdentifications().begin(),
                      c.getPeptideIdentifications().end());
    }

    for (auto & p : infered_peptides)
    {
      // annotate experimental design meta data to merged PeptideIdentifications (the infered peptides)
      // so we can track their origin
      const String& old_identifier = p.getIdentifier();
      vector<ProteinIdentification> & protein_ids = consensus.getProteinIdentifications();
      auto it = std::find_if(protein_ids.begin(), protein_ids.end(), [&old_identifier](const ProteinIdentification & prot_id)
        {
          return prot_id.getIdentifier() == old_identifier;
        });
     
      if (it != protein_ids.end())
      {
        p.setMetaValue("fraction", it->getMetaValue("fraction"));
        p.setMetaValue("fraction_group", it->getMetaValue("fraction_group"));
        StringList spectra_data;
        it->getPrimaryMSRunPath(spectra_data);
        p.setMetaValue("spectra_data", spectra_data[0]);
      }
      else
      {
        // TODO: throw exception
      }

      // We need to set the common identifier for the final inference result.
      // see below and in MzTab.cpp
      p.setIdentifier("Fido"); 
    }

    // reindex peptides to proteins and keep only unique peptides
    if (!in_db.empty())
    {
      PeptideIndexing indexer;
      Param param_pi = indexer.getParameters();
      param_pi.setValue("enzyme:specificity", "none");  // TODO: derive from id files?  
      param_pi.setValue("missing_decoy_action", "silent");
      param_pi.setValue("write_protein_sequence", "true");
      param_pi.setValue("write_protein_description", "true");
      indexer.setParameters(param_pi);

      // stream data in fasta file
      FASTAContainer<TFI_File> fasta_db(in_db);
      PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, infered_protein_groups, infered_peptides);

      if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
          (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
      {
        if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
        {
          return INPUT_FILE_EMPTY;       
        }
        else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
        {
          return UNEXPECTED_RESULT;
        }
        else
        {
          return UNKNOWN_ERROR;
        }
      } 

      // ensure that only one final inference result is generated
      assert(infered_protein_groups.size() == 1);

      // merge (unique) protein identifications
      for (auto & p : infered_protein_groups)
      { 
        // copy hits to protein groups (required for export in mzTab)
        for (auto & h : p.getHits()) 
        { 
          const String & a = h.getAccession();
          ProteinIdentification::ProteinGroup pg;
          pg.accessions.push_back(a);
          infered_protein_groups[0].insertIndistinguishableProteins(pg);
          infered_protein_groups[0].insertProteinGroup(pg);          
        }
        p.setIdentifier("Fido"); 
      }

      // only keep unique peptides (for now)
      IDFilter::keepUniquePeptidesPerProtein(infered_peptides);

      // compute coverage
      infered_protein_groups[0].computeCoverage(infered_peptides);

      // determine observed modifications (exclude fixed mods)
      StringList fixed_mods = {"Carbamidomethyl"};  // TODO: read from search settings, make specific to amino acid (site)
      infered_protein_groups[0].computeModifications(infered_peptides, fixed_mods);

      /////////////////////////////////////////
      // annotate some mzTab related protein statistics

      map<String, map<String, Size >> acc2psms; // map runpath->accession->#PSMs (how many PSMs identify a protein in every run)
      // Note: only helpful if the PSM maps to one protein (or an indistinguishable group) - probably not helpful for shared peptides/PSMs 

      // distinct peptides: different if umodified sequence is different (charge state or modified forms of same peptide are not counted separately) 
      map<String, map<String, set<String>>> acc2distinct_peptides; // map runpath->accession->set of distinct peptides
      // unique peptides (not shared)
      map<String, map<String, set<String>>> acc2unique_peptides;
      for (Size pep_i = 0; pep_i != infered_peptides.size(); ++pep_i)
      {
        // peptide hits
        const PeptideIdentification & peptide_id = infered_peptides[pep_i];
        const String & runpath = peptide_id.getMetaValue("spectra_data");
        const vector<PeptideHit>& peptide_hits = peptide_id.getHits();

        for (Size ph_i = 0; ph_i != peptide_hits.size(); ++ph_i)
        {
          const PeptideHit & peptide_hit = peptide_hits[ph_i];
          const std::vector<PeptideEvidence>& ph_evidences = peptide_hit.getPeptideEvidences();
          if (ph_evidences.empty()) continue; // TODO: check
          const AASequence & aas = peptide_hit.getSequence();
          const String & seq_nomod = aas.toUnmodifiedString();          

          // unique peptide? store peptide sequence in map runpath->protein accession
          if (peptide_hit.extractProteinAccessionsSet().size() == 1) // == unique Note: we need to check the set as the vector may contain several references into the same protein
          {
            // 
            const String acc = *(peptide_hit.extractProteinAccessionsSet().begin());
            acc2unique_peptides[runpath][acc].insert(seq_nomod);
          } 

          for (Size phe_i = 0; phe_i != ph_evidences.size(); ++phe_i)
          {
            const String acc = ph_evidences[phe_i].getProteinAccession();

            // count PSM for the referenced protein
            acc2psms[runpath][acc] += 1; 
            // side note: proteins in indistinguishable groups will have the same numbers of PSMs associated

            // store (later: count) unmodified peptide sequences referencing the protein
            acc2distinct_peptides[runpath][acc].insert(seq_nomod);
          }
        }
      }

      // store run level mzTab statistics in protein hits of inference run (= final result)
      for (auto & p : infered_protein_groups[0].getHits())
      {
        const String acc = p.getAccession();
        
        IntList npsms, ndistinct, nunique;

        // TODO: validate somehow that the order of MS run interation is correct/identical in every part of this tool
        for (auto const ms_files : frac2ms) // for each fraction->ms file(s)
        {
          for (const String & runpath : ms_files.second)
          {
            if (acc2psms.at(runpath).count(acc) > 0)
            {
              npsms.push_back(acc2psms.at(runpath).at(acc));
            }
            else
            {
              npsms.push_back(0);
            }

            if (acc2distinct_peptides.at(runpath).count(acc) > 0)
            {
              auto distinct_peptides = acc2distinct_peptides.at(runpath).at(acc);
              ndistinct.push_back(distinct_peptides.size());
            }
            else
            {
              ndistinct.push_back(0);
            }
            
            if (acc2unique_peptides.at(runpath).count(acc) > 0)
            {
              nunique.push_back(acc2unique_peptides.at(runpath).at(acc).size());
            }
            else
            {
              nunique.push_back(0);
            }
          }
        }

        // will be exported in mzTab PRT section
        p.setMetaValue("num_psms_ms_run", npsms);
        p.setMetaValue("num_peptides_distinct_ms_run", ndistinct);
        p.setMetaValue("num_peptides_unique_ms_run", nunique);
      }
    }

    //-------------------------------------------------------------
    // Peptide quantification
    //-------------------------------------------------------------
    PeptideAndProteinQuant quantifier;
    quantifier.setParameters(pq_param);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides(infered_peptides);

    //-------------------------------------------------------------
    // Protein quantification
    //-------------------------------------------------------------
    // TODO: ProteinQuantifier on (merged?) consensusXML (with 1% FDR?) + inference ids (unfiltered?)? 
    if (infered_protein_groups[0].getIndistinguishableProteins().empty())
    {
      throw Exception::MissingInformation(
       __FILE__, 
       __LINE__, 
       OPENMS_PRETTY_FUNCTION, 
       "No information on indistinguishable protein groups found.");
    }

    quantifier.quantifyProteins(infered_protein_groups[0]);
    //-------------------------------------------------------------
    // Export of MzTab file as final output
    //-------------------------------------------------------------

    // Annotate quants to protein(groups) for easier export in mzTab
    auto const & protein_quants = quantifier.getProteinResults();

    // annotates final quantities to proteins and protein groups in the ID data structure
    PeptideAndProteinQuant::annotateQuantificationsToProteins(protein_quants, infered_protein_groups[0]);
    vector<ProteinIdentification>& proteins = consensus.getProteinIdentifications();
    proteins.insert(proteins.begin(), infered_protein_groups[0]); // insert inference information as first protein identification
    proteins[0].setSearchEngine("Fido");  // Note: currently needed so mzTab Exporter knows how to handle inference data in first prot. ID

    consensus.resolveUniqueIdConflicts(); // TODO: find out why this is needed to get proper UIDs in consensus
    if (debug_level_ >= 666)
    {
      // Note: idXML and consensusXML doesn't support writing quantification at protein groups
      // Note: consensusXML currently doesn't support writing out inference data
      // (they are neverless stored and passed to mzTab for proper export)
      IdXMLFile().store("debug_ids.idXML", proteins, infered_peptides);
      ConsensusXMLFile().store("debug_consensus.consensusXML", consensus);
    }

    // Fill MzTab with meta data and quants annotated in identification data structure
    const bool report_unmapped(true);
    const bool report_unidentified_features(false);
    MzTab m = MzTab::exportConsensusMapToMzTab(
      consensus, 
      String("null"),   
      report_unidentified_features, 
      report_unmapped);
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

