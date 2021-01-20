// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>
#include <QFileInfo>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
  @page TOPP_FLASHDeconv TOPP_FLASHDeconv
  (Need to be modified)

  @brief  @ref
  @code
  @endcode
  @verbinclude
  @htmlinclude
*/
// We do not want this class to show up in the docu:
// NEED to fill this part later


class TOPPFLASHDeconv :
    public TOPPBase
{
public:
  TOPPFLASHDeconv() :
      TOPPBase("FLASHDeconv", "Ultra-fast high-quality deconvolution enables online processing of top-down MS data",
               false)
  {
  }

protected:
  //typedef FLASHDeconvHelperStructs::Parameter Parameter;

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {

    registerInputFile_("in", "<file>", "", "Input file (mzML)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("in_train", "<file>", "", "topPIC result *prsm.csv file for QScore training", false, true);
    setValidFormats_("in_train", ListUtils::create<String>("csv"));

    registerOutputFile_("out", "<file>", "",
                        "output file (tsv) - feature level deconvoluted masses (or ensemble spectrum level deconvluted mass if use_ensemble_spectrum is set to 1) ");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_train", "<file>", "",
                        "train result csv file for QScore training", false, true);
    setValidFormats_("out_train", ListUtils::create<String>("csv"));


    registerOutputFileList_("out_spec", "<file for MS1, file for MS2, ...>", {""},
                            "output files (tsv) - spectrum level deconvoluted masses per ms level", false);
    setValidFormats_("out_spec", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_mzml", "<file>", "",
                        "mzml format output file (mzML) - spectrum level deconvoluted masses per ms level", false);
    setValidFormats_("out_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_promex", "<file>", "",
                        "promex format output file (ms1ft) - only MS1 deconvoluted masses are recorded", false);
    setValidFormats_("out_promex", ListUtils::create<String>("ms1ft"), false);

    registerOutputFileList_("out_topFD", "<file for MS1, file for MS2, ...>", {""},
                            "topFD format output files (msalign) - spectrum level deconvoluted masses per ms level. The file name for MSn should end with msn.msalign to be able to be recognized by TopPIC. "
                            "For example, -out_topFD [name]_ms1.msalign [name]_ms2.msalign", false);
    setValidFormats_("out_topFD", ListUtils::create<String>("msalign"), false);

    registerIntOption_("mzml_mass_charge",
                       "<0:uncharged 1: +1 charged -1: -1 charged>",
                       0,
                       "Charge status of deconvoluted masses in mzml output",
                       false);

    setMinInt_("mzml_mass_charge", -1);
    setMaxInt_("mzml_mass_charge", 1);

    /*registerDoubleList_("tol",
                        "<ms1_tol, ms2_tol, ...>",
                        {10.0, 10.0},
                        "ppm tolerance for MS1, 2, ...  "
                        "(e.g., -tol 10.0 15.0 to specify 10.0 and 15.0 ppm for MS1 and MS2, respectively)",
                        false);
*/
    registerIntOption_("write_detail",
                       "<1:true 0:false>",
                       0,
                       "to write peak info per deconvoluted mass in detail or not in [prefix]_MSn_spec.tsv files. If set to 1, all peak information (m/z, intensity, charge, and isotope index) per mass is reported.",
                       false,
                       false);

    setMinInt_("write_detail", 0);
    setMaxInt_("write_detail", 1);

    registerIntOption_("max_MS_level", "", 2, "maximum MS level (inclusive) for deconvolution", false, true);
    setMinInt_("max_MS_level", 1);

    registerIntOption_("use_ensemble_spectrum",
                       "",
                       0,
                       "if set to 1, all spectra will add up to a single ensemble spectrum (per MS level) for which the deconvolution is performed. out_spec should specify the output spectrum level deconvolution tsv files.",
                       false,
                       false);

    setMinInt_("use_ensemble_spectrum", 0);
    setMaxInt_("use_ensemble_spectrum", 1);

    registerIntOption_("use_RNA_averagine", "", 0, "if set to 1, RNA averagine model is used", false, true);
    setMinInt_("use_RNA_averagine", 0);
    setMaxInt_("use_RNA_averagine", 1);

    Param fd_defaults = FLASHDeconvAlgorithm().getDefaults();
    // overwrite algorithm default so we export everything (important for copying back MSstats results)
    fd_defaults.setValue("tol", DoubleList{10.0, 10.0}, "ppm tolerance");
    fd_defaults.setValue("min_charge", 1);
    fd_defaults.setValue("max_charge", 100);
    fd_defaults.setValue("min_mz", -1.0);
    fd_defaults.addTag("min_mz", "advanced");
    fd_defaults.setValue("max_mz", -1.0);
    fd_defaults.addTag("max_mz", "advanced");
    fd_defaults.setValue("min_rt", -1.0);
    fd_defaults.addTag("min_rt", "advanced");
    fd_defaults.setValue("max_rt", -1.0);
    fd_defaults.addTag("max_rt", "advanced");
    fd_defaults.setValue("min_mass_", 50.0);
    fd_defaults.setValue("max_mass", 100000.0);
    //fd_defaults.addTag("tol", "advanced"); // hide entry
    fd_defaults.setValue("min_peaks", IntList{3, 1});
    fd_defaults.addTag("min_peaks", "advanced");
    fd_defaults.setValue("min_intensity", .0, "intensity threshold");
    fd_defaults.addTag("min_intensity", "advanced");
    fd_defaults.setValue("min_isotope_cosine_",
                         DoubleList{.75, .75},
                         "cosine threshold between avg. and observed isotope pattern for MS1, 2, ... (e.g., -min_isotope_cosine_ 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)");
    //fd_defaults.addTag("min_isotope_cosine_", "advanced");

    fd_defaults.setValue("max_mass_count_",
                         IntList{-1, -1},
                         "maximum mass count per spec for MS1, 2, ... (e.g., -max_mass_count_ 100 50 to specify 100 and 50 for MS1 and MS2, respectively. -1 specifies unlimited)");
    fd_defaults.addTag("max_mass_count_", "advanced");


    fd_defaults.setValue("RT_window", 20.0, "RT window for MS1 deconvolution");
    fd_defaults.addTag("RT_window", "advanced");

    fd_defaults.remove("max_mass_count_");
    fd_defaults.remove("min_mass_count_");

    Param mf_defaults = MassFeatureTrace().getDefaults();
    mf_defaults.setValue("mass_error_da",
                         1.5,
                         "da tolerance for feature tracing. Due to frequent isotope errer, 1.5 Da is recommended.");
    mf_defaults.remove("mass_error_ppm"); // hide entry
    mf_defaults.remove("trace_termination_criterion");
    mf_defaults.remove("reestimate_mt_sd");
    mf_defaults.remove("noise_threshold_int");
    mf_defaults.remove("min_sample_rate");
    mf_defaults.remove("trace_termination_outliers"); // hide entry
    mf_defaults.setValue("min_trace_length", 10.0, "min feature trace length in second");//
    mf_defaults.setValue("quant_method", "area", "");
    mf_defaults.addTag("quant_method", "advanced"); // hide entry
    mf_defaults.setValue("min_isotope_cosine_", -1.0, "if not set, controlled by -Algorithm:min_isotope_cosine_ option");
    mf_defaults.addTag("min_isotope_cosine_", "advanced");

    Param combined;
    combined.insert("Algorithm:", fd_defaults);
    combined.insert("FeatureTracing:", mf_defaults);
    registerFullParam_(combined);
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    OPENMS_LOG_INFO << "Initializing ... " << endl;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String in_file = getStringOption_("in");
    String out_file = getStringOption_("out");
    String in_train_file = getStringOption_("in_train");
    String out_train_file = getStringOption_("out_train");
    auto out_spec_file = getStringList_("out_spec");
    String out_mzml_file = getStringOption_("out_mzml");
    String out_promex_file = getStringOption_("out_promex");
    auto out_topfd_file = getStringList_("out_topFD");

    bool use_RNA_averagine = getIntOption_("use_RNA_averagine") > 0;
    int max_ms_level = getIntOption_("max_MS_level");
    bool ensemble = getIntOption_("use_ensemble_spectrum") > 0;
    bool write_detail = getIntOption_("write_detail") > 0;
    int mzml_charge = getIntOption_("mzml_mass_charge");

    fstream out_stream, out_train_stream, out_promex_stream;
    std::vector<fstream> out_spec_streams, out_topfd_streams;

    fstream fi_out;
    //fi_out.open(in_file+".txt", fstream::out); // TDDO

    out_stream.open(out_file, fstream::out);
    MassFeatureTrace::writeHeader(out_stream);

    if(!out_promex_file.empty()){
      out_promex_stream.open(out_promex_file, fstream::out);
    }
    if(!out_topfd_file.empty()){
      out_topfd_streams = std::vector<fstream>(out_topfd_file.size());
      for(int i=0; i < out_topfd_file.size(); i++){
        out_topfd_streams[i].open(out_topfd_file[i], fstream::out);
      }
    }
    if(!out_spec_file.empty()){
      out_spec_streams = std::vector<fstream>(out_spec_file.size());
      for(int i=0; i < out_spec_file.size(); i++){
        out_spec_streams[i].open(out_spec_file[i], fstream::out);
        DeconvolutedSpectrum::writeDeconvolutedMassesHeader(out_spec_streams[i], i + 1, write_detail);
      }
    }

    std::set<int> train_scan_numbers;
    if (!in_train_file.empty() && !out_train_file.empty())
    {
      out_train_stream.open(out_train_file, fstream::out);
      QScore::writeAttHeader(out_train_stream);
      std::ifstream in_trainstream(in_train_file);
      String line;
      bool start = false;
      while (std::getline(in_trainstream, line))
      {
        if (line.rfind("Data file name", 0) == 0){
          start = true;
          continue;
        }
        if(!start){
          continue;
        }
        vector<String> results;
        stringstream  tmp_stream(line);
        String str;
        while (getline(tmp_stream, str, ',')) {
          results.push_back(str);
        }
        train_scan_numbers.insert(std::stoi(results[4]));
      }
      in_trainstream.close();
    }

    int current_max_ms_level = 0;

    auto spec_cntr = std::vector<int>(max_ms_level, 0);
    // spectrum number with at least one deconvoluted mass per ms level per input file
    auto qspec_cntr = std::vector<int>(max_ms_level, 0);
    // mass number per ms level per input file
    auto mass_cntr = std::vector<int>(max_ms_level, 0);
    // feature number per input file
    int feature_cntr = 0;

    // feature index written in the output file
    int feature_index = 1;

    MSExperiment ensemble_map;
    // generate ensemble spectrum if param.ensemble is set
    if (ensemble)
    {
      for (int i = 0; i < max_ms_level; i++)
      {
        auto spec = MSSpectrum();
        spec.setMSLevel(i + 1);
        std::ostringstream name;
        name << "ensemble MS" << (i + 1) << " spectrum";
        spec.setName(name.str());
        ensemble_map.addSpectrum(spec);
      }
    }

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    MSExperiment map;
    MzMLFile mzml;

    // all for measure elapsed cup wall time
    double elapsed_cpu_secs = 0, elapsed_wall_secs = 0;
    auto elapsed_deconv_cpu_secs = std::vector<double>(max_ms_level, .0);
    auto elapsed_deconv_wall_secs = std::vector<double>(max_ms_level, .0);

    auto begin = clock();
    auto t_start = chrono::high_resolution_clock::now();

    OPENMS_LOG_INFO << "Processing : " << in_file << endl;

    mzml.setLogType(log_type_);
    mzml.load(in_file, map);

    //      double rtDuration = map[map.size() - 1].getRT() - map[0].getRT();
    int ms1_cntr = 0;
    double ms2_cntr = .0; // for debug...
    current_max_ms_level = 0;

    // read input dataset once to count spectra and generate ensemble spectrum if necessary
    for (auto &it : map)
    {
      if (it.empty())
      {
        continue;
      }
      if (it.getMSLevel() > max_ms_level)
      {
        continue;
      }

      int ms_level = it.getMSLevel();
      current_max_ms_level = current_max_ms_level < ms_level ? ms_level : current_max_ms_level;

      if (ensemble)
      {
        auto &espec = ensemble_map[it.getMSLevel() - 1];
        for (auto &p : it)
        {
          espec.push_back(p);
        }
      }

      if (it.getMSLevel() == 1)
      {
        ms1_cntr++;
      }
      if (it.getMSLevel() == 2)
      {
        ms2_cntr++;
      }
    }

    // Max MS Level is adjusted according to the input dataset
    current_max_ms_level = current_max_ms_level > max_ms_level ? max_ms_level : current_max_ms_level;
    // if an ensemble spectrum is analyzed, replace the input dataset with the ensemble one
    if (ensemble)
    {
      for (int i = 0; i < current_max_ms_level; ++i)
      {
        ensemble_map[i].sortByPosition();
      }
      map = ensemble_map;
    }
    // Run FLASHDeconv here

    int scan_number = 0;
    float prev_progress = .0;
    auto last_deconvoluted_spectra = std::unordered_map<UInt, DeconvolutedSpectrum>();
    //auto lastlast_deconvoluted_spectra = std::unordered_map<UInt, DeconvolutedSpectrum>();

    MSExperiment exp;

    auto fd = FLASHDeconvAlgorithm();
    Param fd_param = getParam_().copy("Algorithm:", true);
    //fd_param.setValue("tol", getParam_().getValue("tol"));
    fd.setParameters(fd_param);
    fd.calculateAveragine(use_RNA_averagine);
    auto avg = fd.getAveragine();
    auto mass_tracer = MassFeatureTrace();
    Param mf_param = getParam_().copy("FeatureTracing:", true);
    DoubleList isotope_cosines = fd_param.getValue("min_isotope_cosine_");
    //mf_param.setValue("mass_error_ppm", ms1tol);
    mf_param.setValue("noise_threshold_int", .0, "");
    mf_param.setValue("reestimate_mt_sd", "false", "");
    mf_param.setValue("trace_termination_criterion", "outlier");
    mf_param.setValue("trace_termination_outliers", 20);
    //mf_param.setValue("min_charge_cosine", fd_param.getValue("min_charge_cosine"));
    if (((double)mf_param.getValue("min_isotope_cosine_")) < 0)
    {
      mf_param.setValue("min_isotope_cosine_", isotope_cosines[0]);
    }
    mass_tracer.setParameters(mf_param);
    //std::cout<<mass_tracer.getParameters()<<std::endl;

    OPENMS_LOG_INFO << "Running FLASHDeconv ... " << endl;


    for (auto it = map.begin(); it != map.end(); ++it)
    {
      scan_number = SpectrumLookup::extractScanNumber(it->getNativeID(),
                                                      map.getSourceFiles()[0].getNativeIDTypeAccession());
      if (it->empty())
      {
        continue;
      }

      int ms_level = it->getMSLevel();
      if (ms_level > current_max_ms_level)
      {
        continue;
      }

      if(ms_level == 1){
       // fi_out << "Spec\t"<<it->getRT()<<"\n";
        for(auto &p : *it){
          if(p.getIntensity() <= 0){
            continue;
          }
         // fi_out << p.getMZ() << "\t" << p.getIntensity()<<"\n";
        }
      }

      spec_cntr[ms_level - 1]++;
      auto deconv_begin = clock();
      auto deconv_t_start = chrono::high_resolution_clock::now();

      //auto deconvoluted_spectrum = DeconvolutedSpectrum(*it, scan_number);
      // for MS>1 spectrum, register precursor
      DeconvolutedSpectrum* precursor_spec = nullptr;

      if (ms_level > 1 && last_deconvoluted_spectra.find(ms_level - 1) != last_deconvoluted_spectra.end())
      {
        precursor_spec = &(last_deconvoluted_spectra[ms_level - 1]);
      }
                      // per spec deconvolution
      auto deconvoluted_spectrum = fd.getDeconvolutedSpectrum(*it, precursor_spec, scan_number);

      if (it->getMSLevel() == 2 && !in_train_file.empty() && !out_train_file.empty()
          && !deconvoluted_spectrum.getPrecursorPeakGroup().empty()){
        QScore::writeAttTsv(deconvoluted_spectrum.getOriginalSpectrum().getRT(), deconvoluted_spectrum.getPrecursorPeakGroup(),
                            deconvoluted_spectrum.getPrecursorCharge(),
                            train_scan_numbers.find(scan_number) != train_scan_numbers.end(), avg, out_train_stream);

      }
      if (!out_mzml_file.empty())
      {
        exp.addSpectrum(deconvoluted_spectrum.toSpectrum(mzml_charge));
      }
      elapsed_deconv_cpu_secs[ms_level - 1] += double(clock() - deconv_begin) / CLOCKS_PER_SEC;
      elapsed_deconv_wall_secs[ms_level - 1] += chrono::duration<double>(
          chrono::high_resolution_clock::now() - deconv_t_start).count();

      if (ms_level < current_max_ms_level)
      {
        //if(last_deconvoluted_spectra.find(ms_level) != last_deconvoluted_spectra.end())
        //{
        //  lastlast_deconvoluted_spectra[ms_level] = last_deconvoluted_spectra[ms_level];
        //}
        last_deconvoluted_spectra[ms_level] = deconvoluted_spectrum;
      }

      if (deconvoluted_spectrum.empty())
      {
        continue;
      }
      //if (ms_level < current_max_ms_level)
      //{
      //  last_deconvoluted_spectra[ms_level] = deconvoluted_spectrum_; // to register precursor in the future..
      //}

      if (!ensemble)
      {
        mass_tracer.storeInformationFromDeconvolutedSpectrum(deconvoluted_spectrum);// add deconvoluted mass in mass_tracer
      }

      qspec_cntr[ms_level - 1]++;
      mass_cntr[ms_level - 1] += deconvoluted_spectrum.size();
      if(out_spec_streams.size() > ms_level - 1)
      {
        deconvoluted_spectrum
            .writeDeconvolutedMasses(out_spec_streams[ms_level - 1], in_file, avg, write_detail);
      }
      if (out_topfd_streams.size() > ms_level - 1)
      {
        deconvoluted_spectrum.writeTopFD(out_topfd_streams[ms_level - 1], scan_number, avg);
      }

      //deconvoluted_spectrum_.clearPeakGroupsChargeInfo();
      //deconvoluted_spectrum_.getPrecursorPeakGroup().clearChargeInfo();
      float progress = (float) (it - map.begin()) / map.size();
      if (progress > prev_progress + .01)
      {
        printProgress_(progress);
        prev_progress = progress;
      }
    }

    printProgress_(1); //
    std::cout << std::endl;
    std::unordered_map<UInt, DeconvolutedSpectrum>().swap(last_deconvoluted_spectra); // empty memory
    //std::unordered_map<UInt, DeconvolutedSpectrum>().swap(lastlast_deconvoluted_spectra); // empty memory

    // mass_tracer run
    if (!ensemble)
    {
      mass_tracer
          .findFeatures(in_file, !out_promex_file.empty(), feature_cntr, feature_index, out_stream, out_promex_stream, fd.getAveragine());
    }
    if (!out_mzml_file.empty())
    {
      MzMLFile mzml_file;
      mzml_file.store(out_mzml_file, exp);
    }

    for (int j = 0; j < (int) current_max_ms_level; j++)
    {
      if (spec_cntr[j] == 0)
      {
        continue;
      }

      if (ensemble)
      {
        OPENMS_LOG_INFO << "So far, FLASHDeconv found " << mass_cntr[j] << " masses in the ensemble MS"
                        << (j + 1) << " spectrum" << endl;

      }
      else
      {
        OPENMS_LOG_INFO << "So far, FLASHDeconv found " << mass_cntr[j] << " masses in " << qspec_cntr[j]
                        << " MS" << (j + 1) << " spectra out of "
                        << spec_cntr[j] << endl;
      }
    }
    if (feature_cntr > 0)
    {
      OPENMS_LOG_INFO << "Mass tracer found " << feature_cntr << " features" << endl;
    }

    auto t_end = chrono::high_resolution_clock::now();
    auto end = clock();

    elapsed_cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
    elapsed_wall_secs = chrono::duration<double>(t_end - t_start).count();

    OPENMS_LOG_INFO << "-- done [took " << elapsed_cpu_secs << " s (CPU), " << elapsed_wall_secs
                    << " s (Wall)] --"
                    << endl;

    int total_spec_cntr = 0;
    for (int j = 0; j < (int) current_max_ms_level; j++)
    {
      total_spec_cntr += spec_cntr[j];

      OPENMS_LOG_INFO << "-- deconv per MS" << (j + 1) << " spectrum (except spec loading, feature finding) [took "
                      << 1000.0 * elapsed_deconv_cpu_secs[j] / total_spec_cntr
                      << " ms (CPU), " << 1000.0 * elapsed_deconv_wall_secs[j] / total_spec_cntr << " ms (Wall)] --" << endl;
    }

    //fi_out.close(); //

    out_stream.close();

    if(!out_promex_file.empty()){
      out_promex_stream.close();
    }
    if(!out_topfd_file.empty()){
      for(auto & out_topfd_stream : out_topfd_streams){
        out_topfd_stream.close();
      }
    }
    if(!out_spec_file.empty()){
      for(auto & out_spec_stream : out_spec_streams){
        out_spec_stream.close();
      }
    }

    if (!out_train_file.empty())
    {
      out_train_stream.close();
    }


    return EXECUTION_OK;
  }

  static void printProgress_(float progress)
  {
    float bar_width = 70;
    std::cout << "[";
    int pos = (int) (bar_width * progress);
    for (int i = 0; i < bar_width; ++i)
    {
      if (i < pos)
      {
        std::cout << "=";
      }
      else if (i == pos)
      {
        std::cout << ">";
      }
      else
      {
        std::cout << " ";
      }
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }
};

// the actual main function needed to create an executable
int main(int argc, const char **argv)
{
  TOPPFLASHDeconv tool;
  return tool.main(argc, argv);
}
