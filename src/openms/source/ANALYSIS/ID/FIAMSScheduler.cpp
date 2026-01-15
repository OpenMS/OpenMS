// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Svetlana Kutuzova, Douglas McCloskey $
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FIAMSScheduler.h>

#include <OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <QDir>

namespace OpenMS {
  /// default constructor
  FIAMSScheduler::FIAMSScheduler(
    String filename,
    String base_dir,
    String output_dir,
    bool load_cached
  )
    : 
    filename_(std::move(filename)),
    base_dir_(std::move(base_dir)),
    output_dir_(std::move(output_dir)),
    load_cached_(load_cached),
    samples_()
  {
  loadSamples_();
  }

  void FIAMSScheduler::loadSamples_() {
    CsvFile csv_file(filename_, ',');
    StringList headers;
    csv_file.getRow(0, headers);
    StringList row;
    for (Size i = 1; i < csv_file.rowCount(); ++i) {
      csv_file.getRow(i, row);
      std::map<String, String> mapping;
      for (Size j = 0; j < headers.size(); ++j) {
        mapping[headers[j]] = row[j];
      }
      samples_.push_back(mapping);
    }
  }

  void FIAMSScheduler::run() {
    #pragma omp parallel for
    for (int i = 0; i < (int)samples_.size(); ++i) {
      MSExperiment exp;
      FileHandler().loadExperiment(base_dir_ + samples_[i].at("dir_input") + "/" + samples_[i].at("filename") + ".mzML", exp, {FileTypes::MZML});

      FIAMSDataProcessor fia_processor;
      Param p;
      p.setValue("filename", samples_[i].at("filename"));
      p.setValue("dir_output", output_dir_ + samples_[i].at("dir_output"));
      QDir qd;
      qd.mkpath(p.getValue("dir_output").toString().c_str());
      p.setValue("resolution", std::stof(samples_[i].at("resolution")));
      p.setValue("polarity", samples_[i].at("charge"));
      p.setValue("db:mapping", std::vector<std::string>{base_dir_ + samples_[i].at("db_mapping")});
      p.setValue("db:struct", std::vector<std::string>{base_dir_ + samples_[i].at("db_struct")});
      p.setValue("positive_adducts", base_dir_ + samples_[i].at("positive_adducts"));
      p.setValue("negative_adducts", base_dir_ + samples_[i].at("negative_adducts"));
      fia_processor.setParameters(p);

      String time = samples_[i].at("time");
      std::vector<String> times;
      time.split(";", times);
      for (Size j = 0; j < times.size(); ++j) {
        OPENMS_LOG_INFO << "Started " << samples_[i].at("filename") << " for " << times[j] << " seconds" << std::endl;
        MzTab mztab_output;
        fia_processor.run(exp, std::stof(times[j]), mztab_output, load_cached_);
        OPENMS_LOG_INFO << "Finished " << samples_[i].at("filename") << " for " << times[j] << " seconds" << std::endl;
      }
    }
  }

  const std::vector<std::map<String, String>>& FIAMSScheduler::getSamples() {
    return samples_;
  }

  const String& FIAMSScheduler::getBaseDir() {
    return base_dir_;
  }

  const String& FIAMSScheduler::getOutputDir() {
    return output_dir_;
  }
}
