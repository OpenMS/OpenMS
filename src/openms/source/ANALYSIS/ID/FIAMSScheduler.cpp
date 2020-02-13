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
// $Maintainer: Svetlana Kutuzova, Douglas McCloskey $
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>
#include <OpenMS/ANALYSIS/ID/FIAMSScheduler.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS {
  /// default constructor
  FIAMSScheduler::FIAMSScheduler(
    String filename,
    String base_dir,
    bool load_cached
  )
    : 
    filename_(filename),
    base_dir_(base_dir),
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
      MzMLFile mzml;
      mzml.load(base_dir_ + samples_[i].at("dir_input") + "/" + samples_[i].at("filename") + ".mzML", exp);

      FIAMSDataProcessor fia_processor;
      Param p;
      p.setValue("filename", samples_[i].at("filename"));
      p.setValue("dir_output", base_dir_ + samples_[i].at("dir_output"));
      p.setValue("resolution", std::stof(samples_[i].at("resolution")));
      p.setValue("polarity", samples_[i].at("charge"));
      p.setValue("db:mapping", ListUtils::create<String>(base_dir_ + samples_[i].at("db_mapping")));
      p.setValue("db:struct", ListUtils::create<String>(base_dir_ + samples_[i].at("db_struct")));
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
}
