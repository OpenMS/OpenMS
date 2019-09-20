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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  class ExperimentalSettings;
  namespace Interfaces
  {
    class IMSDataConsumer;
  }
  

  /**
   * @brief File adapter for Swath files.
   *
   * This class can load SWATH files in different storage versions. The most
   * convenient file is a single MzML file which contains one experiment.
   * However, also the loading of a list of files is supported (loadSplit)
   * where it is assumed that each individual file only contains scans from one
   * precursor isolation window (one SWATH). Finally, experimental support for
   * mzXML is available but needs to be selected with a specific compile flag
   * (this is not for everyday use).
   *
   */
  class OPENMS_DLLAPI SwathFile :
    public ProgressLogger
  {
public:

    /// Loads a Swath run from a list of split mzML files
    std::vector<OpenSwath::SwathMap> loadSplit(StringList file_list,
                                               String tmp,
                                               boost::shared_ptr<ExperimentalSettings>& exp_meta, 
                                               String readoptions = "normal");

    /**
      @brief Loads a Swath run from a single mzML file

      Using the @p plugin_consumer, you can provide a custom consumer which will be chained
      into the process of loading the data and making it available (depending on @p readoptions).
      This is useful if you want to modify the data a priori or extract some other information using
      MSDataTransformingConsumer (for example). Make sure it leaves the data intact, such that the 
      returned SwathMaps are actually useful.

      @param [IN] file Input filename
      @param [IN] tmp Temporary directory (for cached data)
      @param [OUT] exp_meta Experimental metadata from mzML file
      @param [IN] readoptions How are spectra accessed after reading - tradeoff between memory usage and time (disk caching)
      @param [IN] plugin_consumer An intermediate custom consumer
      @return Swath maps for MS2 and MS1 (unless readoptions == split, which returns no data)
    */
    std::vector<OpenSwath::SwathMap> loadMzML(const String& file, 
                                              const String& tmp,
                                              boost::shared_ptr<ExperimentalSettings>& exp_meta,
                                              const String& readoptions = "normal",
                                              Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

    /// Loads a Swath run from a single mzXML file
    std::vector<OpenSwath::SwathMap> loadMzXML(String file, 
                                               String tmp,
                                               boost::shared_ptr<ExperimentalSettings>& exp_meta,
                                               String readoptions = "normal");

    /// Loads a Swath run from a single sqMass file
    std::vector<OpenSwath::SwathMap> loadSqMass(String file, boost::shared_ptr<ExperimentalSettings>& /* exp_meta */);

protected:

    /// Cache a file to disk
    OpenSwath::SpectrumAccessPtr doCacheFile_(const String& in, const String& tmp, const String& tmp_fname,
                                              boost::shared_ptr<PeakMap > experiment_metadata);

    /// Only read the meta data from a file and use it to populate exp_meta
    boost::shared_ptr< PeakMap > populateMetaData_(const String& file);

    /// Counts the number of scans in a full Swath file (e.g. concatenated non-split file)
    void countScansInSwath_(const std::vector<MSSpectrum>& exp,
                            std::vector<int>& swath_counter, int& nr_ms1_spectra, 
                            std::vector<OpenSwath::SwathMap>& known_window_boundaries);

  };
}

