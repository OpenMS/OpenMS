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

#ifndef OPENMS_FORMAT_SWATHFILE_H
#define OPENMS_FORMAT_SWATHFILE_H

// Datastructures
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#ifdef OPENMS_FORMAT_SWATHFILE_MZXMLSUPPORT
#include <OpenMS/FORMAT/MzXMLFile.h>
#endif

#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteSwathHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h>

namespace OpenMS
{

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
    std::vector<OpenSwath::SwathMap> loadSplit(StringList file_list, String tmp,
                                               boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions = "normal")
    {
      int progress = 0;
      startProgress(0, file_list.size(), "Loading data");

      std::vector<OpenSwath::SwathMap> swath_maps(file_list.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(file_list.size()); ++i)
      {

#ifdef _OPENMP
#pragma omp critical (OPENMS_SwathFile_loadSplit)
#endif
        {
          std::cout << "Loading file " << i << " with name " << file_list[i] << " using readoptions " << readoptions << std::endl;
        }

        String tmp_fname = "openswath_tmpfile_" + String(i) + ".mzML";

        boost::shared_ptr<PeakMap > exp(new PeakMap);
        OpenSwath::SpectrumAccessPtr spectra_ptr;

        // Populate meta-data
        if (i == 0)
        {
          exp_meta = populateMetaData_(file_list[i]);
        }

        if (readoptions == "normal")
        {
          MzMLFile().load(file_list[i], *exp.get());
          spectra_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
        }
        else if (readoptions == "cache")
        {
          // Cache and load the exp (metadata only) file again
          spectra_ptr = doCacheFile_(file_list[i], tmp, tmp_fname, exp);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "Unknown option " + readoptions);
        }

        OpenSwath::SwathMap swath_map;

        bool ms1 = false;
        double upper = -1, lower = -1;
        if (exp->size() == 0)
        {
          std::cerr << "WARNING: File " << file_list[i] << "\n does not have any scans - I will skip it" << std::endl;
          continue;
        }
        if (exp->getSpectra()[0].getPrecursors().size() == 0)
        {
          std::cout << "NOTE: File " << file_list[i] << "\n does not have any precursors - I will assume it is the MS1 scan." << std::endl;
          ms1 = true;
        }
        else
        {
          // Checks that this is really a SWATH map and extracts upper/lower window
          OpenSwathHelper::checkSwathMap(*exp.get(), lower, upper);
        }

        swath_map.sptr = spectra_ptr;
        swath_map.lower = lower;
        swath_map.upper = upper;
        swath_map.ms1 = ms1;
#ifdef _OPENMP
#pragma omp critical (OPENMS_SwathFile_loadSplit)
#endif
        {
          LOG_DEBUG << "Adding Swath file " << file_list[i] << " with " << swath_map.lower << " to " << swath_map.upper << std::endl;
          swath_maps[i] = swath_map;
          setProgress(progress++);
        }
      }
      endProgress();
      return swath_maps;
    }

    /// Loads a Swath run from a single mzML file
    std::vector<OpenSwath::SwathMap> loadMzML(String file, String tmp,
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions = "normal")
    {
      std::cout << "Loading mzML file " << file << " using readoptions " << readoptions << std::endl;
      String tmp_fname = "openswath_tmpfile";

      startProgress(0, 1, "Loading metadata file " + file);
      boost::shared_ptr<PeakMap> experiment_metadata = populateMetaData_(file);
      exp_meta = experiment_metadata;

      // First pass through the file -> get the meta data
      std::cout << "Will analyze the metadata first to determine the number of SWATH windows and the window sizes." << std::endl;
      std::vector<int> swath_counter;
      int nr_ms1_spectra;
      std::vector<OpenSwath::SwathMap> known_window_boundaries;
      countScansInSwath_(experiment_metadata->getSpectra(), swath_counter, nr_ms1_spectra, known_window_boundaries);
      std::cout << "Determined there to be " << swath_counter.size() <<
        " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
      endProgress();

      FullSwathFileConsumer* dataConsumer;
      boost::shared_ptr<PeakMap> exp(new PeakMap);
      startProgress(0, 1, "Loading data file " + file);
      if (readoptions == "normal")
      {
        dataConsumer = new RegularSwathFileConsumer(known_window_boundaries);
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else if (readoptions == "cache")
      {
        dataConsumer = new CachedSwathFileConsumer(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else if (readoptions == "split")
      {
        dataConsumer = new MzMLSwathFileConsumer(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Unknown or unsupported option " + readoptions);
      }
      LOG_DEBUG << "Finished parsing Swath file " << std::endl; 
      std::vector<OpenSwath::SwathMap> swath_maps;
      dataConsumer->retrieveSwathMaps(swath_maps);
      delete dataConsumer;

      endProgress();
      return swath_maps;
    }

    /// Loads a Swath run from a single mzXML file
    std::vector<OpenSwath::SwathMap> loadMzXML(String file, String tmp,
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions = "normal")
    {
      std::cout << "Loading mzXML file " << file << " using readoptions " << readoptions << std::endl;
      String tmp_fname = "openswath_tmpfile";

      startProgress(0, 1, "Loading metadata file " + file);
      boost::shared_ptr<PeakMap > experiment_metadata(new PeakMap);
      MzXMLFile f;
      f.getOptions().setAlwaysAppendData(true);
      f.getOptions().setFillData(false);
      f.load(file, *experiment_metadata);
      exp_meta = experiment_metadata;

      // First pass through the file -> get the meta data
      std::cout << "Will analyze the metadata first to determine the number of SWATH windows and the window sizes." << std::endl;
      std::vector<int> swath_counter;
      int nr_ms1_spectra;
      std::vector<OpenSwath::SwathMap> known_window_boundaries;
      countScansInSwath_(experiment_metadata->getSpectra(), swath_counter, nr_ms1_spectra, known_window_boundaries);
      std::cout << "Determined there to be " << swath_counter.size() <<
        " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
      endProgress();

      FullSwathFileConsumer* dataConsumer;
      boost::shared_ptr<PeakMap > exp(new PeakMap);
      startProgress(0, 1, "Loading data file " + file);
      if (readoptions == "normal")
      {
        dataConsumer = new RegularSwathFileConsumer(known_window_boundaries);
        MzXMLFile().transform(file, dataConsumer, *exp.get());
      }
      else if (readoptions == "cache")
      {
        dataConsumer = new CachedSwathFileConsumer(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
        MzXMLFile().transform(file, dataConsumer, *exp.get());
      }
      else if (readoptions == "split")
      {
        dataConsumer = new MzMLSwathFileConsumer(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
        MzXMLFile().transform(file, dataConsumer, *exp.get());
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Unknown or unsupported option " + readoptions);
      }
      LOG_DEBUG << "Finished parsing Swath file " << std::endl; 
      std::vector<OpenSwath::SwathMap> swath_maps;
      dataConsumer->retrieveSwathMaps(swath_maps);
      delete dataConsumer;

      endProgress();
      return swath_maps;
    }

    /// Loads a Swath run from a single sqMass file
    std::vector<OpenSwath::SwathMap> loadSqMass(String file, boost::shared_ptr<ExperimentalSettings>& /* exp_meta */)
    {
      startProgress(0, 1, "Loading sqmass data file " + file);

      OpenMS::Internal::MzMLSqliteSwathHandler sql_mass_reader(file);
      std::vector<OpenSwath::SwathMap> swath_maps = sql_mass_reader.readSwathWindows();
      for (Size k = 0; k < swath_maps.size(); k++)
      {
        std::vector<int> indices = sql_mass_reader.readSpectraForWindow(swath_maps [k]);
        OpenMS::Internal::MzMLSqliteHandler handler(file);
        OpenSwath::SpectrumAccessPtr sptr(new OpenMS::SpectrumAccessSqMass(handler, indices));
        swath_maps[k].sptr = sptr;
      }

      // also store the MS1 map
      OpenSwath::SwathMap ms1_map;
      std::vector<int> indices = sql_mass_reader.readMS1Spectra();
      OpenMS::Internal::MzMLSqliteHandler handler(file);
      OpenSwath::SpectrumAccessPtr sptr(new OpenMS::SpectrumAccessSqMass(handler, indices));
      ms1_map.sptr = sptr;
      ms1_map.ms1 = true;
      swath_maps.push_back(ms1_map);
      endProgress();

      std::cout << "Determined there to be " << swath_maps.size() <<
        " SWATH windows and in total " << indices.size() << " MS1 spectra" << std::endl;

      return swath_maps;
    }

protected:

    /// Cache a file to disk
    OpenSwath::SpectrumAccessPtr doCacheFile_(String in, String tmp, String tmp_fname,
                                              boost::shared_ptr<PeakMap > experiment_metadata)
    {
      String cached_file = tmp + tmp_fname + ".cached";
      String meta_file = tmp + tmp_fname;

      // Create new consumer, transform infile, write out metadata
      MSDataCachedConsumer* cachedConsumer = new MSDataCachedConsumer(cached_file, true);
      MzMLFile().transform(in, cachedConsumer, *experiment_metadata.get());
      CachedmzML().writeMetadata(*experiment_metadata.get(), meta_file, true);
      delete cachedConsumer; // ensure that filestream gets closed

      boost::shared_ptr<PeakMap > exp(new PeakMap);
      MzMLFile().load(meta_file, *exp.get());
      return SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
    }

    /// Only read the meta data from a file and use it to populate exp_meta
    boost::shared_ptr< PeakMap > populateMetaData_(String file)
    {
      boost::shared_ptr<PeakMap > experiment_metadata(new PeakMap);
      MzMLFile f;
      f.getOptions().setAlwaysAppendData(true);
      f.getOptions().setFillData(false);
      f.load(file, *experiment_metadata);
      return experiment_metadata;
    }

    /// Counts the number of scans in a full Swath file (e.g. concatenated non-split file)
    void countScansInSwath_(const std::vector<MSSpectrum>& exp,
                            std::vector<int>& swath_counter, int& nr_ms1_spectra, 
                            std::vector<OpenSwath::SwathMap>& known_window_boundaries)
    {
      int ms1_counter = 0;
      for (Size i = 0; i < exp.size(); i++)
      {
        const MSSpectrum& s = exp[i];
        {
          if (s.getMSLevel() == 1)
          {
            ms1_counter++;
          }
          else 
          {
            if (s.getPrecursors().empty())
            {
              throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                  "Found SWATH scan (MS level 2 scan) without a precursor. Cannot determine SWATH window.");
            }
            const std::vector<Precursor> prec = s.getPrecursors();
            double center = prec[0].getMZ();
            bool found = false;
            for (Size j = 0; j < known_window_boundaries.size();  j++)
            {
              // We group by the precursor mz (center of the window) since this
              // should be present
              if (std::fabs(center - known_window_boundaries[j].center) < 1e-6)
              {
                found = true;
                swath_counter[j]++;
              }
            }
            if (!found)
            {
              // we found a new SWATH scan
              swath_counter.push_back(1);
              double lower = prec[0].getMZ() - prec[0].getIsolationWindowLowerOffset();
              double upper = prec[0].getMZ() + prec[0].getIsolationWindowUpperOffset();
              OpenSwath::SwathMap boundary;
              boundary.lower = lower;
              boundary.upper = upper;
              boundary.center = center;
              known_window_boundaries.push_back(boundary);

              LOG_DEBUG << "Adding Swath centered at " << center 
                << " m/z with an isolation window of " << lower << " to " << upper 
                << " m/z." << std::endl;
            }
          }
        }
      }
      nr_ms1_spectra = ms1_counter;

      std::cout << "Determined there to be " << swath_counter.size() <<
        " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
    }

  };
}

#endif
