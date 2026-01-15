// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SwathFile.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataChainingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
#include <OpenMS/FORMAT/FileHandler.h>
//TODO remove MzML after we get transform support for our handlers
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteSwathHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/SYSTEM/File.h>

#include <memory> // for make_shared

namespace OpenMS
{

  using Interfaces::IMSDataConsumer;

  /// Loads a Swath run from a list of split mzML files
  std::vector<OpenSwath::SwathMap> SwathFile::loadSplit(StringList file_list,
        const String& tmp,
    boost::shared_ptr<ExperimentalSettings>& exp_meta,
    const String& readoptions)
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
        FileHandler().loadExperiment(file_list[i], *exp.get(), {FileTypes::MZML});
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
      double upper = -1, lower = -1, center = -1;
      if (exp->empty())
      {
        std::cerr << "WARNING: File " << file_list[i] << "\n does not have any scans - I will skip it" << std::endl;
        continue;
      }
      if (exp->getSpectra()[0].getPrecursors().empty())
      {
        std::cout << "NOTE: File " << file_list[i] << "\n does not have any precursors - I will assume it is the MS1 scan." << std::endl;
        ms1 = true;
      }
      else
      {
        // Checks that this is really a SWATH map and extracts upper/lower window
        OpenSwathHelper::checkSwathMap(*exp.get(), lower, upper, center);
      }

      swath_map.sptr = spectra_ptr;
      swath_map.lower = lower;
      swath_map.upper = upper;
      swath_map.center = center;
      swath_map.ms1 = ms1;
#ifdef _OPENMP
#pragma omp critical (OPENMS_SwathFile_loadSplit)
#endif
      {
        OPENMS_LOG_DEBUG << "Adding Swath file " << file_list[i] << " with " << swath_map.lower << " to " << swath_map.upper << std::endl;
        swath_maps[i] = swath_map;
        setProgress(progress++);
      }
    }
    endProgress();
    return swath_maps;
  }

  /// Loads a Swath run from a single mzML file
  std::vector<OpenSwath::SwathMap> SwathFile::loadMzML(const String& file,
                                                       const String& tmp,
                                                       boost::shared_ptr<ExperimentalSettings>& exp_meta,
                                                       const String& readoptions,
                                                       Interfaces::IMSDataConsumer* plugin_consumer)
  {
    std::cout << "Loading mzML file " << file << " using readoptions " << readoptions << std::endl;
    String tmp_fname = tmp.hasSuffix('/') ? File::getUniqueName() : ""; // use tmp-filename if just a directory was given

    startProgress(0, 1, "Loading metadata file " + file);
    boost::shared_ptr<PeakMap> exp_stripped = populateMetaData_(file);
    exp_meta = exp_stripped;

    // First pass through the file -> get the meta data
    std::cout << "Will analyze the metadata first to determine the number of SWATH windows and the window sizes." << std::endl;
    std::vector<int> swath_counter;
    int nr_ms1_spectra;
    std::vector<OpenSwath::SwathMap> known_window_boundaries;

    countScansInSwath_(exp_stripped->getSpectra(), swath_counter, nr_ms1_spectra, known_window_boundaries);
    std::cout << "Determined there to be " << swath_counter.size()
              << " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
    endProgress();

    std::shared_ptr<FullSwathFileConsumer> dataConsumer;
    startProgress(0, 1, "Loading data file " + file);
    if (readoptions == "normal")
    {
      dataConsumer = std::make_shared<RegularSwathFileConsumer>(known_window_boundaries);
      dataConsumer->setExperimentalSettings(*exp_meta.get());
    }
    else if (readoptions == "cache")
    {
      dataConsumer = std::make_shared<CachedSwathFileConsumer>(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
      dataConsumer->setExperimentalSettings(*exp_meta.get());
    }
    else if (readoptions == "split")
    {
      // WARNING: swath_maps will be empty when querying retrieveSwathMaps()
      dataConsumer = std::make_shared<MzMLSwathFileConsumer>(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
      dataConsumer->setExperimentalSettings(*exp_meta.get());
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Unknown or unsupported option " + readoptions);
    }

    std::vector<Interfaces::IMSDataConsumer *> consumer_list;
    // only use plugin if non-empty
    if (plugin_consumer)
    {
      exp_meta->setMetaValue("nr_ms1_spectra", nr_ms1_spectra); // required for SwathQC::getExpSettingsFunc()
      plugin_consumer->setExperimentalSettings(*exp_meta.get()); // set the meta data
      exp_meta->removeMetaValue("nr_ms1_spectra"); // served its need. remove
      // plugin_consumer->setExpectedSize(nr_ms1_spectra + accumulate(swath_counter)); // not needed currently
      consumer_list.push_back(plugin_consumer);
    }
    consumer_list.push_back(dataConsumer.get());
    MSDataChainingConsumer chaining_consumer(consumer_list);
    MzMLFile().transform(file, &chaining_consumer, false, true); // we do not need to reload metadata, it has already been loaded

    OPENMS_LOG_DEBUG << "Finished parsing Swath file " << std::endl;
    std::vector<OpenSwath::SwathMap> swath_maps;
    dataConsumer->retrieveSwathMaps(swath_maps);
    endProgress();
    return swath_maps;
  }

  /// Loads a Swath run from a single mzXML file
  std::vector<OpenSwath::SwathMap> SwathFile::loadMzXML(const String& file,
    const String& tmp,
    boost::shared_ptr<ExperimentalSettings>& exp_meta,
    const String& readoptions)
  {
    std::cout << "Loading mzXML file " << file << " using readoptions " << readoptions << std::endl;
    String tmp_fname = "openswath_tmpfile";

    startProgress(0, 1, "Loading metadata file " + file);
    boost::shared_ptr<PeakMap > experiment_metadata(new PeakMap);
    FileHandler f;
    f.getOptions().setAlwaysAppendData(true);
    f.getOptions().setFillData(false);
    f.loadExperiment(file, *experiment_metadata, {FileTypes::MZXML});
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
    startProgress(0, 1, "Loading data file " + file);
    if (readoptions == "normal")
    {
      dataConsumer = new RegularSwathFileConsumer(known_window_boundaries);
      MzXMLFile().transform(file, dataConsumer);
    }
    else if (readoptions == "cache")
    {
      dataConsumer = new CachedSwathFileConsumer(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
      MzXMLFile().transform(file, dataConsumer);
    }
    else if (readoptions == "split")
    {
      dataConsumer = new MzMLSwathFileConsumer(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
      MzXMLFile().transform(file, dataConsumer);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Unknown or unsupported option " + readoptions);
    }
    OPENMS_LOG_DEBUG << "Finished parsing Swath file " << std::endl;
    std::vector<OpenSwath::SwathMap> swath_maps;
    dataConsumer->retrieveSwathMaps(swath_maps);
    delete dataConsumer;

    endProgress();
    return swath_maps;
  }

  /// Loads a Swath run from a single sqMass file
  std::vector<OpenSwath::SwathMap> SwathFile::loadSqMass(const String& file, boost::shared_ptr<ExperimentalSettings>& /* exp_meta */)
  {
    startProgress(0, 1, "Loading sqmass data file " + file);

    OpenMS::Internal::MzMLSqliteSwathHandler sql_mass_reader(file);
    std::vector<OpenSwath::SwathMap> swath_maps = sql_mass_reader.readSwathWindows();
    for (Size k = 0; k < swath_maps.size(); k++)
    {
      std::vector<int> indices = sql_mass_reader.readSpectraForWindow(swath_maps[k]);
      OpenMS::Internal::MzMLSqliteHandler handler(file, 0);
      OpenSwath::SpectrumAccessPtr sptr(new OpenMS::SpectrumAccessSqMass(handler, indices));
      swath_maps[k].sptr = sptr;
    }

    // also store the MS1 map
    OpenSwath::SwathMap ms1_map;
    std::vector<int> indices = sql_mass_reader.readMS1Spectra();
    OpenMS::Internal::MzMLSqliteHandler handler(file, 0);
    OpenSwath::SpectrumAccessPtr sptr(new OpenMS::SpectrumAccessSqMass(handler, indices));
    ms1_map.sptr = sptr;
    ms1_map.ms1 = true;
    swath_maps.push_back(ms1_map);
    endProgress();

    std::cout << "Determined there to be " << swath_maps.size() <<
      " SWATH windows and in total " << indices.size() << " MS1 spectra" << std::endl;

    return swath_maps;
  }


  /// Cache a file to disk
  OpenSwath::SpectrumAccessPtr SwathFile::doCacheFile_(const String& in, const String& tmp, const String& tmp_fname,
    const boost::shared_ptr<PeakMap >& experiment_metadata)
  {
    String cached_file = tmp + tmp_fname + ".cached";
    String meta_file = tmp + tmp_fname;

    // Create new consumer, transform infile, write out metadata
    {
      MSDataCachedConsumer cachedConsumer(cached_file, true);
      MzMLFile().transform(in, &cachedConsumer, *experiment_metadata.get());
      Internal::CachedMzMLHandler().writeMetadata(*experiment_metadata.get(), meta_file, true);
    } // ensure that filestream gets closed

    boost::shared_ptr<PeakMap > exp(new PeakMap);
    FileHandler().loadExperiment(meta_file, *exp.get(), {FileTypes::MZML});
    return SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  }

  /// Only read the meta data from a file and use it to populate exp_meta
  boost::shared_ptr< PeakMap > SwathFile::populateMetaData_(const String& file)
  {
    boost::shared_ptr<PeakMap > experiment_metadata(new PeakMap);
    FileHandler f;
    f.getOptions().setAlwaysAppendData(true);
    f.getOptions().setFillData(false);
    f.loadExperiment(file, *experiment_metadata);
    return experiment_metadata;
  }
  /// Counts the number of scans in a full Swath file (e.g. concatenated non-split file)
  void SwathFile::countScansInSwath_(const std::vector<MSSpectrum>& exp,
                                     std::vector<int>& swath_counter, int& nr_ms1_spectra,
                                     std::vector<OpenSwath::SwathMap>& known_window_boundaries,
                                     double TOLERANCE)
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

          // set ion mobility if exists, otherwise will take default value of -1
          double imLower, imUpper;
          if (s.metaValueExists("ion mobility lower limit"))
          {
            imLower = s.getMetaValue("ion mobility lower limit"); // want this to be -1 if no ion mobility
            imUpper = s.getMetaValue("ion mobility upper limit");

          }
          else
          {
            imLower = -1;
            imUpper = -1;
          }

          const OpenSwath::SwathMap boundary(prec[0].getMZ() - prec[0].getIsolationWindowLowerOffset(), 
                                          prec[0].getMZ() + prec[0].getIsolationWindowUpperOffset(), 
                                          prec[0].getMZ(),
                                          imLower,
                                          imUpper,
                                          false);
          bool found = false;
          for (Size j = 0; j < known_window_boundaries.size(); j++)
          {
            // Check if the current scan is within the known window boundaries
            if (known_window_boundaries[j].isEqual(boundary, TOLERANCE))
            {
              found = true;
              swath_counter[j]++;
            }
          }
          if (!found)
          {
            // we found a new SWATH scan
            swath_counter.push_back(1);
            known_window_boundaries.push_back(boundary);

            OPENMS_LOG_DEBUG << "Adding Swath centered at " << boundary.center
              << " m/z with an isolation window of " << boundary.lower << " to " << boundary.upper
              << " m/z and IM start of " << boundary.imLower << " and IM end of " << boundary.imUpper << std::endl;
          }
        }
      }
    }
    nr_ms1_spectra = ms1_counter;

    std::cout << "Determined there to be " << swath_counter.size() <<
      " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
  }
}
