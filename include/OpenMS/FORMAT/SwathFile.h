// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#ifdef OPENMS_FORMAT_SWATHFILE_MZXMLSUPPORT
#include <OpenMS/FORMAT/MzXMLFile.h>
#endif

#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>

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
    std::vector< OpenSwath::SwathMap > loadSplit(StringList file_list, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal")
    {
      int progress = 0;
      startProgress(0, file_list.size(), "Loading data");

      std::vector< OpenSwath::SwathMap > swath_maps;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(file_list.size()); ++i)
      {
        std::cout << "Loading file " << file_list[i] << std::endl;
        String tmp_fname = "openswath_tmpfile_" + String(i) + ".mzML";

        boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
        OpenSwath::SpectrumAccessPtr spectra_ptr;

        // Populate meta-data
        if (i == 0) { populateMetaData_(file_list[i], exp_meta); }

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
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
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
#pragma omp critical (load)
#endif
        {
          swath_maps.push_back( swath_map );
          setProgress(progress++);
        }
      }
      endProgress();
      return swath_maps;
    }

    /// Loads a Swath run from a single mzML file
    std::vector< OpenSwath::SwathMap > loadMzML(String file, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal")
    {
      startProgress(0, 1, "Loading data file " + file);
      std::vector< OpenSwath::SwathMap > swath_maps;
      FullSwathFileConsumer * dataConsumer;
      String tmp_fname = "openswath_tmpfile";

      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
      OpenSwath::SpectrumAccessPtr spectra_ptr;
      OpenSwath::SwathMap swath_map;

      populateMetaData_(file, exp_meta); 

      if (readoptions == "normal")
      {
        dataConsumer = new RegularSwathFileConsumer();
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else if (readoptions == "cache")
      {

        std::cout << "Will analyze the metadata first to determine the number of SWATH windows and the window sizes." << std::endl;
        boost::shared_ptr<MSExperiment<Peak1D> > experiment_metadata(new MSExperiment<Peak1D>);
        // First pass through the file -> get the meta data
        {
          MzMLFile f;
          f.getOptions().setAlwaysAppendData(true);
          f.getOptions().setFillData(false);
          f.load(file, *experiment_metadata);
        }

        std::vector<int> swath_counter;
        int nr_ms1_spectra;
        countScansInSwath_(experiment_metadata->getSpectra(), swath_counter, nr_ms1_spectra);

        std::cout << "Determined there to be " << swath_counter.size() << 
          " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
        dataConsumer = new CachedSwathFileConsumer(tmp, tmp_fname, nr_ms1_spectra, swath_counter);
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Unknown or unsupported option " + readoptions);
      }
      dataConsumer->retrieveSwathMaps(swath_maps);
      delete dataConsumer;

      endProgress();
      return swath_maps;
    }

#ifndef OPENMS_FORMAT_SWATHFILE_MZXMLSUPPORT
    /// Loads a Swath run from a single mzXML file (currently not supported)
    std::vector< OpenSwath::SwathMap > loadMzXML(String, String, 
      boost::shared_ptr<ExperimentalSettings>&, String)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "MzXML not supported");
    }
#else
    std::vector< OpenSwath::SwathMap > loadMzXML(String file, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal")
    {
      startProgress(0, 1, "Loading data file " + file);
      std::vector< OpenSwath::SwathMap > swath_maps;
      boost::shared_ptr<FullSwathFileConsumer> dataConsumer;
      String tmp_fname = "openswath_tmpfile";

      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
      OpenSwath::SpectrumAccessPtr spectra_ptr;
      OpenSwath::SwathMap swath_map;

      if (readoptions == "normal")
      {
        dataConsumer = boost::shared_ptr<RegularSwathFileConsumer>(new RegularSwathFileConsumer());
        Internal::MSMzXMLDataReader<FullSwathFileConsumer> datareader;
        datareader.setConsumer(dataConsumer);
        MzXMLFile().load(file, datareader);
        *exp_meta = datareader;
      }
      else if (readoptions == "cache")
      {

        // First pass through the file -> get the meta data
        std::cout << "Will analyze the metadata first to determine the number of SWATH windows and the window sizes." << std::endl;
        boost::shared_ptr<MSExperiment<Peak1D> > experiment_metadata(new MSExperiment<Peak1D>);
        std::vector<int> swath_counter;
        int nr_ms1_spectra;
        {
          boost::shared_ptr<MSDataTransformingConsumer> noopConsumer = 
            boost::shared_ptr<MSDataTransformingConsumer>( new MSDataTransformingConsumer() ) ; 
          Internal::MSMzXMLDataReader<MSDataTransformingConsumer> datareader;
          datareader.setConsumer(noopConsumer);
          MzXMLFile f;
          f.getOptions().setFillData(false);
          f.load(file, datareader);
          countScansInSwath_(datareader.getRealSpectra(), swath_counter, nr_ms1_spectra);
          *exp_meta = datareader;
        }

        std::cout << "Determined there to be " << swath_counter.size() << 
          " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
        dataConsumer = boost::shared_ptr<CachedSwathFileConsumer>( 
            new CachedSwathFileConsumer(tmp, tmp_fname, nr_ms1_spectra, swath_counter) ) ; 
        Internal::MSMzXMLDataReader<FullSwathFileConsumer> datareader;
        datareader.setConsumer(dataConsumer);
        MzXMLFile().load(file, datareader);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Unknown or unsupported option " + readoptions);
      }
      dataConsumer->retrieveSwathMaps(swath_maps);

      endProgress();
      return swath_maps;
    }
#endif

  protected:

    /// Cache a file to disk
    OpenSwath::SpectrumAccessPtr doCacheFile_(String in, String tmp, String tmp_fname, 
        boost::shared_ptr<MSExperiment<Peak1D> > experiment_metadata )
    {
      String cached_file = tmp + tmp_fname + ".cached";
      String meta_file = tmp + tmp_fname;

      // Create new consumer, transform infile, write out metadata
      CachedMzMLConsumer * cachedConsumer = new CachedMzMLConsumer(cached_file, true);
      MzMLFile().transform(in, cachedConsumer, *experiment_metadata.get());
      CachedmzML().writeMetadata(*experiment_metadata.get(), meta_file, true);
      delete cachedConsumer; // ensure that filestream gets closed

      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
      MzMLFile().load(meta_file, *exp.get());
      return SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  }

    /// Only read the meta data from a file and use it to populate exp_meta
    void populateMetaData_(String file, boost::shared_ptr<ExperimentalSettings>& exp_meta)
    {
      MSExperiment<Peak1D> tmp;
      MSDataTransformingConsumer c;
      MzMLFile().transform(file, &c, tmp);
      *exp_meta = tmp;
    }

    /// Counts the number of scans in a full Swath file (e.g. concatenated non-split file)
    void countScansInSwath_(const std::vector<MSSpectrum<> > exp, 
        std::vector<int> & swath_counter, int & nr_ms1_spectra)
    {
      int ms1_counter = 0;
      int ms2_counter = 0;
      for (Size i = 0; i < exp.size(); i++)
      {
        const MSSpectrum<> & s = exp[i];
        {
          if (s.getMSLevel() == 1)
          {
            ms2_counter = 0;
            ms1_counter++;
          }
          else 
          {
            if (ms2_counter == (int)swath_counter.size())
            {
              swath_counter.push_back(0);
            }
            swath_counter[ms2_counter]++;
            ms2_counter++;
          }
        }
      }
      nr_ms1_spectra = ms1_counter;
    }

  };
}

#endif

