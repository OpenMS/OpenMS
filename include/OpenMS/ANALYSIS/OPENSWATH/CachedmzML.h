// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_CACHEDMZML_H
#define OPENMS_ANALYSIS_OPENSWATH_CACHEDMZML_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h"

#include <fstream>

#define MAGIC_NUMBER 8093

namespace OpenMS 
{
  /**
    @brief An class that uses on-disk caching to read and write spectra and chromatograms

    This class implements the OpenSWATH Spectrum Access interface
    (ISpectrumAccess) using the CachedmzML class which is able to read and
    write a cached mzML file.

  */
  class OPENMS_DLLAPI CachedmzML
    : public ProgressLogger
  {
    int int_field;
    double dbl_field;

  public:

    typedef MSExperiment<Peak1D> MapType;
    typedef MSSpectrum<Peak1D> SpectrumType;
    typedef MSChromatogram<ChromatogramPeak> ChromatogramType;
#if 1
    // double means twice the file size
    typedef double DatumSingleton;
#else
    // float means half the file size
    typedef float DatumSingleton;
#endif
    typedef std::vector< DatumSingleton> Datavector;

    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor
    CachedmzML()
    {
    }

    /// Default destructor
    ~CachedmzML()
    {
    }

    /// Assignment operator
    CachedmzML & operator=(const CachedmzML & rhs)
    {
      if (&rhs == this)
        return *this;

      spectra_index = rhs.spectra_index;
      chrom_index = rhs.chrom_index;

      return *this;
    }
    //@}

    /** @name Read / Write an MSExperiment 
    */
    //@{
    /// Write complete spectra as a dump to the disk
    void writeMemdump(MapType& exp, String out)
    {
      std::ofstream ofs(out.c_str(), std::ios::binary);
      Size exp_size = exp.size();
      Size chrom_size = exp.getChromatograms().size();
      int magic_number = MAGIC_NUMBER;
      ofs.write( (char *)&magic_number, sizeof(magic_number));
      ofs.write( (char *)&exp_size , sizeof(exp_size));
      ofs.write( (char *)&chrom_size , sizeof(chrom_size));

      startProgress(0, exp.size() + exp.getChromatograms().size(), "storing binary spectra");
      for(Size i =0; i<exp.size(); i++)
      {
        setProgress(i);
        write_spectrum(exp[i], ofs);
      }

      for(Size i =0; i<exp.getChromatograms().size(); i++)
      {
        setProgress(i);
        write_chromatogram(exp.getChromatograms()[i], ofs);
      }

      ofs.close();
      endProgress();
    }

    /// Read all spectra from a dump from the disk
    void readMemdump(MapType& exp_reading, String filename) const
    {
      std::ifstream ifs(filename.c_str(), std::ios::binary);
      Size exp_size, chrom_size;
      Peak1D current_peak;
      std::string identifier;

      int magic_number;
      ifs.read( (char *)&magic_number , sizeof(magic_number));
      if(magic_number != MAGIC_NUMBER)
      {
        throw "wrong file, does not start with MAGIC_NUMBER";
      }

      ifs.read( (char *)&exp_size , sizeof(exp_size));
      ifs.read( (char *)&chrom_size , sizeof(chrom_size));

      exp_reading.reserve(exp_size);
      startProgress(0, exp_size + chrom_size, "reading binary spectra");
      for(Size i =0; i<exp_size; i++)
      {
        setProgress(i);
        SpectrumType spectrum;
        _read_Spectrum(spectrum, ifs);
        exp_reading.push_back(spectrum);
      }
       std::vector < ChromatogramType > chromatograms;
      for(Size i =0; i<chrom_size; i++)
      {
        setProgress(i);
        ChromatogramType chromatogram;
        _read_Chromatogram(chromatogram, ifs);
        chromatograms.push_back(chromatogram);
      }
      exp_reading.setChromatograms(chromatograms);

      ifs.close();
      endProgress();
    }
    //@}

    /** @name Read a single MSSpectrum
    */
    //@{
    /// Read a single spectrum from the given filename
    void readSingleSpectrum(MSSpectrum<Peak1D>& spectrum, const String & filename, const Size & idx) const
    {
      // open stream, read
      std::ifstream ifs(filename.c_str(), std::ios::binary);
      readSingleSpectrum(spectrum, ifs, idx);
    }

    // Read a single spectrum from the given filestream
    void readSingleSpectrum(MSSpectrum<Peak1D>& spectrum, std::ifstream& ifs, const Size & idx) const
    {
      // go to the specified index
      ifs.seekg(idx);
      _read_Spectrum(spectrum, ifs);
    }
    //@}

    /** @name Access to the binary indices
    */
    //@{
    const std::vector< Size >& getSpectraIndex() const
    {
      return spectra_index;
    }

    const std::vector< Size >& getChromatogramIndex() const
    {
      return chrom_index;
    }
    //@}

    /// Create an index on the location of all the spectra and chromatograms
    void createMemdumpIndex(String filename) 
    {
      std::ifstream ifs(filename.c_str(), std::ios::binary);
      Size exp_size, chrom_size;
      Peak1D current_peak;

      std::string identifier;

      spectra_index.clear();
      chrom_index.clear();
      int magic_number;
      int extra_offset = sizeof(dbl_field) + sizeof(int_field);
      int chrom_offset = 0;

      ifs.read( (char *)&magic_number , sizeof(magic_number));
      if(magic_number != MAGIC_NUMBER)
      {
        throw "wrong file, does not start with MAGIC_NUMBER";
      }

      // For spectra and chromatograms go through file, read the size of the
      // spectrum/chromatogram and record the starting index of the element, then
      // skip ahead to the next spectrum/chromatogram.
      ifs.read( (char *)&exp_size , sizeof(exp_size));
      ifs.read( (char *)&chrom_size , sizeof(chrom_size));
      startProgress(0, exp_size + chrom_size, "Creating index for binary spectra");
      for(Size i =0; i<exp_size; i++)
      {
        setProgress(i);

        Size spec_size;
        spectra_index.push_back(ifs.tellg());
        ifs.read( (char *)&spec_size , sizeof(spec_size));
        ifs.seekg( (int)ifs.tellg() + extra_offset + (sizeof(DatumSingleton))*2*(spec_size));

      }

      for(Size i =0; i<chrom_size; i++)
      {
        setProgress(i);

        Size chrom_size;
        chrom_index.push_back(ifs.tellg());
        ifs.read( (char *)&chrom_size , sizeof(chrom_size));
        ifs.seekg( (int)ifs.tellg() + chrom_offset + (sizeof(DatumSingleton))*2*(chrom_size));

      }

      ifs.close();
      endProgress();
    }

    /// Write only the meta data of an MSExperiment
    void writeMetadata(MapType exp, String out_meta)
    {
      // delete the actual data for all spectra and chromatograms, leave only metadata
      std::vector < MSChromatogram < ChromatogramPeak> > chromatograms = exp.getChromatograms(); // copy
      for(Size i =0; i<exp.size(); i++)
      {
        exp[i].clear(false);
      }
      for(Size i =0; i<exp.getChromatograms().size(); i++)
      {
        // delete the actual data, leave only metadata
        //exp.getChromatograms()[i].clear(false);
        chromatograms[i].clear(false);
      }
      exp.setChromatograms(chromatograms);

      // store the meta data that is left in out_meta file
      MzMLFile f;
      f.store(out_meta,exp);
    }

    /// fast access without copying
    static inline void readSpectrumFast(OpenSwath::BinaryDataArrayPtr data1,
        OpenSwath::BinaryDataArrayPtr data2, std::ifstream & ifs, int ms_level,
        double rt)
    {
      Size spec_size = -1;
      ifs.read((char *) &spec_size, sizeof(spec_size));
      ifs.read((char *) &ms_level, sizeof(ms_level));
      ifs.read((char *) &rt, sizeof(rt));

      data1->data.resize(spec_size);
      data2->data.resize(spec_size);
      ifs.read((char *) &(data1->data)[0], spec_size * sizeof(double));
      ifs.read((char *) &(data2->data)[0], spec_size * sizeof(double));
    }

    /// fast access without copying
    static inline void readChromatogramFast(OpenSwath::BinaryDataArrayPtr data1,
        OpenSwath::BinaryDataArrayPtr data2, std::ifstream & ifs)
    {
      Size spec_size = -1;
      ifs.read((char *) &spec_size, sizeof(spec_size));
      data1->data.resize(spec_size);

      data2->data.resize(spec_size);
      ifs.read((char *) &(data1->data)[0], spec_size * sizeof(double));
      ifs.read((char *) &(data2->data)[0], spec_size * sizeof(double));
    }

  private:

    // read a single spectrum directly into a datavector (assuming file is already at the correct position)
    void _read_Spectrum(Datavector & data1, Datavector & data2, std::ifstream & ifs, int & ms_level, double & rt) const
    {
        Size spec_size = -1;
        ifs.read( (char *)&spec_size , sizeof(spec_size));
        ifs.read( (char *)&ms_level , sizeof(ms_level));
        ifs.read( (char *)&rt, sizeof(rt));

        data1.resize(spec_size);
        data2.resize(spec_size);
        ifs.read( (char *)&data1[0], spec_size * sizeof(DatumSingleton));
        ifs.read( (char *)&data2[0], spec_size * sizeof(DatumSingleton));
    }

    // read a single chromatogram directly into a datavector (assuming file is already at the correct position)
    void _read_Chromatogram(Datavector & data1, Datavector & data2, std::ifstream & ifs) const
    {
        Size spec_size = -1;
        ifs.read( (char *)&spec_size , sizeof(spec_size));
        data1.resize(spec_size);
        data2.resize(spec_size);
        ifs.read( (char *)&data1[0], spec_size * sizeof(DatumSingleton));
        ifs.read( (char *)&data2[0], spec_size * sizeof(DatumSingleton));
    }

    // read a single spectrum directly into an OpenMS MSSpectrum (assuming file is already at the correct position)
    void _read_Spectrum(SpectrumType & spectrum, std::ifstream & ifs) const
    {
        Datavector mz_data;
        Datavector int_data;

        int ms_level;
        double rt;
        _read_Spectrum(mz_data, int_data, ifs, ms_level, rt);
        spectrum.reserve(mz_data.size());
        spectrum.setMSLevel(ms_level);
        spectrum.setRT(rt);

        for(Size j=0; j<mz_data.size(); j++)
        {
          Peak1D p;
          p.setMZ(mz_data[j]);
          p.setIntensity(int_data[j]);
          spectrum.push_back(p);
        }

    }

    // read a single chromatogram directly into an OpenMS MSChromatograms (assuming file is already at the correct position)
    void _read_Chromatogram(ChromatogramType & chromatogram, std::ifstream & ifs) const
    {
        Datavector rt_data;
        Datavector int_data;
        _read_Chromatogram(rt_data, int_data, ifs);
        chromatogram.reserve(rt_data.size());

        for(Size j=0; j<rt_data.size(); j++)
        {
          ChromatogramPeak p;
          p.setRT(rt_data[j]);
          p.setIntensity(int_data[j]);
          chromatogram.push_back(p);
        }

    }

    // write a single spectrum to filestream
    void write_spectrum(SpectrumType& spectrum, std::ofstream& ofs) 
    {
        Size exp_size = spectrum.size();
        ofs.write( (char *)&exp_size , sizeof(exp_size));
        int_field = spectrum.getMSLevel();
        ofs.write( (char *)&int_field , sizeof(int_field));
        dbl_field = spectrum.getRT();
        ofs.write( (char *)&dbl_field , sizeof(dbl_field));

#if 0
        ofs.write( (char *)&exp[i].front() , exp[i].size()*sizeof(exp[i].front()));
        std::cout << " storing spectrum " << i << " with size " << exp[i].size() << std::endl;
#else
        Datavector mz_data;
        Datavector int_data;
        for(Size j=0; j<spectrum.size(); j++)
        {
          mz_data.push_back(spectrum[j].getMZ());
          int_data.push_back(spectrum[j].getIntensity());
        }
        ofs.write( (char *)&mz_data.front() , mz_data.size()*sizeof(mz_data.front()));
        ofs.write( (char *)&int_data.front() , int_data.size()*sizeof(int_data.front()));
#endif
        //std::cout << exp[i] << std::endl;
    }

    // write a single chromatogram to filestream
    void write_chromatogram(const ChromatogramType& chromatogram, std::ofstream& ofs)
    {
        Size exp_size = chromatogram.size();
        ofs.write( (char *)&exp_size , sizeof(exp_size));
        Datavector rt_data;
        Datavector int_data;
        for(Size j=0; j< chromatogram.size(); j++)
        {
          rt_data.push_back(chromatogram[j].getRT());
          int_data.push_back(chromatogram[j].getIntensity());
        }
        ofs.write( (char *)&rt_data.front() , rt_data.size()*sizeof(rt_data.front()));
        ofs.write( (char *)&int_data.front() , int_data.size()*sizeof(int_data.front()));
    }

    std::vector< Size > spectra_index;
    std::vector< Size > chrom_index;

  };
}
#endif
