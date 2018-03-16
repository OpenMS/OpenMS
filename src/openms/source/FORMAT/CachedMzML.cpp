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

#include <OpenMS/FORMAT/CachedMzML.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

namespace OpenMS
{

  CachedmzML::CachedmzML()
  {
  }

  CachedmzML::~CachedmzML()
  {
  }

  CachedmzML& CachedmzML::operator=(const CachedmzML& rhs)
  {
    if (&rhs == this)
      return *this;

    spectra_index_ = rhs.spectra_index_;
    chrom_index_ = rhs.chrom_index_;

    return *this;
  }

  void CachedmzML::writeMemdump(MapType& exp, String out)
  {
    std::ofstream ofs(out.c_str(), std::ios::binary);
    Size exp_size = exp.size();
    Size chrom_size = exp.getChromatograms().size();
    int file_identifier = CACHED_MZML_FILE_IDENTIFIER;
    ofs.write((char*)&file_identifier, sizeof(file_identifier));

    startProgress(0, exp.size() + exp.getChromatograms().size(), "storing binary data");
    for (Size i = 0; i < exp.size(); i++)
    {
      setProgress(i);
      writeSpectrum_(exp[i], ofs);
    }

    for (Size i = 0; i < exp.getChromatograms().size(); i++)
    {
      setProgress(i);
      writeChromatogram_(exp.getChromatograms()[i], ofs);
    }

    ofs.write((char*)&exp_size, sizeof(exp_size));
    ofs.write((char*)&chrom_size, sizeof(chrom_size));
    ofs.close();
    endProgress();
  }

  void CachedmzML::readMemdump(MapType& exp_reading, String filename) const
  {
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if (ifs.fail())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    Size exp_size, chrom_size;
    Peak1D current_peak;

    int file_identifier;
    ifs.read((char*)&file_identifier, sizeof(file_identifier));
    if (file_identifier != CACHED_MZML_FILE_IDENTIFIER)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "File might not be a cached mzML file (wrong file magic number). Aborting!", filename);
    }

    ifs.seekg(0, ifs.end); // set file pointer to end
    ifs.seekg(ifs.tellg(), ifs.beg); // set file pointer to end, in forward direction
    ifs.seekg(- static_cast<int>(sizeof(exp_size) + sizeof(chrom_size)), ifs.cur); // move two fields to the left, start reading
    ifs.read((char*)&exp_size, sizeof(exp_size));
    ifs.read((char*)&chrom_size, sizeof(chrom_size));
    ifs.seekg(sizeof(file_identifier), ifs.beg); // set file pointer to beginning (after identifier), start reading

    exp_reading.reserve(exp_size);
    startProgress(0, exp_size + chrom_size, "reading binary data");
    for (Size i = 0; i < exp_size; i++)
    {
      setProgress(i);
      SpectrumType spectrum;
      readSpectrum_(spectrum, ifs);
      exp_reading.addSpectrum(spectrum);
    }
    std::vector<ChromatogramType> chromatograms;
    for (Size i = 0; i < chrom_size; i++)
    {
      setProgress(i);
      ChromatogramType chromatogram;
      readChromatogram_(chromatogram, ifs);
      chromatograms.push_back(chromatogram);
    }
    exp_reading.setChromatograms(chromatograms);

    ifs.close();
    endProgress();
  }

  const std::vector<std::streampos>& CachedmzML::getSpectraIndex() const
  {
    return spectra_index_;
  }

  const std::vector<std::streampos>& CachedmzML::getChromatogramIndex() const
  {
    return chrom_index_;
  }

  void CachedmzML::createMemdumpIndex(String filename)
  {
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if (ifs.fail())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    Size exp_size, chrom_size;
    Peak1D current_peak;

    ifs.seekg(0, ifs.beg); // set file pointer to beginning, start reading
    spectra_index_.clear();
    chrom_index_.clear();
    int file_identifier;
    int extra_offset = sizeof(dbl_field_) + sizeof(int_field_);
    int chrom_offset = 0;

    ifs.read((char*)&file_identifier, sizeof(file_identifier));
    if (file_identifier != CACHED_MZML_FILE_IDENTIFIER)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "File might not be a cached mzML file (wrong file magic number). Aborting!", filename);
    }

    // For spectra and chromatograms go through file, read the size of the
    // spectrum/chromatogram and record the starting index of the element, then
    // skip ahead to the next spectrum/chromatogram.

    ifs.seekg(0, ifs.end); // set file pointer to end
    ifs.seekg(ifs.tellg(), ifs.beg); // set file pointer to end, in forward direction
    ifs.seekg(- static_cast<int>(sizeof(exp_size) + sizeof(chrom_size)), ifs.cur); // move two fields to the left, start reading
    ifs.read((char*)&exp_size, sizeof(exp_size));
    ifs.read((char*)&chrom_size, sizeof(chrom_size));
    ifs.seekg(sizeof(file_identifier), ifs.beg); // set file pointer to beginning (after identifier), start reading

    startProgress(0, exp_size + chrom_size, "Creating index for binary spectra");
    for (Size i = 0; i < exp_size; i++)
    {
      setProgress(i);

      Size spec_size;
      Size float_arr;
      spectra_index_.push_back(ifs.tellg());
      ifs.read((char*)&spec_size, sizeof(spec_size));
      ifs.read((char*)&float_arr, sizeof(float_arr));
      ifs.seekg(extra_offset + (sizeof(DatumSingleton)) * 2 * (spec_size), ifs.cur);

      // Read the extra data arrays
      for (Size k = 0; k < float_arr; k++)
      {
        Size len, len_name;
        ifs.read((char*)&len, sizeof(len));
        ifs.read((char*)&len_name, sizeof(len_name));
        ifs.seekg(len_name * sizeof(char), ifs.cur);
        ifs.seekg(sizeof(DatumSingleton) * len, ifs.cur);
      }
    }

    for (Size i = 0; i < chrom_size; i++)
    {
      setProgress(i);

      Size ch_size;
      Size float_arr;
      chrom_index_.push_back(ifs.tellg());
      ifs.read((char*)&ch_size, sizeof(ch_size));
      ifs.read((char*)&float_arr, sizeof(float_arr));
      ifs.seekg(chrom_offset + (sizeof(DatumSingleton)) * 2 * (ch_size), ifs.cur);

      // Read the extra data arrays
      for (Size k = 0; k < float_arr; k++)
      {
        Size len, len_name;
        ifs.read((char*)&len, sizeof(len));
        ifs.read((char*)&len_name, sizeof(len_name));
        ifs.seekg(len_name * sizeof(char), ifs.cur);
        ifs.seekg(sizeof(DatumSingleton) * len, ifs.cur);
      }
    }

    ifs.close();
    endProgress();
  }

  void CachedmzML::writeMetadata(MapType exp, String out_meta, bool addCacheMetaValue)
  {
    // delete the actual data for all spectra and chromatograms, leave only metadata
    // TODO : remove copy
    std::vector<MSChromatogram > chromatograms = exp.getChromatograms(); // copy
    for (Size i = 0; i < exp.size(); i++)
    {
      exp[i].clear(false);
    }
    for (Size i = 0; i < exp.getChromatograms().size(); i++)
    {
      chromatograms[i].clear(false);
    }
    exp.setChromatograms(chromatograms);

    if (addCacheMetaValue)
    {
      // set dataprocessing on each spectrum/chromatogram
      boost::shared_ptr< DataProcessing > dp = boost::shared_ptr< DataProcessing >(new DataProcessing);
      std::set<DataProcessing::ProcessingAction> actions;
      actions.insert(DataProcessing::FORMAT_CONVERSION);
      dp->setProcessingActions(actions);
      dp->setMetaValue("cached_data", "true");
      for (Size i=0; i<exp.size(); ++i)
      {
        exp[i].getDataProcessing().push_back(dp);
      }
      std::vector<MSChromatogram > l_chromatograms = exp.getChromatograms();
      for (Size i=0; i<l_chromatograms.size(); ++i)
      {
        l_chromatograms[i].getDataProcessing().push_back(dp);
      }
      exp.setChromatograms(l_chromatograms);
    }

    // store the meta data using the regular MzMLFile
    MzMLFile().store(out_meta, exp);
  }

  void CachedmzML::readSpectrum_(std::vector<Datavector>& data, std::ifstream& ifs, int& ms_level, double& rt) const
  {
    Size spec_size = -1;
    Size extra_arrays = -1;
    ifs.read((char*)&spec_size, sizeof(spec_size));
    ifs.read((char*)&extra_arrays, sizeof(extra_arrays));
    ifs.read((char*)&ms_level, sizeof(ms_level));
    ifs.read((char*)&rt, sizeof(rt));

    data.resize(2);
    data[0].resize(spec_size);
    data[1].resize(spec_size);

    if (spec_size > 0)
    {
      ifs.read((char*)&data[0][0], spec_size * sizeof(DatumSingleton));
      ifs.read((char*)&data[1][0], spec_size * sizeof(DatumSingleton));
    }

    for (Size k = 0; k < extra_arrays; k++)
    {
      std::cout << " read spectrum_  : skip over " << std::endl;
      Size len, len_name;
      ifs.read((char*)&len, sizeof(len));
      ifs.read((char*)&len_name, sizeof(len_name));
      ifs.seekg(len_name * sizeof(char), ifs.cur);
      data.resize( data.size() + 1);
      data.back().resize(spec_size);
      ifs.read((char*)&data.back()[0], spec_size * sizeof(DatumSingleton));
    }
  }

  void CachedmzML::readChromatogram_(std::vector<Datavector>& data, std::ifstream& ifs) const
  {
    Size chrom_size = -1;
    Size extra_arrays = -1;
    ifs.read((char*)&chrom_size, sizeof(chrom_size));
    ifs.read((char*)&extra_arrays, sizeof(extra_arrays));
    data.resize(2);
    data[0].resize(chrom_size);
    data[1].resize(chrom_size);

    if (chrom_size > 0)
    {
      ifs.read((char*)&data[0][0], chrom_size * sizeof(DatumSingleton));
      ifs.read((char*)&data[1][0], chrom_size * sizeof(DatumSingleton));
    }

    for (Size k = 0; k < extra_arrays; k++)
    {
      data.resize( data.size() + 1);
      data.back().resize(chrom_size);
      ifs.read((char*)&data.back()[0], chrom_size * sizeof(DatumSingleton));
    }
  }

  void CachedmzML::readSpectrum_(SpectrumType& spectrum, std::ifstream& ifs) const
  {
    std::vector<Datavector> data;
    int ms_level;
    double rt;
    readSpectrum_(data, ifs, ms_level, rt);
    spectrum.reserve(data[0].size());
    spectrum.setMSLevel(ms_level);
    spectrum.setRT(rt);

    for (Size j = 0; j < data[0].size(); j++)
    {
      Peak1D p;
      p.setMZ(data[0][j]);
      p.setIntensity(data[1][j]);
      spectrum.push_back(p);
    }

  }

  void CachedmzML::readChromatogram_(ChromatogramType& chromatogram, std::ifstream& ifs) const
  {
    std::vector<Datavector> data;
    readChromatogram_(data, ifs);
    chromatogram.reserve(data[0].size());

    for (Size j = 0; j < data[0].size(); j++)
    {
      ChromatogramPeak p;
      p.setRT(data[0][j]);
      p.setIntensity(data[1][j]);
      chromatogram.push_back(p);
    }

  }

  void CachedmzML::writeSpectrum_(const SpectrumType& spectrum, std::ofstream& ofs)
  {
    Size exp_size = spectrum.size();
    ofs.write((char*)&exp_size, sizeof(exp_size));
    Size arr_s = spectrum.getFloatDataArrays().size() + spectrum.getIntegerDataArrays().size();
    ofs.write((char*)&arr_s, sizeof(arr_s));
    int_field_ = spectrum.getMSLevel();
    ofs.write((char*)&int_field_, sizeof(int_field_));
    dbl_field_ = spectrum.getRT();
    ofs.write((char*)&dbl_field_, sizeof(dbl_field_));

    // Catch empty spectrum: we do not write any data and since the "size" we
    // just wrote is zero, no data will be read
    if (spectrum.empty())
    {
      return;
    }

    Datavector mz_data;
    Datavector int_data;
    mz_data.reserve(spectrum.size());
    int_data.reserve(spectrum.size());
    for (Size j = 0; j < spectrum.size(); j++)
    {
      mz_data.push_back(spectrum[j].getMZ());
      int_data.push_back(static_cast<double>(spectrum[j].getIntensity()));
    }

    ofs.write((char*)&mz_data.front(), mz_data.size() * sizeof(mz_data.front()));
    ofs.write((char*)&int_data.front(), int_data.size() * sizeof(int_data.front()));

    for (const auto& fda : spectrum.getFloatDataArrays() )
    {
      std::cout << " got float array " << fda.getName() << std::endl;
      Size len = fda.size();
      ofs.write((char*)&len, sizeof(len));
      Size len_name = fda.getName().size();
      ofs.write((char*)&len_name, sizeof(len_name));
      ofs.write((char*)&fda.getName().front(), len_name * sizeof(fda.getName().front()));
      // now go to the actual data
      Datavector tmp;
      tmp.reserve(fda.size());
      for (const auto& val : fda) {tmp.push_back(val);}
      ofs.write((char*)&tmp.front(), tmp.size() * sizeof(tmp.front()));
    }
    for (const auto& ida : spectrum.getIntegerDataArrays() )
    {
      Size len = ida.size();
      ofs.write((char*)&len, sizeof(len));
      Size len_name = ida.getName().size();
      ofs.write((char*)&len_name, sizeof(len_name));
      ofs.write((char*)&ida.getName().front(), len_name * sizeof(ida.getName().front()));
      // now go to the actual data
      Datavector tmp;
      tmp.reserve(ida.size());
      for (const auto& val : ida) {tmp.push_back(val);}
      ofs.write((char*)&tmp.front(), tmp.size() * sizeof(tmp.front()));
    }
  }

  void CachedmzML::writeChromatogram_(const ChromatogramType& chromatogram, std::ofstream& ofs)
  {
    Size exp_size = chromatogram.size();
    ofs.write((char*)&exp_size, sizeof(exp_size));
    Size arr_s = chromatogram.getFloatDataArrays().size() + chromatogram.getIntegerDataArrays().size();
    ofs.write((char*)&arr_s, sizeof(arr_s));

    // Catch empty chromatogram: we do not write any data and since the "size" we
    // just wrote is zero, no data will be read
    if (chromatogram.empty())
    {
      return;
    }

    Datavector rt_data;
    Datavector int_data;
    rt_data.reserve(chromatogram.size());
    int_data.reserve(chromatogram.size());
    for (Size j = 0; j < chromatogram.size(); j++)
    {
      rt_data.push_back(chromatogram[j].getRT());
      int_data.push_back(chromatogram[j].getIntensity());
    }
    ofs.write((char*)&rt_data.front(), rt_data.size() * sizeof(rt_data.front()));
    ofs.write((char*)&int_data.front(), int_data.size() * sizeof(int_data.front()));

    for (const auto& fda : chromatogram.getFloatDataArrays() )
    {
      std::cout << " got float array " << fda.getName() << std::endl;
      Size len = fda.size();
      ofs.write((char*)&len, sizeof(len));
      Size len_name = fda.getName().size();
      ofs.write((char*)&len_name, sizeof(len_name));
      ofs.write((char*)&fda.getName().front(), len_name * sizeof(fda.getName().front()));
      std::cout << " written  " << len_name * sizeof(fda.getName().front()) << std::endl;
      Datavector tmp;
      tmp.reserve(fda.size());
      for (const auto& val : fda) {tmp.push_back(val);}
      ofs.write((char*)&tmp.front(), tmp.size() * sizeof(tmp.front()));
    }
    for (const auto& ida : chromatogram.getIntegerDataArrays() )
    {
      Size len = ida.size();
      ofs.write((char*)&len, sizeof(len));
      Size len_name = ida.getName().size();
      ofs.write((char*)&len_name, sizeof(len_name));
      ofs.write((char*)&ida.getName().front(), len_name * sizeof(ida.getName().front()));
      // now go to the actual data
      Datavector tmp;
      tmp.reserve(ida.size());
      for (const auto& val : ida) {tmp.push_back(val);}
      ofs.write((char*)&tmp.front(), tmp.size() * sizeof(tmp.front()));
    }
  }

}

