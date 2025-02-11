// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

namespace OpenMS::Internal
{
  CachedMzMLHandler::CachedMzMLHandler() = default;

  CachedMzMLHandler::~CachedMzMLHandler() = default;

  CachedMzMLHandler& CachedMzMLHandler::operator=(const CachedMzMLHandler& rhs)
  {
    if (&rhs == this)
    {
      return *this;
    }
    spectra_index_ = rhs.spectra_index_;
    chrom_index_ = rhs.chrom_index_;

    return *this;
  }

  void CachedMzMLHandler::writeMemdump(const MapType& exp, const String& out) const
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

  void CachedMzMLHandler::readMemdump(MapType& exp_reading, const String& filename) const
  {
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if (ifs.fail())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    Size exp_size, chrom_size;

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
      readSpectrum(spectrum, ifs);
      exp_reading.addSpectrum(spectrum);
    }
    std::vector<ChromatogramType> chromatograms;
    for (Size i = 0; i < chrom_size; i++)
    {
      setProgress(i);
      ChromatogramType chromatogram;
      readChromatogram(chromatogram, ifs);
      chromatograms.push_back(chromatogram);
    }
    exp_reading.setChromatograms(chromatograms);

    ifs.close();
    endProgress();
  }

  const std::vector<std::streampos>& CachedMzMLHandler::getSpectraIndex() const
  {
    return spectra_index_;
  }

  const std::vector<std::streampos>& CachedMzMLHandler::getChromatogramIndex() const
  {
    return chrom_index_;
  }

  void CachedMzMLHandler::createMemdumpIndex(const String& filename)
  {
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if (ifs.fail())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    Size exp_size, chrom_size;

    ifs.seekg(0, ifs.beg); // set file pointer to beginning, start reading
    spectra_index_.clear();
    chrom_index_.clear();
    int file_identifier;
    int extra_offset = sizeof(DoubleType) + sizeof(IntType);
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

  void CachedMzMLHandler::writeMetadata(MapType exp, const String& out_meta, bool addCacheMetaValue)
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

  void CachedMzMLHandler::writeMetadata_x(const MapType& exp, const String& out_meta, const bool addCacheMetaValue)
  {
    // delete the actual data for all spectra and chromatograms, leave only metadata
    // TODO : remove copy
    const ExperimentalSettings& qq = exp;
    MSExperiment out_exp;
    out_exp = qq;
    // std::vector<MSChromatogram > chromatograms = exp.getChromatograms(); // copy
    for (const auto& s: exp)
    {
      out_exp.addSpectrum(s);
      out_exp.getSpectra().back().clear(false);
    }
    for (const auto& c: exp.getChromatograms())
    {
      out_exp.addChromatogram(c);
      out_exp.getChromatograms().back().clear(false);
    }

    if (addCacheMetaValue)
    {
      // set dataprocessing on each spectrum/chromatogram
      boost::shared_ptr< DataProcessing > dp = boost::shared_ptr< DataProcessing >(new DataProcessing);
      std::set<DataProcessing::ProcessingAction> actions;
      actions.insert(DataProcessing::FORMAT_CONVERSION);
      dp->setProcessingActions(actions);
      dp->setMetaValue("cached_data", "true");
      for (Size i=0; i<out_exp.size(); ++i)
      {
        out_exp[i].getDataProcessing().push_back(dp);
      }
      std::vector<MSChromatogram > l_chromatograms = out_exp.getChromatograms();
      for (Size i=0; i<l_chromatograms.size(); ++i)
      {
        l_chromatograms[i].getDataProcessing().push_back(dp);
      }
      out_exp.setChromatograms(l_chromatograms);
    }

    // store the meta data using the regular MzMLFile
    MzMLFile().store(out_meta, out_exp);
  }

  std::vector<OpenSwath::BinaryDataArrayPtr> CachedMzMLHandler::readSpectrumFast(std::ifstream& ifs, int& ms_level, double& rt)
  {
    std::vector<OpenSwath::BinaryDataArrayPtr> data;
    data.push_back(OpenSwath::BinaryDataArrayPtr(new OpenSwath::BinaryDataArray));
    data.push_back(OpenSwath::BinaryDataArrayPtr(new OpenSwath::BinaryDataArray));

    Size spec_size = -1;
    Size nr_float_arrays = -1;
    ifs.read((char*) &spec_size, sizeof(spec_size));
    ifs.read((char*) &nr_float_arrays, sizeof(nr_float_arrays));
    ifs.read((char*) &ms_level, sizeof(ms_level));
    ifs.read((char*) &rt, sizeof(rt));

    if (static_cast<int>(spec_size) < 0)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Read an invalid spectrum length, something is wrong here. Aborting.", "filestream");
    }

    readDataFast_(ifs, data, spec_size, nr_float_arrays);
    return data;
  }

  void CachedMzMLHandler::readDataFast_(std::ifstream& ifs,
                                        std::vector<OpenSwath::BinaryDataArrayPtr>& data,
                                        const Size& data_size,
                                        const Size& nr_float_arrays)
  {
    OPENMS_PRECONDITION(data.size() == 2, "Input data needs to have 2 slots.")

    data[0]->data.resize(data_size);
    data[1]->data.resize(data_size);

    if (data_size > 0)
    {
      ifs.read((char*) &(data[0]->data)[0], data_size * sizeof(DatumSingleton));
      ifs.read((char*) &(data[1]->data)[0], data_size * sizeof(DatumSingleton));
    }
    if (nr_float_arrays == 0)
    {
      return;
    }
    char* buffer = new(std::nothrow) char[1024];
    for (Size k = 0; k < nr_float_arrays; k++)
    {
      data.push_back(OpenSwath::BinaryDataArrayPtr(new OpenSwath::BinaryDataArray));
      Size len, len_name;
      ifs.read((char*)&len, sizeof(len));
      ifs.read((char*)&len_name, sizeof(len_name));

      // We will not read data longer than 1024 bytes as this will not fit into
      // our buffer (and is user-generated input data)
      if (len_name > 1023)
      {
        ifs.seekg(len_name * sizeof(char), ifs.cur);
      }
      else
      {
        ifs.read(buffer, len_name);
        buffer[len_name] = '\0';
      }
      data.back()->data.resize(len);
      data.back()->description = buffer;
      ifs.read((char*)&(data.back()->data)[0], len * sizeof(DatumSingleton));
    }
    delete[] buffer;
    return;
  }

  std::vector<OpenSwath::BinaryDataArrayPtr> CachedMzMLHandler::readChromatogramFast(std::ifstream& ifs)
  {
    std::vector<OpenSwath::BinaryDataArrayPtr> data;
    data.push_back(OpenSwath::BinaryDataArrayPtr(new OpenSwath::BinaryDataArray));
    data.push_back(OpenSwath::BinaryDataArrayPtr(new OpenSwath::BinaryDataArray));

    Size chrom_size = -1;
    Size nr_float_arrays = -1;
    ifs.read((char*) &chrom_size, sizeof(chrom_size));
    ifs.read((char*) &nr_float_arrays, sizeof(nr_float_arrays));

    if (static_cast<int>(chrom_size) < 0)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Read an invalid chromatogram length, something is wrong here. Aborting.", "filestream");
    }

    readDataFast_(ifs, data, chrom_size, nr_float_arrays);
    return data;
  }

  void CachedMzMLHandler::readSpectrum(SpectrumType& spectrum, std::ifstream& ifs)
  {
    int ms_level;
    double rt;
    std::vector<OpenSwath::BinaryDataArrayPtr> data = readSpectrumFast(ifs, ms_level, rt);
    spectrum.reserve(data[0]->data.size());
    spectrum.setMSLevel(ms_level);
    spectrum.setRT(rt);

    for (Size j = 0; j < data[0]->data.size(); j++)
    {
      Peak1D p;
      p.setMZ(data[0]->data[j]);
      p.setIntensity(data[1]->data[j]);
      spectrum.push_back(p);
    }

    for (Size j = 2; j < data.size(); j++)
    {
      spectrum.getFloatDataArrays().push_back(MSSpectrum::FloatDataArray());
      spectrum.getFloatDataArrays().back().reserve(data[j]->data.size());
      spectrum.getFloatDataArrays().back().setName(data[j]->description);
      for (const auto& k : data[j]->data) spectrum.getFloatDataArrays().back().push_back(k);
    }
  }

  void CachedMzMLHandler::readChromatogram(ChromatogramType& chromatogram, std::ifstream& ifs)
  {
    std::vector<OpenSwath::BinaryDataArrayPtr> data = readChromatogramFast(ifs);
    chromatogram.reserve(data[0]->data.size());

    for (Size j = 0; j < data[0]->data.size(); j++)
    {
      ChromatogramPeak p;
      p.setRT(data[0]->data[j]);
      p.setIntensity(data[1]->data[j]);
      chromatogram.push_back(p);
    }

    MSChromatogram::FloatDataArrays fdas;
    for (Size j = 2; j < data.size(); j++)
    {
      MSChromatogram::FloatDataArray fda;
      fda.reserve(data[j]->data.size());
      for (const auto& k : fda) fda.push_back(k);
      fda.setName(data[j]->description);
      fdas.push_back(fda);
    }
    chromatogram.setFloatDataArrays(fdas);
  }

  void CachedMzMLHandler::writeSpectrum_(const SpectrumType& spectrum, std::ofstream& ofs) const
  {
    Size exp_size = spectrum.size();
    ofs.write((char*)&exp_size, sizeof(exp_size));
    Size arr_s = spectrum.getFloatDataArrays().size() + spectrum.getIntegerDataArrays().size();
    ofs.write((char*)&arr_s, sizeof(arr_s));
    IntType int_field_ = spectrum.getMSLevel();
    ofs.write((char*)&int_field_, sizeof(int_field_));
    DoubleType dbl_field_ = spectrum.getRT();
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

    Datavector tmp;
    for (const auto& fda : spectrum.getFloatDataArrays() )
    {
      Size len = fda.size();
      ofs.write((char*)&len, sizeof(len));
      Size len_name = fda.getName().size();
      ofs.write((char*)&len_name, sizeof(len_name));
      ofs.write((char*)&fda.getName().front(), len_name * sizeof(fda.getName().front()));
      // now go to the actual data
      tmp.clear();
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
      tmp.clear();
      tmp.reserve(ida.size());
      for (const auto& val : ida) {tmp.push_back(val);}
      ofs.write((char*)&tmp.front(), tmp.size() * sizeof(tmp.front()));
    }
  }

  void CachedMzMLHandler::writeChromatogram_(const ChromatogramType& chromatogram, std::ofstream& ofs) const
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

    Datavector tmp;
    for (const auto& fda : chromatogram.getFloatDataArrays() )
    {
      Size len = fda.size();
      ofs.write((char*)&len, sizeof(len));
      Size len_name = fda.getName().size();
      ofs.write((char*)&len_name, sizeof(len_name));
      ofs.write((char*)&fda.getName().front(), len_name * sizeof(fda.getName().front()));
      // now go to the actual data
      tmp.clear();
      tmp.reserve(fda.size());
      for (const auto& val : fda)
      {
        tmp.push_back(val);
      }
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
      tmp.clear();
      tmp.reserve(ida.size());
      for (const auto& val : ida)
      {
        tmp.push_back(val);
      }
      ofs.write((char*)&tmp.front(), tmp.size() * sizeof(tmp.front()));
    }
  }

}//namespace OpenMS  //namespace Internal

