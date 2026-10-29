// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/DATAACCESS/HTTPInputSource.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/FileNameUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <fstream>
#include <exception>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  MzMLFile::MzMLFile() :
    XMLFile("/SCHEMAS/mzML_1_10.xsd", "1.1.0"),
    indexed_schema_location_("/SCHEMAS/mzML_idx_1_10.xsd")
  {
  }

  MzMLFile::~MzMLFile() = default;

  namespace
  {
    // ponytail: fixed cap; per-chunk header reparse and memory bandwidth stop paying off above ~32 threads (EPYC 9654).
    constexpr int PARALLEL_LOAD_MAX_THREADS = 32;
    /// Contiguous spectrum ranges per thread, so uneven spectra (MS1 frames vs. MS2) still balance.
    constexpr int PARALLEL_LOAD_CHUNKS_PER_THREAD = 4;
    /// Smaller chunks cost more in header reparsing and thread start-up than they save (small files).
    constexpr std::streamoff PARALLEL_LOAD_MIN_CHUNK_BYTES = 1 << 20;

    /// Appends bytes [begin, end) of @p in to @p out (reusing its capacity).
    void appendRange(std::ifstream& in, std::streamoff begin, std::streamoff end, std::string& out)
    {
      const Size old_size = out.size();
      out.resize(old_size + static_cast<Size>(end - begin));
      in.seekg(begin);
      in.read(out.data() + old_size, static_cast<std::streamsize>(end - begin));
      if (in.gcount() != static_cast<std::streamsize>(end - begin))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Short read in indexed mzML");
      }
    }

    std::string readRange(std::ifstream& in, std::streamoff begin, std::streamoff end)
    {
      std::string bytes;
      appendRange(in, begin, end, bytes);
      return bytes;
    }
  }

  bool MzMLFile::loadParallelIndexed_(const std::string& filename, PeakMap& map)
  {
#ifdef _OPENMP
    const int threads = std::min(omp_get_max_threads(), PARALLEL_LOAD_MAX_THREADS);
    // Only the plain full load is parallelized; everything else keeps the serial streaming parse.
    if (threads < 2 || FileNameUtils::compressionType(filename) != FileTypes::UNKNOWN || options_.hasFilters() ||
        !options_.getFillData() || options_.getMetadataOnly() || options_.getSkipChromatograms())
    {
      return false;
    }
    IndexedMzMLDecoder decoder;
    const std::streampos index_offset = decoder.findIndexListOffset(filename);
    IndexedMzMLDecoder::OffsetVector spectrum_offsets, chromatogram_offsets;
    if (index_offset == std::streampos(-1) ||
        decoder.parseOffsets(filename, index_offset, spectrum_offsets, chromatogram_offsets) != 0 ||
        spectrum_offsets.size() < 2)
    {
      return false;
    }

    // Split the document: everything before the first spectrum is the header (prolog, mzML/run metadata and the
    // spectrumList start tag); each chunk is header + a contiguous run of <spectrum> elements + closing tags, so a
    // normal MzMLHandler parses metadata and peaks exactly as in the serial load.
    std::ifstream in(filename, std::ios::binary);
    const std::streamoff first = spectrum_offsets.front().second;
    const std::string header = readRange(in, 0, first);
    const bool wrapped = header.find("<indexedmzML") != std::string::npos;
    const std::string closing = std::string("</spectrumList></run></mzML>") + (wrapped ? "</indexedmzML>" : "");
    const std::streamoff tail_begin = spectrum_offsets.back().second;
    const std::string tail = readRange(in, tail_begin, index_offset);
    const Size spectra_end_in_tail = tail.find("</spectrumList>");
    if (spectra_end_in_tail == std::string::npos)
    {
      return false;
    }
    const std::streamoff spectra_end = tail_begin + static_cast<std::streamoff>(spectra_end_in_tail);
    for (Size i = 1; i < spectrum_offsets.size(); ++i)
    {
      if (spectrum_offsets[i].second <= spectrum_offsets[i - 1].second) return false; // not in document order
    }

    // Chunk boundaries balanced by bytes, always at spectrum starts.
    const Size n_chunks = std::max<Size>(1, std::min<Size>({spectrum_offsets.size(), static_cast<Size>(threads) * PARALLEL_LOAD_CHUNKS_PER_THREAD,
                                                             static_cast<Size>((spectra_end - first) / PARALLEL_LOAD_MIN_CHUNK_BYTES)}));
    if (n_chunks < 2)
    {
      return false; // too small to gain from splitting
    }
    std::vector<std::streamoff> bounds{first};
    for (Size c = 1; c < n_chunks; ++c)
    {
      const std::streamoff target = first + (spectra_end - first) * static_cast<std::streamoff>(c) / static_cast<std::streamoff>(n_chunks);
      auto it = std::lower_bound(spectrum_offsets.begin(), spectrum_offsets.end(), target,
                                 [](const auto& entry, std::streamoff value) { return std::streamoff(entry.second) < value; });
      if (it != spectrum_offsets.end() && std::streamoff(it->second) > bounds.back()) bounds.push_back(it->second);
    }
    bounds.push_back(spectra_end);
    const Size chunks = bounds.size() - 1;

    // Parts: spectrum chunks, then (if any) everything from </spectrumList> to the index: the chromatogram list.
    const bool has_chromatograms = !chromatogram_offsets.empty();
    std::vector<PeakMap> parts(chunks + (has_chromatograms ? 1 : 0));
    std::exception_ptr error;
    #pragma omp parallel num_threads(std::min<int>(threads, static_cast<int>(parts.size())))
    {
      std::ifstream part_in(filename, std::ios::binary);
      // One buffer per thread, reused across chunks: freeing large buffers in a process with many threads
      // costs a TLB shootdown on every core, which dominated with several concurrent loads.
      std::string buffer;
      #pragma omp for schedule(dynamic, 1)
      for (SignedSize p = 0; p < static_cast<SignedSize>(parts.size()); ++p)
      {
        try
        {
          buffer.assign(header);
          if (static_cast<Size>(p) < chunks)
          {
            appendRange(part_in, bounds[p], bounds[p + 1], buffer);
            buffer += closing;
          }
          else
          {
            appendRange(part_in, spectra_end, index_offset, buffer);
            if (wrapped) buffer += "</indexedmzML>";
          }
          MzMLFile part_file;
          part_file.setOptions(options_);
          Internal::MzMLHandler handler(parts[p], filename, getVersion(), part_file);
          handler.setOptions(options_);
          part_file.parseBuffer_(buffer, &handler);
        }
        catch (...)
        {
          #pragma omp critical(MzMLFileParallelLoad)
          if (!error) error = std::current_exception();
        }
      }
    }
    if (error) std::rethrow_exception(error);

    // Settings come from the header (identical in every part); spectra are concatenated in document order.
    Size n_spectra = 0;
    for (Size p = 0; p < chunks; ++p) n_spectra += parts[p].size();
    const Size n_chromatograms = has_chromatograms ? parts.back().getChromatograms().size() : 0;
    if (n_spectra != spectrum_offsets.size() || n_chromatograms != chromatogram_offsets.size())
    {
      return false; // index disagrees with the document; the caller reparses serially
    }
    std::vector<MSSpectrum> spectra;
    spectra.reserve(n_spectra);
    for (Size p = 0; p < chunks; ++p)
    {
      for (auto& spectrum : parts[p].getSpectra()) spectra.push_back(std::move(spectrum));
    }
    std::vector<MSChromatogram> chromatograms;
    if (has_chromatograms) chromatograms = std::move(parts.back().getChromatograms());
    map = std::move(parts.front());
    map.setSpectra(std::move(spectra));
    map.setChromatograms(std::move(chromatograms));
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);
    return true;
#else
    (void)filename;
    (void)map;
    return false;
#endif
  }

  PeakFileOptions& MzMLFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions& MzMLFile::getOptions() const
  {
    return options_;
  }

  void MzMLFile::setOptions(const PeakFileOptions& options)
  {
    options_ = options;
  }

  bool MzMLFile::hasIndex(const std::string& filename)
  {
    const std::streampos NOT_FOUND {-1};
    return NOT_FOUND != IndexedMzMLDecoder().findIndexListOffset(filename);
  }

  // reimplemented in order to handle index MzML
  bool MzMLFile::isValid(const std::string& filename, std::ostream& os)
  {
    //determine if this is indexed mzML or not
    bool indexed = false;
    TextFile file(filename, true, 4);
    std::string s = StringUtils::concatenate(file);
    if (StringUtils::hasSubstring(s, "<indexedmzML"))
    {
      indexed = true;
    }
    // find the corresponding schema
    std::string current_location;
    if (indexed)
    {
      current_location = File::find(indexed_schema_location_);
    }
    else
    {
      current_location = File::find(schema_location_);
    }

    return XMLValidator().isValid(filename, current_location, os);
  }

  bool MzMLFile::isSemanticallyValid(const std::string& filename, StringList& errors, StringList& warnings)
  {
    // load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/ms-mapping.xml"), mapping);

    // validate
    Internal::MzMLValidator v(mapping, ControlledVocabulary::getPSIMSCV());
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

  void MzMLFile::loadSize(const std::string& filename, Size& scount, Size& ccount)
  {
    PeakMap dummy;
    Internal::MzMLHandler handler(dummy, filename, getVersion(), *this);
    handler.setOptions(options_);
    if (options_.hasFilters())
    {
      handler.setLoadDetail(Internal::XMLHandler::LD_COUNTS_WITHOPTIONS);
    }
    else
    { // no filters where specified. Just take the 'counts' attributes from the mzML file and end parsing
      handler.setLoadDetail(Internal::XMLHandler::LD_RAWCOUNTS);
    }
    
    safeParse_(filename, &handler);
    handler.getCounts(scount, ccount);
  }

  void MzMLFile::safeParse_(const std::string& filename, Internal::XMLHandler* handler)
  {
    try
    {
      parse_(filename, handler);
    }
    catch (Exception::ParseError&)
    {
      throw; // names the file and the problem (e.g. XMLHandler::fatalError())
    }
    catch (Exception::BaseException& e)
    {
      // the problem and the file, not the source location of the exception (debug output, see TOPPBase); a URL without
      // its credentials
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, HTTPInputSource::redactURL(filename),
                                  std::string(e.getName()) + ": " + e.what());
    }
  }

  void MzMLFile::loadBuffer(const std::string& buffer, PeakMap& map)
  {
    map.reset();

    Internal::MzMLHandler handler(map, "memory", getVersion(), *this);
    handler.setOptions(options_);
    parseBuffer_(buffer, &handler);
    map.updateRanges();
  }

  void MzMLFile::load(const std::string& filename, PeakMap& map)
  {
    const bool remote = StringUtils::hasPrefix(filename, "http://") || StringUtils::hasPrefix(filename, "https://")
                        || StringUtils::hasPrefix(filename, "ftp://") || StringUtils::hasPrefix(filename, "s3://");
    // the parallel load reads byte ranges of a local file; remote documents are streamed
    if (remote || !loadParallelIndexed_(filename, map))
    {
      map.reset();

      //set DocumentIdentifier
      if (remote)
      {
        map.setLoadedFileType(FileTypes::MZML);
        map.setLoadedNetFilePath(filename);
      }
      else
      {
        map.setLoadedFileType(filename);
        map.setLoadedFilePath(filename);
      }

      Internal::MzMLHandler handler(map, filename, getVersion(), *this);
      handler.setOptions(options_);
      safeParse_(filename, &handler);
    }
    map.updateRanges();
  }

  void MzMLFile::store(const std::string& filename, const PeakMap& map) const
  {
    Internal::MzMLHandler handler(map, filename, getVersion(), *this);
    handler.setOptions(options_);
    save_(filename, &handler);
  }

  void MzMLFile::storeBuffer(std::string& output, const PeakMap& map) const
  {
    Internal::MzMLHandler handler(map, "dummy", getVersion(), *this);
    handler.setOptions(options_);
    {
      std::stringstream os;

      //set high precision for writing of floating point numbers
      os.precision(writtenDigits(double()));

      // write data and close stream
      handler.writeTo(os);
      output = os.str();
    }
  }

  void MzMLFile::transform(const std::string& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count, bool skip_first_pass)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    if (!skip_first_pass) transformFirstPass_(filename_in, consumer, skip_full_count);

    // Second pass through the data, now read the spectra!
    {
      PeakMap dummy;
      Internal::MzMLHandler handler(dummy, filename_in, getVersion(), *this);
      handler.setOptions(options_);
      handler.setMSDataConsumer(consumer);
      safeParse_(filename_in, &handler);
    }
  }

  void MzMLFile::transform(const std::string& filename_in, Interfaces::IMSDataConsumer* consumer, PeakMap& map, bool skip_full_count, bool skip_first_pass)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    if (!skip_first_pass)
    {
      transformFirstPass_(filename_in, consumer, skip_full_count);
    }
    // Second pass through the data, now read the spectra!
    {
      PeakFileOptions tmp_options(options_);
      Internal::MzMLHandler handler(map, filename_in, getVersion(), *this);
      tmp_options.setAlwaysAppendData(true);
      handler.setOptions(tmp_options);
      handler.setMSDataConsumer(consumer);

      safeParse_(filename_in, &handler);
    }
  }

  void MzMLFile::transformFirstPass_(const std::string& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count)
  {
    // Create temporary objects and counters
    PeakFileOptions tmp_options(options_);
    Size scount = 0, ccount = 0;
    PeakMap experimental_settings;
    Internal::MzMLHandler handler(experimental_settings, filename_in, getVersion(), *this);

    // set temporary options for handler
    tmp_options.setMetadataOnly( skip_full_count );
    handler.setOptions(tmp_options);
    handler.setLoadDetail(Internal::XMLHandler::LD_RAWCOUNTS);

    safeParse_(filename_in, &handler);

    // After parsing, collect information
    handler.getCounts(scount, ccount);
    consumer->setExpectedSize(scount, ccount);
    consumer->setExperimentalSettings(experimental_settings);
  }

  std::map<UInt, MzMLFile::SpecInfo> MzMLFile::getCentroidInfo(const std::string& filename, const Size first_n_spectra_only)
  {
    bool oldoption = options_.getFillData();
    options_.setFillData(true); // we want the data as well (to allow estimation from data if metadata is missing)
    std::map<UInt, SpecInfo> ret;
    Size first_n_spectra_only_remaining = first_n_spectra_only;
    auto f = [&ret, &first_n_spectra_only_remaining](const MSSpectrum& s)
    {
        UInt lvl = s.getMSLevel();
        switch (s.getType(true))
        {
          case (MSSpectrum::SpectrumType::CENTROID):
            ++ret[lvl].count_centroided;
            --first_n_spectra_only_remaining;
            break;
          case (MSSpectrum::SpectrumType::PROFILE):
            ++ret[lvl].count_profile;
            --first_n_spectra_only_remaining;
            break;
          case (MSSpectrum::SpectrumType::UNKNOWN):  // this can only happen for spectra with very few peaks (or completely empty spectra)
            ++ret[lvl].count_unknown;
            break;
          default:
            // make sure we did not forget a case
            throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }

        // stop parsing after 10 or so spectra
        if (first_n_spectra_only_remaining == 0)
        {
          throw Internal::XMLHandler::EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
    };
    MSDataTransformingConsumer c;
    c.setSpectraProcessingFunc(f);
    transform(filename, &c, true, true); // no first pass

    // restore old state
    options_.setFillData(oldoption);

    return ret;
  }

} // namespace OpenMS
