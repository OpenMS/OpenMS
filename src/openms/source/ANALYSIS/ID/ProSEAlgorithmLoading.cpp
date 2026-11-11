// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include "ProSEAlgorithmImpl.h"

#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <iterator>
#include <string_view>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

  namespace
  {
    // The positions 0 .. n - 1 of @p keys (no NaN) ordered by (key, position) ascending, as MSExperiment::sortSpectra()
    // orders spectra by retention time: a strict total order, so the order does not depend on how it is sorted. Parts
    // sorted by the OpenMP threads, then merged pairwise.
    std::vector<Size> orderByKeyInParallel(const std::vector<double>& keys)
    {
      const Size n = keys.size();
      std::vector<Size> order(n);
      std::iota(order.begin(), order.end(), Size(0));
      const auto less = [&keys](Size a, Size b) { return keys[a] < keys[b] || (!(keys[b] < keys[a]) && a < b); };
      int parts = 1;
#ifdef _OPENMP
      parts = static_cast<int>(std::max<Size>(1, std::min<Size>(static_cast<Size>(omp_get_max_threads()), n / 4096)));
#endif
      std::vector<Size> bounds(static_cast<Size>(parts) + 1);
      for (int p = 0; p <= parts; ++p) { bounds[p] = n / static_cast<Size>(parts) * static_cast<Size>(p) + std::min(n % static_cast<Size>(parts), static_cast<Size>(p)); }
#pragma omp parallel for schedule(static, 1)
      for (int p = 0; p < parts; ++p) { std::sort(order.begin() + bounds[p], order.begin() + bounds[p + 1], less); }
      std::vector<Size> merged(n);
      for (int width = 1; width < parts; width *= 2)
      {
        const int pairs = (parts + 2 * width - 1) / (2 * width);
#pragma omp parallel for schedule(dynamic, 1)
        for (int q = 0; q < pairs; ++q)
        {
          const int p = 2 * width * q;
          const Size lo = bounds[p], mid = bounds[std::min(p + width, parts)], hi = bounds[std::min(p + 2 * width, parts)];
          std::merge(order.begin() + lo, order.begin() + mid, order.begin() + mid, order.begin() + hi, merged.begin() + lo, less);
        }
        order.swap(merged);
      }
      return order;
    }

    // FASTAFile::load() with the OpenMP threads: the file is cut into pieces at starts of entries,
    // the threads read the pieces with FASTAFile::readNext(), and the pieces are joined in file
    // order. Same entries as load(): readNext() reads an entry from its '>' up to a line break
    // followed by '>'. A piece starts at a '>' after a line break whose line holds no '>' (not a
    // header line, so usually the last line of a sequence), where the reader of load() starts an
    // entry as well, and the reader of the piece before stops when it arrives there. Whenever this
    // does not work out (small file, no such place, a reader that fails or does not arrive at the
    // next piece, e.g. as it skipped an entry without sequence there), load() reads the file and
    // reports errors as before.
    void readFASTA(const std::string& filename, std::vector<FASTAFile::FASTAEntry>& entries)
    {
      std::vector<std::streamoff> pieces{0}; // where each piece starts; the first one where readStart() starts
#ifdef _OPENMP
      if (!omp_in_parallel() && omp_get_max_threads() > 1)
      {
        std::ifstream in(filename, std::ios::binary);
        in.seekg(0, std::ios::end);
        const std::streamoff size = in.tellg();
        // several pieces per thread: the threads do not read equally fast
        const std::streamoff count = std::min<std::streamoff>(4 * omp_get_max_threads(), size / (std::streamoff(1) << 18));
        std::string window(Size(1) << 16, '\0');
        for (std::streamoff i = 1; i < count && in.good(); ++i)
        {
          const std::streamoff from = size / count * i;
          in.seekg(from);
          in.read(window.data(), static_cast<std::streamsize>(window.size()));
          const std::string_view text(window.data(), static_cast<Size>(in.gcount()));
          in.clear(); // reading up to the end of the file is fine
          // lines of the window that are complete: between two line breaks
          for (Size line = text.find('\n'); line != std::string_view::npos && line + 1 < text.size();)
          {
            const Size next = text.find('\n', line + 1);
            if (next == std::string_view::npos || next + 1 >= text.size()) { break; }
            if (text[next + 1] == '>' && text.substr(line + 1, next - line - 1).find('>') == std::string_view::npos)
            {
              if (from + static_cast<std::streamoff>(next) + 1 > pieces.back()) { pieces.push_back(from + static_cast<std::streamoff>(next) + 1); }
              break;
            }
            line = next;
          }
        }
      }
#endif
      if (pieces.size() < 2)
      {
        FASTAFile().load(filename, entries);
        return;
      }

      std::vector<std::vector<FASTAFile::FASTAEntry>> read(pieces.size());
      bool complete = true;
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize i = 0; i < static_cast<SignedSize>(pieces.size()); ++i)
      {
        bool ok = false;
        try // exceptions must not leave the parallel region
        {
          FASTAFile file;
          FASTAFile::FASTAEntry entry;
          file.readStart(filename);
          ok = i == 0 || file.setPosition(pieces[i]);
          if (ok && i + 1 == static_cast<SignedSize>(pieces.size()))
          { // the last piece: up to the end of the file, as load()
            while (file.readNext(entry)) { read[i].push_back(std::move(entry)); }
          }
          else if (ok)
          {
            const std::streampos end = pieces[i + 1];
            while (ok && file.position() < end)
            {
              ok = file.readNext(entry);
              if (ok) { read[i].push_back(std::move(entry)); }
            }
            ok = ok && file.position() == end;
          }
        }
        catch (...)
        {
          ok = false;
        }
        if (!ok)
        {
#pragma omp critical (ProSEAlgorithm_loadFASTA)
          complete = false;
        }
      }
      if (!complete)
      {
        FASTAFile().load(filename, entries);
        return;
      }

      Size total = 0;
      for (const auto& piece : read) { total += piece.size(); }
      entries.clear();
      entries.reserve(total);
      for (auto& piece : read)
      {
        entries.insert(entries.end(), std::make_move_iterator(piece.begin()), std::make_move_iterator(piece.end()));
      }
    }
  }

  namespace Internal::ProSEAlgorithmHelpers
  {
    // The database of a search: readFASTA(), without repeated entries (ProSEAlgorithm::removeDuplicateProteins(), which
    // names the file in its warning or error), before the database is used.
    void loadFASTA(const std::string& filename, std::vector<FASTAFile::FASTAEntry>& entries)
    {
      readFASTA(filename, entries);
      ProSEAlgorithm::removeDuplicateProteins(entries, filename);
    }
  } // namespace Internal::ProSEAlgorithmHelpers

  namespace
  {
    // Removes the spectra of @p spectra that @p keep does not accept (none without @p keep)
    void removeUnkept(std::vector<MSSpectrum>& spectra, const std::function<bool(const MSSpectrum&)>& keep)
    {
      if (!keep) { return; }
      spectra.erase(std::remove_if(spectra.begin(), spectra.end(), [&keep](const MSSpectrum& s) { return !keep(s); }), spectra.end());
    }

    // MzMLFile::load() with @p threads OpenMP threads (at least 2). The spectra are cut into chunks at
    // their <spectrum> start tags; every chunk is parsed together with the file's header and closing
    // tags by an MzMLFile of its own (MzMLFile::loadBuffer), and the spectra are joined in file order.
    // The last chunk runs up to </run>, so it also holds the chromatograms. Each thread reads its chunks
    // from the file itself, so the data held beyond the spectra is at most about 9 MB per thread,
    // whatever the size of the file. The header is parsed first on this thread: it gives the experimental
    // settings, registers its meta value names in file order and initialises the XML parser.
    // Returns false if the file does not have the layout the chunks rely on (compressed; a comment,
    // CDATA section or DOCTYPE before the end of the run; not exactly one spectrum list; a count
    // attribute other than the number of spectra found; few spectra) or a chunk fails to parse. The
    // caller then reads the file with MzMLFile::load(), which also reports its errors.
    // With @p keep only the spectra it accepts are kept (by each chunk, right after parsing it).
    bool loadMzMLChunked(const std::string& filename, const PeakFileOptions& options, int threads, PeakMap& exp,
                         ChunkHelp* help = nullptr, int chunk_threads = 0,
                         const std::function<bool(const MSSpectrum&)>& keep = {})
    {
      std::ifstream in(filename, std::ios::binary);
      if (!in) { return false; }
      in.seekg(0, std::ios::end);
      const Size size = static_cast<Size>(in.tellg());
      char magic[2] = {0, 0};
      in.seekg(0);
      in.read(magic, 2);
      if (!in || (magic[0] == 'B' && magic[1] == 'Z') || (magic[0] == '\x1f' && magic[1] == '\x8b') || (magic[0] == 'P' && magic[1] == 'K'))
      {
        return false; // compressed (bzip2, gzip, zip) or too short
      }

      // 1. Where the markup of interest starts, found by the threads in blocks of the file.
      enum class Mark : unsigned char { spectrum, list_start, list_end, run_end, declaration };
      const auto is_tag = [](std::string_view text, std::string_view name)
      {
        return text.size() > name.size() && text.substr(0, name.size()) == name
               && std::string_view(" \t\r\n>/").find(text[name.size()]) != std::string_view::npos;
      };
      constexpr Size block = Size(1) << 20;
      constexpr Size lookahead = 16; // longer than the longest name compared below, "</spectrumList" plus one
      const Size n_blocks = (size + block - 1) / block;
      std::vector<std::vector<std::pair<Size, Mark>>> marks(n_blocks);
      bool ok = true;
#pragma omp parallel num_threads(threads)
      {
        std::ifstream file; // opened by the thread's first block, in the try below
        std::string buffer;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize sb = 0; sb < static_cast<SignedSize>(n_blocks); ++sb)
        {
          bool read = false;
          try // exceptions (e.g. std::bad_alloc of the buffers) must not leave the parallel region
          {
            if (!file.is_open()) { file.open(filename, std::ios::binary); }
            const Size b = static_cast<Size>(sb);
            const Size from = b * block;
            const Size end = std::min(block, size - from); // marks starting before end belong to this block
            buffer.resize(std::min(block + lookahead, size - from));
            file.seekg(static_cast<std::streamoff>(from));
            file.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
            read = static_cast<bool>(file);
            for (Size i = 0; read && i < end; ++i)
            {
              const char* p = static_cast<const char*>(std::memchr(buffer.data() + i, '<', end - i));
              if (p == nullptr) { break; }
              i = static_cast<Size>(p - buffer.data());
              const std::string_view text(p, std::min(lookahead, buffer.size() - i));
              if (is_tag(text, "<spectrumList")) { marks[b].emplace_back(from + i, Mark::list_start); }
              else if (is_tag(text, "<spectrum")) { marks[b].emplace_back(from + i, Mark::spectrum); }
              else if (is_tag(text, "</spectrumList")) { marks[b].emplace_back(from + i, Mark::list_end); }
              else if (is_tag(text, "</run")) { marks[b].emplace_back(from + i, Mark::run_end); }
              else if (text.substr(0, 2) == "<!") { marks[b].emplace_back(from + i, Mark::declaration); }
            }
          }
          catch (...)
          {
            read = false;
          }
          if (!read)
          {
#pragma omp critical (ProSEAlgorithm_loadMzMLChunked)
            ok = false;
            file.clear();
          }
        }
      }
      constexpr Size none = std::numeric_limits<Size>::max();
      Size list_start = none, list_end = none, run_end = none;
      std::vector<Size> starts; // of the spectra
      for (const auto& block_marks : marks)
      {
        for (const auto& [pos, mark] : block_marks)
        {
          switch (mark)
          {
            case Mark::list_start: ok = ok && list_start == none; list_start = pos; break;
            case Mark::spectrum: ok = ok && list_start != none && list_end == none; starts.push_back(pos); break;
            case Mark::list_end: ok = ok && list_start != none && list_end == none; list_end = pos; break;
            case Mark::run_end: if (run_end == none) { ok = ok && list_end != none; run_end = pos; } break;
            case Mark::declaration: ok = ok && run_end != none; break;
          }
        }
      }
      marks = {};
      if (!ok || run_end == none || starts.size() < 64) { return false; }

      // 2. The header (everything before the first spectrum) and the count attribute of its
      //    <spectrumList>, which each chunk sets to its own number of spectra.
      std::string header(starts[0], '\0');
      in.seekg(0);
      in.read(header.data(), static_cast<std::streamsize>(header.size()));
      if (!in) { return false; }
      Size count_begin = none, count_end = none;
      for (Size p = list_start + std::string_view("<spectrumList").size();;)
      {
        p = header.find_first_not_of(" \t\r\n", p);
        if (p == std::string::npos) { return false; }
        if (header[p] == '>') { break; }
        const Size name_end = header.find_first_of("= \t\r\n>", p);
        Size q = header.find_first_not_of(" \t\r\n", name_end);
        if (name_end == std::string::npos || q == std::string::npos || header[q] != '=') { return false; }
        q = header.find_first_not_of(" \t\r\n", q + 1);
        if (q == std::string::npos || (header[q] != '"' && header[q] != '\'')) { return false; }
        const Size value_end = header.find(header[q], q + 1);
        if (value_end == std::string::npos) { return false; }
        if (std::string_view(header).substr(p, name_end - p) == "count") { count_begin = q + 1; count_end = value_end; }
        p = value_end + 1;
      }
      if (count_begin == none || header.compare(count_begin, count_end - count_begin, std::to_string(starts.size())) != 0)
      {
        return false;
      }
      const std::string end_of_run = std::string("</run></mzML>") + (header.find("<indexedmzML") != std::string::npos ? "</indexedmzML>" : "");
      const auto chunk_text = [&](Size n_spectra, std::string& text)
      {
        text.assign(header, 0, count_begin);
        text += std::to_string(n_spectra);
        text.append(header, count_end, std::string::npos);
      };

      // 3. Chunks of about equal size: two per thread (the threads do not parse equally fast; every chunk
      //    costs the set-up of a parser and handler), of at most 8 MB, and of at least 16 spectra. A read in the
      //    background counts chunk_threads instead (see chunkThreads()).
      constexpr Size max_chunk_bytes = Size(8) << 20;
      const Size n = starts.size();
      const Size bytes = list_end - starts[0];
      const Size n_chunks = std::min(n / 16, std::max(2 * static_cast<Size>(std::max(threads, chunk_threads)), (bytes + max_chunk_bytes - 1) / max_chunk_bytes));
      std::vector<Size> first{0}; // the first spectrum of each chunk
      for (Size c = 1; c < n_chunks; ++c)
      {
        const Size f = static_cast<Size>(std::lower_bound(starts.begin(), starts.end(), starts[0] + bytes / n_chunks * c) - starts.begin());
        if (f > first.back() && f < n) { first.push_back(f); }
      }
      first.push_back(n);

      // A chunk decodes its spectra on its thread; in batches of the default size, so that it does not
      // hold the encoded data of all its spectra at once.
      PeakFileOptions chunk_options = options;
      chunk_options.setMaxDataPoolSize(PeakFileOptions().getMaxDataPoolSize());
      try
      {
        std::string text;
        chunk_text(0, text);
        text += "</spectrumList>" + end_of_run;
        MzMLFile mzml;
        mzml.getOptions() = chunk_options;
        mzml.loadBuffer(text, exp);
      }
      catch (...)
      {
        return false;
      }
      // The spectra that each chunk keeps, until all are parsed and put in their places (below)
      const Size n_chunks_used = first.size() - 1;
      std::vector<std::vector<MSSpectrum>> parts(n_chunks_used);
      std::vector<Size> kept(n_chunks_used, 0);
      std::vector<MSChromatogram> chromatograms;
      // The chunks are handed out one by one, to the threads of this reader and to those of a helper (another
      // contention group: no OpenMP critical between the two).
      std::atomic<Size> next_chunk{0};
      std::atomic<bool> chunks_ok{true};
      const auto parse_chunks = [&]([[maybe_unused]] const int team)
      {
#pragma omp parallel num_threads(team)
      {
        std::ifstream file; // opened by the thread's first chunk, in the try below
        std::string chunk;
        for (Size c = next_chunk++; c < n_chunks_used; c = next_chunk++)
        {
          bool parsed = false;
          try // exceptions (a parse error, std::bad_alloc of the chunk text, ...) must not leave the parallel region
          {
            if (!file.is_open()) { file.open(filename, std::ios::binary); }
            const bool last = c + 1 == n_chunks_used;
            const Size from = starts[first[c]];
            const Size to = last ? run_end : starts[first[c + 1]];
            chunk_text(first[c + 1] - first[c], chunk);
            const Size at = chunk.size();
            chunk.resize(at + (to - from));
            file.seekg(static_cast<std::streamoff>(from));
            file.read(chunk.data() + at, static_cast<std::streamsize>(to - from));
            parsed = static_cast<bool>(file);
            if (parsed)
            {
              chunk += last ? end_of_run : "</spectrumList>" + end_of_run;
              MzMLFile mzml;
              mzml.getOptions() = chunk_options;
              PeakMap part;
              mzml.loadBuffer(chunk, part);
              parsed = part.size() <= first[c + 1] - first[c];
              if (parsed) { removeUnkept(part.getSpectra(), keep); }
              if (parsed)
              {
                for (MSSpectrum& spectrum : part) { spectrum.sortByPosition(); } // see loadMS2Spectra()
                kept[c] = part.size();
                if (last) { chromatograms = std::move(part.getChromatograms()); }
                parts[c] = std::move(part.getSpectra());
              }
            }
          }
          catch (...)
          {
            parsed = false;
          }
          if (!parsed)
          {
            chunks_ok = false;
            file.clear();
          }
        }
      }
      };
      {
        // closed (waiting for the chunks of the helper) before what the chunks write to goes, also on an exception
        const OnScopeExit end_help([help]() { if (help != nullptr) { help->close(); } });
#ifdef _OPENMP
        // At least 32 parsers together (each holds a chunk of up to 8 MB): the 32 chunks of a file of up to 256 MB in
        // one round. A reader alone has at most 16 threads (no faster with more), but a thread that waits for the
        // spectra has nothing else to do (at 64 threads, where the reader has 16 threads of its own and the index build
        // waits for it: the chunks of astral_B1 parsed in 91 ms instead of 136 ms, those of hfx_B1 in 62 ms as before).
        // Up to all the threads of the search (joinThreads(): the reader's and the waiting thread's, as many as it was
        // given), as a large file has many more chunks than parsers: the 600 chunks of a 4.9 GB timsTOF HeLa run were
        // parsed in 2.2 s by 32 parsers and in 1.4 s by 64 (once StringManager::parseInt() no longer serialised them).
        constexpr int min_parsers = 32;
        if (help != nullptr)
        {
          help->open([&]()
          {
            const int team = std::min(omp_get_max_threads(), std::max(min_parsers, help->joinThreads(threads)) - threads);
            if (team >= 1) { parse_chunks(team); }
          });
        }
#endif
        parse_chunks(threads);
      }
      if (!ok || !chunks_ok) { return false; }

      // The rest in parallel, also with the threads of a thread that waits for the spectra: the spectra are put into
      // their places in the list, in retention time order as MSExperiment::sortSpectra(false) puts them (equal ones in
      // file order), unless a retention time is NaN (then in file order: the order of sortSpectra() depends on its sort,
      // which loadMS2Spectra() then runs), and MSExperiment::updateRanges(). (A list in file order, sorted afterwards,
      // took 0.29 s on one thread for the 619,000 spectra of a timsTOF HeLa run, which are not in retention time order:
      // 0.2 s to move them along the cycles of the permutation.)
#ifdef _OPENMP
      // Not on a runtime with one number of threads for the process (SysInfo::ompThreadCountIsGlobal()): set
      // here, it would be the number of the thread that waits for the spectra, too. The regions below then run with
      // the process's number next to that thread, for the fraction of a second they take.
      const int own_threads = SysInfo::ompThreadCountIsGlobal() ? 0 : omp_get_max_threads();
      if (own_threads > 0) { omp_set_num_threads(help != nullptr ? help->joinThreads(threads) : threads); }
      const OnScopeExit restore([own_threads]() { if (own_threads > 0) { omp_set_num_threads(own_threads); } });
#endif
      std::vector<Size> offset(n_chunks_used + 1, 0); // the first spectrum of each chunk among the spectra kept
      for (Size c = 0; c < n_chunks_used; ++c) { offset[c + 1] = offset[c] + kept[c]; }
      const Size kept_total = offset.back();
      std::vector<double> rts(kept_total);
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize c = 0; c < static_cast<SignedSize>(n_chunks_used); ++c)
      {
        for (Size j = 0; j < kept[c]; ++j) { rts[offset[c] + j] = parts[c][j].getRT(); }
      }
      const bool any_nan = std::any_of(rts.begin(), rts.end(), [](double rt) { return std::isnan(rt); });
      const bool in_order = std::is_sorted(rts.begin(), rts.end());
      std::vector<Size> place; // the place of the spectrum with file position i; empty: i
      if (!any_nan && !in_order)
      {
        const std::vector<Size> order = orderByKeyInParallel(rts);
        place.resize(kept_total);
#pragma omp parallel for schedule(static)
        for (SignedSize k = 0; k < static_cast<SignedSize>(kept_total); ++k) { place[order[k]] = static_cast<Size>(k); }
      }
      std::vector<MSSpectrum> spectra;
      resizeWithParallelFaults(spectra, kept_total);
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize c = 0; c < static_cast<SignedSize>(n_chunks_used); ++c)
      {
        for (Size j = 0; j < kept[c]; ++j)
        {
          const Size i = offset[c] + j;
          spectra[place.empty() ? i : place[i]] = std::move(parts[c][j]);
        }
        std::vector<MSSpectrum>().swap(parts[c]);
      }
      exp.setSpectra(std::move(spectra));
      exp.setChromatograms(std::move(chromatograms));
      exp.setLoadedFileType(filename);
      exp.setLoadedFilePath(filename);
      exp.updateRanges();
      return true;
    }
  }

  namespace Internal::ProSEAlgorithmHelpers
  {
    // The MS2 spectra (with @p survey_scans also the MS1 spectra) of a spectrum file (mzML, mzPeak, Bruker .d or Thermo .raw), sorted by RT, read with
    // @p threads OpenMP threads: an mzML file in chunks parsed in parallel (loadMzMLChunked) with two
    // or more (or @p chunk_threads of two or more, with a helper: @p help, see ChunkHelp), otherwise with FileHandler.
    // With @p keep only the spectra it accepts (see ProSEAlgorithm::loadSpectra()).
    PeakMap loadMS2Spectra(const std::string& filename, int threads, ChunkHelp* help, int chunk_threads,
                           bool survey_scans, const std::function<bool(const MSSpectrum&)>& keep)
    {
      PeakMap spectra;
      // (of a remote location, e.g. s3://bucket/run.mzML?region=..., by its name without the query; such a file is
      // streamed by MzMLFile below, as loadMzMLChunked() reads local files only)
      const FileTypes::Type type = FileHandler::getTypeByFileName(FileHandler::displayName(filename));
      FileHandler f;
      f.getOptions().clearMSLevels();
      f.getOptions().addMSLevel(2);
      if (survey_scans) { f.getOptions().addMSLevel(1); } // chimeric:enabled (see ProSEAlgorithm::search(spectra, ctx, ...))
      // MzMLHandler decodes the spectra in a parallel region every maxDataPoolSize spectra (100 by
      // default): with several threads, fewer and larger regions
      if (threads > 1) { f.getOptions().setMaxDataPoolSize(1000); }
      const bool chunked = (threads > 1 || (help != nullptr && chunk_threads > 1)) && type == FileTypes::MZML && !FileHandler::isRemote(filename)
                           && loadMzMLChunked(filename, f.getOptions(), threads, spectra, help, chunk_threads, keep);
      if (chunked)
      {
        ChromatogramTools().convertSpectraToChromatograms<PeakMap>(spectra, true); // as FileHandler::loadExperiment()
        // No malloc_trim(0) here: the pages that the reader threads freed in their arenas (109-165 MB with 16 readers)
        // stay with the process. Returning them cost 5-6 % wall time on 20 small files at 64 threads (int14: 0.744 vs
        // 0.700 of int12-c2, 3 rounds) for 0.20 GiB less max RSS.
      }
      else if (keep && (type == FileTypes::MZML || type == FileTypes::MZPEAK))
      {
        // streamed, so that only the kept spectra are held
        spectra.clear(true);
        MSDataTransformingConsumer consumer;
        // (moved: the reader does not use a spectrum once its consumer has it)
        consumer.setSpectraProcessingFunc([&spectra, &keep](MSSpectrum& spectrum) { if (keep(spectrum)) { spectra.addSpectrum(std::move(spectrum)); } });
        f.transformExperiment(filename, &consumer, {FileTypes::MZML, FileTypes::MZPEAK});
      }
      else
      {
        // An mzPeak file is read whole here (MzPeakFile decodes its spectra with the OpenMP threads of this thread,
        // and skips the MS1 spectra before decoding them); not in the background like an mzML file (see search()).
        f.loadExperiment(filename, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW, FileTypes::MZPEAK});
        removeUnkept(spectra.getSpectra(), keep); // (Bruker .d, Thermo .raw, mzPeak)
      }
      // by retention time, and the peaks of each spectrum by m/z: the chunks of loadMzMLChunked() sorted the peaks of
      // their spectra on their threads (11 ms on one thread for the 16 million peaks of an Astral run); spectra listed in
      // retention time order (as mzML files usually list them) are not sorted again (see inRetentionTimeOrder())
      if (!chunked || !inRetentionTimeOrder(spectra)) { spectra.sortSpectra(!chunked); }
      return spectra;
    }

    // The threads that a BackgroundConsumer (below) lends from this thread while it runs: the read policies below count
    // them as this thread's (the threads this thread was given), so that a read decides as without the consumer.
    thread_local int consumer_lent_threads = 0;
    thread_local int reader_lent_threads = 0;   // lent to a background reader (BackgroundSpectra::start()), until taken back

    bool helpersOversubscribe()
    {
#ifdef _OPENMP
      static const int cpus = omp_get_num_procs(); // the logical CPUs available to the process
      return threadsGiven() >= cpus;
#else
      return false;
#endif
    }
    thread_local BackgroundSpectra* lending_read = nullptr;  // the read that lends threads of this thread, if any
    thread_local BackgroundSpectra* deferred_read = nullptr; // the read that waits for startDeferredRead(), if any
    thread_local BackgroundConsumer* lending_consumer = nullptr; // the consumer that lends threads of this thread, if any

    void takeBackThreadsOfEndedRead()
    {
#ifdef _OPENMP
      if (lending_read != nullptr && parallelLevel() == 0) { lending_read->takeBackThreadsIfRead(); } // see start()
      if (lending_consumer != nullptr && parallelLevel() == 0) { lending_consumer->takeBackThreadsIfDone(); } // see BackgroundConsumer
#endif
    }

    void startDeferredRead()
    {
      if (deferred_read != nullptr) { deferred_read->startNow(); }
    }

    bool readLendsThreads() { return lending_read != nullptr; }
  } // namespace Internal::ProSEAlgorithmHelpers

  PeakMap ProSEAlgorithm::loadSpectra(const std::string& filename, const std::function<bool(const MSSpectrum&)>& keep)
  {
    return loadMS2Spectra(filename, readerThreads(false), nullptr, 0, false, keep);
  }
} // namespace OpenMS
