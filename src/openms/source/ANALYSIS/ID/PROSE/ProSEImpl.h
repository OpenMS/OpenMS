// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

// Private implementation header of ProSEAlgorithm: the declarations, types and inline helpers that its source
// files (ProSEAlgorithm*.cpp) share. Not installed.

#pragma once

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEAlgorithm.h>

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <algorithm>
#if defined(_OPENMP) && defined(__GLIBCXX__) // libstdc++'s parallel mode (not libc++)
#include <parallel/algorithm> // __gnu_parallel::nth_element (medianAndMADBySelection)
#define PROSE_GNU_PARALLEL
#endif
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <exception>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

#if defined(__linux__)
  #include <cstdint>
  #include <sys/mman.h> // madvise (resizeWithParallelFaults)
  #include <unistd.h>   // sysconf
#endif

namespace OpenMS
{
  namespace Internal::ProSEAlgorithmHelpers
  {
#ifdef _OPENMP
/// omp_get_level(), which the OpenMP 2.0 of MSVC lacks: there an inactive enclosing region (if(false), one thread)
/// counts as none, which only makes the thread lending of BackgroundSpectra a no-op inside such a region.
inline int parallelLevel()
{
  #if _OPENMP >= 200805
  return omp_get_level();
  #else
  return omp_in_parallel() ? 1 : 0;
  #endif
}
#endif

    // Destroys the elements of a large container with the OpenMP threads, then the container itself. Destroying the
    // spectra, hits and identifications of a search one after the other is a serial step of 10 - 100 ms each at the
    // end of a search; the elements own their memory, so they can be destroyed in any order.
    // In equal static parts per thread: glibc returns a block to the arena of the thread that allocated it, under that
    // arena's lock. The peak lists were made by the static loop of preprocessSpectra_ (the same parts of the same
    // team), so each thread frees into its own arena; in dynamic chunks most frees waited for the lock of another
    // thread's arena (futex system calls: HeLa 25pc_50_25 at 64 threads 0.44 s for the query spectra instead of 0.06).
    // (Only where the memory goes back depends on this: with another team size or a container changed in between,
    // some frees cross arenas, with the same result.)
    template <typename Container>
    void freeInParallel(Container& container)
    {
      const SignedSize n = static_cast<SignedSize>(container.size());
#pragma omp parallel for schedule(static) if (n >= 256)
      for (SignedSize i = 0; i < n; ++i)
      {
        [[maybe_unused]] const auto released = std::move(container[i]);
      }
      Container().swap(container);
    }

    // Destroys the elements of a large container on the OpenMP threads that allocated them (owner(i): the thread
    // number of element i in a team of the same size; other numbers are spread over the team), then the container:
    // as freeInParallel(), for elements that a dynamically scheduled loop made. Each thread looks at all owners and
    // frees its own elements, into its own arena.
    template <typename Container, typename Owner>
    void freeOnOwnerThreads(Container& container, const Owner& owner)
    {
      const SignedSize n = static_cast<SignedSize>(container.size());
#pragma omp parallel if (n >= 256)
      {
#ifdef _OPENMP
        const int team = omp_get_num_threads();
        const int me = omp_get_thread_num();
#else
        const int team = 1;
        const int me = 0;
#endif
        for (SignedSize i = 0; i < n; ++i)
        {
          const int o = owner(i);
          if ((o >= 0 ? o % team : 0) == me) { [[maybe_unused]] const auto released = std::move(container[i]); }
        }
      }
      Container().swap(container);
    }

    inline void freeInParallel(PeakMap& spectra)
    {
      freeInParallel(spectra.getSpectra());
      spectra.clear(true);
    }

    // Resizes @p v to @p n elements as resize() does; on Linux, the pages of the storage that it allocates for them are
    // first faulted in by the OpenMP threads (madvise(MADV_POPULATE_WRITE) on parts of the range, Linux 5.14 or later; no
    // access to the storage itself), so that the serial construction of the elements does not take the page faults of
    // the whole block one after the other (619,000 empty spectra, 465 MB: 0.21 s, then 0.015 s for the faults in parallel
    // and 0.07 s for the construction at 64 threads). Elsewhere, or if the call fails, as resize().
    template <typename T>
    void resizeWithParallelFaults(std::vector<T>& v, const Size n)
    {
#if defined(__linux__)
      if (n > v.capacity())
      {
        v.reserve(n);
        const uintptr_t page = static_cast<uintptr_t>(sysconf(_SC_PAGESIZE));
        const uintptr_t begin = (reinterpret_cast<uintptr_t>(v.data() + v.size()) + page - 1) / page * page;
        const uintptr_t end = reinterpret_cast<uintptr_t>(v.data() + n) / page * page; // whole pages of the new elements
        const SignedSize pages = end > begin ? static_cast<SignedSize>((end - begin) / page) : 0;
        if (pages >= 1024)
        {
          constexpr int populate_write = 23; // MADV_POPULATE_WRITE (EINVAL before Linux 5.14: then the faults come later)
          const SignedSize parts = 256;
#pragma omp parallel for schedule(dynamic, 1)
          for (SignedSize part = 0; part < parts; ++part)
          {
            const uintptr_t from = begin + static_cast<uintptr_t>(pages * part / parts) * page;
            const uintptr_t to = begin + static_cast<uintptr_t>(pages * (part + 1) / parts) * page;
            if (to > from) { madvise(reinterpret_cast<void*>(from), to - from, populate_write); }
          }
        }
      }
#endif
      v.resize(n);
    }

    // Whether MSExperiment::sortSpectra(false) would leave @p spectra as they are: their retention times ascend (equal ones
    // keep their order in its sort) and none is NaN (which would make the order depend on the sort). An O(n) check that
    // saves the O(n log n) sort (0.1-0.2 s for the 250,000 spectra of a timsTOF run).
    inline bool inRetentionTimeOrder(const PeakMap& spectra)
    {
      const MSSpectrum::RTLess rt_less;
      return std::none_of(spectra.begin(), spectra.end(), [](const MSSpectrum& spectrum) { return std::isnan(spectrum.getRT()); })
             && std::adjacent_find(spectra.begin(), spectra.end(),
                                   [&rt_less](const MSSpectrum& a, const MSSpectrum& b) { return rt_less(b, a); }) == spectra.end();
    }
    /// The m/z of the monoisotopic precursor peak of an observed precursor m/z that FragmentIndex matched at the isotope
    /// offset @p isotope_offset (13C spacings, theoretical minus observed: it searches shifted_mass = precursor_mass +
    /// isotope_offset * C13C12, so M_theo ~ M_obs + isotope_offset * C13C12): observed + offset * C13C12 / charge. The
    /// precursor errors of the PSMs (PRECURSOR_ERROR_PPM_USERPARAM), of the calibration and of the run statistics are
    /// those of this m/z.
    inline double monoisotopicPrecursorMZ(double observed_mz, int isotope_offset, double charge)
    {
      return observed_mz + static_cast<double>(isotope_offset) * Constants::C13C12_MASSDIFF_U / charge;
    }
    // Calls f when the scope ends, also when an exception leaves it (f must not throw then).
    template <typename F>
    class OnScopeExit
    {
    public:
      explicit OnScopeExit(F f) : f_(std::move(f)) {}
      OnScopeExit(const OnScopeExit&) = delete;
      OnScopeExit& operator=(const OnScopeExit&) = delete;
      ~OnScopeExit() { f_(); }

    private:
      F f_;
    };

    // Background reads of spectra that lend OpenMP threads of this thread (BackgroundSpectra, below):
    // takeBackThreadsOfEndedRead() gives this thread its threads back if the read has ended (called between the phases
    // of a search, outside of parallel regions); startDeferredRead() starts a read that waits for the rescoring of the
    // search before it; readLendsThreads(): whether a read lends threads of this thread now.
    void takeBackThreadsOfEndedRead();
    void startDeferredRead();
    bool readLendsThreads();
    // helpersOversubscribe(): whether helper threads next to the OpenMP threads of this thread (a reader of its own,
    // a release of dropped hits or spectra) would exceed the logical CPUs of the process. With as many threads as CPUs,
    // the threads that spin at a barrier take the CPUs of the threads that work: at 128 threads on 64 cores with SMT
    // the 8,000-spectrum runs of the benchmark took twice the wall of 64 threads (1.9 s with OMP_WAIT_POLICY=passive,
    // 3.1 s with the default spinning, 1.5 s at 64 threads). The helpers then take threads of this thread instead.
    bool helpersOversubscribe();
    // The first exception that a parallel loop caught (lowest index; others are dropped), rethrown after the loop: an
    // exception must not leave an OpenMP region
    struct FirstError
    {
      std::exception_ptr error;
      SignedSize index = std::numeric_limits<SignedSize>::max();
      void set(SignedSize i)
      {
#pragma omp critical (ProSEAlgorithm_firstError)
        if (i < index)
        {
          index = i;
          error = std::current_exception();
        }
      }
      void rethrow() const
      {
        if (error) std::rethrow_exception(error);
      }
    };

    // A cleavage rule of the form "(?<=[KR])" with an optional "(?!P)" (Trypsin, Trypsin/P, Lys-C, Arg-C, ...; as
    // FragmentIndex's regex-free digest): a site after a residue of the first set unless the next one blocks it
    struct SegmentRule
    {
      std::array<bool, 256> after{};
      std::array<bool, 256> not_before{};
    };

    // True (and @p rule filled) if @p regex is exactly "(?<=[" + letters A-Z + "])", optionally followed by "(?!" + one
    // letter + ")"; every other rule is left to the regex
    inline bool parseSegmentRule(const std::string& regex, SegmentRule& rule)
    {
      const auto is_letter = [](char c) { return c >= 'A' && c <= 'Z'; };
      if (regex.compare(0, 5, "(?<=[") != 0) return false;
      Size i = 5;
      for (; i < regex.size() && is_letter(regex[i]); ++i) rule.after[static_cast<unsigned char>(regex[i])] = true;
      if (i == 5 || regex.compare(i, 2, "])") != 0) return false;
      i += 2;
      if (i == regex.size()) return true;
      if (regex.size() != i + 5 || regex.compare(i, 3, "(?!") != 0 || !is_letter(regex[i + 3]) || regex[i + 4] != ')') return false;
      rule.not_before[static_cast<unsigned char>(regex[i + 3])] = true;
      return true;
    }
    // The per-spectrum containers of scoreSpectraAgainstIndex_(), one set per thread for the whole run: the candidates
    // of the query and the theoretical spectrum of the scorer. Constructed and destroyed once per spectrum, they
    // returned their buffers to the allocator 8,000 times per thread and took fresh ones; kept per thread and cleared
    // between spectra, they keep the capacity the largest spectrum of the thread needed and allocate nothing else. The
    // query appends to the candidates, so they are cleared before every query. (A function-local thread_local, not one
    // in the loop body, so that the OpenMP data-sharing clauses of the loop do not see it.)
    struct ScoringScratch_
    {
      FragmentIndex::SpectrumMatchesTopN top_sms;
      PeakSpectrum theo_spectrum;
      std::vector<uint8_t> theo_ion_series;
    };
    inline ScoringScratch_& scoringScratch_()
    {
      thread_local ScoringScratch_ scratch;
      return scratch;
    }

    // A helper thread that releases spectra, joined when this goes, after the release held before (see
    // releaseSpectraLater_())
    struct SpectraRelease
    {
      std::shared_ptr<void> before;
      std::thread thread;

      ~SpectraRelease()
      {
        if (thread.joinable()) { thread.join(); }
      }
    };

    Param rescoringPercolatorParameters(Param pp, int threads, const ProSEAlgorithm::PercolatorTraining& training,
                                        const Param& extra);

    std::string rescoringEnzyme(const std::string& enzyme, const ProSEAlgorithm::PercolatorTraining& training);

    void annotateStratumQValues(PeptideIdentificationList& peptide_ids);

    // Lets the thread that waits for a spectra reader (BackgroundSpectra::get()) parse chunks of its loadMzMLChunked()
    // with its own OpenMP threads, which would idle otherwise, instead of only waiting for the end of the read.
    class ChunkHelp
    {
    public:
      /// By the reader: from now on, help() runs @p work (which parses chunks until none is left) ...
      void open(std::function<void()> work)
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!done_) { work_ = std::move(work); }
        changed_.notify_all();
      }

      /// ... until the reader closes: waits for the helpers that run it. Also by the reader once it reads no more chunks.
      void close() noexcept
      {
        std::unique_lock<std::mutex> lock(mutex_);
        work_ = nullptr;
        done_ = true;
        changed_.notify_all();
        changed_.wait(lock, [this]() { return helpers_ == 0; });
      }

      /// By the waiting thread: runs the work (with the OpenMP threads of this thread) if the reader opens it, or returns
      /// once the reader closes. The thread then waits for the spectra.
      void help()
      {
        std::function<void()> work;
        {
          std::unique_lock<std::mutex> lock(mutex_);
#ifdef _OPENMP
          waiting_threads_ = omp_get_max_threads();
#else
          waiting_threads_ = 1;
#endif
          changed_.wait(lock, [this]() { return work_ || done_; });
          if (!work_) { return; }
          work = work_;
          ++helpers_;
        }
        const OnScopeExit leave([this]()
        {
          std::lock_guard<std::mutex> lock(mutex_);
          --helpers_;
          changed_.notify_all();
        });
        work();
      }

      /// By the reader: whether the reader took its threads from the thread that waits (see readBorrowsThreads())
      void setLent(const bool lent) { lent_ = lent; }

      /// By the reader, after the chunks: the OpenMP threads for its remaining parallel work, its own @p reader_threads
      /// and, once the thread that waits for the spectra does (help() was called; it does nothing else until it has
      /// them), those of that thread, too (if they were lent: the threads of both; otherwise as many as the larger team)
      int joinThreads(const int reader_threads)
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (waiting_threads_ == 0) { return reader_threads; }
        return lent_ ? reader_threads + waiting_threads_ : std::max(reader_threads, waiting_threads_);
      }

    private:
      std::mutex mutex_;
      std::condition_variable changed_;
      std::function<void()> work_;
      int helpers_ = 0;
      bool done_ = false;
      int waiting_threads_ = 0; // OpenMP threads of the thread that waits for the spectra; 0: it does not wait yet
      bool lent_ = false;
    };

    void loadFASTA(const std::string& filename, std::vector<FASTAFile::FASTAEntry>& entries);

    PeakMap loadMS2Spectra(const std::string& filename, int threads, ChunkHelp* help = nullptr, int chunk_threads = 0,
                           bool survey_scans = false, const std::function<bool(const MSSpectrum&)>& keep = {});
    extern thread_local int consumer_lent_threads;
    extern thread_local int reader_lent_threads;

    // The OpenMP threads of this thread, with those that a BackgroundConsumer has lent
    [[maybe_unused]] inline int threadsGiven()
    {
#ifdef _OPENMP
      return omp_get_max_threads() + consumer_lent_threads + reader_lent_threads; // (the budget, loans included)
#else
      return 1;
#endif
    }

    // The OpenMP threads a spectrum file is read with, of those of this thread: at most 16 (the chunked read gets no
    // faster with more). A read in the foreground gets all of them. A read in the background (BackgroundSpectra) takes
    // its threads from this thread, which helps with the rest of the read when it takes the spectra: next to the index
    // build an eighth of them, at least one (measured at 16 and 64 threads on 8,000-spectrum files of six instrument
    // types, against a reader with a quarter or all of the threads next to all threads of this thread: no more busy
    // threads than these, 12-17% less CPU time, 3% more wall time at 16 threads and 4-5% less at 64; a fixed share
    // without the help: 7-10% more wall time at 16 threads, as either side idles its share once it is done); next to the
    // search of another file a quarter, at least one, which this thread takes back between the phases of the search
    // once the read has ended.
    inline int readerThreads([[maybe_unused]] bool background, [[maybe_unused]] bool next_to_index_build = true)
    {
#ifdef _OPENMP
      const int threads = threadsGiven();
      if (!background) { return std::min(threads, 16); }
      return std::clamp(threads / (next_to_index_build ? 8 : 4), 1, 16);
#else
      return 1;
#endif
    }

    // The threads whose number sets the chunks of a read in the background (see loadMzMLChunked()): as many as such a
    // read had before it took its threads from those of this thread (a quarter, at least two; or all, at most 16, if
    // the index build waits for it), so that a file is cut into the same chunks as before. The order in which the
    // parsers of the chunks register the names of the meta values they meet first depends on the chunks.
    inline int chunkThreads([[maybe_unused]] bool awaited)
    {
#ifdef _OPENMP
      const int threads = threadsGiven();
      return awaited ? std::min(threads, 16) : std::clamp(threads / 4, std::min(threads, 2), 16);
#else
      return 1;
#endif
    }

    // The spectra of a file that a helper thread reads (loadMS2Spectra()) with some of the OpenMP threads of this thread
    // (see readerThreads()), while this thread runs its parallel regions with the others: together they keep to the
    // threads they were given (-threads). With all of its threads next to the reader, this thread kept up to twice as
    // many busy (29 cores at 16 threads while the FASTA file was read and the peptides generated and sorted), which
    // costs time wherever the process may not use more cores than it was given (a CPU quota, a shared machine). This
    // thread gets its threads back when it takes the spectra (get()) or no longer waits for them (destruction, start()
    // of the next read). The helper sets its OpenMP threads itself: a new thread starts with the default of the process
    // (OMP_NUM_THREADS or all cores), not with what was set for this one (e.g. by TOPPBase from -threads).
    // Whether a read in the background takes its threads from this thread (see BackgroundSpectra): only if this thread
    // has fewer threads than the process has physical cores (SysInfo::getAvailablePhysicalCores()). With as many threads
    // as cores, or more, the reader runs next to them with threads of its own (chunkThreads()), as before the budget
    // (measured on 64 cores, 20 files: at 64 and 128 threads 1-2% more wall time and 1-7% more CPU time than with the
    // budget; at 16 threads unchanged). Never on a runtime that keeps one number of threads for the whole process
    // (SysInfo::ompThreadCountIsGlobal(), MSVC's vcomp): the loans are per-thread numbers, and a helper thread
    // that set its own would set this thread's (a search of three HeLa files with -threads 12 ran on 3 threads: the
    // reader's quarter, which nothing set back; 128 s instead of about 45 s).
    inline bool readBorrowsThreads()
    {
#ifdef _OPENMP
      if (SysInfo::ompThreadCountIsGlobal()) { return false; }
      static const int cores = static_cast<int>(SysInfo::getAvailablePhysicalCores());
      return threadsGiven() < cores || helpersOversubscribe(); // (or the reader's own threads would exceed the CPUs)
#else
      return false;
#endif
    }

    class BackgroundSpectra;
    extern thread_local BackgroundSpectra* lending_read;
    extern thread_local BackgroundSpectra* deferred_read;

    class BackgroundSpectra
    {
    public:
      BackgroundSpectra() = default;

      /// @p survey_scans: the MS1 spectra too (chimeric:enabled, see loadMS2Spectra())
      explicit BackgroundSpectra(const std::string& filename, const bool awaited = false, const bool survey_scans = false)
      {
        start(filename, awaited, 0, survey_scans);
      }

      BackgroundSpectra(const BackgroundSpectra&) = delete;
      BackgroundSpectra& operator=(const BackgroundSpectra&) = delete;

      ~BackgroundSpectra()
      {
        if (deferred_read == this) { deferred_read = nullptr; }
        restoreThreads_();
        join_(); // waits for the helper
      }

      /// Starts to read @p filename, with @p reader_threads (default: readerThreads(true)) of the threads of this thread,
      /// after the read before (if any) has ended and this thread has its threads back. @p awaited: the index build waits
      /// for the spectra (see chunkThreads()). Used by the thread that created it only.
      void start(const std::string& filename, const bool awaited = false, const int reader_threads = 0,
                 const bool survey_scans = false)
      {
        if (deferred_read == this) { deferred_read = nullptr; }
        deferred_ = nullptr;
        restoreThreads_();
        const int chunk_threads = chunkThreads(awaited);
        const bool borrow = readBorrowsThreads();
        const int threads = !borrow ? chunk_threads : reader_threads > 0 ? reader_threads : readerThreads(true);
        join_(); // waits for the helper of a read before
#if defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)
        // WebAssembly without threads has no helper thread: load when the result is collected (get()), without a
        // ChunkHelp, whose help() would wait for a reader that has not started
        constexpr std::launch policy = std::launch::deferred;
        help_ = nullptr;
#else
        constexpr std::launch policy = std::launch::async;
        help_ = std::make_shared<ChunkHelp>();
        help_->setLent(borrow);
#endif
        auto read = [&filename, threads, chunk_threads, survey_scans, help = help_]()
        {
          const OnScopeExit done([&help]() { if (help) { help->close(); } }); // get() waits no longer for chunks to parse
#ifdef _OPENMP
          // the regions of the read take their threads from their clauses; a global number stays the process's
          if (!SysInfo::ompThreadCountIsGlobal()) { omp_set_num_threads(threads); }
#endif
          return loadMS2Spectra(filename, threads, help.get(), chunk_threads, survey_scans);
        };
        if (policy == std::launch::deferred)
        {
          spectra_ = std::async(std::launch::deferred, std::move(read));
        }
        else
        {
          // Not std::async: libc++'s std::async detaches its thread, and in WebAssembly (Emscripten) a detached thread
          // that ends gives its Web Worker back to the pool only when the event loop of the main thread runs, which it
          // does not while a search runs: every file read in the background kept a Worker, and with 4 threads the 7th
          // file of a search waited for one forever. A join (join_()) gives the Worker back at once.
          std::packaged_task<PeakMap()> task(std::move(read));
          spectra_ = task.get_future();
          reader_ = std::thread(std::move(task));
        }
#ifdef _OPENMP
        if (!borrow) { return; }
        own_threads_ = omp_get_max_threads();
        omp_set_num_threads(std::max(1, own_threads_ - threads));
        reader_lent_threads = own_threads_ - std::max(1, own_threads_ - threads);
        // taken back early (takeBackThreadsOfEndedRead()) only outside of parallel regions, where this thread runs a
        // single task, whose ICV this is
        if (parallelLevel() == 0) { lending_read = this; }
#endif
      }

      /// Like start(@p filename, false, readerThreads(true, false)), but only when this thread calls startDeferredRead()
      /// (or get()): during the rescoring of the search before, whose Percolator training uses at most 9 threads.
      void startDeferred(const std::string& filename, const bool survey_scans = false)
      {
#ifdef _OPENMP
        if (parallelLevel() != 0 || !readBorrowsThreads()) { start(filename, false, readerThreads(true, false), survey_scans); return; } // see start()
#endif
        restoreThreads_();
        join_(); // waits for the helper of a read before
        help_ = nullptr;
        deferred_ = &filename;
        deferred_survey_scans_ = survey_scans;
        deferred_read = this;
      }

      /// Starts a deferred read (startDeferred())
      void startNow()
      {
        if (deferred_ != nullptr) { start(*deferred_, false, readerThreads(true, false), deferred_survey_scans_); }
      }

      bool valid() const { return spectra_.valid() || deferred_ != nullptr; }

      /// Whether the spectra are read (valid() only)
      bool ready() const { return deferred_ == nullptr && spectra_.wait_for(std::chrono::seconds(0)) == std::future_status::ready; }

      /// Gives this thread its threads back if the spectra are read
      void takeBackThreadsIfRead()
      {
        if (own_threads_ > 0 && ready()) { restoreThreads_(); }
      }

      /// The spectra (valid() only; waits for them; rethrows what the reader threw). Meanwhile, the threads of this
      /// thread that the reader leaves parse chunks of the file with it (see ChunkHelp).
      PeakMap get()
      {
        startNow();
        {
          const OnScopeExit restore([this]() { restoreThreads_(); }); // also if help() throws
          if (help_) { help_->help(); }
        }
        const OnScopeExit join([this]() { join_(); });
        return spectra_.get();
      }

    private:
      /// Drops the spectra of a read before (if any), after its helper has ended
      void join_() noexcept
      {
        if (reader_.joinable()) { reader_.join(); }
        spectra_ = std::future<PeakMap>();
      }

      void restoreThreads_() noexcept
      {
#ifdef _OPENMP
        if (own_threads_ > 0) { omp_set_num_threads(own_threads_); reader_lent_threads = 0; }
#endif
        own_threads_ = 0;
        if (lending_read == this) { lending_read = nullptr; }
      }

      std::shared_ptr<ChunkHelp> help_;
      std::future<PeakMap> spectra_;
      std::thread reader_; // the helper that reads (none for a deferred read)
      int own_threads_ = 0; // of this thread, while the helper reads; 0: not reduced
      const std::string* deferred_ = nullptr; // the file of a deferred read
      bool deferred_survey_scans_ = false;    // ... read with its MS1 spectra
    };

    // The consumer of the multi-file search's results (ProSEAlgorithm::setFileResultConsumer()), one result at a time, on a
    // helper thread with two of the OpenMP threads of this thread (one if it has fewer than 8), which lends them (lowers
    // its own) until it takes them back where it starts a phase after the helper has ended
    // (takeBackThreadsOfEndedRead()), when it starts the next consumer, or at the end (finish(), destruction). Only if a read in the background never takes threads from
    // this thread (readBorrowsThreads() false): this thread's number of threads is then changed here only, so it can be
    // set back as it was. Otherwise, or with a single thread, or outside of the top level, the consumer runs here. On a
    // runtime that keeps one number of threads for the process (SysInfo::ompThreadCountIsGlobal()) the helper
    // does not set its number with omp_set_num_threads() (that would set this thread's) but with
    // SysInfo::setOmpThreadsOfThisThread(), which the idXML writer and the QPX export size their teams by.
    class BackgroundConsumer;
    extern thread_local BackgroundConsumer* lending_consumer;

    class BackgroundConsumer
    {
    public:
      BackgroundConsumer() = default;
      BackgroundConsumer(const BackgroundConsumer&) = delete;
      BackgroundConsumer& operator=(const BackgroundConsumer&) = delete;
      ~BackgroundConsumer() { finish(); }

      /// Calls @p consumer(@p index, @p result) on a helper thread or here (see above), after the consumer before.
      void start(const ProSEAlgorithm::FileResultConsumer& consumer, const Size index, const ProSEAlgorithm::SearchResult& result)
      {
        finish();
        int helper_threads = 0;
#ifdef _OPENMP
        const int threads = omp_get_max_threads();
        // Two threads, not a quarter: the consumer (the idXML output of ProSE) has the search of the next file to run
        // next to, and with glibc's malloc (default arena settings) each further thread of its team can populate an
        // arena of its own, whose freed memory may remain resident. Search of 6 timsTOF HeLa runs in one process at
        // 64 threads (kim, 64 cores; 16 / 2 helper threads):
        // max RSS 25.9 / 24.5 GiB, wall time 38.2 / 37.2 s, identical output (the idXML of the largest run written in
        // about 2.5 s instead of 1.5 s, still during the search of the next one). The next consumer waits for the one
        // before, so searches shorter than that write wait for it: that run followed by two small ones (0.4 GB) took
        // 17.4 instead of 15.4 s (4 / 8 helper threads: 15.9 / 15.6 s, and less of the memory saved).
        bool lend = true; // the helper's threads are taken from this thread's team
        if (parallelLevel() == 0 && threads > 1 && !readBorrowsThreads()) { helper_threads = std::min(2, std::max(1, threads / 4)); }
        // With as many threads as CPUs the reader borrows threads of this thread (see helpersOversubscribe()), and the
        // two loans cannot be juggled at once: the consumer then gets one thread of its own (the idXML writer is serial),
        // one more than the CPUs, so that the write still overlaps the search of the next file instead of running here.
        else if (parallelLevel() == 0 && threads > 1 && helpersOversubscribe()) { helper_threads = 1; lend = false; }
#endif
        if (helper_threads > 0)
        {
          try
          {
            done_ = false;
            helper_ = std::thread([this, &consumer, index, &result, helper_threads]()
            {
#ifdef _OPENMP
              SysInfo::setOmpThreadsOfThisThread(helper_threads); // (see above; also where the runtime keeps a number per thread)
              if (!SysInfo::ompThreadCountIsGlobal()) { omp_set_num_threads(helper_threads); }
#endif
              try { consumer(index, result); } catch (...) {} // (the consumer handles its errors)
              done_.store(true, std::memory_order_release);
            });
#ifdef _OPENMP
            if (lend)
            {
              own_threads_ = threads;
              omp_set_num_threads(threads - helper_threads);
              consumer_lent_threads = helper_threads;
              lending_consumer = this;
            }
#endif
            return;
          }
          catch (...) {} // no thread: here
        }
        try { consumer(index, result); } catch (...) {}
      }

      /// Gives this thread its threads back if the consumer has ended (at the top level only)
      void takeBackThreadsIfDone()
      {
        if (helper_.joinable() && done_.load(std::memory_order_acquire)) { finish(); }
      }

      /// Waits for the consumer and gives this thread its threads back
      void finish() noexcept
      {
        if (helper_.joinable()) { helper_.join(); }
#ifdef _OPENMP
        if (own_threads_ > 0) { omp_set_num_threads(own_threads_); consumer_lent_threads = 0; }
#endif
        own_threads_ = 0;
        if (lending_consumer == this) { lending_consumer = nullptr; }
      }

    private:
      std::thread helper_;
      std::atomic<bool> done_{false};
      int own_threads_ = 0; // of this thread, while the helper runs; 0: not lowered
    };
  } // namespace Internal::ProSEAlgorithmHelpers

  // static
  struct ProSEAlgorithm::SpectrumGenerators_
  {
    TheoreticalSpectrumGenerator standard; ///< configured ion series
    TheoreticalSpectrumGenerator electron; ///< configured ion series plus c and z+1 ions
    bool by_activation = false;

    /// True if @p spectrum is matched against, and scored with, c and z+1 ions as well: its
    /// precursor was activated by electrons
    bool electronIons(const MSSpectrum& spectrum) const
    {
      return by_activation && isElectronActivated_(spectrum);
    }

    /// The generator for @p spectrum, chosen by the activation of its precursor
    const TheoreticalSpectrumGenerator& forSpectrum(const MSSpectrum& spectrum) const
    {
      return electronIons(spectrum) ? electron : standard;
    }
  };
} // namespace OpenMS
