// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIABoundedQueue.h"
#include "OpenDIARecords.h"
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <latch>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>

using namespace OpenDIA;
namespace
{
  constexpr ByteCount MIB = 1024 * 1024;
  const RecordContext ctx{0, 71, 29, 503};
  void check(bool condition, const char* message)
  {
    if (!condition) { throw std::runtime_error(message); }
  }
  template <typename F> void rejects(F&& action)
  {
    bool rejected = false;
    try { action(); }
    catch (const OpenMS::Exception::InvalidValue&) { rejected = true; }
    check(rejected, "Expected explicit invalid/oversized accounting failure");
  }
  MemoryLease charge(ByteBudget& budget, ByteCount bytes)
  {
    auto reservation = budget.tryReserve(bytes);
    check(reservation.has_value(), "Unexpected admission failure");
    return reservation->commit(bytes);
  }
  void ownership()
  {
    ByteBudget budget(MIB, 0);
    std::vector<RawPoint> points(2);
    points.reserve(8);
    points[0].coordinates.native_id.assign(1024, 'x');
    const auto bytes = RawTile::requiredBytes(points);
    check(bytes >= points.capacity() * sizeof(RawPoint) + 1025, "Capacity/string charge missing");
    std::optional<RawTile> retained;
    {
      RawTile original(ctx, {2, 2}, points, charge(budget, bytes));
      retained = original;
      auto another = original;
      check(budget.snapshot().committed_bytes == bytes, "Shared records double charged");
      another = RawTile(ctx, {}, {});
    }
    check(budget.snapshot().committed_bytes == bytes && retained->size() == 2, "Accounting escaped shared storage");
    BoundedQueue<RawTile> queue(1);
    check(queue.push(std::move(*retained)), "Queue rejected record");
    retained.reset();
    auto last = queue.pop();
    queue.cancel();
    check(budget.snapshot().committed_bytes == bytes, "Dequeue released retained record accounting");
    last.reset();
    check(budget.snapshot().committed_bytes == 0, "Last record leaked accounting");
    auto undercharged = points;
    const auto copied_bytes = RawTile::requiredBytes(undercharged);
    rejects([&] { RawTile invalid(ctx, {2, 2}, std::move(undercharged), charge(budget, copied_bytes - 1)); });
    auto invalid_points = points;
    invalid_points[0].intensity = std::numeric_limits<double>::infinity();
    rejects([&] { RawTile invalid(ctx, {2, 2}, invalid_points, charge(budget, RawTile::requiredBytes(invalid_points))); });
    check(points.size() == 2 && points[0].intensity == 0 && budget.snapshot().committed_bytes == 0,
          "Failure changed input population or leaked reservation");
    {
      std::vector<Hypothesis> hypotheses(1);
      hypotheses[0].provenance.library.assign(1024, 'l');
      hypotheses[0].provenance.prediction = std::string(2048, 'p');
      std::vector<OpenMS::Size> offsets{0, 1};
      std::vector<Fragment> fragments(1);
      const auto page_bytes = HypothesisPage::requiredBytes(hypotheses, offsets, fragments);
      std::optional<HypothesisPage> kept;
      {
        HypothesisPage page(ctx, {1, 1}, hypotheses, offsets, fragments, charge(budget, page_bytes));
        kept = page;
      }
      check(budget.snapshot().committed_bytes == page_bytes && kept->fragmentCount(0) == 1,
            "Hypothesis storage lost accounting");
      kept.reset();
      auto undercharged_hypotheses = hypotheses;
      auto undercharged_offsets = offsets;
      auto undercharged_fragments = fragments;
      const auto copied_page_bytes = HypothesisPage::requiredBytes(undercharged_hypotheses, undercharged_offsets, undercharged_fragments);
      rejects([&] { HypothesisPage bad(ctx, {1, 1}, std::move(undercharged_hypotheses), std::move(undercharged_offsets),
                                      std::move(undercharged_fragments), charge(budget, copied_page_bytes - 1)); });
      rejects([&] { HypothesisPage bad(ctx, {1, 1}, hypotheses, {0, 2}, fragments, charge(budget, page_bytes)); });
      check(budget.snapshot().committed_bytes == 0, "Hypothesis exception leaked bytes");
    }
    {
      std::vector<EvidencePoint> evidence(2);
      evidence[0].value = 0;
      evidence[1].value.reset();
      const auto evidence_bytes = EvidenceBlock::requiredBytes(evidence);
      std::optional<EvidenceBlock> kept;
      {
        EvidenceBlock block(ctx, {2, 2}, evidence, charge(budget, evidence_bytes));
        kept = block;
      }
      check(budget.snapshot().committed_bytes == evidence_bytes && kept->point(0).value == 0 && !kept->point(1).value,
            "Evidence missingness/accounting changed");
      kept.reset();
      auto undercharged_evidence = evidence;
      const auto copied_evidence_bytes = EvidenceBlock::requiredBytes(undercharged_evidence);
      rejects([&] { EvidenceBlock bad(ctx, {2, 2}, std::move(undercharged_evidence), charge(budget, copied_evidence_bytes - 1)); });
      rejects([&] { EvidenceBlock bad(ctx, {1, 1}, evidence, charge(budget, evidence_bytes)); });
    }
    check(budget.snapshot().committed_bytes == 0 && budget.snapshot().reserved_bytes == 0, "Record accounting leaked");
  }
  MobilityUnit mobilityUnit(OpenMS::DriftTimeUnit unit)
  {
    check(unit == OpenMS::DriftTimeUnit::MILLISECOND || unit == OpenMS::DriftTimeUnit::VSSC, "Unknown fixture IM unit");
    return unit == OpenMS::DriftTimeUnit::MILLISECOND ? MobilityUnit::MILLISECOND : MobilityUnit::INVERSE_REDUCED;
  }
  std::vector<RawPoint> flatten(const OpenMS::MSExperiment& experiment)
  {
    std::vector<RawPoint> result;
    OpenMS::Size count = 0;
    for (const auto& spectrum : experiment) { count = checkedByteAdd(count, spectrum.size()); }
    result.reserve(count);
    for (OpenMS::Size scan = 0; scan < experiment.size(); ++scan)
    {
      const auto& spectrum = experiment[scan];
      check(spectrum.isSorted() && spectrum.getPrecursors().size() <= 1, "Unsupported fixture shape");
      Coordinates c;
      c.rt_seconds = spectrum.getRT(); c.ms_level = spectrum.getMSLevel();
      c.native_id = spectrum.getNativeID(); c.scan_id = scan + 2000; c.frame_id = scan + 1000;
      if (!spectrum.getPrecursors().empty())
      {
        const auto& p = spectrum.getPrecursors().front();
        c.isolation_lower_mz = p.getMZ() - p.getIsolationWindowLowerOffset();
        c.isolation_upper_mz = p.getMZ() + p.getIsolationWindowUpperOffset();
      }
      const auto unit = spectrum.getDriftTimeUnit();
      if (unit == OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE) { c.faims_voltage = spectrum.getDriftTime(); }
      else if (unit != OpenMS::DriftTimeUnit::NONE) { c.mobility = Mobility{spectrum.getDriftTime(), mobilityUnit(unit)}; }
      for (OpenMS::Size i = 0; i < spectrum.size(); ++i)
      {
        c.mz = spectrum[i].getMZ();
        if (spectrum.containsIMData())
        {
          const auto [index, im_unit] = spectrum.getIMData();
          c.mobility = Mobility{spectrum.getFloatDataArrays().at(index).at(i), mobilityUnit(im_unit)};
        }
        result.push_back({c, spectrum[i].getIntensity()});
      }
    }
    return result;
  }
  void equal(const Coordinates& a, const Coordinates& b)
  {
    check(a.rt_seconds == b.rt_seconds && a.mz == b.mz && a.frame_id == b.frame_id && a.scan_id == b.scan_id &&
          a.ms_level == b.ms_level && a.isolation_lower_mz == b.isolation_lower_mz && a.isolation_upper_mz == b.isolation_upper_mz &&
          a.faims_voltage == b.faims_voltage && a.native_id == b.native_id && a.mobility.has_value() == b.mobility.has_value(),
          "Pipeline changed acquisition metadata");
    if (a.mobility) { check(a.mobility->value == b.mobility->value && a.mobility->unit == b.mobility->unit, "Pipeline changed IM"); }
  }
  struct RawJob { ByteReservation remaining; RawTile tile; OpenMS::Size offset; };
  struct EvidenceJob { EvidenceBlock block; OpenMS::Size offset; };

  void pipeline(const std::string& relative, bool cancel_stalled)
  {
    constexpr OpenMS::Size CHUNK = 64;
    // Caller-declared fixture decoder/native allowances, not a bounded mzML
    // decoder implementation. External RSS is measured by the test runner.
    constexpr ByteCount NATIVE = 128 * MIB, DECODER = 64 * MIB, WRITER = 65536;
    // A conservative fixture-only native-ID allowance permits admission before
    // decoding. This is not a production decoder or an allocator RSS limit.
    constexpr ByteCount strings = 4096 + 256;
    const auto raw_ceiling = checkedByteAdd(4096, checkedByteMultiply(CHUNK, checkedByteAdd(sizeof(RawPoint), strings)));
    const auto evidence_ceiling = checkedByteAdd(4096, checkedByteMultiply(CHUNK, checkedByteAdd(sizeof(EvidencePoint), strings)));
    const auto scratch_ceiling = checkedByteAdd(4096, checkedByteAdd(checkedByteMultiply(CHUNK, 2 * sizeof(double)), strings));
    const auto job_ceiling = checkedByteAdd(checkedByteAdd(raw_ceiling, evidence_ceiling), scratch_ceiling);
    const auto queue_bytes = 2 * (sizeof(std::optional<RawJob>) + sizeof(std::optional<EvidenceJob>)) + 4096;
    const auto fixed_bytes = DECODER + WRITER + queue_bytes;
    ByteBudget budget(checkedByteAdd(NATIVE, checkedByteAdd(fixed_bytes, checkedByteMultiply(3, job_ceiling))), NATIVE);
    {
      auto fixed = charge(budget, fixed_bytes);
      OpenMS::MSExperiment source;
      OpenMS::MzMLFile().load(std::string(OPENDIA_TEST_DATA_DIR) + "/" + relative, source);
      const auto reference = flatten(source);
      check(reference.size() > CHUNK * 8, "Fixture too small for stalled writer stress");
      for (const auto& point : reference)
      {
        check(point.coordinates.native_id.capacity() <= 4096, "Fixture exceeds declared native-ID allowance");
      }
      // Reserve bounded writer/drain and queue nodes before admitting producers.
      std::vector<char> writer_buffer(WRITER);
      BoundedQueue<RawJob> raw(2);
      BoundedQueue<EvidenceJob> output(2);
      std::stop_source stop;
      std::latch writer_gate(1);
      std::atomic<bool> gate_open = false;
      auto open_gate = [&] { if (!gate_open.exchange(true)) { writer_gate.count_down(); } };
      auto cancel = [&] { stop.request_stop(); budget.cancel(); raw.cancel(); output.cancel(); open_gate(); };
      std::mutex failure_mutex;
      std::exception_ptr failure;
      auto guarded = [&](auto action)
      {
        try { action(); }
        catch (...)
        {
          { std::lock_guard lock(failure_mutex); if (!failure) { failure = std::current_exception(); } }
          cancel();
        }
      };
      const auto filename = std::filesystem::current_path() / ("OpenDIAResourceInterop_" + OpenMS::File::getUniqueName(false) + ".tsv");
      struct RemoveOutput
      {
        std::filesystem::path path;
        ~RemoveOutput() { std::error_code error; std::filesystem::remove(path, error); }
      } cleanup{filename};
      OpenMS::Size written = 0;
      std::jthread writer([&] { guarded([&]
      {
        std::ofstream file;
        file.rdbuf()->pubsetbuf(writer_buffer.data(), writer_buffer.size());
        file.open(filename);
        check(file.good(), "Could not open test writer output");
        file << std::setprecision(std::numeric_limits<double>::max_digits10);
        writer_gate.wait();
        while (auto job = output.pop(stop.get_token()))
        {
          check(job->offset == written, "Writer lost/reordered records");
          for (OpenMS::Size i = 0; i < job->block.size(); ++i)
          {
            const auto row = job->block.point(i);
            const auto& expected = reference.at(written);
            equal(row.coordinates, expected.coordinates);
            check(row.value == expected.intensity, "Writer changed measured zero/missing/value");
            file << written << '\t' << row.coordinates.mz << '\t' << *row.value << '\n';
            ++written;
          }
        }
        file.close();
        check(file.good(), "Writer failed");
      }); });
      std::jthread extractor([&] { guarded([&]
      {
        while (auto job = raw.pop(stop.get_token()))
        {
          auto scratch = job->remaining.split(scratch_ceiling).commit(scratch_ceiling);
          auto evidence_memory = job->remaining.split(evidence_ceiling).commit(evidence_ceiling);
          std::vector<double> mz, intensity;
          std::vector<EvidencePoint> evidence;
          const auto count = job->tile.size();
          mz.reserve(count); intensity.reserve(count); evidence.reserve(count);
          for (OpenMS::Size i = 0; i < count; ++i)
          {
            const auto point = job->tile.point(i);
            mz.push_back(point.coordinates.mz); intensity.push_back(point.intensity);
            evidence.push_back({0, point.coordinates.scan_id, i, point.coordinates, point.intensity});
          }
          // Exercise the actual public contiguous-array extraction boundary.
          auto begin = mz.cbegin(), end = mz.cend(), cursor = begin, values = intensity.cbegin();
          double actual = 0, expected = 0;
          OpenMS::ChromatogramExtractorAlgorithm().extract_value_tophat(begin, cursor, end, values, mz.front(), actual, 0.05, false);
          for (OpenMS::Size i = 0; i < count; ++i)
          {
            if (mz[i] > mz.front() - 0.025 && mz[i] < mz.front() + 0.025) { expected += intensity[i]; }
          }
          check(std::abs(actual - expected) <= 1e-10 * std::max(1.0, std::abs(expected)), "Public extraction changed");
          EvidenceBlock block(ctx, {count, count}, std::move(evidence), std::move(evidence_memory));
          EvidenceJob next{std::move(block), job->offset};
          if (!output.push(std::move(next), stop.get_token())) { break; }
        }
        output.close();
      }); });
      std::jthread producer([&] { guarded([&]
      {
        for (OpenMS::Size offset = 0; offset < reference.size();)
        {
          auto reservation = budget.reserve(job_ceiling, stop.get_token());
          if (!reservation) { break; }
          auto memory = reservation->split(raw_ceiling).commit(raw_ceiling);
          std::vector<RawPoint> points;
          points.reserve(CHUNK);
          auto end = offset;
          while (end < reference.size() && end - offset < CHUNK && reference[end].coordinates.scan_id == reference[offset].coordinates.scan_id)
          { points.push_back(reference[end++]); }
          const auto count = points.size();
          RawTile tile(ctx, {count, count}, std::move(points), std::move(memory));
          RawJob job{std::move(*reservation), std::move(tile), offset};
          if (!raw.push(std::move(job), stop.get_token())) { break; }
          offset = end;
        }
        raw.close();
      }); });
      // Stalling is deliberate, and synchronized by observed full drain plus a
      // blocked byte admission, not a timing-only sleep. Always unwind all threads.
      const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
      bool stalled = false;
      while (!stalled)
      {
        // The gated writer cannot shrink its full queue. Admission waiters can
        // disappear as the extractor releases a raw job, so retain the observed
        // condition instead of checking a second, potentially different snapshot.
        stalled = output.snapshot().size == 2 && budget.snapshot().waiters > 0;
        if (stalled || stop.stop_requested() || std::chrono::steady_clock::now() >= deadline) { break; }
        std::this_thread::yield();
      }
      if (cancel_stalled || !stalled) { cancel(); } else { open_gate(); }
      producer.join(); extractor.join(); writer.join();
      if (failure) { std::rethrow_exception(failure); }
      check(stalled, "Did not exercise full queue and blocked byte admission");
      const auto accounting = budget.snapshot();
      check(accounting.reserved_bytes == 0 && accounting.committed_bytes == fixed_bytes && accounting.peak_total_bytes <= accounting.usable_bytes,
            "Pipeline leaked/exceeded accounted bytes");
      check(raw.snapshot().peak_size <= 2 && output.snapshot().peak_size == 2, "Queue exceeded item cap");
      // Oversized admission fails explicitly. It must not mutate the source set.
      const auto count_before = reference.size();
      rejects([&] { budget.reserve(checkedByteAdd(accounting.usable_bytes, 1)); });
      check(reference.size() == count_before, "Oversized rejection truncated input population");
      if (!cancel_stalled)
      {
        check(written == reference.size(), "Successful drain omitted scientific rows");
        std::ifstream file(filename);
        OpenMS::Size index = 0, row = 0;
        double mz = 0, intensity = 0;
        while (file >> row >> mz >> intensity)
        {
          check(row == index && mz == reference.at(index).coordinates.mz && intensity == reference.at(index).intensity,
                "Reopened writer output changed point population/value");
          ++index;
        }
        check(file.eof() && index == reference.size(), "Writer output was truncated or malformed");
      }
      else { check(written == 0 && accounting.cancelled, "Cancelled stalled writer processed unexpected rows"); }
      OPENMS_LOG_INFO << "WP02 " << relative << " cancelled=" << cancel_stalled << " input_rows=" << reference.size()
                      << " written_rows=" << written << " reserved_peak=" << accounting.peak_reserved_bytes
                      << " committed_peak=" << accounting.peak_committed_bytes << " accounted_peak=" << accounting.peak_total_bytes
                      << " native_allowance=" << NATIVE << " usable_bytes=" << accounting.usable_bytes
                      << " declared_total_bytes=" << accounting.total_bytes << " raw_queue_peak=" << raw.snapshot().peak_size
                      << " writer_queue_peak=" << output.snapshot().peak_size << '\n';
    }
    check(budget.snapshot().committed_bytes == 0 && budget.snapshot().reserved_bytes == 0, "Final pipeline accounting not empty");
  }
}
int main()
{
  try
  {
    ownership();
    for (const auto& fixture : {"topp/OpenSwathWorkflow_1_input.mzML", "topp/OpenSwathWorkflow_23_input.mzML", "class_tests/openms/data/IM_FAIMS_test.mzML"})
    {
      pipeline(fixture, false);
    }
    pipeline("topp/OpenSwathWorkflow_23_input.mzML", true);
  }
  catch (const std::exception& error)
  {
    OPENMS_LOG_ERROR << "Resource integration failed: " << error.what() << '\n';
    return 1;
  }
}
