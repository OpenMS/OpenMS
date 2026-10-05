// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IdentificationDataInference.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace
{
using ID = OpenMS::IdentificationData;
using IO = OpenMS::IdentificationDataFile;
using Clock = std::chrono::steady_clock;

void report(const std::string& operation, Clock::time_point start, OpenMS::UInt64 rows, double checksum)
{
  const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
  std::cout << operation << '\t' << rows << '\t' << std::setprecision(12) << seconds << '\t' << checksum << '\n';
}

OpenMS::UInt64 count(const std::vector<IO::RunDescriptor>& runs)
{
  OpenMS::UInt64 rows = 0;
  for (const auto& run : runs)
    rows += run.match_count;
  return rows;
}

struct ScanResult
{
  OpenMS::UInt64 queries = 0;
  OpenMS::UInt64 matches = 0;
  double checksum = 0;
};

ScanResult scan(const std::string& path, bool full, const std::string& operation)
{
  IO::ScanOptions options;
  if (! full)
  {
    options.projection.molecule = false;
    options.projection.evidence = false;
    options.projection.annotations = false;
    options.projection.metadata = false;
    // Keep all numeric scores so heterogeneous real runs, including scoreless
    // empty runs, need no invented common score index.
  }
  ScanResult result;
  const auto start = Clock::now();
  const auto statistics = IO::scan(path, options, {}, [&](const auto&, const auto& batch) {
    result.matches += batch.size();
    for (const auto& match : batch)
    {
      if (! match.scores.empty() && match.scores[0]) result.checksum += *match.scores[0];
      if (full && match.data.representation.empty()) throw std::runtime_error("Missing projected molecule");
    }
  });
  result.queries = statistics.queries;
  report(operation, start, result.matches, result.checksum);
  if (result.matches != statistics.matches || result.matches != count(IO::inspect(path))) throw std::runtime_error("Scanned row count differs");
  std::cout << "# " << operation << " resident_manifest_bytes=" << statistics.descriptor_bytes << '\n';
  return result;
}

void benchmarkIdXML(const std::string& input, const std::string& output)
{
  ScanResult expected;
  {
    std::vector<OpenMS::ProteinIdentification> proteins;
    OpenMS::PeptideIdentificationList peptides;
    auto start = Clock::now();
    OpenMS::IdXMLFile().load(input, proteins, peptides);
    const auto loaded = Clock::now();
    for (const auto& query : peptides)
      expected.matches += query.getHits().size();
    // Counting is outside the XML parser timing.
    std::cout << "idxml-load\t" << expected.matches << '\t' << std::setprecision(12) << std::chrono::duration<double>(loaded - start).count()
              << "\t0\n";
    expected.queries = peptides.size();

    start = Clock::now();
    auto data = OpenMS::IdentificationDataAdapter::fromLegacy(proteins, peptides);
    report("import", start, expected.matches, 0);
    OpenMS::UInt64 imported_queries = 0, imported_matches = 0, sources = 0;
    for (const auto& run : data.getRuns())
    {
      imported_queries += run.getNumberOfIdentifications();
      imported_matches += run.getNumberOfMatches();
      sources += run.getSourceBlocks().size();
      for (const auto& source : run.getSourceBlocks())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches())
            if (! match.getScoreValues().empty() && ! std::isnan(match.getScoreValues()[0])) expected.checksum += match.getScoreValues()[0];
    }
    if (imported_queries != expected.queries || imported_matches != expected.matches)
      throw std::runtime_error("Import changed query or match counts");
    std::cout << "# runs=" << data.getRuns().size() << " sources=" << sources << " queries=" << expected.queries << " matches=" << expected.matches
              << " inference_results=" << data.getInferenceResults().size() << '\n';

    // Release the legacy input before storing; owning import peak RSS still
    // includes both representations and must not be called streaming memory.
    proteins.clear();
    proteins.shrink_to_fit();
    peptides.clear();
    peptides.shrink_to_fit();
    start = Clock::now();
    IO::store(output, data);
    report("store", start, expected.matches, 0);
  }

  OpenMS::UInt64 bytes = 0, parquet_files = 0;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(output))
    if (entry.is_regular_file())
    {
      bytes += entry.file_size();
      parquet_files += entry.path().extension() == ".parquet";
    }
  std::cout << "# input_bytes=" << std::filesystem::file_size(input) << " native_bytes=" << bytes << " parquet_files=" << parquet_files << '\n';

  // These scans are deliberately warm after writing. Run the standalone scan
  // modes in fresh processes to measure their own peak RSS.
  for (bool full : {false, true})
  {
    const auto observed = scan(output, full, full ? "full-scan" : "scan");
    if (observed.queries != expected.queries || observed.matches != expected.matches || observed.checksum != expected.checksum)
      throw std::runtime_error("Native scan changed query/match counts or the score checksum");
  }
}

ID synthetic(OpenMS::UInt64 rows, OpenMS::UInt64 run_count, bool inference = false)
{
  if (! run_count || rows < run_count) throw std::invalid_argument("Require matches >= runs > 0");
  ID data;
  for (OpenMS::UInt64 r = 0; r < run_count; ++r)
  {
    auto& run = data.addRun("search-" + std::to_string(r));
    ID::SourceFile source;
    source.identifier = "raw-" + std::to_string(r);
    source.path = "/data/run-" + std::to_string(r) + ".mzML";
    const auto source_id = run.addSource(source);
    ID::ScoreDefinition raw;
    raw.name = "search score";
    raw.software = "synthetic benchmark";
    const auto raw_id = run.addScore(raw);
    ID::ScoreDefinition pep;
    pep.name = "Posterior Error Probability";
    pep.higher_better = false;
    pep.calibration = "synthetic";
    run.addScore(pep);
    const auto local_rows = rows / run_count + (r < rows % run_count);
    for (OpenMS::UInt64 i = 0; i < local_rows; ++i)
    {
      ID::Observation query;
      query.data_id = "scan=" + std::to_string(i + 1);
      query.rt = i * 0.2;
      query.mz = 500.0 + (i % 200);
      const auto query_id = run.addIdentification(source_id, query);
      ID::MatchData match;
      match.representation = i % 2 ? "PEPTIDEK" : "PEPTIDER";
      match.charge = 2 + i % 2;
      match.target_decoy = i % 10 ? ID::TargetDecoy::TARGET : ID::TargetDecoy::DECOY;
      match.parent_evidence.push_back({{"synthetic.fasta", "P" + std::to_string(i % 2000)}, 1, 8, "K", "A"});
      if (inference)
      {
        // Repeated peptidoforms retain the same mapping across every input run.
        constexpr char alphabet[] = "ACDEFGHIKLMNPQRSTVWY";
        auto key = i % 10000;
        match.representation = "PEPTIDE";
        for (unsigned position = 0; position < 4; ++position)
        {
          match.representation += alphabet[key % 20];
          key /= 20;
        }
        match.representation += 'K';
        match.parent_evidence.front().parent.accession = "P" + std::to_string((i % 10000) / 5);
      }
      match.setMetaValue("rank", 1);
      match.setMetaValue("engine", "synthetic");
      run.addMatch(query_id, match, {10.0 + (i % 100), (i % 100) / 1000.0});
    }
    run.setPrimaryScore(raw_id);
  }
  return data;
}
} // namespace

int main(int argc, char** argv)
{
  try
  {
    if (argc < 3)
    {
      std::cerr << "Usage: IdentificationDataBenchmark generate DATASET MATCHES RUNS\n"
                   "       IdentificationDataBenchmark idxml INPUT.idXML OUTPUT.idparquet\n"
                   "       IdentificationDataBenchmark scan|full-scan|load DATASET\n"
                   "       IdentificationDataBenchmark filter DATASET OUTPUT\n"
                   "       IdentificationDataBenchmark inference MATCHES RUNS\n";
      return 2;
    }
    const std::string mode = argv[1], path = argv[2];
    if (mode == "inference")
    {
      if (argc != 4) throw std::invalid_argument("inference needs MATCHES RUNS");
      const auto rows = std::stoull(path), runs = std::stoull(argv[3]);
      auto start = Clock::now();
      auto data = synthetic(rows, runs, true);
      report("construct", start, rows, 0);
      std::vector<OpenMS::IdentificationDataInference::Input> inputs;
      for (const auto& run : data.getRuns())
        inputs.push_back({run.getUuid(), run.getScoreId(1)});
      start = Clock::now();
      const auto result = OpenMS::IdentificationDataInference::infer(data, inputs, "benchmark");
      report("inference", start, rows, result.proteins.getHits().size());
      OpenMS::UInt64 members = 0;
      for (const auto& input : result.inputs)
        members += input.matches.size();
      if (members != rows || result.assignments.size() != rows) throw std::runtime_error("Inference provenance count differs");
    }
    else if (mode == "idxml")
    {
      if (argc != 4) throw std::invalid_argument("idxml needs INPUT OUTPUT");
      benchmarkIdXML(path, argv[3]);
    }
    else if (mode == "generate")
    {
      if (argc != 5) throw std::invalid_argument("generate needs MATCHES RUNS");
      const auto rows = std::stoull(argv[3]), runs = std::stoull(argv[4]);
      auto start = Clock::now();
      auto data = synthetic(rows, runs);
      report("construct", start, rows, 0);
      start = Clock::now();
      IO::store(path, data);
      report("store", start, rows, 0);
      if (count(IO::inspect(path)) != rows) throw std::runtime_error("Stored row count differs");
    }
    else if (mode == "scan" || mode == "full-scan")
    {
      if (argc != 3) throw std::invalid_argument(mode + " needs DATASET");
      scan(path, mode == "full-scan", mode);
    }
    else if (mode == "load")
    {
      const auto start = Clock::now();
      ID data;
      IO::load(path, data);
      OpenMS::UInt64 rows = 0;
      for (const auto& run : data.getRuns())
        rows += run.getNumberOfMatches();
      report(mode, start, rows, 0);
      if (rows != count(IO::inspect(path))) throw std::runtime_error("Loaded row count differs");
    }
    else if (mode == "filter")
    {
      if (argc != 4) throw std::invalid_argument("filter needs OUTPUT");
      const auto start = Clock::now();
      OpenMS::UInt64 retained = 0;
      IO::filter(
        path, argv[3],
        [&](const auto&, const auto& match) {
          const bool keep = match.match_id % 2 == 0;
          retained += keep;
          return keep;
        },
        ID::InferencePolicy::DISCARD);
      report(mode, start, retained, 0);
      if (retained != count(IO::inspect(argv[3]))) throw std::runtime_error("Filtered row count differs");
    }
    else
      throw std::invalid_argument("Unknown operation");
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
