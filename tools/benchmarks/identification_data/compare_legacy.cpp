// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <bit>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace OpenMS;
using Clock = std::chrono::steady_clock;
using Native = IdentificationData;

struct Digest
{
  UInt64 rows = 0, hash = 0;
  double scores = 0;
  void add(const std::string& sequence, const std::string& observation, Int charge, double score, double rt, double mz)
  {
    UInt64 value = 14695981039346656037ULL;
    for (const auto& text : {sequence, observation})
    {
      for (unsigned char c : text)
      {
        value ^= c;
        value *= 1099511628211ULL;
      }
      value ^= 255;
      value *= 1099511628211ULL;
    }
    for (UInt64 n : {UInt64(charge), std::bit_cast<UInt64>(score), std::bit_cast<UInt64>(rt), std::bit_cast<UInt64>(mz)})
    {
      value ^= n;
      value *= 1099511628211ULL;
    }
    ++rows;
    hash += value;
    scores += score;
  }
};

void timing(const std::string& phase, Clock::time_point start)
{ std::cout << "TIME\t" << phase << '\t' << std::setprecision(12) << std::chrono::duration<double>(Clock::now() - start).count() << '\n'; }
void emit(const Digest& d)
{ std::cout << "DIGEST\t" << d.rows << '\t' << d.hash << '\t' << std::setprecision(17) << d.scores << '\n'; }
Digest digest(const PeptideIdentificationList& peptides)
{
  Digest d;
  for (const auto& q : peptides)
    for (const auto& h : q.getHits())
      d.add(h.getSequence().toString(), q.getSpectrumReference(), h.getCharge(), h.getScore(), q.getRT(), q.getMZ());
  return d;
}
Digest digest(const Native& data)
{
  Digest d;
  for (const auto& r : data.getRuns())
    for (const auto& s : r.getSources())
      for (const auto& q : s.identifications)
        for (const auto& h : q.getMatches())
          d.add(h.representation, q.data_id, h.charge, h.getScoreValues()[r.getPrimaryScore()->value], *q.rt, *q.mz);
  return d;
}

void generate(UInt64 rows, Size runs, std::vector<ProteinIdentification>& proteins, PeptideIdentificationList& peptides)
{
  if (! runs || rows < runs) throw std::runtime_error("Require rows >= runs > 0");
  std::vector<AASequence> sequences;
  for (UInt64 i = 0; i < 10000; ++i)
  {
    UInt64 n = i;
    std::string sequence = "PEP";
    constexpr char alphabet[] = "ACDEFGHIKLMNPQRSTVWY";
    for (unsigned k = 0; k < 4; ++k)
    {
      sequence += alphabet[n % 20];
      n /= 20;
    }
    sequence += i % 10 == 0 ? "M(Oxidation)CK" : (i % 10 == 1 ? "MC(Carbamidomethyl)K" : "MCK");
    sequences.push_back(AASequence::fromString(sequence));
  }
  peptides.reserve(rows);
  UInt64 global = 0;
  for (Size r = 0; r < runs; ++r)
  {
    ProteinIdentification protein;
    protein.setIdentifier("search-" + std::to_string(r));
    DateTime date;
    date.set("2026-10-05 12:00:00");
    protein.setDateTime(date);
    protein.setSearchEngine("comparison");
    protein.setSearchEngineVersion("1");
    protein.setScoreType("search score");
    protein.setHigherScoreBetter(true);
    protein.setPrimaryMSRunPath({"/data/raw-" + std::to_string(r) + ".mzML"});
    ProteinIdentification::SearchParameters params;
    params.db = "comparison.fasta";
    params.charges = "2:3";

    params.variable_modifications = {"Oxidation (M)", "Carbamidomethyl (C)"};
    protein.setSearchParameters(params);
    const auto local_rows = rows / runs + (r < rows % runs);
    for (UInt64 p = 0; p < std::min<UInt64>(200, local_rows); ++p)
    {
      ProteinHit hit;
      hit.setAccession("P" + std::to_string((global + p) % 200));
      hit.setScore(0);
      hit.setMetaValue("target_decoy", "target");
      protein.insertHit(hit);
    }
    for (UInt64 i = 0; i < local_rows; ++i, ++global)
    {
      PeptideIdentification q;
      q.setIdentifier(protein.getIdentifier());
      q.setScoreType("search score");
      q.setHigherScoreBetter(true);
      q.setRT(global * 0.25);
      q.setMZ(400.0 + (global % 100) / 4.0);
      q.setSpectrumReference("scan=" + std::to_string(global + 1));
      PeptideHit hit;
      hit.setSequence(sequences[global % sequences.size()]);
      hit.setCharge(2 + global % 2);
      hit.setScore(10 + global % 100);
      hit.setRank(1);
      hit.setMetaValue("engine", "comparison");
      hit.setMetaValue("PEP", (global % 100) / 128.0);
      hit.setMetaValue("target_decoy", "target");
      hit.addPeptideEvidence(PeptideEvidence("P" + std::to_string(global % 200), 0, 9, 'K', 'A'));
      q.insertHit(hit);
      peptides.push_back(std::move(q));
    }
    proteins.push_back(std::move(protein));
  }
}

int main(int argc, char** argv)
{
  try
  {
    if (argc < 4)
      throw std::runtime_error("Usage: IdentificationDataLegacyBenchmark write FORMAT PATH ROWS RUNS [THREADS] | read FORMAT PATH [THREADS] | "
                               "convert FORMAT PATH INPUT.idXML [THREADS]");
    const std::string mode = argv[1], format = argv[2], path = argv[3];
    const int required = mode == "write" ? 6 : (mode == "convert" ? 5 : 4);
    if (argc != required && argc != required + 1) throw std::runtime_error("Bad argument count");
    IdentificationDataFile::Options native_options;
    if (argc == required + 1)
    {
      if (format != "native") throw std::runtime_error("Thread argument is only supported for native I/O");
      native_options.threads = std::stoull(argv[required]);
      if (! native_options.threads || native_options.threads > static_cast<Size>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Thread count must be positive and fit in int");
    }
    if (mode == "write" || mode == "convert")
    {
      if (std::filesystem::exists(path)) throw std::runtime_error("Existing output");
      std::vector<ProteinIdentification> proteins;
      PeptideIdentificationList peptides;
      auto start = Clock::now();
      if (mode == "write") generate(std::stoull(argv[4]), std::stoull(argv[5]), proteins, peptides);
      else
        IdXMLFile().load(argv[4], proteins, peptides);
      timing("prepare", start);
      emit(digest(peptides));
      if (format == "native")
      {
        start = Clock::now();
        auto data = IdentificationDataAdapter::fromLegacy(proteins, peptides);
        timing("convert", start);
        emit(digest(data));
        proteins.clear();
        peptides.clear();
        start = Clock::now();
        IdentificationDataFile::store(path, data, native_options);
        timing("write", start);
      }
      else if (format == "idxml")
      {
        start = Clock::now();
        IdXMLFile().store(path, proteins, peptides);
        timing("write", start);
      }
      else
        throw std::runtime_error("Unknown format");
    }
    else if (mode == "read")
    {
      auto start = Clock::now();
      if (format == "native")
      {
        Native data;
        IdentificationDataFile::load(path, data, native_options);
        timing("read", start);
        emit(digest(data));
      }
      else
      {
        std::vector<ProteinIdentification> proteins;
        PeptideIdentificationList peptides;
        if (format == "idxml") IdXMLFile().load(path, proteins, peptides);
        else
          throw std::runtime_error("Unknown format");
        timing("read", start);
        emit(digest(peptides));
      }
    }
    else
      throw std::runtime_error("Unknown operation");
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
