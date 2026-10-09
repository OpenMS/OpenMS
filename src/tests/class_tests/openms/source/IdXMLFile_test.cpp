// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/TestFileValidation.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <locale>
#include <new>
#include <sstream>
#include <streambuf>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __linux__
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

// The fault-injection sections below (an allocation failure, a change of the working directory, or a file renamed into
// the place of the output while store() runs) replace the global operator new of this test program; the replacement
// must also receive the allocations of libOpenMS. That is the case on Linux (and checked at run time, see
// allocFaultReachesLibOpenMS()); elsewhere, or if OPENMS_TEST_NO_FAULT_INJECTION is defined, they are skipped with a
// message.
#if defined(__linux__) && !defined(OPENMS_TEST_NO_FAULT_INJECTION)
#define FAULT_INJECTION_TESTS 1
#else
#define FAULT_INJECTION_TESTS 0
#endif

#if FAULT_INJECTION_TESTS
// Fault injection for the allocation tests below. Armed on one thread, the n-th allocation of that thread after arming
// runs action() or, without an action, throws std::bad_alloc; with then_throw, it runs action() and then throws
// std::bad_alloc. Allocations made inside action() are not counted. Over-aligned allocations (operator new with
// std::align_val_t) are not replaced and not counted.
namespace AllocFault
{
  thread_local long countdown = 0;
  thread_local bool fired = false;
  thread_local bool in_action = false;
  thread_local bool then_throw = false;
  thread_local void (*action)() = nullptr;
  void arm(long n, void (*on_fire)() = nullptr, bool throw_after_action = false)
  {
    fired = false;
    action = on_fire;
    then_throw = throw_after_action;
    countdown = n;
  }
  void disarm() { countdown = 0; action = nullptr; then_throw = false; }
}
void* operator new(std::size_t size)
{
  if (AllocFault::countdown > 0 && !AllocFault::in_action && --AllocFault::countdown == 0)
  {
    AllocFault::fired = true;
    if (AllocFault::action == nullptr) throw std::bad_alloc();
    AllocFault::in_action = true;
    AllocFault::action();
    AllocFault::in_action = false;
    if (AllocFault::then_throw) throw std::bad_alloc();
  }
  if (void* p = std::malloc(size == 0 ? 1 : size)) return p;
  throw std::bad_alloc();
}
void* operator new[](std::size_t size) { return ::operator new(size); }
void* operator new(std::size_t size, const std::nothrow_t&) noexcept
{
  try { return ::operator new(size); } catch (...) { return nullptr; }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept
{
  try { return ::operator new(size); } catch (...) { return nullptr; }
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
void operator delete(void* p, const std::nothrow_t&) noexcept { std::free(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { std::free(p); }
#endif

///////////////////////////

using namespace OpenMS;

namespace
{
  // registers a tool-defined modification; ModificationsDB is process-wide, so every section uses its own name
  const ResidueModification* defineMod4b(const std::string& id, char origin, const std::string& formula)
  {
    ResidueModification d;
    d.setId(id);
    d.setOrigin(origin);
    d.setTermSpecificity(ResidueModification::ANYWHERE);
    d.setFullId();
    d.setDiffFormula(EmpiricalFormula(formula));
    d.setDiffMonoMass(EmpiricalFormula(formula).getMonoWeight());
    return ModificationsDB::getInstance()->registerDefinition(d);
  }

  // a definition record for a name that is NOT registered in this process
  std::string freshRecord4b(const std::string& id, char origin, const std::string& formula)
  {
    ResidueModification d;
    d.setId(id);
    d.setOrigin(origin);
    d.setTermSpecificity(ResidueModification::ANYWHERE);
    d.setFullId();
    d.setDiffFormula(EmpiricalFormula(formula));
    d.setDiffMonoMass(EmpiricalFormula(formula).getMonoWeight());
    return d.toDefinitionString();
  }

  std::string slurp4b(const std::string& path)
  {
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  bool fileContains4b(const std::string& path, const std::string& needle)
  {
    return slurp4b(path).find(needle) != std::string::npos;
  }

  // exposes the block formatter of the parallel writer
  struct IdXMLFileBlockWriter : public IdXMLFile
  {
    using IdXMLFile::formatBlock_;
  };

  // whether the message of an error of store() names the file and says that a partial file may remain (store() does not
  // remove its output)
  bool saysPartialFileRemains(const std::string& message, const std::string& file)
  {
    return message.find("writing '" + file + "' did not complete: a partial file may remain") != std::string::npos;
  }

  // the message of the exception of type E that f() throws; empty if it throws none or another one
  template <typename E, typename F>
  std::string messageOf(F&& f)
  {
    try
    {
      f();
    }
    catch (const E& e)
    {
      return e.what();
    }
    catch (...)
    {
    }
    return std::string();
  }

  // whether @p file holds a partial idXML: written to, but not to its end
  bool isPartialIdXML(const std::string& file)
  {
    if (!File::exists(file)) return false;
    const std::string content = slurp4b(file);
    return !content.empty() && content.find("</IdXML>") == std::string::npos;
  }

  // a stream buffer that cannot grow, as a std::stringbuf whose allocation fails
  struct BadAllocStreamBuf : public std::streambuf
  {
  protected:
    int_type overflow(int_type) override { throw std::bad_alloc(); }
  };

  // a stream buffer that fails without throwing
  struct FailingStreamBuf : public std::streambuf
  {
  protected:
    int_type overflow(int_type) override { return traits_type::eof(); }
  };

#if FAULT_INJECTION_TESTS
  void noAction() {}

  // whether the replacement of operator new above receives an allocation made inside libOpenMS (it does not, e.g., with
  // a statically linked C++ runtime)
  bool allocFaultReachesLibOpenMS()
  {
    AllocFault::arm(1, &noAction);
    const std::string name = File::getUniqueName(false);
    const bool fired = AllocFault::fired;
    AllocFault::disarm();
    return fired && !name.empty();
  }

  // the working directory test of store(): the action of AllocFault changes to directory B once the output in A exists
  std::filesystem::path cwd_test_dir_b;
  std::filesystem::path cwd_test_out_a;
  bool cwd_test_changed = false;
  void changeToDirBOnceOutputExists()
  {
    std::error_code ec;
    if (!cwd_test_changed && std::filesystem::exists(cwd_test_out_a, ec))
    {
      std::filesystem::current_path(cwd_test_dir_b, ec);
      cwd_test_changed = !ec;
    }
  }

  // the open-failure test of store(): the action of AllocFault notes whether the output exists
  std::filesystem::path open_test_file;
  bool open_test_file_seen = false;
  void noteWhetherOpenTestFileExists()
  {
    std::error_code ec;
    open_test_file_seen = std::filesystem::exists(open_test_file, ec);
  }

  // the replacement test of store(): once store() has written to its output (its size is no longer 0, so opening it is
  // over), the action of AllocFault moves the output aside and renames the sentinel into its place, as another process
  // could
  std::filesystem::path swap_test_out;
  std::filesystem::path swap_test_aside;
  std::filesystem::path swap_test_sentinel;
  bool swap_test_swapped = false;
  void swapInSentinelOnceOutputWritten()
  {
    std::error_code ec;
    const std::uintmax_t size = std::filesystem::file_size(swap_test_out, ec);
    if (swap_test_swapped || ec || size == 0) return;
    std::filesystem::rename(swap_test_out, swap_test_aside, ec);
    if (ec) return;
    std::filesystem::rename(swap_test_sentinel, swap_test_out, ec);
    swap_test_swapped = !ec;
  }

  // the replacement test of store() while it opens its output: as above, but while the output exists and is still
  // empty, i.e. also within opening it
  void swapInSentinelWhileOutputEmpty()
  {
    std::error_code ec;
    const std::uintmax_t size = std::filesystem::file_size(swap_test_out, ec);
    if (swap_test_swapped || ec || size != 0) return;
    std::filesystem::rename(swap_test_out, swap_test_aside, ec);
    if (ec) return;
    std::filesystem::rename(swap_test_sentinel, swap_test_out, ec);
    swap_test_swapped = !ec;
  }

#ifdef __GLIBCXX__
  // The reporting test of store(): a code conversion facet whose unshift() throws makes the final close() of store()
  // throw (libstdc++). Just before it throws, it arms AllocFault, so that the n-th allocation after the throw fails: an
  // allocation made to report the error.
  struct CloseTestError
  {
  };
  std::exception_ptr close_test_error; // what the facet throws (rethrown: throwing it allocates nothing); CloseTestError if null
  long close_test_alloc_n = 0;         // 0: no allocation fails
  struct ArmAndThrowOnUnshift : public std::codecvt<char, char, std::mbstate_t>
  {
  protected:
    bool do_always_noconv() const noexcept override { return false; }
    result do_out(state_type&, const char* from, const char* from_end, const char*& from_next, char* to, char* to_end,
                  char*& to_next) const override
    {
      const std::ptrdiff_t n = std::min(from_end - from, to_end - to);
      std::copy(from, from + n, to);
      from_next = from + n;
      to_next = to + n;
      return from_next == from_end ? ok : partial;
    }
    result do_unshift(state_type&, char*, char*, char*&) const override
    {
      if (close_test_alloc_n > 0) AllocFault::arm(close_test_alloc_n);
      if (close_test_error) std::rethrow_exception(close_test_error);
      throw CloseTestError();
    }
  };

  // Runs in a child process (an allocation failure may end it): store() whose close() throws the error of @p kind (0:
  // Exception::ConversionError, 1: std::bad_alloc, 2: CloseTestError, which is not a std::exception), and the n-th
  // allocation after that fails (none if n is 0). Returns 0 if store() raised that error with its type (the OpenMS
  // exception with the note, unless the failed allocation was one of the note's), 1 otherwise; plus 16 if the
  // allocation failure was injected.
  int closeTestChild(int kind, long n, const std::string& file, const std::vector<ProteinIdentification>& prots,
                     const PeptideIdentificationList& peps)
  {
    struct rlimit no_core = {0, 0};
    ::setrlimit(RLIMIT_CORE, &no_core); // a child that is ended leaves no core file
    if (kind == 0) close_test_error = std::make_exception_ptr(Exception::ConversionError(__FILE__, __LINE__, "closeTestChild", "an error while closing"));
    if (kind == 1) close_test_error = std::make_exception_ptr(std::bad_alloc());
    close_test_alloc_n = n;
    std::locale::global(std::locale(std::locale(), new ArmAndThrowOnUnshift));
    IdXMLFile writer; // constructed before store(): only allocations after the throw fail
    bool with_type = false;
    try
    {
      writer.store(file, prots, peps);
      AllocFault::disarm();
    }
    catch (const Exception::ConversionError& e)
    {
      AllocFault::disarm();
      with_type = kind == 0 && (AllocFault::fired || saysPartialFileRemains(e.what(), file));
    }
    catch (const std::bad_alloc&)
    {
      AllocFault::disarm();
      with_type = kind == 1;
    }
    catch (const CloseTestError&)
    {
      AllocFault::disarm();
      with_type = kind == 2;
    }
    catch (...)
    {
      AllocFault::disarm();
    }
    return (with_type ? 0 : 1) | (AllocFault::fired ? 16 : 0);
  }
#endif
#endif

  [[maybe_unused]] const char* const fault_injection_unavailable = "SKIPPED: fault injection (a replacement of the global operator new "
      "that also receives the allocations of libOpenMS) is available on Linux only, and not if "
      "OPENMS_TEST_NO_FAULT_INJECTION is defined";
  const char* const fault_injection_unreached = "SKIPPED: the replacement of operator new in this test program does "
      "not receive the allocations of libOpenMS (e.g. a statically linked C++ runtime)";

  // first occurrence only; returns false when @p from is absent
  bool replaceInFile4b(const std::string& path, const std::string& from, const std::string& to)
  {
    std::string s = slurp4b(path);
    const std::size_t pos = s.find(from);
    if (pos == std::string::npos) return false;
    s.replace(pos, from.size(), to);
    std::ofstream out(path);
    out << s;
    return true;
  }
}

START_TEST(IdXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((IdXMLFile()))
  IdXMLFile* ptr = nullptr;
  IdXMLFile* nullPointer = nullptr;
  ptr = new IdXMLFile();
  TEST_NOT_EQUAL(ptr,nullPointer)
  delete ptr;
END_SECTION

START_SECTION(void load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids) )
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(peptide_ids.size(),3)
END_SECTION


START_SECTION(void load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids, std::string& document_id) )
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids, document_id);

  TEST_STRING_EQUAL(document_id,"LSID1234")
  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(peptide_ids.size(),3)

  /////////////// protein id 1 //////////////////
  TEST_EQUAL(protein_ids[0].getScoreType(),"MOWSE")
  TEST_EQUAL(protein_ids[0].isHigherScoreBetter(),true)
  TEST_EQUAL(protein_ids[0].getSearchEngine(),"Mascot")
  TEST_EQUAL(protein_ids[0].getSearchEngineVersion(),"2.1.0")
  TEST_EQUAL(protein_ids[0].getDateTime().getDate(),"2006-01-12")
  TEST_EQUAL(protein_ids[0].getDateTime().getTime(),"12:13:14")
  TEST_EQUAL(StringUtils::hasPrefix(protein_ids[0].getIdentifier(), "Mascot_2006-01-12T12:13:14"), true)
  TEST_EQUAL(protein_ids[0].getSearchParameters().db,"MSDB")
  TEST_EQUAL(protein_ids[0].getSearchParameters().db_version,"1.0")
  TEST_EQUAL(protein_ids[0].getSearchParameters().charges,"1, 2")
  TEST_EQUAL(protein_ids[0].getSearchParameters().mass_type,ProteinIdentification::PeakMassType::AVERAGE)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().fragment_mass_tolerance,0.3)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_mass_tolerance,1.0)
  TEST_EQUAL(std::string(protein_ids[0].getMetaValue("name")),"ProteinIdentification")

  TEST_EQUAL(protein_ids[0].getProteinGroups().size(), 1);
  TEST_EQUAL(protein_ids[0].getProteinGroups()[0].probability, 0.88);
  TEST_EQUAL(protein_ids[0].getProteinGroups()[0].accessions.size(), 2);
  TEST_EQUAL(protein_ids[0].getProteinGroups()[0].accessions[0], "PROT1");
  TEST_EQUAL(protein_ids[0].getProteinGroups()[0].accessions[1], "PROT2");

  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins().size(), 1);
  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins()[0].accessions.size(),
             2);
  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins()[0].accessions[0],
             "PROT1");
  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins()[0].accessions[1],
             "PROT2");

  TEST_EQUAL(protein_ids[0].getHits().size(),2)
  //protein hit 1
  TEST_REAL_SIMILAR(protein_ids[0].getHits()[0].getScore(),34.4)
  TEST_EQUAL(protein_ids[0].getHits()[0].getAccession(),"PROT1")
  TEST_EQUAL(protein_ids[0].getHits()[0].getSequence(),"ABCDEFG")
  TEST_EQUAL(std::string(protein_ids[0].getHits()[0].getMetaValue("name")),"ProteinHit")
  //protein hit 2
  TEST_REAL_SIMILAR(protein_ids[0].getHits()[1].getScore(),24.4)
  TEST_EQUAL(protein_ids[0].getHits()[1].getAccession(),"PROT2")
  TEST_EQUAL(protein_ids[0].getHits()[1].getSequence(),"ABCDEFG")

  //peptide id 1
  TEST_EQUAL(peptide_ids[0].getScoreType(),"MOWSE")
  TEST_EQUAL(peptide_ids[0].isHigherScoreBetter(),false)
  TEST_EQUAL(StringUtils::hasPrefix(peptide_ids[0].getIdentifier(), "Mascot_2006-01-12T12:13:14"), true)
  TEST_REAL_SIMILAR(peptide_ids[0].getMZ(),675.9)
  TEST_REAL_SIMILAR(peptide_ids[0].getRT(),1234.5)
  TEST_EQUAL((peptide_ids[0].getSpectrumReference()),"17")
  TEST_EQUAL(std::string(peptide_ids[0].getMetaValue("name")),"PeptideIdentification")
  TEST_EQUAL(peptide_ids[0].getHits().size(),2)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),0.9)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence(), AASequence::fromString("PEPTIDER"))
  TEST_EQUAL(peptide_ids[0].getHits()[0].getCharge(),1)
  vector<PeptideEvidence> pes0 = peptide_ids[0].getHits()[0].getPeptideEvidences();
  TEST_EQUAL(pes0.size(),2)
  TEST_EQUAL(pes0[0].getProteinAccession(),"PROT1")
  TEST_EQUAL(pes0[1].getProteinAccession(),"PROT2")
  TEST_EQUAL(pes0[0].getAABefore(),'A')
  TEST_EQUAL(pes0[0].getAAAfter(),'B')
  TEST_EQUAL(std::string(peptide_ids[0].getHits()[0].getMetaValue("name")),"PeptideHit")
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[1].getScore(),1.4)
  vector<PeptideEvidence> pes1 = peptide_ids[0].getHits()[1].getPeptideEvidences();
  TEST_EQUAL(peptide_ids[0].getHits()[1].getSequence(), AASequence::fromString("PEPTIDERR"))
  TEST_EQUAL(peptide_ids[0].getHits()[1].getCharge(),1)
  TEST_EQUAL(pes1.size(),0)
  //peptide id 2
  TEST_EQUAL(peptide_ids[1].getScoreType(),"MOWSE")
  TEST_EQUAL(peptide_ids[1].isHigherScoreBetter(),true)
  TEST_EQUAL(StringUtils::hasPrefix(peptide_ids[1].getIdentifier(), "Mascot_2006-01-12T12:13:14"), true)
  TEST_EQUAL(peptide_ids[1].getHits().size(),2)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),44.4)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence(), AASequence::fromString("PEPTIDERRR"))
  TEST_EQUAL(peptide_ids[1].getHits()[0].getCharge(),2)
  vector<PeptideEvidence> pes2 = peptide_ids[1].getHits()[0].getPeptideEvidences();
  TEST_EQUAL(pes2.size(),0)
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),33.3)
  TEST_EQUAL(peptide_ids[1].getHits()[1].getSequence(), AASequence::fromString("PEPTIDERRRR"))
  TEST_EQUAL(peptide_ids[1].getHits()[1].getCharge(),2)
  vector<PeptideEvidence> pes3 = peptide_ids[1].getHits()[1].getPeptideEvidences();
  TEST_EQUAL(pes3.size(),0)

  /////////////// protein id 2 //////////////////
  TEST_EQUAL(protein_ids[1].getScoreType(),"MOWSE")
  TEST_EQUAL(protein_ids[1].isHigherScoreBetter(),true)
  TEST_EQUAL(protein_ids[1].getSearchEngine(),"Mascot")
  TEST_EQUAL(protein_ids[1].getSearchEngineVersion(),"2.1.1")
  TEST_EQUAL(protein_ids[1].getDateTime().getDate(),"2007-01-12")
  TEST_EQUAL(protein_ids[1].getDateTime().getTime(),"12:13:14")
  TEST_EQUAL(StringUtils::hasPrefix(protein_ids[1].getIdentifier(), "Mascot_2007-01-12T12:13:14"), true)
  TEST_EQUAL(protein_ids[1].getSearchParameters().db,"MSDB")
  TEST_EQUAL(protein_ids[1].getSearchParameters().db_version,"1.1")
  TEST_EQUAL(protein_ids[1].getSearchParameters().charges,"1, 2, 3")
  TEST_EQUAL(protein_ids[1].getSearchParameters().mass_type,ProteinIdentification::PeakMassType::MONOISOTOPIC)
  TEST_REAL_SIMILAR(protein_ids[1].getSearchParameters().fragment_mass_tolerance,0.3)
  TEST_REAL_SIMILAR(protein_ids[1].getSearchParameters().precursor_mass_tolerance,1.0)
  TEST_EQUAL(protein_ids[1].getSearchParameters().fixed_modifications.size(),2)
  TEST_EQUAL(protein_ids[1].getSearchParameters().fixed_modifications[0],"Fixed")
  TEST_EQUAL(protein_ids[1].getSearchParameters().fixed_modifications[1],"Fixed2")
  TEST_EQUAL(protein_ids[1].getSearchParameters().variable_modifications.size(),2)
  TEST_EQUAL(protein_ids[1].getSearchParameters().variable_modifications[0],"Variable")
  TEST_EQUAL(protein_ids[1].getSearchParameters().variable_modifications[1],"Variable2")
  TEST_EQUAL(protein_ids[1].getHits().size(),1)
  //protein hit 1
  TEST_REAL_SIMILAR(protein_ids[1].getHits()[0].getScore(),100.0)
  TEST_EQUAL(protein_ids[1].getHits()[0].getAccession(),"PROT3")
  TEST_EQUAL(protein_ids[1].getHits()[0].getSequence(),"")
  //peptide id 3
  TEST_EQUAL(peptide_ids[2].getScoreType(),"MOWSE")
  TEST_EQUAL(peptide_ids[2].isHigherScoreBetter(),true)
  TEST_EQUAL(StringUtils::hasPrefix(peptide_ids[2].getIdentifier(), "Mascot_2007-01-12T12:13:14"), true)
  TEST_EQUAL(peptide_ids[2].getHits().size(),1)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[2].getHits()[0].getScore(),1.4)
  TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence(), AASequence::fromString("PEPTIDERRRRR"))
  TEST_EQUAL(peptide_ids[2].getHits()[0].getCharge(),1)
  vector<PeptideEvidence> pes4 = peptide_ids[2].getHits()[0].getPeptideEvidences();
  TEST_EQUAL(pes4.size(),1)
  TEST_EQUAL(pes4[0].getProteinAccession(),"PROT3")
  TEST_EQUAL(pes4[0].getAABefore(), PeptideEvidence::UNKNOWN_AA)
  TEST_EQUAL(pes4[0].getAAAfter(), PeptideEvidence::UNKNOWN_AA)
END_SECTION

START_SECTION(void store(std::string filename, const std::vector<ProteinIdentification>& protein_ids, const PeptideIdentificationList& peptide_ids, const std::string& document_id="") )

  // load, store, and reload data
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  PeptideIdentificationList peptide_ids, peptide_ids2;
  std::string document_id, document_id2;
  std::string target_file = OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML");
  IdXMLFile().load(target_file, protein_ids2, peptide_ids2, document_id2);

  std::string actual_file;
  NEW_TMP_FILE(actual_file)
  IdXMLFile().store(actual_file, protein_ids2, peptide_ids2, document_id2);

  FuzzyStringComparator fuzzy;
  fuzzy.setWhitelist(ListUtils::create<std::string>("<?xml-stylesheet"));
  fuzzy.setAcceptableAbsolute(0.0001);
  bool result = fuzzy.compareFiles(actual_file, target_file);
  TEST_EQUAL(result, true);
END_SECTION


START_SECTION([EXTRA] static bool isValid(const std::string& filename))
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  PeptideIdentificationList peptide_ids, peptide_ids2;
  std::string filename;
  IdXMLFile f;

  //test if empty file is valid
  NEW_TMP_FILE(filename)
  f.store(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(f.isValid(filename, std::cerr),true);

  //test if full file is valid
  NEW_TMP_FILE(filename);
  std::string document_id;
  f.load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids2, peptide_ids2, document_id);
  protein_ids2[0].setMetaValue("stringvalue",std::string("bla"));
  protein_ids2[0].setMetaValue("intvalue",4711);
  protein_ids2[0].setMetaValue("floatvalue",5.3);
  f.store(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(f.isValid(filename, std::cerr),true);

  //check if meta information can be loaded
  f.load(filename, protein_ids2, peptide_ids2, document_id);
END_SECTION

START_SECTION(([EXTRA] No protein identification bug))
  IdXMLFile id_xmlfile;
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  id_xmlfile.load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_no_proteinhits.idXML"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(), 1)
  TEST_EQUAL(protein_ids[0].getHits().size(), 0)
  TEST_EQUAL(peptide_ids.size(), 10)
  TEST_EQUAL(peptide_ids[0].getHits().size(), 1)

  std::string filename;
  NEW_TMP_FILE(filename)
  id_xmlfile.store(filename , protein_ids, peptide_ids);

  vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;
  id_xmlfile.load(filename, protein_ids2, peptide_ids2);

  // identifiers contain a random number when loaded to avoid ambiguities when merging ProtIDs; make them equal for our purposes here
  protein_ids2[0].setIdentifier(protein_ids[0].getIdentifier());
  for (auto& pep : peptide_ids2)
  {
    pep.setIdentifier(peptide_ids[0].getIdentifier());
  }

  TEST_TRUE(protein_ids == protein_ids2)
  TEST_TRUE(peptide_ids == peptide_ids2)

END_SECTION

START_SECTION(([EXTRA] XLMS data labeled cross-linker))
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  std::string input_file= OPENMS_GET_TEST_DATA_PATH("IdXML_XLMS_labelled.idXML");
  IdXMLFile().load(input_file, protein_ids, peptide_ids);

  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[1].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[8].annotation, "[alpha|xi$b8]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[20].annotation, "[alpha|xi$b9]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[25].charge, 3)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[25].annotation, "[alpha|xi$y8]")

END_SECTION

START_SECTION([EXTRA] Compressed file writing - gzip round-trip)
  // Load reference data
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids);

  // Store as gzip-compressed file
  std::string tmp_gz;
  NEW_TMP_FILE(tmp_gz);
  tmp_gz += ".gz";
  IdXMLFile().store(tmp_gz, protein_ids, peptide_ids);

  // Load back from compressed file
  std::vector<ProteinIdentification> protein_ids_gz;
  PeptideIdentificationList peptide_ids_gz;
  IdXMLFile().load(tmp_gz, protein_ids_gz, peptide_ids_gz);

  // Verify round-trip integrity
  TEST_EQUAL(protein_ids_gz.size(), protein_ids.size())
  TEST_EQUAL(peptide_ids_gz.size(), peptide_ids.size())
  TEST_EQUAL(protein_ids_gz[0].getHits().size(), protein_ids[0].getHits().size())
  TEST_EQUAL(protein_ids_gz[0].getHits()[0].getAccession(), protein_ids[0].getHits()[0].getAccession())
END_SECTION

START_SECTION([EXTRA] Compressed file writing - bzip2 round-trip)
  // Load reference data
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids);

  // Store as bzip2-compressed file
  std::string tmp_bz2;
  NEW_TMP_FILE(tmp_bz2);
  tmp_bz2 += ".bz2";
  IdXMLFile().store(tmp_bz2, protein_ids, peptide_ids);

  // Load back from compressed file
  std::vector<ProteinIdentification> protein_ids_bz2;
  PeptideIdentificationList peptide_ids_bz2;
  IdXMLFile().load(tmp_bz2, protein_ids_bz2, peptide_ids_bz2);

  // Verify round-trip integrity
  TEST_EQUAL(protein_ids_bz2.size(), protein_ids.size())
  TEST_EQUAL(peptide_ids_bz2.size(), peptide_ids.size())
  TEST_EQUAL(protein_ids_bz2[0].getHits().size(), protein_ids[0].getHits().size())
  TEST_EQUAL(protein_ids_bz2[0].getHits()[0].getAccession(), protein_ids[0].getHits()[0].getAccession())
END_SECTION

START_SECTION([EXTRA] store/load - a tool-defined modification travels with its definition)
{
  TEST_TRUE(defineMod4b("TestIdXML:Adduct", 'K', "C9H11N2O8P") != nullptr)
  ProteinIdentification prot;
  prot.setIdentifier("run4b");
  prot.setDateTime(DateTime::now());
  PeptideIdentification pep;
  pep.setIdentifier("run4b");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("AEADNLDDK(TestIdXML:Adduct)K"));
  pep.insertHit(hit);
  std::vector<ProteinIdentification> prots(1, prot);
  PeptideIdentificationList peps;
  peps.push_back(pep);

  std::string filename;
  NEW_TMP_FILE(filename)
  IdXMLFile().store(filename, prots, peps);
  TEST_TRUE(fileContains4b(filename, "name=\"modification_definitions\""))
  TEST_TRUE(fileContains4b(filename, "1|TestIdXML:Adduct|TestIdXML:Adduct (K)|"))

  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(filename, prots_in, peps_in);
  TEST_EQUAL(prots_in.size(), 1)
  TEST_EQUAL(peps_in.size(), 1)
  if (prots_in.size() == 1 && peps_in.size() == 1 && !peps_in[0].getHits().empty())
  {
    TEST_TRUE(prots_in[0].getSearchParameters().metaValueExists(Constants::UserParam::MODIFICATION_DEFINITIONS))
    const AASequence& back = peps_in[0].getHits()[0].getSequence();
    TEST_EQUAL(back.toString(), "AEADNLDDK(TestIdXML:Adduct)K")
    TEST_EQUAL(back.getFormula(Residue::Full, 0).toString(), "C54H86N15O28P1")
  }
}
END_SECTION

START_SECTION([EXTRA] store - runs with equal parameters share one block that carries both definition sets)
{
  TEST_TRUE(defineMod4b("TestIdXML:RunA", 'K', "C2H2O") != nullptr)
  TEST_TRUE(defineMod4b("TestIdXML:RunB", 'R', "CH2") != nullptr)
  std::vector<ProteinIdentification> prots(2);
  prots[0].setIdentifier("runA");
  prots[0].setDateTime(DateTime::now());
  prots[1].setIdentifier("runB");
  prots[1].setDateTime(DateTime::now());
  PeptideIdentificationList peps(2);
  peps[0].setIdentifier("runA");
  peps[1].setIdentifier("runB");
  PeptideHit ha, hb;
  ha.setSequence(AASequence::fromString("PEPK(TestIdXML:RunA)IDE"));
  hb.setSequence(AASequence::fromString("PEPR(TestIdXML:RunB)IDE"));
  peps[0].insertHit(ha);
  peps[1].insertHit(hb);

  std::string filename;
  NEW_TMP_FILE(filename)
  IdXMLFile().store(filename, prots, peps);
  const std::string text = slurp4b(filename);
  Size blocks = 0;
  for (std::size_t pos = text.find("<SearchParameters "); pos != std::string::npos; pos = text.find("<SearchParameters ", pos + 1)) ++blocks;
  TEST_EQUAL(blocks, 1)
  TEST_TRUE(text.find("1|TestIdXML:RunA|TestIdXML:RunA (K)|") != std::string::npos)
  TEST_TRUE(text.find("1|TestIdXML:RunB|TestIdXML:RunB (R)|") != std::string::npos)

  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(filename, prots_in, peps_in);
  TEST_EQUAL(prots_in.size(), 2)
  TEST_EQUAL(peps_in.size(), 2)
  if (peps_in.size() == 2 && !peps_in[0].getHits().empty() && !peps_in[1].getHits().empty())
  {
    TEST_EQUAL(peps_in[0].getHits()[0].getSequence().toString(), "PEPK(TestIdXML:RunA)IDE")
    TEST_EQUAL(peps_in[1].getHits()[0].getSequence().toString(), "PEPR(TestIdXML:RunB)IDE")
  }
}
END_SECTION

START_SECTION([EXTRA] store - a defined variable modification on no hit keeps its definition)
{
  TEST_TRUE(defineMod4b("TestIdXML:VarOnly", 'S', "HPO3") != nullptr)
  ProteinIdentification prot;
  prot.setIdentifier("run4b_var");
  prot.setDateTime(DateTime::now());
  ProteinIdentification::SearchParameters sp;
  sp.variable_modifications.push_back("TestIdXML:VarOnly (S)");
  prot.setSearchParameters(sp);
  PeptideIdentification pep;
  pep.setIdentifier("run4b_var");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  pep.insertHit(hit);
  std::vector<ProteinIdentification> prots(1, prot);
  PeptideIdentificationList peps;
  peps.push_back(pep);

  std::string filename;
  NEW_TMP_FILE(filename)
  IdXMLFile().store(filename, prots, peps);
  TEST_TRUE(fileContains4b(filename, "1|TestIdXML:VarOnly|TestIdXML:VarOnly (S)|"))
}
END_SECTION

START_SECTION([EXTRA] load - definitions are registered before the sequences are parsed)
{
  const ModificationsDB* db = ModificationsDB::getInstance();
  TEST_FALSE(db->hasDefinedModification("TestIdXML:Fresh"))
  ProteinIdentification prot;
  prot.setIdentifier("run4b_fresh");
  prot.setDateTime(DateTime::now());
  ProteinIdentification::SearchParameters sp;
  sp.setMetaValue(Constants::UserParam::MODIFICATION_DEFINITIONS, freshRecord4b("TestIdXML:Fresh", 'K', "C2H2O"));
  prot.setSearchParameters(sp);
  PeptideIdentification pep;
  pep.setIdentifier("run4b_fresh");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTKIDE"));
  pep.insertHit(hit);
  std::vector<ProteinIdentification> prots(1, prot);
  PeptideIdentificationList peps;
  peps.push_back(pep);

  std::string filename;
  NEW_TMP_FILE(filename)
  IdXMLFile().store(filename, prots, peps);
  TEST_FALSE(db->hasDefinedModification("TestIdXML:Fresh")) // storing registers nothing
  // a hit using the not-yet-registered name, as a file from another process would carry it
  TEST_TRUE(replaceInFile4b(filename, "sequence=\"PEPTKIDE\"", "sequence=\"PEPTK(TestIdXML:Fresh)IDE\""))

  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(filename, prots_in, peps_in);
  TEST_TRUE(db->hasDefinedModification("TestIdXML:Fresh"))
  TEST_EQUAL(peps_in.size(), 1)
  if (peps_in.size() == 1 && !peps_in[0].getHits().empty())
  {
    const AASequence& back = peps_in[0].getHits()[0].getSequence();
    TEST_EQUAL(back.toString(), "PEPTK(TestIdXML:Fresh)IDE")
    TEST_REAL_SIMILAR(back.getMonoWeight(), AASequence::fromString("PEPTKIDE").getMonoWeight() + EmpiricalFormula("C2H2O").getMonoWeight())
  }
}
END_SECTION

START_SECTION([EXTRA] store - many peptide identifications are written in input order with any number of threads and the first error is reported)
{
  // more peptide identifications than one block of the parallel writer, in two runs and some of no run; hits not in
  // score order
  std::vector<ProteinIdentification> prots(2);
  prots[0].setIdentifier("runPar1");
  prots[1].setIdentifier("runPar2");
  prots[0].insertHit(ProteinHit(0.0, 1, "ACC1", ""));
  prots[1].insertHit(ProteinHit(0.0, 1, "ACC2", ""));
  for (ProteinIdentification& prot : prots) prot.setDateTime(DateTime::now());

  const Size n = 400;
  const std::vector<std::string> sequences = {"PEPTIDEA", "PEPTIDEB", "PEPTIDEC"};
  const std::vector<double> scores = {1.0, 3.0, 3.0}; // B and C tie: they keep their order
  PeptideIdentificationList peps(n);
  for (Size l = 0; l < n; ++l)
  {
    PeptideIdentification& pep = peps[l];
    pep.setIdentifier(l % 50 == 13 ? "runNone" : l % 2 == 0 ? "runPar1" : "runPar2"); // no run: not written
    pep.setScoreType("score");
    pep.setHigherScoreBetter(true);
    pep.setRT(double(l));
    pep.setMZ(500.0 + double(l));
    pep.setSpectrumReference("scan=" + StringUtils::toStr(l));
    pep.setMetaValue("par_test_index", int(l));
    if (l % 50 == 7) continue; // no hits: not written
    for (Size h = 0; h < sequences.size(); ++h)
    {
      PeptideHit hit(scores[h], 0, 2, AASequence::fromString(sequences[h]));
      hit.addPeptideEvidence(PeptideEvidence(l % 2 == 0 ? "ACC1" : "ACC2", 8 * int(h), 8 * int(h) + 7, '-', 'P'));
      hit.setMetaValue("par_test_int", int(3 * l + h));
      hit.setMetaValue("par_test_double", double(l) + 0.25 * double(h));
      hit.setMetaValue("par_test_string", std::string("a<b&\"c\""));
      hit.setMetaValue("par_test_strings", StringList{"x,y", "z"});
      hit.setMetaValue("par_test_ints", IntList{1, int(l)});
      hit.setMetaValue("par_test_doubles", DoubleList{0.5, double(l)});
      PeptideHit::PeakAnnotation annotation;
      annotation.annotation = "y1";
      annotation.charge = 1;
      annotation.mz = 100.0 + double(l);
      annotation.intensity = 1.0;
      hit.setPeakAnnotations({annotation});
      pep.insertHit(hit);
    }
  }

  // the same file with one and with several threads
  std::string file_1, file_n;
  NEW_TMP_FILE(file_1)
  NEW_TMP_FILE(file_n)
#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  omp_set_num_threads(1);
#endif
  IdXMLFile().store(file_1, prots, peps);
#ifdef _OPENMP
  omp_set_num_threads(std::max(max_threads, 4));
#endif
  IdXMLFile().store(file_n, prots, peps);
  TEST_EQUAL(slurp4b(file_1) == slurp4b(file_n), true)

  // run by run, in input order
  std::vector<Size> expected;
  for (Size l = 0; l < n; l += 2) expected.push_back(l);
  for (Size l = 1; l < n; l += 2) if (l % 50 != 7 && l % 50 != 13) expected.push_back(l);
  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(file_n, prots_in, peps_in);
  TEST_EQUAL(prots_in.size(), 2)
  TEST_EQUAL(peps_in.size(), expected.size())
  ABORT_IF(peps_in.size() != expected.size())
  bool in_order = true;
  for (Size k = 0; k < expected.size(); ++k)
  {
    const PeptideIdentification& pep = peps_in[k];
    const Size l = expected[k];
    in_order &= pep.getIdentifier() == prots_in[l % 2].getIdentifier(); // idXML does not keep the identifiers
    in_order &= pep.getRT() == double(l) && pep.getSpectrumReference() == "scan=" + StringUtils::toStr(l);
    in_order &= int(pep.getMetaValue("par_test_index")) == int(l);
    in_order &= pep.getHits().size() == 3 && pep.getHits()[0].getSequence().toString() == "PEPTIDEB" &&
                pep.getHits()[1].getSequence().toString() == "PEPTIDEC" && pep.getHits()[2].getSequence().toString() == "PEPTIDEA";
  }
  TEST_EQUAL(in_order, true)
  for (const Size k : {Size(0), Size(17), expected.size() - 1})
  {
    const Size l = expected[k];
    const PeptideHit& hit = peps_in[k].getHits()[0];
    TEST_EQUAL(int(hit.getMetaValue("par_test_int")), int(3 * l + 1))
    TEST_EQUAL(double(hit.getMetaValue("par_test_double")), double(l) + 0.25)
    TEST_EQUAL(hit.getMetaValue("par_test_string").toString(), "a<b&\"c\"")
    TEST_EQUAL(ListUtils::concatenate(hit.getMetaValue("par_test_strings").toStringList(), ";"), "x,y;z")
    TEST_EQUAL(ListUtils::concatenate(hit.getMetaValue("par_test_ints").toIntList(), ";"), "1;" + StringUtils::toStr(l))
    TEST_EQUAL(ListUtils::concatenate(hit.getMetaValue("par_test_doubles").toDoubleList(), ";"), "0.5;" + StringUtils::toStr(double(l)))
    TEST_EQUAL(hit.getPeakAnnotations().size(), 1)
    TEST_EQUAL(hit.getPeakAnnotations()[0].mz, 100.0 + double(l))
    TEST_EQUAL(hit.getPeptideEvidences().size(), 1)
    TEST_EQUAL(hit.getPeptideEvidences()[0].getProteinAccession(), l % 2 == 0 ? "ACC1" : "ACC2")
  }

  // unknown accessions in two peptide identifications of the first run: the one first in input order is reported
  PeptideIdentificationList bad = peps;
  for (PeptideHit& hit : bad[36].getHits()) hit.setPeptideEvidences({PeptideEvidence("UNKNOWN_A", 0, 7, '-', 'P')});
  for (PeptideHit& hit : bad[300].getHits()) hit.setPeptideEvidences({PeptideEvidence("UNKNOWN_B", 0, 7, '-', 'P')});
  std::string file_bad;
  NEW_TMP_FILE(file_bad)
  for (int repeat = 0; repeat < 5; ++repeat)
  {
    std::string message;
    try
    {
      IdXMLFile().store(file_bad, prots, bad);
    }
    catch (const Exception::ElementNotFound& e)
    {
      message = e.what();
    }
    TEST_EQUAL(message.find("No accession UNKNOWN_A found in run 'runPar1'") != std::string::npos, true)
  }
  // all of them unknown: every block fails, the first one is reported
  for (Size l = 0; l < n; ++l)
  {
    for (PeptideHit& hit : bad[l].getHits()) hit.setPeptideEvidences({PeptideEvidence("UNKNOWN_" + StringUtils::toStr(l), 0, 7, '-', 'P')});
  }
  for (int repeat = 0; repeat < 5; ++repeat)
  {
    std::string message;
    try
    {
      IdXMLFile().store(file_bad, prots, bad);
    }
    catch (const Exception::ElementNotFound& e)
    {
      message = e.what();
    }
    TEST_EQUAL(message.find("No accession UNKNOWN_0 found in run 'runPar1'") != std::string::npos, true)
  }
  // an EMPTY meta value cannot be written (see XMLHandler::writeUserParamValue_()); on every hit, all threads fail at
  // once and the ConversionError arrives intact
  PeptideIdentificationList empty_values = peps;
  for (PeptideIdentification& pep : empty_values)
  {
    for (PeptideHit& hit : pep.getHits()) hit.setMetaValue("par_test_empty", DataValue());
  }
  for (int repeat = 0; repeat < 5; ++repeat)
  {
    TEST_EXCEPTION(Exception::ConversionError, IdXMLFile().store(file_bad, prots, empty_values))
  }
  std::remove(file_bad.c_str()); // incomplete, not to be validated
#ifdef _OPENMP
  omp_set_num_threads(max_threads);
#endif
}
END_SECTION

START_SECTION([EXTRA] store - hits with NaN scores are written in the order of PeptideIdentification::sort() with any number of threads)
{
  // NaN scores are not a strict weak order; the hits must still be written in the order that sorting them gives
  std::vector<ProteinIdentification> prots(3);
  for (Size r = 0; r < prots.size(); ++r)
  {
    prots[r].setIdentifier("runNaN" + StringUtils::toStr(r));
    prots[r].setDateTime(DateTime::now());
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<std::vector<double>> patterns = {{2.0, 1.0, nan, 5.0, nan, nan}, {nan, 1.0}, {1.0, nan, 2.0},
                                                     {nan, nan, 3.0, 1.0, 2.0, 4.0, nan, 0.5}, {3.0, 1.0, 2.0}};
  const std::string residues = "ACDEFGHI";
  const Size n = 120;
  PeptideIdentificationList peps(n);
  for (Size l = 0; l < n; ++l)
  {
    PeptideIdentification& pep = peps[l];
    pep.setIdentifier("runNaN" + StringUtils::toStr(l % prots.size()));
    pep.setScoreType("score");
    pep.setHigherScoreBetter(l % 7 < 4);
    const std::vector<double>& scores = patterns[l % patterns.size()];
    for (Size h = 0; h < scores.size(); ++h)
    {
      pep.insertHit(PeptideHit(scores[h], 0, 2, AASequence::fromString(std::string("SAMPLE") + residues[h])));
    }
  }

  std::string file_1, file_n;
  NEW_TMP_FILE(file_1)
  NEW_TMP_FILE(file_n)
#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  omp_set_num_threads(1);
#endif
  IdXMLFile().store(file_1, prots, peps);
#ifdef _OPENMP
  omp_set_num_threads(std::max(max_threads, 4));
#endif
  IdXMLFile().store(file_n, prots, peps);
#ifdef _OPENMP
  omp_set_num_threads(max_threads);
#endif
  TEST_EQUAL(slurp4b(file_1) == slurp4b(file_n), true)

  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(file_n, prots_in, peps_in);
  TEST_EQUAL(peps_in.size(), n)
  ABORT_IF(peps_in.size() != n)
  bool same_order = true;
  Size k = 0;
  for (Size r = 0; r < prots.size(); ++r)
  {
    for (Size l = r; l < n; l += prots.size(), ++k)
    {
      PeptideIdentification sorted = peps[l];
      sorted.sort();
      const std::vector<PeptideHit>& hits_in = peps_in[k].getHits();
      same_order &= hits_in.size() == sorted.getHits().size();
      for (Size h = 0; same_order && h < hits_in.size(); ++h)
      {
        same_order &= hits_in[h].getSequence() == sorted.getHits()[h].getSequence();
        same_order &= std::isnan(hits_in[h].getScore()) == std::isnan(sorted.getHits()[h].getScore());
      }
    }
  }
  TEST_EQUAL(same_order, true)
}
END_SECTION

START_SECTION([EXTRA] store - a failing store() throws and says so; the partial file is left in place)
{
  // The block formatter of the parallel writer: a failure of the block's stream buffer propagates. Without exceptions on
  // the stream, the std::bad_alloc would only set badbit and the partial block would be written as if complete.
  const auto write_block = [](std::ostream& out) { out << "\t\t<PeptideIdentification score_type=\"score\" >\n" << 1.25 << '\n'; };
  {
    std::ostringstream formatted, direct;
    IdXMLFileBlockWriter::formatBlock_(formatted, write_block);
    write_block(direct);
    TEST_EQUAL(formatted.str(), direct.str())
  }
  {
    BadAllocStreamBuf buffer;
    std::ostream out(&buffer);
    TEST_EXCEPTION(std::bad_alloc, IdXMLFileBlockWriter::formatBlock_(out, write_block))
  }
  {
    FailingStreamBuf buffer;
    std::ostream out(&buffer);
    TEST_EXCEPTION(std::ios_base::failure, IdXMLFileBlockWriter::formatBlock_(out, write_block))
  }

  std::vector<ProteinIdentification> prots(1);
  prots[0].setIdentifier("runFail");
  prots[0].setDateTime(DateTime::now());
  prots[0].insertHit(ProteinHit(0.0, 1, "ACC_FAIL", ""));
  const auto make_peps = [](Size n)
  {
    PeptideIdentificationList peps(n);
    for (Size l = 0; l < n; ++l)
    {
      peps[l].setIdentifier("runFail");
      peps[l].setScoreType("score");
      peps[l].setHigherScoreBetter(true);
      peps[l].setRT(double(l));
      peps[l].setMZ(500.0 + double(l));
      peps[l].setSpectrumReference("scan=" + StringUtils::toStr(l));
      for (const std::string sequence : {"PEPTIDEA", "PEPTIDEB"})
      {
        PeptideHit hit(double(l), 0, 2, AASequence::fromString(sequence));
        hit.addPeptideEvidence(PeptideEvidence("ACC_FAIL", 0, 7, '-', 'P'));
        hit.setMetaValue("fail_test_string", std::string(64, 'x'));
        peps[l].insertHit(hit);
      }
    }
    return peps;
  };

  // A block that throws: the error keeps its type, its message names the file and says that a partial file may remain,
  // and the partial file is left in place (store() does not remove its output).
  {
    PeptideIdentificationList bad = make_peps(100);
    for (PeptideHit& hit : bad[40].getHits()) hit.setMetaValue("fail_test_empty", DataValue()); // cannot be written
    std::string file_bad;
    NEW_TMP_FILE(file_bad)
    const std::string message = messageOf<Exception::ConversionError>([&] { IdXMLFile().store(file_bad, prots, bad); });
    TEST_TRUE(saysPartialFileRemains(message, file_bad))
    TEST_TRUE(isPartialIdXML(file_bad))
    std::remove(file_bad.c_str()); // partial, not to be validated
  }

  // An exception from any other part of the file, here the protein section (a meta value without a value cannot be
  // written), as well: for a file that store() created, and for a regular file that existed before (its content is
  // replaced when store() opens it), the partial file is left in place and the message says so.
  std::vector<ProteinIdentification> bad_prots = prots;
  {
    std::vector<ProteinHit> hits = bad_prots[0].getHits();
    hits[0].setMetaValue("fail_test_empty", DataValue());
    bad_prots[0].setHits(hits);
  }
  for (const bool existed : {false, true})
  {
    std::string file_bad;
    NEW_TMP_FILE(file_bad)
    std::remove(file_bad.c_str());
    if (existed) { std::ofstream(file_bad) << "previous content\n"; }
    const std::string message = messageOf<Exception::ConversionError>([&] { IdXMLFile().store(file_bad, bad_prots, make_peps(10)); });
    TEST_TRUE(saysPartialFileRemains(message, file_bad))
    TEST_TRUE(isPartialIdXML(file_bad))
    std::remove(file_bad.c_str()); // partial, not to be validated
  }

#ifdef __linux__
  // Names that store() did not create and that are not a regular file of their own are left in place when store()
  // fails, as every output is (regression guard: an earlier version removed some outputs of a failed store()).
  // The names are unique (File::getUniqueName), so concurrent runs of this test do not interfere, and they are not
  // NEW_TMP_FILEs: the end-of-test validation would read them.
  PeptideIdentificationList bad_peps = make_peps(100);
  for (PeptideHit& hit : bad_peps[40].getHits()) hit.setMetaValue("fail_test_empty", DataValue());
  for (const bool fail_in_block : {true, false})
  {
    // a file with a second hard link: neither name is removed (failure in a block of peptide identifications, or in
    // the protein section)
    const std::string original = "IdXMLFile_test_" + File::getUniqueName(false) + "_original.idXML";
    const std::string link = "IdXMLFile_test_" + File::getUniqueName(false) + "_hardlink.idXML";
    { std::ofstream(original) << "previous content\n"; }
    std::error_code ec;
    std::filesystem::create_hard_link(original, link, ec);
    ABORT_IF(bool(ec))
    if (fail_in_block) { TEST_EXCEPTION(Exception::ConversionError, IdXMLFile().store(link, prots, bad_peps)) }
    else { TEST_EXCEPTION(Exception::ConversionError, IdXMLFile().store(link, bad_prots, make_peps(10))) }
    TEST_EQUAL(File::exists(link), true)
    TEST_EQUAL(File::exists(original), true)
    std::remove(link.c_str());
    std::remove(original.c_str());
  }

  // Disk full: an idXML symlinked to /dev/full. One identification fits into the buffer of the file, so the error shows
  // only when the stream is flushed on close; 400 identifications overflow it while the blocks are written. Either way,
  // store() throws and says that a partial file may remain; the link stays, /dev/full itself is untouched.
  struct stat device;
  if (::stat("/dev/full", &device) == 0 && S_ISCHR(device.st_mode))
  {
#ifdef _OPENMP
    const int max_threads = omp_get_max_threads();
#endif
    for (const Size n : {Size(1), Size(400)})
    {
      for ([[maybe_unused]] const int threads : {1, 4})
      {
#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif
        const PeptideIdentificationList peps = make_peps(n);
        const std::string link = "IdXMLFile_test_" + File::getUniqueName(false) + "_dev_full.idXML";
        ABORT_IF(::symlink("/dev/full", link.c_str()) != 0)
        const std::string message = messageOf<Exception::UnableToCreateFile>([&] { IdXMLFile().store(link, prots, peps); });
        TEST_TRUE(saysPartialFileRemains(message, link))
        struct stat link_stat;
        TEST_EQUAL(::lstat(link.c_str(), &link_stat) == 0 && S_ISLNK(link_stat.st_mode), true)
        std::remove(link.c_str());
      }
    }
#ifdef _OPENMP
    omp_set_num_threads(max_threads);
#endif
    TEST_EQUAL(::stat("/dev/full", &device) == 0 && S_ISCHR(device.st_mode), true)
  }
#else
  STATUS("SKIPPED (Linux only): a file with a second hard link, and an idXML symlinked to /dev/full")
#endif

  // the same identifications into a regular file: written completely
  std::string file_ok;
  NEW_TMP_FILE(file_ok)
  IdXMLFile().store(file_ok, prots, make_peps(400));
  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(file_ok, prots_in, peps_in);
  TEST_EQUAL(peps_in.size(), 400)
}
END_SECTION

START_SECTION([EXTRA] store - an allocation failure anywhere in store() is raised and leaves the output in place)
{
#if FAULT_INJECTION_TESTS
  if (!allocFaultReachesLibOpenMS())
  {
    STATUS(fault_injection_unreached)
  }
  else
  {
    // Every allocation of the calling thread in store() fails once, one after the other. store() must raise the failure
    // (the std::bad_alloc unchanged; an OpenMS exception with a message that names the file and says that a partial file
    // may remain), and must not remove its output: once the output exists when the allocation fails (it existed before,
    // or opening created it), it is left in place. A store() that completes (the failed allocation was handled, e.g. a
    // nothrow allocation) writes the complete file.
    std::vector<ProteinIdentification> prots(1);
    prots[0].setIdentifier("runAlloc");
    prots[0].setDateTime(DateTime::now());
    prots[0].insertHit(ProteinHit(0.0, 1, "ACC_ALLOC", ""));
    PeptideIdentificationList peps(3); // one block, formatted by the calling thread
    for (Size l = 0; l < peps.size(); ++l)
    {
      peps[l].setIdentifier("runAlloc");
      peps[l].setScoreType("score");
      peps[l].setHigherScoreBetter(true);
      peps[l].setRT(double(l));
      peps[l].setMZ(500.0 + double(l));
      PeptideHit hit(double(l), 0, 2, AASequence::fromString("PEPTIDER"));
      hit.addPeptideEvidence(PeptideEvidence("ACC_ALLOC", 0, 7, '-', '-'));
      hit.setMetaValue("alloc_test", std::string(64, 'x'));
      peps[l].insertHit(hit);
    }
    // longer than any short-string buffer, so that copying the name allocates. Not a NEW_TMP_FILE: the sweep leaves
    // partial files, and VALIDATE_TMP_FILES at the end would check it.
    const std::string file = "IdXMLFile_test_" + File::getUniqueName(false) + "_" + std::string(100, 'n') + ".idXML";
    IdXMLFile().store(file, prots, peps);
    const std::string complete = slurp4b(file);
    const std::string previous = "previous content\n";
    open_test_file = file;
    for (const bool existed : {false, true})
    {
      long injected = 0, failed_stores = 0, raised_bad_alloc = 0, left_in_place = 0, wrong = 0, first_wrong = 0;
      for (long n = 1; n < 10000000; ++n)
      {
        std::remove(file.c_str());
        if (existed) { std::ofstream(file) << previous; }
        open_test_file_seen = false;
        bool threw = false, bad_alloc = false, openms_error = false;
        std::string message;
        IdXMLFile writer; // constructed before arming: only the allocations of store() fail
        AllocFault::arm(n, &noteWhetherOpenTestFileExists, true); // notes whether the output exists, then fails
        try
        {
          writer.store(file, prots, peps);
        }
        catch (const Exception::BaseException& e)
        {
          threw = openms_error = true;
          message = e.what();
        }
        catch (const std::bad_alloc&)
        {
          threw = bad_alloc = true;
        }
        catch (...)
        {
          threw = true;
        }
        const bool fired = AllocFault::fired;
        AllocFault::disarm();
        const bool exists = File::exists(file);
        const bool ok = threw ? (bad_alloc || (openms_error && saysPartialFileRemains(message, file)))
                                  && (exists || !open_test_file_seen)
                              : exists && slurp4b(file) == complete;
        if (!ok && wrong++ == 0) first_wrong = n;
        if (!fired) break; // n is past the last allocation of store()
        ++injected;
        if (threw) ++failed_stores;
        if (bad_alloc) ++raised_bad_alloc;
        if (threw && open_test_file_seen && exists) ++left_in_place;
      }
      STATUS("file existed before: " << existed << "; " << injected << " allocations failed one at a time, " << failed_stores
             << " stores threw (" << raised_bad_alloc << " the std::bad_alloc), " << left_in_place
             << " left the output in place, " << wrong
             << " raised another error, an OpenMS error without the note, or removed the output (first at allocation "
             << first_wrong << ")")
      TEST_TRUE(failed_stores > 0)
      TEST_TRUE(left_in_place > 0)
      TEST_EQUAL(wrong, 0)
    }

#ifdef __GLIBCXX__
    // Opening can create or truncate the file and then throw: libstdc++ allocates the stream buffer after it opened the
    // file. That allocation is the first one at which a file that did not exist before exists. If it fails, store()
    // raises the std::bad_alloc unchanged and leaves the empty file in place, whether it created it or the file was
    // already empty before.
    long n_open = 0;
    for (long n = 1; n < 10000000 && n_open == 0; ++n)
    {
      std::remove(file.c_str());
      open_test_file_seen = false;
      AllocFault::arm(n, &noteWhetherOpenTestFileExists);
      try
      {
        IdXMLFile().store(file, prots, peps);
      }
      catch (...)
      {
      }
      const bool fired = AllocFault::fired;
      AllocFault::disarm();
      if (!fired) break;
      if (open_test_file_seen) n_open = n;
    }
    STATUS("first allocation at which the output exists: " << n_open)
    TEST_TRUE(n_open > 0)
    for (const bool existed_empty : {false, true})
    {
      std::remove(file.c_str());
      if (existed_empty) { std::ofstream create(file); }
      bool bad_alloc = false;
      AllocFault::arm(n_open);
      try
      {
        IdXMLFile().store(file, prots, peps);
      }
      catch (const std::bad_alloc&)
      {
        bad_alloc = true;
      }
      catch (...)
      {
      }
      AllocFault::disarm();
      TEST_TRUE(bad_alloc)
      TEST_TRUE(File::exists(file) && slurp4b(file).empty())
    }
#endif
    std::remove(file.c_str());
  }
#else
  STATUS(fault_injection_unavailable)
#endif
}
END_SECTION

START_SECTION([EXTRA] store - a change of the working directory while store() runs does not make it touch another file)
{
#if FAULT_INJECTION_TESTS
  if (!allocFaultReachesLibOpenMS())
  {
    STATUS(fault_injection_unreached)
  }
  else
  {
    // store("out.idXML") in directory A writes A/out.idXML. If the working directory changes to B (e.g. in another
    // thread) before store() fails, B/out.idXML, a file this call did not write, must stay unchanged (regression guard:
    // an earlier version removed the output of a failed store() by its name), and the partial A/out.idXML is left in
    // place. The change happens at each allocation of the calling thread after the file was created, one after the
    // other; store() fails in the protein section (a meta value without a value).
    namespace fs = std::filesystem;
    // not NEW_TMP_FILE: the directories are removed at the end of the section
    const fs::path base = fs::absolute("IdXMLFile_test_" + File::getUniqueName(false) + "_cwd");
    const fs::path dir_a = base / "A";
    cwd_test_dir_b = base / "B";
    fs::create_directories(dir_a);
    fs::create_directories(cwd_test_dir_b);
    cwd_test_out_a = dir_a / "out.idXML";
    const fs::path out_b = cwd_test_dir_b / "out.idXML";
    const std::string unrelated = "a file in B that store() does not write\n";

    std::vector<ProteinIdentification> prots(1);
    prots[0].setIdentifier("runCwd");
    prots[0].setDateTime(DateTime::now());
    ProteinHit protein_hit(0.0, 1, "ACC_CWD", "");
    protein_hit.setMetaValue("fail_test_empty", DataValue()); // cannot be written
    prots[0].insertHit(protein_hit);
    PeptideIdentificationList peps(2);
    for (Size l = 0; l < peps.size(); ++l)
    {
      peps[l].setIdentifier("runCwd");
      peps[l].setScoreType("score");
      PeptideHit hit(double(l), 0, 2, AASequence::fromString("PEPTIDER"));
      peps[l].insertHit(hit);
    }

    struct RestoreWorkingDirectory
    {
      fs::path saved = fs::current_path();
      ~RestoreWorkingDirectory() { std::error_code ec; fs::current_path(saved, ec); }
    } restore_working_directory;
    long changed = 0, b_damaged = 0, a_left = 0, first_b_damaged = 0;
    for (long n = 1; n < 10000000; ++n)
    {
      { std::ofstream(out_b.string()) << unrelated; }
      fs::remove(cwd_test_out_a);
      fs::current_path(dir_a);
      cwd_test_changed = false;
      AllocFault::arm(n, &changeToDirBOnceOutputExists);
      try
      {
        IdXMLFile().store("out.idXML", prots, peps);
      }
      catch (...)
      {
      }
      const bool fired = AllocFault::fired;
      AllocFault::disarm();
      fs::current_path(restore_working_directory.saved);
      if (!fired) break; // n is past the last allocation of store()
      if (!cwd_test_changed) continue; // the change came before the file was created: B/out.idXML was the output
      ++changed;
      if (!fs::exists(out_b) || slurp4b(out_b.string()) != unrelated)
      {
        if (b_damaged++ == 0) first_b_damaged = n;
      }
      if (fs::exists(cwd_test_out_a)) ++a_left;
    }
    STATUS(changed << " stores changed the working directory after the file was created; B/out.idXML removed or changed in "
           << b_damaged << " (first at allocation " << first_b_damaged << "), A/out.idXML left in " << a_left)
    TEST_TRUE(changed > 0)
    TEST_EQUAL(b_damaged, 0)
    TEST_EQUAL(a_left, changed)

    // a working directory that no longer exists: a relative name cannot be resolved, store() reports that it cannot
    // create the file (as when opening it fails)
    const fs::path gone = base / "gone";
    fs::create_directories(gone);
    fs::current_path(gone);
    fs::remove(gone);
    TEST_EXCEPTION(Exception::UnableToCreateFile, IdXMLFile().store("out.idXML", prots, peps))
    fs::current_path(restore_working_directory.saved);
    std::error_code ec;
    fs::remove_all(base, ec);
  }
#else
  STATUS(fault_injection_unavailable)
#endif
}
END_SECTION

START_SECTION([EXTRA] store - a relative name is written where the working directory can be used but not its absolute path)
{
#ifdef __linux__
  // Creating a file by a relative name needs permissions on the working directory only; its absolute path also needs
  // search permission on every ancestor. Here the parent of the working directory has no permissions: store() of a
  // relative name must write the file, as any other writer can (a privileged user is not restricted; then the case is
  // not reached and only reported). Not NEW_TMP_FILE: the directories are removed at the end of the section.
  namespace fs = std::filesystem;
  const fs::path base = fs::absolute("IdXMLFile_test_" + File::getUniqueName(false) + "_ancestor");
  const fs::path parent = base / "parent";
  const fs::path work = parent / "work";
  fs::create_directories(work);
  // restores the permissions and the working directory and removes the directories, also if the section stops early
  struct RestoreDirectories
  {
    fs::path saved = fs::current_path();
    fs::path locked;
    fs::path remove_tree;
    ~RestoreDirectories()
    {
      std::error_code ec;
      if (!locked.empty()) fs::permissions(locked, fs::perms::owner_all, ec);
      fs::current_path(saved, ec);
      if (!remove_tree.empty()) fs::remove_all(remove_tree, ec);
    }
  } restore;
  restore.remove_tree = base;

  std::vector<ProteinIdentification> prots(1);
  prots[0].setIdentifier("runAncestor");
  prots[0].setDateTime(DateTime::now());
  prots[0].insertHit(ProteinHit(0.0, 1, "ACC_ANCESTOR", ""));
  PeptideIdentificationList peps(2);
  for (Size l = 0; l < peps.size(); ++l)
  {
    peps[l].setIdentifier("runAncestor");
    peps[l].setScoreType("score");
    PeptideHit hit(double(l), 0, 2, AASequence::fromString("PEPTIDER"));
    hit.addPeptideEvidence(PeptideEvidence("ACC_ANCESTOR", 0, 7, '-', '-'));
    peps[l].insertHit(hit);
  }
  const fs::path reference = base / "reference.idXML";
  IdXMLFile().store(reference.string(), prots, peps);
  const std::string complete = slurp4b(reference.string());

  fs::current_path(work);
  std::error_code lock_error; // a file system without permissions: the case is not reached
  fs::permissions(parent, fs::perms::none, lock_error);
  if (!lock_error) restore.locked = parent;
  const bool reached = !lock_error && ::access((work / "control").c_str(), F_OK) != 0 && errno == EACCES;
  { std::ofstream("control") << "x"; }
  const bool relative_name_usable = File::exists("control");
  STATUS("absolute path of the working directory not searchable: " << reached << "; a relative name can be created: "
         << relative_name_usable)
  if (reached && relative_name_usable)
  {
    bool stored = false;
    try
    {
      IdXMLFile().store("out.idXML", prots, peps);
      stored = true;
    }
    catch (const Exception::BaseException& e)
    {
      STATUS("store() threw: " << e.what());
    }
    TEST_TRUE(stored)
    TEST_TRUE(File::exists("out.idXML") && slurp4b("out.idXML") == complete)
  }

  // a pre-existing empty file that cannot be opened for writing stays (opening it changed nothing)
  std::error_code ec;
  fs::permissions(parent, fs::perms::owner_all, ec);
  restore.locked.clear();
  fs::current_path(restore.saved);
  const fs::path read_only = work / "read_only_empty.idXML";
  { std::ofstream create(read_only.string()); }
  std::error_code read_only_error;
  fs::permissions(read_only, fs::perms::owner_read, read_only_error);
  if (!read_only_error && ::access(read_only.c_str(), W_OK) != 0)
  {
    TEST_EXCEPTION(Exception::UnableToCreateFile, IdXMLFile().store(read_only.string(), prots, peps))
    TEST_TRUE(fs::exists(read_only))
  }
  fs::permissions(read_only, fs::perms::owner_all, ec);
#else
  STATUS("SKIPPED (Linux only): POSIX permissions on the ancestors of the working directory")
#endif
}
END_SECTION

START_SECTION([EXTRA] store - a file renamed into the place of the output while store() runs is not removed)
{
#if FAULT_INJECTION_TESTS
  if (!allocFaultReachesLibOpenMS())
  {
    STATUS(fault_injection_unreached)
  }
  else
  {
    // store() fails in the protein section (a meta value without a value), after the proteins before it have overflowed
    // the stream buffer, so the output has been written to. At each allocation of the calling thread after that, one
    // after the other, the output is moved aside and another file is renamed into its place, as another process could.
    // store() removes nothing: that file must stay, with its content, and the output moved aside is left in place
    // (regression guard: an earlier version removed the output of a failed store() by its name). Not NEW_TMP_FILE: the
    // directory is removed at the end of the section.
    namespace fs = std::filesystem;
    struct RemoveTree
    {
      fs::path tree;
      ~RemoveTree() { std::error_code ec; fs::remove_all(tree, ec); }
    } remove_tree{fs::absolute("IdXMLFile_test_" + File::getUniqueName(false) + "_replaced")};
    fs::create_directories(remove_tree.tree);
    swap_test_out = remove_tree.tree / "out.idXML";
    swap_test_aside = remove_tree.tree / "aside.idXML";
    swap_test_sentinel = remove_tree.tree / "sentinel.idXML";
    const std::string sentinel = "a file that was renamed into the place of the output\n";

    std::vector<ProteinIdentification> prots(1);
    prots[0].setIdentifier("runReplaced");
    prots[0].setDateTime(DateTime::now());
    for (Size i = 0; i < 120; ++i)
    {
      ProteinHit hit(0.0, 1, "ACC_REPLACED_" + StringUtils::toStr(i), "");
      hit.setMetaValue("replaced_test", std::string(64, 'x'));
      prots[0].insertHit(hit);
    }
    ProteinHit bad_hit(0.0, 1, "ACC_REPLACED_BAD", "");
    bad_hit.setMetaValue("fail_test_empty", DataValue()); // cannot be written
    prots[0].insertHit(bad_hit);
    PeptideIdentificationList peps(2);
    for (Size l = 0; l < peps.size(); ++l)
    {
      peps[l].setIdentifier("runReplaced");
      peps[l].setScoreType("score");
      PeptideHit hit(double(l), 0, 2, AASequence::fromString("PEPTIDER"));
      peps[l].insertHit(hit);
    }

    long swapped = 0, sentinel_lost = 0, first_lost = 0, other_outcome = 0, aside_lost = 0;
    for (long n = 1; n < 10000000; ++n)
    {
      std::error_code ec;
      fs::remove(swap_test_out, ec);
      fs::remove(swap_test_aside, ec);
      { std::ofstream(swap_test_sentinel.string()) << sentinel; }
      swap_test_swapped = false;
      bool conversion_error = false;
      AllocFault::arm(n, &swapInSentinelOnceOutputWritten);
      try
      {
        IdXMLFile().store(swap_test_out.string(), prots, peps);
      }
      catch (const Exception::ConversionError&)
      {
        conversion_error = true;
      }
      catch (...)
      {
      }
      const bool fired = AllocFault::fired;
      AllocFault::disarm();
      if (!fired) break; // n is past the last allocation of store()
      if (!swap_test_swapped) continue; // the output had not been written to yet
      ++swapped;
      if (!conversion_error) ++other_outcome;
      if (!fs::exists(swap_test_out) || slurp4b(swap_test_out.string()) != sentinel)
      {
        if (sentinel_lost++ == 0) first_lost = n;
      }
      if (!fs::exists(swap_test_aside)) ++aside_lost;
    }
    STATUS(swapped << " stores had another file renamed into the place of their output; it was removed or changed in "
           << sentinel_lost << " (first at allocation " << first_lost << "); stores that did not fail as expected: "
           << other_outcome << "; the output moved aside was removed in " << aside_lost)
    TEST_TRUE(swapped > 0)
    TEST_EQUAL(sentinel_lost, 0)
    TEST_EQUAL(other_outcome, 0)
    TEST_EQUAL(aside_lost, 0)
  }
#else
  STATUS(fault_injection_unavailable)
#endif
}
END_SECTION

START_SECTION([EXTRA] store - a file renamed into the place of the output while store() opens it is not removed)
{
#if FAULT_INJECTION_TESTS
  if (!allocFaultReachesLibOpenMS())
  {
    STATUS(fault_injection_unreached)
  }
  else
  {
    // As in the section above, but the output is moved aside and another file is renamed into its place while the
    // output exists and is still empty: at each such allocation of the calling thread, one after the other. With
    // libstdc++, the first is the allocation of the stream buffer within opening the file. Two cases:
    // - a file with content is renamed into place and store() continues; it fails in the protein section (a meta
    //   value without a value);
    // - an empty file is renamed into place and the allocation fails (std::bad_alloc); within opening, opening fails.
    //   The data have no other error here: an allocation that fails while an OpenMS exception is constructed
    //   terminates the program (the constructors are noexcept).
    // store() removes nothing: the file renamed into place must stay, the same file with its content, and the output
    // moved aside is left in place. Not NEW_TMP_FILE: the directory is removed at the end of the section.
    namespace fs = std::filesystem;
    struct RemoveTree
    {
      fs::path tree;
      ~RemoveTree() { std::error_code ec; fs::remove_all(tree, ec); }
    } remove_tree{fs::absolute("IdXMLFile_test_" + File::getUniqueName(false) + "_replaced_open")};
    fs::create_directories(remove_tree.tree);
    swap_test_out = remove_tree.tree / "out.idXML";
    swap_test_aside = remove_tree.tree / "aside.idXML";
    swap_test_sentinel = remove_tree.tree / "sentinel.idXML";

    std::vector<ProteinIdentification> prots(1);
    prots[0].setIdentifier("runReplacedOpen");
    prots[0].setDateTime(DateTime::now());
    prots[0].insertHit(ProteinHit(0.0, 1, "ACC_REPLACED_OPEN", ""));
    const std::vector<ProteinIdentification> prots_valid = prots;
    ProteinHit bad_hit(0.0, 1, "ACC_REPLACED_OPEN_BAD", "");
    bad_hit.setMetaValue("fail_test_empty", DataValue()); // cannot be written
    prots[0].insertHit(bad_hit);
    PeptideIdentificationList peps(2);
    for (Size l = 0; l < peps.size(); ++l)
    {
      peps[l].setIdentifier("runReplacedOpen");
      peps[l].setScoreType("score");
      PeptideHit hit(double(l), 0, 2, AASequence::fromString("PEPTIDER"));
      peps[l].insertHit(hit);
    }

    for (const bool empty_and_fail : {false, true})
    {
      const std::string sentinel = empty_and_fail ? std::string()
                                                  : std::string("a file that was renamed into the place of the output\n");
      long swapped = 0, sentinel_lost = 0, first_lost = 0, other_outcome = 0, aside_lost = 0;
      for (long n = 1; n < 10000000; ++n)
      {
        std::error_code ec;
        fs::remove(swap_test_out, ec);
        fs::remove(swap_test_aside, ec);
        { std::ofstream(swap_test_sentinel.string()) << sentinel; }
        struct stat sentinel_stat;
        ABORT_IF(::lstat(swap_test_sentinel.c_str(), &sentinel_stat) != 0)
        swap_test_swapped = false;
        // with content: store() must fail in the protein section; empty: it fails or, after a failed allocation that it
        // handled (a nothrow allocation), completes into the output moved aside
        bool expected_error = empty_and_fail;
        AllocFault::arm(n, &swapInSentinelWhileOutputEmpty, empty_and_fail);
        try
        {
          IdXMLFile().store(swap_test_out.string(), empty_and_fail ? prots_valid : prots, peps);
        }
        catch (const Exception::ConversionError&)
        {
          expected_error = !empty_and_fail;
        }
        catch (...)
        {
          expected_error = empty_and_fail; // the std::bad_alloc, or UnableToCreateFile if opening failed
        }
        const bool fired = AllocFault::fired;
        AllocFault::disarm();
        if (!fired) break; // n is past the last allocation of store()
        if (!swap_test_swapped) continue; // the output did not exist yet, or had been written to
        ++swapped;
        if (!expected_error) ++other_outcome;
        struct stat out_stat;
        const bool same_file = ::lstat(swap_test_out.c_str(), &out_stat) == 0 && out_stat.st_dev == sentinel_stat.st_dev
                               && out_stat.st_ino == sentinel_stat.st_ino;
        if (!same_file || slurp4b(swap_test_out.string()) != sentinel)
        {
          if (sentinel_lost++ == 0) first_lost = n;
        }
        if (!fs::exists(swap_test_aside)) ++aside_lost;
      }
      STATUS((empty_and_fail ? "an empty file renamed into place, then the allocation fails: "
                             : "a file with content renamed into place: ")
             << swapped << " stores; the file renamed into place was removed or changed in " << sentinel_lost
             << " (first at allocation " << first_lost << "); stores that did not fail as expected: " << other_outcome
             << "; the output moved aside was removed in " << aside_lost)
      TEST_TRUE(swapped > 0)
      TEST_EQUAL(sentinel_lost, 0)
      TEST_EQUAL(other_outcome, 0)
      TEST_EQUAL(aside_lost, 0)
    }
  }
#else
  STATUS(fault_injection_unavailable)
#endif
}
END_SECTION

START_SECTION([EXTRA] store - an exception from closing the file does not replace the error that store() reports)
{
#ifdef __GLIBCXX__
  // A code conversion facet whose unshift() throws makes closing a written file stream throw. The stream of store() is
  // closed when store() fails (here: a meta value of a peptide hit that cannot be written, in a block of the parallel
  // writer); that must not replace the error, nor its note on the partial file. Without another error, the exception
  // of close() (not a std::exception) is raised unchanged, and the file stays. libstdc++ closes the file and then
  // rethrows from close(); other standard libraries leave a stream whose close() threw in a state that its destructor
  // cannot handle, so the section is libstdc++-only.
  struct CloseError
  {
  };
  struct ThrowOnUnshift : public std::codecvt<char, char, std::mbstate_t>
  {
  protected:
    bool do_always_noconv() const noexcept override { return false; }
    result do_out(state_type&, const char* from, const char* from_end, const char*& from_next, char* to, char* to_end,
                  char*& to_next) const override
    {
      const std::ptrdiff_t n = std::min(from_end - from, to_end - to);
      std::copy(from, from + n, to);
      from_next = from + n;
      to_next = to + n;
      return from_next == from_end ? ok : partial;
    }
    result do_unshift(state_type&, char*, char*, char*&) const override { throw CloseError(); }
  };
  struct RestoreGlobalLocale
  {
    std::locale saved;
    ~RestoreGlobalLocale() { std::locale::global(saved); }
  };

  std::vector<ProteinIdentification> prots(1);
  prots[0].setIdentifier("runClose");
  prots[0].setDateTime(DateTime::now());
  prots[0].insertHit(ProteinHit(0.0, 1, "ACC_CLOSE", ""));
  PeptideIdentificationList peps_valid(100);
  for (Size l = 0; l < peps_valid.size(); ++l)
  {
    peps_valid[l].setIdentifier("runClose");
    peps_valid[l].setScoreType("score");
    PeptideHit hit(double(l), 0, 2, AASequence::fromString("PEPTIDER"));
    hit.addPeptideEvidence(PeptideEvidence("ACC_CLOSE", 0, 7, '-', '-'));
    peps_valid[l].insertHit(hit);
  }
  PeptideIdentificationList peps = peps_valid;
  for (PeptideHit& hit : peps[40].getHits()) hit.setMetaValue("fail_test_empty", DataValue()); // cannot be written
  std::string file;
  NEW_TMP_FILE(file)
  {
    // the facet throws as intended: closing a written stream throws CloseError
    RestoreGlobalLocale restore{std::locale::global(std::locale(std::locale(), new ThrowOnUnshift))};
    std::ofstream probe(file);
    probe << "x";
    bool close_threw = false;
    try
    {
      probe.close();
    }
    catch (const CloseError&)
    {
      close_threw = true;
    }
    TEST_TRUE(close_threw)
  }
  std::remove(file.c_str());
  {
    RestoreGlobalLocale restore{std::locale::global(std::locale(std::locale(), new ThrowOnUnshift))};
    const std::string message = messageOf<Exception::ConversionError>([&] { IdXMLFile().store(file, prots, peps); });
    TEST_TRUE(saysPartialFileRemains(message, file))
  }
  TEST_TRUE(isPartialIdXML(file))
  std::remove(file.c_str()); // partial, not to be validated
  {
    // only close() fails
    RestoreGlobalLocale restore{std::locale::global(std::locale(std::locale(), new ThrowOnUnshift))};
    bool close_error = false;
    try
    {
      IdXMLFile().store(file, prots, peps_valid);
    }
    catch (const CloseError&)
    {
      close_error = true;
    }
    catch (...)
    {
    }
    TEST_TRUE(close_error)
  }
  TEST_TRUE(File::exists(file))
  std::remove(file.c_str());
#else
  STATUS("SKIPPED (libstdc++ only): a stream whose close() threw")
#endif
}
END_SECTION

START_SECTION([EXTRA] store - an error keeps its type even if memory runs out while store() reports it)
{
#if FAULT_INJECTION_TESTS && defined(__GLIBCXX__)
  if (!allocFaultReachesLibOpenMS())
  {
    STATUS(fault_injection_unreached)
  }
  else
  {
    // The final close() of store() throws (see closeTestChild()): an OpenMS exception, a std::bad_alloc, or an exception
    // that is not a std::exception. Then every allocation after the throw fails once, one after the other, each in a
    // child process: an allocation failure while an OpenMS exception is constructed ends the program (the constructors
    // are noexcept). store() must raise the error with its type (an OpenMS exception with the note, unless the failed
    // allocation was one of the note's; the others unchanged), and the output stays.
    std::vector<ProteinIdentification> prots(1);
    prots[0].setIdentifier("runReport");
    prots[0].setDateTime(DateTime::now());
    prots[0].insertHit(ProteinHit(0.0, 1, "ACC_REPORT", ""));
    PeptideIdentificationList peps(3); // one block, formatted by the calling thread
    for (Size l = 0; l < peps.size(); ++l)
    {
      peps[l].setIdentifier("runReport");
      peps[l].setScoreType("score");
      PeptideHit hit(double(l), 0, 2, AASequence::fromString("PEPTIDER"));
      hit.addPeptideEvidence(PeptideEvidence("ACC_REPORT", 0, 7, '-', '-'));
      peps[l].insertHit(hit);
    }
    // longer than any short-string buffer, so that copying the name allocates; not a NEW_TMP_FILE (removed below)
    const std::string file = "IdXMLFile_test_" + File::getUniqueName(false) + "_" + std::string(100, 'r') + ".idXML";
    for (const int kind : {0, 1, 2})
    {
      long runs = 0, ended = 0, wrong = 0, first_wrong = -1;
      for (long n = 0; n < 10000; ++n)
      {
        std::remove(file.c_str());
        std::cout.flush();
        std::cerr.flush();
        const pid_t pid = ::fork();
        ABORT_IF(pid < 0)
        if (pid == 0) ::_exit(closeTestChild(kind, n, file, prots, peps));
        int status = 0;
        while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
        ++runs;
        const bool exited = WIFEXITED(status);
        if (!exited) ++ended;
        const bool ok = exited && (WEXITSTATUS(status) & 1) == 0 && File::exists(file);
        if (!ok && wrong++ == 0) first_wrong = n;
        if (n > 0 && exited && (WEXITSTATUS(status) & 16) == 0) break; // n is past the last allocation after the throw
      }
      STATUS((kind == 0 ? "Exception::ConversionError" : kind == 1 ? "std::bad_alloc" : "not a std::exception")
             << " from close(): " << runs << " runs (the first without a failed allocation); ended by a signal: " << ended
             << "; raised another error, an OpenMS error without the note, or removed the output: " << wrong
             << " (first at allocation " << first_wrong << ")")
      TEST_EQUAL(wrong, 0)
    }
    std::remove(file.c_str());
  }
#elif FAULT_INJECTION_TESTS
  STATUS("SKIPPED (libstdc++ only): a stream whose close() throws")
#else
  STATUS(fault_injection_unavailable)
#endif
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// check the temporary files written above against their XML schema (types without a validator are skipped)
VALIDATE_TMP_FILES

END_TEST
