# RibonucleotideDB Provider Separation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Decouple file I/O (JSON + TSV parsing) from RibonucleotideDB indexing logic using the same data provider pattern established for ModificationsDB.

**Architecture:** Extract `readFromJSON_()` and `readFromFile_()` into concrete `RibonucleotideDataProvider` implementations. RibonucleotideDB gets a public provider-based constructor for DI/testing, while `getInstance()` remains unchanged. The ambiguity map is populated by the DB during indexing (not by providers) since it requires cross-referencing already-loaded ribonucleotides.

**Tech Stack:** C++20, nlohmann::json, Qt (QFile/QTextStream for Unicode), OpenMS test framework

---

### Task 1: Create RibonucleotideDataProvider Interface

**Files:**
- Create: `src/openms/include/OpenMS/CHEMISTRY/RibonucleotideDataProvider.h`

**Context:** This mirrors `ModificationDataProvider.h` but for ribonucleotides. The key difference is that ribonucleotide providers must also communicate ambiguity information (some mods map to pairs of alternatives). We handle this with a `RibonucleotideEntry` struct.

**Step 1: Write the interface header**

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <memory>
#include <vector>

namespace OpenMS
{
  /// A single ribonucleotide entry as returned by a RibonucleotideDataProvider.
  struct OPENMS_DLLAPI RibonucleotideEntry
  {
    std::unique_ptr<Ribonucleotide> ribonucleotide;
    /// If non-empty, this is an ambiguous code mapping to two alternative codes
    String alternative_code_1;
    String alternative_code_2;
    bool isAmbiguous() const { return !alternative_code_1.empty(); }
  };

  /**
    @brief Interface for providing Ribonucleotide data to RibonucleotideDB.

    Implementations abstract the source of ribonucleotide data,
    enabling dependency injection. File-based providers (ModomicsJSONDataProvider,
    RibonucleotideTSVDataProvider) handle I/O; InMemoryRibonucleotideDataProvider
    supports testing.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI RibonucleotideDataProvider
  {
  public:
    virtual ~RibonucleotideDataProvider() = default;

    /**
      @brief Load ribonucleotides from whatever source this provider wraps.
      @return Vector of ribonucleotide entries with ownership transferred to caller.
      @note Providers may only be called once. Subsequent calls may return empty results.
    */
    virtual std::vector<RibonucleotideEntry> loadRibonucleotides() = 0;
  };

  /**
    @brief Data provider that serves pre-built ribonucleotides from memory.

    Useful for unit testing RibonucleotideDB without requiring files on disk.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI InMemoryRibonucleotideDataProvider : public RibonucleotideDataProvider
  {
  public:
    explicit InMemoryRibonucleotideDataProvider(std::vector<RibonucleotideEntry> entries)
      : entries_(std::move(entries))
    {
    }

    InMemoryRibonucleotideDataProvider(const InMemoryRibonucleotideDataProvider&) = delete;
    InMemoryRibonucleotideDataProvider& operator=(const InMemoryRibonucleotideDataProvider&) = delete;
    InMemoryRibonucleotideDataProvider(InMemoryRibonucleotideDataProvider&&) = default;
    InMemoryRibonucleotideDataProvider& operator=(InMemoryRibonucleotideDataProvider&&) = default;

    std::vector<RibonucleotideEntry> loadRibonucleotides() override
    {
      return std::move(entries_);
    }

  private:
    std::vector<RibonucleotideEntry> entries_;
  };

} // namespace OpenMS
```

**Step 2: Commit**
```bash
git add src/openms/include/OpenMS/CHEMISTRY/RibonucleotideDataProvider.h
git commit -m "feat: add RibonucleotideDataProvider interface and InMemoryRibonucleotideDataProvider"
```

---

### Task 2: Create ModomicsJSONDataProvider

**Files:**
- Create: `src/openms/include/OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h`
- Create: `src/openms/source/CHEMISTRY/ModomicsJSONDataProvider.cpp`

**Context:** This extracts `readFromJSON_()` from RibonucleotideDB.cpp (lines 229-283) along with the helper functions `entryIsWellFormed_()` (lines 78-104), `getBaseLossFormula_()` (lines 107-139), and `parseEntry_()` (lines 142-226). These are currently free functions in the `OpenMS` namespace in RibonucleotideDB.cpp. Move them into the provider's `.cpp` as file-local helpers (anonymous namespace).

**Key difference from readFromJSON_:** The original method directly populates `ribonucleotides_`, `code_map_`, `max_code_length_`, and `ambiguity_map_`. The provider instead returns `vector<RibonucleotideEntry>` — the DB handles indexing later.

The `ambiguity_map_` population in the original calls `getRibonucleotide()` (line 268) to resolve alternative codes to pointers. Since the provider doesn't have access to the DB, it stores the alternative codes as strings in `RibonucleotideEntry` and lets the DB resolve them during indexing.

**Step 1: Write the header**

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>

namespace OpenMS
{
  /**
    @brief Provides ribonucleotide data from Modomics JSON files.

    Wraps the Modomics JSON format (from https://www.genesilico.pl/modomics/api/modifications).

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ModomicsJSONDataProvider : public RibonucleotideDataProvider
  {
  public:
    explicit ModomicsJSONDataProvider(const String& filename);
    std::vector<RibonucleotideEntry> loadRibonucleotides() override;

  private:
    String filename_;
  };

} // namespace OpenMS
```

**Step 2: Write the implementation**

Move `readFromJSON_()` body, `entryIsWellFormed_()`, `getBaseLossFormula_()`, `parseEntry_()` into `ModomicsJSONDataProvider.cpp`. Key changes:

- `entryIsWellFormed_()`, `getBaseLossFormula_()`, `parseEntry_()` go into anonymous namespace (they're already free functions in the original)
- The `ParsedEntry_` struct stays as a file-local helper
- Instead of populating `code_map_`/`ribonucleotides_`/`ambiguity_map_`, build and return `vector<RibonucleotideEntry>`
- Instead of calling `getRibonucleotide(entry.alternative_1)` (which requires the DB), store alternative codes as strings

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <nlohmann/json.hpp>

using namespace std;

// Specialize nlohmann::adl_serializer for OpenMS::EmpiricalFormula
namespace nlohmann {
    template <>
    struct adl_serializer<OpenMS::EmpiricalFormula>
    {
        static void from_json(const json& j, OpenMS::EmpiricalFormula& ef)
        {
            std::string formula_string = j.get<std::string>();
            ef = OpenMS::EmpiricalFormula(formula_string);
        }

        static void to_json(json& j, const OpenMS::EmpiricalFormula& ef)
        {
            j = ef.toString();
        }
    };
}

namespace OpenMS
{
  namespace // anonymous namespace for file-local helpers
  {
    // [COPY entryIsWellFormed_() from RibonucleotideDB.cpp lines 78-104 verbatim]
    // [COPY getBaseLossFormula_() from RibonucleotideDB.cpp lines 107-139 verbatim]

    // Same as ParsedEntry_ in RibonucleotideDB.cpp
    struct ParsedEntry_
    {
      unique_ptr<Ribonucleotide> ribo;
      String alternative_1;
      String alternative_2;
      bool isAmbiguous() { return !alternative_1.empty(); }
    };

    // [COPY parseEntry_() from RibonucleotideDB.cpp lines 142-226 verbatim]
  } // anonymous namespace

  ModomicsJSONDataProvider::ModomicsJSONDataProvider(const String& filename)
    : filename_(filename)
  {
  }

  std::vector<RibonucleotideEntry> ModomicsJSONDataProvider::loadRibonucleotides()
  {
    using json = nlohmann::json;

    String full_path = File::find(filename_);

    QFile file(full_path.toQString());
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, full_path);
    }

    QTextStream source(&file);
    source.setAutoDetectUnicode(true);
    Size line_count = 0;
    json mod_obj;
    try
    {
      mod_obj = json::parse(String(source.readAll()));
    }
    catch (Exception::ParseError& e)
    {
      OPENMS_LOG_ERROR << "Error: Failed to parse Modomics JSON. Reason:\n" << e.getName() << " - " << e.what() << endl;
      throw;
    }

    std::vector<RibonucleotideEntry> result;
    for (auto& element : mod_obj)
    {
      line_count++;
      try
      {
        entryIsWellFormed_(element);
        ParsedEntry_ parsed = parseEntry_(element);

        if (parsed.ribo->getCode() != "")
        {
          RibonucleotideEntry entry;
          if (parsed.isAmbiguous())
          {
            entry.alternative_code_1 = parsed.alternative_1;
            entry.alternative_code_2 = parsed.alternative_2;
          }
          entry.ribonucleotide = std::move(parsed.ribo);
          result.push_back(std::move(entry));
        }
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Error: Failed to parse input element " << line_count << ". Reason:\n" << e.getName() << " - " << e.what() << "\nSkipping this line." << endl;
      }
    }
    return result;
  }

} // namespace OpenMS
```

**Step 3: Commit**
```bash
git add src/openms/include/OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h \
        src/openms/source/CHEMISTRY/ModomicsJSONDataProvider.cpp
git commit -m "feat: add ModomicsJSONDataProvider for JSON-based ribonucleotide loading"
```

---

### Task 3: Create RibonucleotideTSVDataProvider

**Files:**
- Create: `src/openms/include/OpenMS/CHEMISTRY/RibonucleotideTSVDataProvider.h`
- Create: `src/openms/source/CHEMISTRY/RibonucleotideTSVDataProvider.cpp`

**Context:** Extracts `readFromFile_()` (lines 286-334) and `parseRow_()` (lines 337-429) from RibonucleotideDB.cpp. The TSV parsing handles Unicode prime characters and ambiguity codes.

**Key difference from readFromFile_:** Instead of populating DB members directly, returns `vector<RibonucleotideEntry>`. The `parseRow_()` helper that references `ambiguity_map_` and `getRibonucleotide()` (line 421) now stores alternative codes as strings instead.

**Step 1: Write the header**

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>

namespace OpenMS
{
  /**
    @brief Provides ribonucleotide data from tab-separated text files.

    Reads the OpenMS custom RNA modifications TSV format with columns:
    name, short_name, new_nomenclature, originating_base, rnamods_abbrev,
    html_abbrev, formula, monoisotopic_mass, average_mass.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI RibonucleotideTSVDataProvider : public RibonucleotideDataProvider
  {
  public:
    explicit RibonucleotideTSVDataProvider(const String& filename);
    std::vector<RibonucleotideEntry> loadRibonucleotides() override;

  private:
    String filename_;
  };

} // namespace OpenMS
```

**Step 2: Write the implementation**

Move `readFromFile_()` and `parseRow_()` into the provider. Key changes:
- `parseRow_()` becomes a file-local helper function
- Ambiguity handling: instead of calling `getRibonucleotide(code1)` / `getRibonucleotide(code2)` to get pointers, store codes as strings in `RibonucleotideEntry`
- `parseRow_()` returns a `RibonucleotideEntry` instead of `unique_ptr<Ribonucleotide>`

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RibonucleotideTSVDataProvider.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtCore/QFile>
#include <QtCore/QTextStream>

using namespace std;

namespace OpenMS
{
  namespace // anonymous namespace for file-local helpers
  {
    // Adapted from RibonucleotideDB::parseRow_ — returns RibonucleotideEntry instead of unique_ptr
    RibonucleotideEntry parseRow_(const std::string& row, Size line_count)
    {
      // [COPY parseRow_() body from RibonucleotideDB.cpp lines 339-428]
      // Key changes:
      // 1. Return type is RibonucleotideEntry (not unique_ptr<Ribonucleotide>)
      // 2. Lines 413-421 (ambiguity handling): instead of:
      //      ambiguity_map_[parts[1]] = make_pair(getRibonucleotide(code1), getRibonucleotide(code2));
      //    do:
      //      entry.alternative_code_1 = code1;
      //      entry.alternative_code_2 = code2;
      // 3. Wrap in a RibonucleotideEntry at the end
    }
  } // anonymous namespace

  RibonucleotideTSVDataProvider::RibonucleotideTSVDataProvider(const String& filename)
    : filename_(filename)
  {
  }

  std::vector<RibonucleotideEntry> RibonucleotideTSVDataProvider::loadRibonucleotides()
  {
    // [COPY readFromFile_() body from RibonucleotideDB.cpp lines 288-333]
    // Key changes:
    // 1. Build and return vector<RibonucleotideEntry> instead of populating DB members
    // 2. Call parseRow_() (free function) and push_back into result vector
    // 3. No code_map_ or max_code_length_ updates (DB handles indexing)
  }

} // namespace OpenMS
```

**Step 3: Commit**
```bash
git add src/openms/include/OpenMS/CHEMISTRY/RibonucleotideTSVDataProvider.h \
        src/openms/source/CHEMISTRY/RibonucleotideTSVDataProvider.cpp
git commit -m "feat: add RibonucleotideTSVDataProvider for TSV-based ribonucleotide loading"
```

---

### Task 4: Add Provider-Based Constructor to RibonucleotideDB

**Files:**
- Modify: `src/openms/include/OpenMS/CHEMISTRY/RibonucleotideDB.h`
- Modify: `src/openms/source/CHEMISTRY/RibonucleotideDB.cpp`

**Context:** Add a public constructor that accepts providers, add `loadFromProviders_()`, and rewire the default constructor/singleton to use providers internally. Remove `readFromJSON_()`, `readFromFile_()`, `parseRow_()`.

**Step 1: Modify the header**

Changes to `RibonucleotideDB.h`:

1. Add `#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>`
2. Add public constructor: `explicit RibonucleotideDB(std::vector<std::unique_ptr<RibonucleotideDataProvider>> providers);`
3. Remove protected methods: `readFromFile_()`, `readFromJSON_()`, `parseRow_()`
4. Add private: `void loadFromProviders_(std::vector<std::unique_ptr<RibonucleotideDataProvider>>& providers);`

Result:
```cpp
#pragma once

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <memory>
#include <unordered_map>

namespace OpenMS
{
  class OPENMS_DLLAPI RibonucleotideDB final
  {
  public:
    using ConstRibonucleotidePtr = const Ribonucleotide *;
    typedef std::vector<std::unique_ptr<Ribonucleotide>>::const_iterator ConstIterator;

    static RibonucleotideDB* getInstance();

    /**
      @brief Construct from data providers (non-singleton, for DI/testing).
    */
    explicit RibonucleotideDB(std::vector<std::unique_ptr<RibonucleotideDataProvider>> providers);

    ~RibonucleotideDB() = default;

    RibonucleotideDB(const RibonucleotideDB& other) = delete;
    RibonucleotideDB& operator=(const RibonucleotideDB& other) = delete;

    // ... all public accessors unchanged (begin, end, getRibonucleotide, etc.) ...

  protected:
    RibonucleotideDB();  // used by singleton path

    // readFromFile_, readFromJSON_, parseRow_ REMOVED

    std::vector<std::unique_ptr<Ribonucleotide>> ribonucleotides_;
    std::unordered_map<std::string, Size> code_map_;
    std::map<std::string, std::pair<ConstRibonucleotidePtr, ConstRibonucleotidePtr>> ambiguity_map_;
    Size max_code_length_;

  private:
    void loadFromProviders_(std::vector<std::unique_ptr<RibonucleotideDataProvider>>& providers);
  };
}
```

**Step 2: Modify the implementation**

Changes to `RibonucleotideDB.cpp`:

1. Add includes for `ModomicsJSONDataProvider.h` and `RibonucleotideTSVDataProvider.h`
2. Remove the `nlohmann` serializer, `ParsedEntry_`, `entryIsWellFormed_()`, `getBaseLossFormula_()`, `parseEntry_()`, `readFromJSON_()`, `readFromFile_()`, `parseRow_()` — all moved to providers
3. Remove now-unnecessary includes: `nlohmann/json.hpp`, `QtCore/QFile`, `QtCore/QTextStream`
4. Update default constructor to create providers:

```cpp
RibonucleotideDB::RibonucleotideDB() : max_code_length_(0)
{
  std::vector<std::unique_ptr<RibonucleotideDataProvider>> providers;
  providers.push_back(std::make_unique<ModomicsJSONDataProvider>("CHEMISTRY/Modomics.json"));
  providers.push_back(std::make_unique<RibonucleotideTSVDataProvider>("CHEMISTRY/Custom_RNA_modifications.tsv"));
  loadFromProviders_(providers);
}
```

5. Add provider-based constructor:

```cpp
RibonucleotideDB::RibonucleotideDB(std::vector<std::unique_ptr<RibonucleotideDataProvider>> providers)
  : max_code_length_(0)
{
  loadFromProviders_(providers);
}
```

6. Add `loadFromProviders_()`:

```cpp
void RibonucleotideDB::loadFromProviders_(std::vector<std::unique_ptr<RibonucleotideDataProvider>>& providers)
{
  // Collect deferred ambiguity entries (need all ribonucleotides loaded first)
  std::vector<std::tuple<std::string, String, String>> deferred_ambiguities;

  for (auto& provider : providers)
  {
    auto entries = provider->loadRibonucleotides();
    for (auto& entry : entries)
    {
      if (entry.isAmbiguous())
      {
        // Defer ambiguity resolution until all ribonucleotides are indexed
        deferred_ambiguities.emplace_back(
          entry.ribonucleotide->getCode(),
          entry.alternative_code_1,
          entry.alternative_code_2);
      }
      code_map_[entry.ribonucleotide->getCode()] = ribonucleotides_.size();
      max_code_length_ = std::max(max_code_length_, entry.ribonucleotide->getCode().size());
      ribonucleotides_.push_back(std::move(entry.ribonucleotide));
    }
  }

  // Now resolve ambiguities (all codes are indexed)
  for (const auto& [code, alt1, alt2] : deferred_ambiguities)
  {
    try
    {
      ambiguity_map_[code] = std::make_pair(getRibonucleotide(alt1), getRibonucleotide(alt2));
    }
    catch (Exception::ElementNotFound& e)
    {
      OPENMS_LOG_ERROR << "Error resolving ambiguity for " << code << ": " << e.what() << std::endl;
    }
  }
}
```

**Important note on ambiguity resolution:** In the original `readFromJSON_()`, ambiguities are resolved inline (line 268) by calling `getRibonucleotide()` on already-loaded entries. This works because the Modomics JSON lists non-ambiguous entries before ambiguous ones. However, with the provider pattern, we defer all ambiguity resolution to after ALL providers have loaded, which is more robust. This is a behavioral improvement.

**Step 3: Commit**
```bash
git add src/openms/include/OpenMS/CHEMISTRY/RibonucleotideDB.h \
        src/openms/source/CHEMISTRY/RibonucleotideDB.cpp
git commit -m "refactor: add provider-based constructor to RibonucleotideDB, remove file I/O methods"
```

---

### Task 5: Update CMake and Tests

**Files:**
- Modify: `src/openms/include/OpenMS/CHEMISTRY/sources.cmake`
- Modify: `src/openms/source/CHEMISTRY/sources.cmake`
- Create: `src/tests/class_tests/openms/source/RibonucleotideDataProvider_test.cpp`
- Modify: `src/tests/class_tests/openms/executables.cmake`
- Modify: `src/tests/class_tests/openms/source/RibonucleotideDB_test.cpp` (remove stale NOT_TESTABLE sections)

**Step 1: Update header sources.cmake**

Add after `RibonucleotideDB.h` (line 40):
```cmake
RibonucleotideDataProvider.h
ModomicsJSONDataProvider.h
RibonucleotideTSVDataProvider.h
```

Actually, add alphabetically:
- `ModomicsJSONDataProvider.h` after `ModifiedPeptideGenerator.h` (line 23)
- `RibonucleotideDataProvider.h` after `RibonucleotideDB.h` (line 40)
- `RibonucleotideTSVDataProvider.h` after `RibonucleotideDataProvider.h`

**Step 2: Update source sources.cmake**

Add alphabetically:
- `ModomicsJSONDataProvider.cpp` after `ModifiedPeptideGenerator.cpp` (line 24)
- `RibonucleotideTSVDataProvider.cpp` after `RibonucleotideDB.cpp` (line 38)

**Step 3: Write unit tests**

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>
#include <OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h>
#include <OpenMS/CHEMISTRY/RibonucleotideTSVDataProvider.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>

using namespace OpenMS;
using namespace std;

START_TEST(RibonucleotideDataProvider, "$Id$")

START_SECTION(InMemoryRibonucleotideDataProvider)
{
  // Create a test ribonucleotide
  vector<RibonucleotideEntry> entries;
  RibonucleotideEntry e;
  e.ribonucleotide = make_unique<Ribonucleotide>();
  e.ribonucleotide->setName("test-adenosine");
  e.ribonucleotide->setCode("tA");
  e.ribonucleotide->setOrigin('A');
  entries.push_back(std::move(e));

  InMemoryRibonucleotideDataProvider provider(std::move(entries));
  auto result = provider.loadRibonucleotides();
  TEST_EQUAL(result.size(), 1)
  TEST_STRING_EQUAL(result[0].ribonucleotide->getCode(), "tA")
  TEST_STRING_EQUAL(result[0].ribonucleotide->getName(), "test-adenosine")
  TEST_EQUAL(result[0].isAmbiguous(), false)
}
END_SECTION

START_SECTION(ModomicsJSONDataProvider)
{
  ModomicsJSONDataProvider provider("CHEMISTRY/Modomics.json");
  auto entries = provider.loadRibonucleotides();
  TEST_EQUAL(entries.size() > 100, true) // Modomics has >100 entries

  // Verify a known entry exists
  bool found_am = false;
  for (const auto& e : entries)
  {
    if (e.ribonucleotide->getCode() == "Am")
    {
      found_am = true;
      TEST_STRING_EQUAL(e.ribonucleotide->getName(), "2'-O-methyladenosine")
      break;
    }
  }
  TEST_EQUAL(found_am, true)
}
END_SECTION

START_SECTION(RibonucleotideTSVDataProvider)
{
  RibonucleotideTSVDataProvider provider("CHEMISTRY/Custom_RNA_modifications.tsv");
  auto entries = provider.loadRibonucleotides();
  TEST_EQUAL(entries.size() > 0, true)

  // Verify a known custom entry exists
  bool found = false;
  for (const auto& e : entries)
  {
    if (e.ribonucleotide->getCode() == "msU?")
    {
      found = true;
      TEST_EQUAL(e.isAmbiguous(), true)
      break;
    }
  }
  TEST_EQUAL(found, true)
}
END_SECTION

START_SECTION([EXTRA] RibonucleotideDB with InMemoryRibonucleotideDataProvider)
{
  // Create minimal test data
  vector<RibonucleotideEntry> entries;

  RibonucleotideEntry e1;
  e1.ribonucleotide = make_unique<Ribonucleotide>();
  e1.ribonucleotide->setName("test-adenosine");
  e1.ribonucleotide->setCode("tA");
  e1.ribonucleotide->setOrigin('A');
  entries.push_back(std::move(e1));

  RibonucleotideEntry e2;
  e2.ribonucleotide = make_unique<Ribonucleotide>();
  e2.ribonucleotide->setName("test-cytidine");
  e2.ribonucleotide->setCode("tC");
  e2.ribonucleotide->setOrigin('C');
  entries.push_back(std::move(e2));

  vector<unique_ptr<RibonucleotideDataProvider>> providers;
  providers.push_back(make_unique<InMemoryRibonucleotideDataProvider>(std::move(entries)));

  RibonucleotideDB db(std::move(providers));

  // Test code lookup
  auto* r = db.getRibonucleotide("tA");
  TEST_STRING_EQUAL(r->getName(), "test-adenosine")

  auto* r2 = db.getRibonucleotide("tC");
  TEST_STRING_EQUAL(r2->getName(), "test-cytidine")

  // Test iterator
  size_t count = 0;
  for (auto it = db.begin(); it != db.end(); ++it) { ++count; }
  TEST_EQUAL(count, 2)

  // Test not found
  TEST_EXCEPTION(Exception::ElementNotFound, db.getRibonucleotide("nonexistent"))
}
END_SECTION

END_TEST
```

**Step 4: Update executables.cmake**

Add `RibonucleotideDataProvider_test` after `RibonucleotideDB_test` (line 435):

```cmake
  RibonucleotideDataProvider_test
```

**Step 5: Remove stale NOT_TESTABLE sections from RibonucleotideDB_test.cpp**

Remove sections for `readFromJSON_` and `readFromFile_` (lines 39-49) since these methods no longer exist.

**Step 6: Commit**
```bash
git add src/openms/include/OpenMS/CHEMISTRY/sources.cmake \
        src/openms/source/CHEMISTRY/sources.cmake \
        src/tests/class_tests/openms/source/RibonucleotideDataProvider_test.cpp \
        src/tests/class_tests/openms/executables.cmake \
        src/tests/class_tests/openms/source/RibonucleotideDB_test.cpp
git commit -m "test: add RibonucleotideDataProvider tests, update cmake, clean stale sections"
```

---

### Task 6: Build and Integration Verification

**Step 1: Build**
```bash
cmake --build OpenMS-build -j$(nproc)
```

**Step 2: Run all related tests**
```bash
ctest --test-dir OpenMS-build -R "RibonucleotideDB_test|RibonucleotideDataProvider_test|NASequence_test|RibonucleotideDB_test" --output-on-failure
```

Expected: All tests pass.

**Step 3: Run the full existing test to verify no regressions**
```bash
ctest --test-dir OpenMS-build -R "RibonucleotideDB_test" -V
```

Verify the existing singleton tests still pass — `getInstance()`, `begin()`, `getRibonucleotide("Am")`, `getRibonucleotide("msU?")`, `getRibonucleotideAlternatives("msU?")`, `getRibonucleotidePrefix("m1AmCGU")`, `getBaselossFormula()`.

**Step 4: Commit verification results**

No commit needed — this is verification only.

---

### Task 7: Final Cleanup

**Step 1: Verify pyOpenMS bindings**

Check `src/pyOpenMS/pxds/RibonucleotideDB.pxd` — it only exposes `getInstance()`, `getRibonucleotide()`, `getRibonucleotidePrefix()`, and `getRibonucleotideAlternatives()`. The removed methods (`readFromFile_`, `readFromJSON_`, `parseRow_`) are protected and NOT exposed in the `.pxd`. No pyOpenMS changes needed.

**Step 2: Verify no remaining references to removed methods**

```bash
grep -r "readFromJSON_\|readFromFile_\|parseRow_" src/openms/ --include="*.h" --include="*.cpp"
```

Expected: No results (all moved to providers or removed).
