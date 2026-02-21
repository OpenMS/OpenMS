# ModificationsDB Separation of Concerns Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Decouple file I/O from ModificationsDB/CrossLinksDB by introducing a ModificationDataProvider interface, enabling dependency injection for testability, flexibility, and lazy loading.

**Architecture:** A polymorphic `ModificationDataProvider` interface abstracts data sourcing. Concrete providers handle file parsing (UnimodXML, OBO) or in-memory data (testing). ModificationsDB receives providers via constructor injection. The existing `getInstance()` singleton API is preserved by internally creating file-based providers. CrossLinksDB passes its own OBO provider to the parent constructor.

**Tech Stack:** C++17, OpenMS build system (CMake with sources.cmake), OpenMS test framework (ClassTest.h)

---

### Task 1: Create ModificationDataProvider interface header

**Files:**
- Create: `src/openms/include/OpenMS/CHEMISTRY/ModificationDataProvider.h`
- Modify: `src/openms/include/OpenMS/CHEMISTRY/sources.cmake`

**Step 1: Write the interface header**

```cpp
// src/openms/include/OpenMS/CHEMISTRY/ModificationDataProvider.h
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/OpenMSConfig.h>

#include <memory>
#include <vector>

namespace OpenMS
{
  /**
    @brief Interface for providing ResidueModification data to ModificationsDB.

    Implementations of this interface abstract the source of modification data,
    enabling dependency injection. File-based providers (UnimodXMLDataProvider,
    OBODataProvider) handle I/O; InMemoryDataProvider supports testing.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ModificationDataProvider
  {
  public:
    virtual ~ModificationDataProvider() = default;

    /**
      @brief Load modifications from whatever source this provider wraps.
      @return Vector of modifications with ownership transferred to caller.
    */
    virtual std::vector<std::unique_ptr<ResidueModification>> loadModifications() = 0;
  };

  /**
    @brief Data provider that serves pre-built modifications from memory.

    Useful for unit testing ModificationsDB without requiring files on disk.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI InMemoryDataProvider : public ModificationDataProvider
  {
  public:
    explicit InMemoryDataProvider(std::vector<std::unique_ptr<ResidueModification>> mods)
      : mods_(std::move(mods))
    {
    }

    std::vector<std::unique_ptr<ResidueModification>> loadModifications() override
    {
      return std::move(mods_);
    }

  private:
    std::vector<std::unique_ptr<ResidueModification>> mods_;
  };

} // namespace OpenMS
```

**Step 2: Register header in sources.cmake**

In `src/openms/include/OpenMS/CHEMISTRY/sources.cmake`, add `ModificationDataProvider.h` to `sources_list_h` (alphabetically after `ModificationsDB.h`):

```cmake
ModificationsDB.h
ModificationDataProvider.h
ModifiedNASequenceGenerator.h
```

**Step 3: Build to verify header compiles**

Run: `cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -20`
Expected: Build succeeds (header-only, no new .cpp yet)

**Step 4: Commit**

```bash
git add src/openms/include/OpenMS/CHEMISTRY/ModificationDataProvider.h \
        src/openms/include/OpenMS/CHEMISTRY/sources.cmake
git commit -m "feat: add ModificationDataProvider interface and InMemoryDataProvider"
```

---

### Task 2: Create UnimodXMLDataProvider

**Files:**
- Create: `src/openms/include/OpenMS/CHEMISTRY/UnimodXMLDataProvider.h`
- Create: `src/openms/source/CHEMISTRY/UnimodXMLDataProvider.cpp`
- Modify: `src/openms/include/OpenMS/CHEMISTRY/sources.cmake`
- Modify: `src/openms/source/CHEMISTRY/sources.cmake`

**Step 1: Write the header**

```cpp
// src/openms/include/OpenMS/CHEMISTRY/UnimodXMLDataProvider.h
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
    @brief Loads ResidueModification data from a Unimod XML file.

    Wraps UnimodXMLFile::load() behind the ModificationDataProvider interface.
    Each modification gets its FullId set after parsing.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI UnimodXMLDataProvider : public ModificationDataProvider
  {
  public:
    /// @param filename Path to a Unimod XML file (resolved via File::find)
    explicit UnimodXMLDataProvider(const String& filename);

    std::vector<std::unique_ptr<ResidueModification>> loadModifications() override;

  private:
    String filename_;
  };

} // namespace OpenMS
```

**Step 2: Write the implementation**

Extract the logic from `ModificationsDB::readFromUnimodXMLFile()` (ModificationsDB.cpp:470-493).
The provider does parsing + setFullId(). It does NOT do indexing (no modification_names_ access).

```cpp
// src/openms/source/CHEMISTRY/UnimodXMLDataProvider.cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/UnimodXMLDataProvider.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>

using namespace std;

namespace OpenMS
{

  UnimodXMLDataProvider::UnimodXMLDataProvider(const String& filename)
    : filename_(filename)
  {
  }

  std::vector<std::unique_ptr<ResidueModification>> UnimodXMLDataProvider::loadModifications()
  {
    // UnimodXMLFile::load returns raw pointers; wrap them in unique_ptr
    vector<ResidueModification*> raw_mods;
    UnimodXMLFile().load(filename_, raw_mods);

    vector<unique_ptr<ResidueModification>> result;
    result.reserve(raw_mods.size());
    for (auto* m : raw_mods)
    {
      m->setFullId();
      result.push_back(unique_ptr<ResidueModification>(m));
    }
    return result;
  }

} // namespace OpenMS
```

**Step 3: Register in both sources.cmake files**

In `src/openms/include/OpenMS/CHEMISTRY/sources.cmake`, add `UnimodXMLDataProvider.h` to `sources_list_h` (alphabetically near end, before `Tagger.h`):

```cmake
TheoreticalSpectrumGenerator.h
UnimodXMLDataProvider.h
TheoreticalSpectrumGeneratorXLMS.h
```

Wait — alphabetical order places it between `Tagger.h` and `TheoreticalSpectrumGenerator.h`. Check actual alphabetical position. `U` comes after `T`, so add after `TheoreticalSpectrumGeneratorXLMS.h`:

```cmake
TheoreticalSpectrumGeneratorXLMS.h
UnimodXMLDataProvider.h
)
```

In `src/openms/source/CHEMISTRY/sources.cmake`, add `UnimodXMLDataProvider.cpp` to `sources_list`:

```cmake
TheoreticalSpectrumGeneratorXLMS.cpp
UnimodXMLDataProvider.cpp
)
```

**Step 4: Build to verify**

Run: `cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -20`
Expected: Build succeeds

**Step 5: Commit**

```bash
git add src/openms/include/OpenMS/CHEMISTRY/UnimodXMLDataProvider.h \
        src/openms/source/CHEMISTRY/UnimodXMLDataProvider.cpp \
        src/openms/include/OpenMS/CHEMISTRY/sources.cmake \
        src/openms/source/CHEMISTRY/sources.cmake
git commit -m "feat: add UnimodXMLDataProvider for file-based UniMod loading"
```

---

### Task 3: Create OBODataProvider

**Files:**
- Create: `src/openms/include/OpenMS/CHEMISTRY/OBODataProvider.h`
- Create: `src/openms/source/CHEMISTRY/OBODataProvider.cpp`
- Modify: `src/openms/include/OpenMS/CHEMISTRY/sources.cmake`
- Modify: `src/openms/source/CHEMISTRY/sources.cmake`

**Step 1: Write the header**

```cpp
// src/openms/include/OpenMS/CHEMISTRY/OBODataProvider.h
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
    @brief Loads ResidueModification data from an OBO ontology file (PSI-MOD or XLMOD).

    Parses OBO format files and returns modifications as unique_ptr.
    The `cross_links_only` flag controls filtering:
    - false (default): returns all modifications EXCEPT cross-linkers (reactionSites==2).
      Used by ModificationsDB for PSI-MOD.obo and XLMOD.obo (mono-links only).
    - true: returns only cross-linkers (reactionSites==2), skipping mono-links.
      Used by CrossLinksDB for XLMOD.obo.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI OBODataProvider : public ModificationDataProvider
  {
  public:
    /**
      @param filename Path to an OBO file (resolved via File::find)
      @param cross_links_only If true, return only cross-linkers; if false, exclude them
    */
    explicit OBODataProvider(const String& filename, bool cross_links_only = false);

    std::vector<std::unique_ptr<ResidueModification>> loadModifications() override;

  private:
    String filename_;
    bool cross_links_only_;
  };

} // namespace OpenMS
```

**Step 2: Write the implementation**

Move the OBO parsing logic from `ModificationsDB::readFromOBOFile()` (ModificationsDB.cpp:559-858) and
`CrossLinksDB::readFromOBOFile()` (CrossLinksDB.cpp:37-333) into this provider.

The two existing implementations differ only in their filtering logic:
- ModificationsDB: `reading_cross_link = true` when reactionSites==2, skips those
- CrossLinksDB: `reading_mono_link = true` when reactionSites==1, skips those

Unify into one implementation with the `cross_links_only_` flag.

The provider handles:
1. File parsing (lines 559-811 of ModificationsDB.cpp) — produces `multimap<String, ResidueModification>`
2. Modification creation and synonym setup (lines 813-858) — converts parsed entries into owned ResidueModification objects

The provider does NOT handle:
- Alias resolution for existing UniMod entries (this requires access to `modification_names_` in ModificationsDB)
- Registering lookup keys in `modification_names_`

For modifications with `getUniModRecordId() > 0`, the provider still returns them. ModificationsDB's `loadFromProviders_()` will detect these and handle the alias-to-existing-UniMod logic (adding PSI-MOD accession as an alias for the existing UniMod entry, then discarding the duplicate mod).

```cpp
// src/openms/source/CHEMISTRY/OBODataProvider.cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/OBODataProvider.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

namespace OpenMS
{

  OBODataProvider::OBODataProvider(const String& filename, bool cross_links_only)
    : filename_(filename),
      cross_links_only_(cross_links_only)
  {
  }

  std::vector<std::unique_ptr<ResidueModification>> OBODataProvider::loadModifications()
  {
    ResidueModification mod;
    multimap<String, ResidueModification> all_mods;

    ifstream is(File::find(filename_).c_str());
    String line, line_wo_spaces, id;
    String origin = "";

    bool skip_current_term = false;

    // Phase 1: Parse the OBO file into all_mods multimap
    while (getline(is, line, '\n'))
    {
      line.trim();
      line_wo_spaces = line;
      line_wo_spaces.removeWhitespaces();

      if (line.empty() || line[0] == '!')
      {
        continue;
      }

      if (line_wo_spaces == "[Term]")
      {
        if (!id.empty() && !skip_current_term)
        {
          vector<String> origins;
          origin.split(",", origins);

          std::sort(origins.begin(), origins.end());
          auto unique_end = unique(origins.begin(), origins.end());
          origins.resize(distance(origins.begin(), unique_end));

          for (auto orig_it = origins.begin(); orig_it != origins.end(); ++orig_it)
          {
            if ((orig_it->size() == 1) && (*orig_it != "B") && (*orig_it != "J") && (*orig_it != "Z"))
            {
              mod.setOrigin((*orig_it)[0]);
              all_mods.insert(make_pair(id, mod));
            }
          }

          if (origin.hasSubstring("ProteinN-term"))
          {
            mod.setTermSpecificity(cross_links_only_ ? ResidueModification::N_TERM : ResidueModification::PROTEIN_N_TERM);
            mod.setOrigin('X');
            all_mods.insert(make_pair(id, mod));
          }
          if (origin.hasSubstring("ProteinC-term"))
          {
            mod.setTermSpecificity(cross_links_only_ ? ResidueModification::C_TERM : ResidueModification::PROTEIN_C_TERM);
            mod.setOrigin('X');
            all_mods.insert(make_pair(id, mod));
          }

          id = "";
          origin = "";
          mod = ResidueModification();
        }
        else if (skip_current_term)
        {
          id = "";
          origin = "";
          mod = ResidueModification();
          skip_current_term = false;
        }
      }
      else if (line_wo_spaces.hasPrefix("id:"))
      {
        id = line.substr(line.find(':') + 1).trim();
        mod.setId(id);
        mod.setPSIMODAccession(id);
      }
      else if (line_wo_spaces.hasPrefix("name:"))
      {
        String name = line.substr(line.find(':') + 1).trim();
        mod.setFullName(name);
        if (mod.getId().hasSubstring("XLMOD"))
        {
          mod.setName(name);
          mod.setId(name);
          mod.setFullName(name);
        }
      }
      else if (line_wo_spaces.hasPrefix("is_a:"))
      {
        // TODO
      }
      else if (line_wo_spaces.hasPrefix("def:"))
      {
        line.remove('[');
        line.remove(']');
        line.remove(',');
        vector<String> split;
        line.split(' ', split);
        for (Size i = 0; i != split.size(); ++i)
        {
          if (split[i].hasPrefix("UniMod:"))
          {
            String identifier = split[i].substr(7, split[i].size());
            mod.setUniModRecordId(identifier.toInt());
          }
        }
      }
      else if (line_wo_spaces.hasPrefix("comment:"))
      {
        // TODO
      }
      else if (line_wo_spaces.hasPrefix("synonym:"))
      {
        vector<String> val_split;
        line.split('"', val_split);
        if (val_split.size() < 3)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "missing \" characters to enclose argument!");
        }
        mod.addSynonym(val_split[1]);

        if (line_wo_spaces.hasSubstring("PSI-MOD-label"))
        {
          mod.setName(val_split[1]);
        }
      }
      else if (line_wo_spaces.hasPrefix("property_value:"))
      {
        String val = line_wo_spaces.substr(15, line_wo_spaces.size() - 15);
        val.trim();

        if (val.hasSubstring("\"none\""))
        {
          continue;
        }

        vector<String> val_split;
        val.split('"', val_split);
        if (val_split.size() != 3)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "missing \" characters to enclose argument!");
        }
        if (val.hasPrefix("DiffAvg:"))
        {
          mod.setDiffAverageMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("DiffFormula:"))
        {
          vector<String> tmp_split;
          line.split('"', tmp_split);
          tmp_split[1].removeWhitespaces();
          mod.setDiffFormula(EmpiricalFormula(tmp_split[1]));
        }
        else if (val.hasPrefix("DiffMono:"))
        {
          mod.setDiffMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("Formula:"))
        {
          mod.setFormula(val_split[1]);
        }
        else if (val.hasPrefix("MassAvg:"))
        {
          mod.setAverageMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("MassMono:"))
        {
          mod.setMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("Origin:"))
        {
          origin = val_split[1];
        }
        else if (val.hasPrefix("Source:"))
        {
          mod.setSourceClassification(val_split[1]);
        }
        else if (val.hasPrefix("TermSpec:"))
        {
          mod.setTermSpecificity(val_split[1]);
        }
        else if (val.hasPrefix("reactionSites:"))
        {
          if (val_split[1] == "2" && !cross_links_only_)
          {
            // Cross-linker: skip when loading for ModificationsDB
            skip_current_term = true;
          }
          else if (val_split[1] == "1" && cross_links_only_)
          {
            // Mono-link: skip when loading for CrossLinksDB
            skip_current_term = true;
          }
        }
        else if (val.hasPrefix("monoisotopicMass:"))
        {
          mod.setDiffMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("specificities:"))
        {
          origin = val_split[1];
          origin.remove('(');
          origin.remove(')');
          origin.substitute("&", ",");
        }
      }
    }

    // Handle last term
    if (!id.empty() && !skip_current_term)
    {
      vector<String> origins;
      origin.split(",", origins);

      std::sort(origins.begin(), origins.end());
      auto unique_end = unique(origins.begin(), origins.end());
      origins.resize(distance(origins.begin(), unique_end));

      for (auto orig_it = origins.begin(); orig_it != origins.end(); ++orig_it)
      {
        if ((orig_it->size() == 1) && (*orig_it != "B") && (*orig_it != "J") && (*orig_it != "Z"))
        {
          mod.setOrigin((*orig_it)[0]);
          all_mods.insert(make_pair(id, mod));
        }
      }

      if (origin.hasSubstring("ProteinN-term"))
      {
        mod.setTermSpecificity(cross_links_only_ ? ResidueModification::N_TERM : ResidueModification::PROTEIN_N_TERM);
        mod.setOrigin('X');
        all_mods.insert(make_pair(id, mod));
      }
      if (origin.hasSubstring("ProteinC-term"))
      {
        mod.setTermSpecificity(cross_links_only_ ? ResidueModification::C_TERM : ResidueModification::PROTEIN_C_TERM);
        mod.setOrigin('X');
        all_mods.insert(make_pair(id, mod));
      }
    }

    // Phase 2: Convert parsed entries into owned ResidueModification objects.
    // Mods with UniMod record IDs are still returned — ModificationsDB handles
    // alias resolution (mapping PSI-MOD accessions to existing UniMod entries).
    vector<unique_ptr<ResidueModification>> result;

    for (auto it = all_mods.begin(); it != all_mods.end(); ++it)
    {
      if (it->second.getUniModRecordId() > 0)
      {
        // Return a minimal mod carrying the PSI-MOD <-> UniMod mapping.
        // ModificationsDB will use this for alias resolution only.
        auto alias_mod = make_unique<ResidueModification>(it->second);
        result.push_back(std::move(alias_mod));
      }
      else
      {
        if ((it->second.getOrigin() != 'X') ||
           ((it->second.getTermSpecificity() != ResidueModification::ANYWHERE) &&
            (it->second.getDiffMonoMass() != 0)))
        {
          auto new_mod = make_unique<ResidueModification>(it->second);

          // Collect all synonyms for this modification
          set<String> synonyms = it->second.getSynonyms();
          synonyms.insert(it->first);
          synonyms.insert(it->second.getFullName());
          synonyms.insert(it->second.getPSIMODAccession());

          // Set the full ID (based on name, not short ID)
          new_mod->setId(it->second.getFullName());
          new_mod->setFullId();
          new_mod->setId(it->second.getId());

          // Re-add synonyms that were collected plus the generated fullId
          for (const auto& syn : synonyms)
          {
            new_mod->addSynonym(syn);
          }
          // Also add the generated fullId as synonym for indexing
          new_mod->addSynonym(new_mod->getFullId());

          result.push_back(std::move(new_mod));
        }
      }
    }

    return result;
  }

} // namespace OpenMS
```

**Important note for implementer:** The OBO parsing phase has duplicated code for handling the last term (after the while loop exits). This is intentional — it matches the existing code structure. A future cleanup could extract this into a helper lambda, but that is out of scope.

**Step 3: Register in sources.cmake files**

In `src/openms/include/OpenMS/CHEMISTRY/sources.cmake`, add `OBODataProvider.h` to `sources_list_h` (alphabetically, between `NASequence.h` and `NucleicAcidSpectrumGenerator.h`):

Wait — `O` comes after `N`, so add between `NASequence.h` and `NucleicAcidSpectrumGenerator.h`:

```cmake
NASequence.h
OBODataProvider.h
NucleicAcidSpectrumGenerator.h
```

In `src/openms/source/CHEMISTRY/sources.cmake`, add `OBODataProvider.cpp`:

```cmake
NASequence.cpp
OBODataProvider.cpp
NucleicAcidSpectrumGenerator.cpp
```

**Step 4: Build to verify**

Run: `cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -20`
Expected: Build succeeds

**Step 5: Commit**

```bash
git add src/openms/include/OpenMS/CHEMISTRY/OBODataProvider.h \
        src/openms/source/CHEMISTRY/OBODataProvider.cpp \
        src/openms/include/OpenMS/CHEMISTRY/sources.cmake \
        src/openms/source/CHEMISTRY/sources.cmake
git commit -m "feat: add OBODataProvider for OBO-format modification loading"
```

---

### Task 4: Add provider-based constructor to ModificationsDB

**Files:**
- Modify: `src/openms/include/OpenMS/CHEMISTRY/ModificationsDB.h`
- Modify: `src/openms/source/CHEMISTRY/ModificationsDB.cpp`

**Step 1: Modify the header**

Add a forward declaration and new constructor/method to `ModificationsDB.h`.

Changes to make:

a) Add `#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>` after the existing includes (line 16).

b) Add the provider-based constructor in the `public` section after `isInstantiated()` (after line 60):

```cpp
    /**
      @brief Construct from data providers (no file I/O performed by this constructor).

      Use this constructor for dependency injection, e.g., with InMemoryDataProvider
      for testing or custom providers for alternative data sources.

      @param providers Vector of data providers; each will be called once to load modifications.
    */
    explicit ModificationsDB(std::vector<std::unique_ptr<ModificationDataProvider>> providers);
```

c) Add `loadFromProviders_` as a private helper (after `readFromUnimodXMLFile` declaration, around line 282):

```cpp
    /// Loads and indexes modifications from the given providers
    void loadFromProviders_(std::vector<std::unique_ptr<ModificationDataProvider>>& providers);
```

**Step 2: Implement the provider constructor and loadFromProviders_**

In `ModificationsDB.cpp`:

a) Add include at top:

```cpp
#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>
```

b) Add the new constructor (after the existing file-path constructor, around line 96):

```cpp
ModificationsDB::ModificationsDB(std::vector<std::unique_ptr<ModificationDataProvider>> providers)
{
  loadFromProviders_(providers);
  is_instantiated_ = true;
}
```

c) Add `loadFromProviders_` implementation. This method calls each provider, then indexes the returned modifications. For OBO-sourced mods with UniMod record IDs, it performs alias resolution against existing entries:

```cpp
void ModificationsDB::loadFromProviders_(std::vector<std::unique_ptr<ModificationDataProvider>>& providers)
{
  for (auto& provider : providers)
  {
    auto mods = provider->loadModifications();
    for (auto& m : mods)
    {
      // Check if this is an OBO alias mapping (has UniMod record ID and PSI-MOD accession).
      // These mods map PSI-MOD accessions to existing UniMod entries.
      if (m->getUniModRecordId() > 0 && !m->getPSIMODAccession().empty())
      {
        #pragma omp critical(OpenMS_ModificationsDB)
        {
          auto existing = modification_names_.find(m->getUniModAccession());
          if (existing != modification_names_.end())
          {
            // Add PSI-MOD accession as alias for existing UniMod entries
            for (const auto* existing_mod : existing->second)
            {
              modification_names_[m->getPSIMODAccession()].insert(existing_mod);
            }
            // Don't add this mod to mods_ — it's just an alias mapping
            continue;
          }
        }
      }

      // Regular modification: add to database and register lookup keys
      #pragma omp critical(OpenMS_ModificationsDB)
      {
        ResidueModification* raw = m.release();
        mods_.push_back(raw);

        // Register standard lookup keys
        modification_names_[raw->getFullId()].insert(raw);
        modification_names_[raw->getId()].insert(raw);
        modification_names_[raw->getFullName()].insert(raw);
        if (!raw->getUniModAccession().empty())
        {
          modification_names_[raw->getUniModAccession()].insert(raw);
        }

        // Register additional synonyms (OBO providers attach these)
        for (const auto& syn : raw->getSynonyms())
        {
          modification_names_[syn].insert(raw);
        }
      }
    }
  }
}
```

d) Modify `initializeModificationsDB` (line 62-72) to create providers and use the new constructor:

```cpp
ModificationsDB* ModificationsDB::initializeModificationsDB(OpenMS::String unimod_file, OpenMS::String custommod_file, OpenMS::String psimod_file, OpenMS::String xlmod_file)
{
  std::vector<std::unique_ptr<ModificationDataProvider>> providers;
  if (!unimod_file.empty())
  {
    providers.push_back(std::make_unique<UnimodXMLDataProvider>(std::move(unimod_file)));
  }
  if (!custommod_file.empty())
  {
    providers.push_back(std::make_unique<UnimodXMLDataProvider>(std::move(custommod_file)));
  }
  if (!psimod_file.empty())
  {
    providers.push_back(std::make_unique<OBODataProvider>(std::move(psimod_file)));
  }
  if (!xlmod_file.empty())
  {
    providers.push_back(std::make_unique<OBODataProvider>(std::move(xlmod_file)));
  }

  static ModificationsDB* db_ = new ModificationsDB(std::move(providers));
  return db_;
}
```

Add the required includes at top of ModificationsDB.cpp:

```cpp
#include <OpenMS/CHEMISTRY/UnimodXMLDataProvider.h>
#include <OpenMS/CHEMISTRY/OBODataProvider.h>
```

e) Remove the old file-path constructor body (lines 74-96), replacing it with a delegating call using providers. Since the old constructor is private and only called by `initializeModificationsDB` (which we just changed to use providers), the old constructor is now dead code. Remove it from both the header and implementation:

In `ModificationsDB.h`, remove line 251:
```cpp
explicit ModificationsDB(const OpenMS::String& unimod_file = "CHEMISTRY/unimod.xml", ...);
```

In `ModificationsDB.cpp`, remove lines 74-96 (the old constructor body).

f) Mark `readFromUnimodXMLFile` and `readFromOBOFile` as deprecated but keep the implementations for now (they are still used by CrossLinksDB until Task 5). No changes needed to these methods yet.

**Step 3: Build and run existing ModificationsDB test**

Run: `cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -30`
Expected: Build succeeds

Run: `ctest --test-dir OpenMS-build -R ModificationsDB_test --output-on-failure`
Expected: All existing tests pass (getInstance() path unchanged)

**Step 4: Commit**

```bash
git add src/openms/include/OpenMS/CHEMISTRY/ModificationsDB.h \
        src/openms/source/CHEMISTRY/ModificationsDB.cpp
git commit -m "feat: add provider-based constructor to ModificationsDB

initializeModificationsDB now creates providers internally.
Removes the file-path constructor; all loading goes through providers."
```

---

### Task 5: Adapt CrossLinksDB to use providers

**Files:**
- Modify: `src/openms/include/OpenMS/CHEMISTRY/CrossLinksDB.h`
- Modify: `src/openms/source/CHEMISTRY/CrossLinksDB.cpp`

**Step 1: Modify CrossLinksDB header**

a) Add a private static helper method and remove `readFromOBOFile`:

Replace the current `readFromOBOFile` public declaration (line 34) and the private constructor (line 45) with:

```cpp
  private:
    CrossLinksDB();
    CrossLinksDB(const CrossLinksDB& residue_db);
    ~CrossLinksDB() override;
    CrossLinksDB& operator=(const CrossLinksDB& aa);

    /// Creates the OBO provider for cross-linker loading
    static std::vector<std::unique_ptr<ModificationDataProvider>> makeCrossLinkProviders_();
```

Keep `getAllSearchModifications` as public (it's an override).
Remove the `readFromOBOFile` public declaration — it moves into the provider.

b) Add `#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>` to the header.

**Step 2: Modify CrossLinksDB implementation**

Replace the constructor (lines 19-25) with:

```cpp
CrossLinksDB::CrossLinksDB()
  : ModificationsDB(makeCrossLinkProviders_())
{
}

std::vector<std::unique_ptr<ModificationDataProvider>> CrossLinksDB::makeCrossLinkProviders_()
{
  std::vector<std::unique_ptr<ModificationDataProvider>> providers;
  providers.push_back(std::make_unique<OBODataProvider>("CHEMISTRY/XLMOD.obo", /*cross_links_only=*/true));
  return providers;
}
```

Remove the `CrossLinksDB::readFromOBOFile` method entirely (~250 lines, 37-333). This logic is now in `OBODataProvider::loadModifications()`.

Update includes:

```cpp
#include <OpenMS/CHEMISTRY/CrossLinksDB.h>
#include <OpenMS/CHEMISTRY/OBODataProvider.h>
```

Remove includes that are no longer needed: `<OpenMS/SYSTEM/File.h>`, `<fstream>`.

The destructor stays as-is (the parent ModificationsDB destructor handles cleanup).

Actually, looking at CrossLinksDB's destructor — it manually deletes mods_ entries. But the parent destructor does the same thing. This means mods get double-deleted. Check whether the parent destructor is virtual (it is: `virtual ~ModificationsDB()`). If CrossLinksDB's destructor runs first and clears mods_, the parent destructor would iterate an empty vector, which is safe. Keep CrossLinksDB's destructor as-is to match existing behavior.

**Step 3: Build and run tests**

Run: `cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -30`
Expected: Build succeeds

Run: `ctest --test-dir OpenMS-build -R CrossLinksDB_test --output-on-failure`
Expected: All existing tests pass

Run: `ctest --test-dir OpenMS-build -R ModificationsDB_test --output-on-failure`
Expected: All existing tests still pass

**Step 4: Commit**

```bash
git add src/openms/include/OpenMS/CHEMISTRY/CrossLinksDB.h \
        src/openms/source/CHEMISTRY/CrossLinksDB.cpp
git commit -m "refactor: adapt CrossLinksDB to use OBODataProvider

CrossLinksDB now passes an OBODataProvider with cross_links_only=true
to the parent ModificationsDB constructor. Removes duplicated OBO
parsing code (~250 lines)."
```

---

### Task 6: Clean up deprecated methods in ModificationsDB

**Files:**
- Modify: `src/openms/source/CHEMISTRY/ModificationsDB.cpp`
- Modify: `src/openms/include/OpenMS/CHEMISTRY/ModificationsDB.h`

**Step 1: Remove readFromOBOFile and readFromUnimodXMLFile**

Now that both ModificationsDB and CrossLinksDB use providers, the old loading methods are no longer called by any internal code.

In `ModificationsDB.h`, remove the private declarations (lines 279-282):
```cpp
void readFromOBOFile(const String& filename);
void readFromUnimodXMLFile(const String& filename);
```

In `ModificationsDB.cpp`, remove:
- `readFromUnimodXMLFile` implementation (lines 470-493)
- `readFromOBOFile` implementation (lines 559-859)

Also remove includes that are now only needed for those methods:
- `#include <OpenMS/FORMAT/UnimodXMLFile.h>` (now in UnimodXMLDataProvider.cpp)
- `#include <fstream>` (now in OBODataProvider.cpp)

**Step 2: Build and run full test suite**

Run: `cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -30`
Expected: Build succeeds (no external code called these private methods)

Run: `ctest --test-dir OpenMS-build -R "ModificationsDB_test|CrossLinksDB_test" --output-on-failure`
Expected: All tests pass

**Step 3: Commit**

```bash
git add src/openms/include/OpenMS/CHEMISTRY/ModificationsDB.h \
        src/openms/source/CHEMISTRY/ModificationsDB.cpp
git commit -m "refactor: remove deprecated readFromOBOFile and readFromUnimodXMLFile

These methods are replaced by OBODataProvider and UnimodXMLDataProvider.
All loading now goes through the provider-based constructor."
```

---

### Task 7: Write provider unit tests

**Files:**
- Create: `src/tests/class_tests/openms/source/ModificationDataProvider_test.cpp`
- Modify: `src/tests/class_tests/openms/executables.cmake`

**Step 1: Write the test file**

```cpp
// src/tests/class_tests/openms/source/ModificationDataProvider_test.cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>
#include <OpenMS/CHEMISTRY/UnimodXMLDataProvider.h>
#include <OpenMS/CHEMISTRY/OBODataProvider.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

using namespace OpenMS;
using namespace std;

START_TEST(ModificationDataProvider, "$Id$")

/////////////////////////////////////////////////////////////
// InMemoryDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(InMemoryDataProvider::loadModifications())
{
  // Create test modifications
  vector<unique_ptr<ResidueModification>> test_mods;
  auto mod1 = make_unique<ResidueModification>();
  mod1->setId("TestMod1");
  mod1->setFullName("Test Modification 1");
  mod1->setOrigin('C');
  mod1->setDiffMonoMass(57.02);
  mod1->setFullId();
  test_mods.push_back(std::move(mod1));

  auto mod2 = make_unique<ResidueModification>();
  mod2->setId("TestMod2");
  mod2->setFullName("Test Modification 2");
  mod2->setOrigin('M');
  mod2->setDiffMonoMass(15.99);
  mod2->setFullId();
  test_mods.push_back(std::move(mod2));

  InMemoryDataProvider provider(std::move(test_mods));
  auto result = provider.loadModifications();

  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(result[0]->getId(), "TestMod1")
  TEST_EQUAL(result[1]->getId(), "TestMod2")
  TEST_REAL_SIMILAR(result[0]->getDiffMonoMass(), 57.02)
}
END_SECTION

/////////////////////////////////////////////////////////////
// UnimodXMLDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(UnimodXMLDataProvider::loadModifications())
{
  UnimodXMLDataProvider provider("CHEMISTRY/unimod.xml");
  auto mods = provider.loadModifications();

  // unimod.xml should contain many modifications
  TEST_EQUAL(mods.size() > 100, true)

  // Check that FullId is set on all modifications
  bool all_have_fullid = true;
  for (const auto& m : mods)
  {
    if (m->getFullId().empty())
    {
      all_have_fullid = false;
      break;
    }
  }
  TEST_EQUAL(all_have_fullid, true)

  // Check that a well-known modification exists
  bool found_oxidation = false;
  for (const auto& m : mods)
  {
    if (m->getFullId() == "Oxidation (M)")
    {
      found_oxidation = true;
      break;
    }
  }
  TEST_EQUAL(found_oxidation, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// OBODataProvider tests (non-crosslink mode)
/////////////////////////////////////////////////////////////

START_SECTION(OBODataProvider::loadModifications() non-crosslink)
{
  OBODataProvider provider("CHEMISTRY/PSI-MOD.obo", false);
  auto mods = provider.loadModifications();

  // PSI-MOD should contain modifications
  TEST_EQUAL(mods.size() > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// OBODataProvider tests (crosslink mode)
/////////////////////////////////////////////////////////////

START_SECTION(OBODataProvider::loadModifications() crosslink)
{
  OBODataProvider provider("CHEMISTRY/XLMOD.obo", true);
  auto mods = provider.loadModifications();

  // XLMOD in cross-link mode should contain cross-linkers
  TEST_EQUAL(mods.size() > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ModificationsDB with InMemoryDataProvider
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] ModificationsDB with InMemoryDataProvider)
{
  // Create a ModificationsDB from in-memory data (no files)
  vector<unique_ptr<ResidueModification>> test_mods;
  auto mod = make_unique<ResidueModification>();
  mod->setId("Phospho");
  mod->setFullName("Phosphorylation");
  mod->setOrigin('S');
  mod->setDiffMonoMass(79.966);
  mod->setFullId();
  test_mods.push_back(std::move(mod));

  vector<unique_ptr<ModificationDataProvider>> providers;
  providers.push_back(make_unique<InMemoryDataProvider>(std::move(test_mods)));

  ModificationsDB db(std::move(providers));

  TEST_EQUAL(db.getNumberOfModifications(), 1)
  TEST_EQUAL(db.has("Phospho (S)"), true)
  TEST_EQUAL(db.has("Phospho"), true)
  TEST_EQUAL(db.has("Nonexistent"), false)
  TEST_REAL_SIMILAR(db.getModification(0)->getDiffMonoMass(), 79.966)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
```

**Step 2: Register test in executables.cmake**

In `src/tests/class_tests/openms/executables.cmake`, add to `chemistry_executables_list` (alphabetically after `ModificationsDB_test`):

```cmake
  ModificationsDB_test
  ModificationDataProvider_test
  ModifiedNASequenceGenerator_test
```

**Step 3: Build and run the new test**

Run: `cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -30`
Expected: Build succeeds

Run: `ctest --test-dir OpenMS-build -R ModificationDataProvider_test --output-on-failure`
Expected: All tests pass

**Step 4: Run full existing test suite for regression**

Run: `ctest --test-dir OpenMS-build -R "ModificationsDB_test|CrossLinksDB_test" --output-on-failure`
Expected: All existing tests still pass

**Step 5: Commit**

```bash
git add src/tests/class_tests/openms/source/ModificationDataProvider_test.cpp \
        src/tests/class_tests/openms/executables.cmake
git commit -m "test: add unit tests for ModificationDataProvider implementations

Tests InMemoryDataProvider, UnimodXMLDataProvider, OBODataProvider,
and ModificationsDB construction from providers."
```

---

### Task 8: Final integration verification

**Files:** None (verification only)

**Step 1: Run full chemistry test suite**

Run: `ctest --test-dir OpenMS-build -R "chemistry" --output-on-failure`
Expected: All chemistry tests pass

**Step 2: Run broader test suite to catch indirect breakage**

Run: `ctest --test-dir OpenMS-build -j$(nproc) --output-on-failure 2>&1 | tail -50`
Expected: No new failures (check against baseline)

**Step 3: Verify no files reference removed methods**

Search for any remaining references to `readFromOBOFile` or `readFromUnimodXMLFile` outside of the provider code:

```bash
grep -rn "readFromOBOFile\|readFromUnimodXMLFile" src/ --include="*.cpp" --include="*.h" | grep -v "OBODataProvider\|UnimodXMLDataProvider\|_test.cpp"
```

Expected: No results (all references are in provider code or tests)

**Step 4: Commit if any fixes were needed**

Only commit if fixes were required during integration testing.
