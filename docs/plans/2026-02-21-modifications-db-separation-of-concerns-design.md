# ModificationsDB Separation of Concerns Design

## Problem

ModificationsDB (and its subclass CrossLinksDB) load data files directly in their constructors, mixing I/O with domain logic. This makes unit testing slow (requires real files), prevents injecting alternative data sources, and couples the database to specific file formats.

## Approach: Data Provider Interface

Introduce a `ModificationDataProvider` interface that abstracts data sourcing. ModificationsDB receives providers via constructor injection. The existing `getInstance()` API is preserved by internally creating file-based providers.

## Core Interface

```cpp
class OPENMS_DLLAPI ModificationDataProvider
{
public:
  virtual ~ModificationDataProvider() = default;
  virtual std::vector<std::unique_ptr<ResidueModification>> loadModifications() = 0;
};
```

Returns `vector<unique_ptr<ResidueModification>>` for clear ownership transfer.

## Concrete Providers

### UnimodXMLDataProvider

Wraps the existing `UnimodXMLFile` parser. Takes a filename, returns parsed modifications.

```cpp
class UnimodXMLDataProvider : public ModificationDataProvider {
public:
  explicit UnimodXMLDataProvider(const String& filename);
  std::vector<std::unique_ptr<ResidueModification>> loadModifications() override;
private:
  String filename_;
};
```

### OBODataProvider

Wraps the OBO parsing logic currently in `ModificationsDB::readFromOBOFile()`. A `cross_links_only` flag controls whether it returns only cross-linkers (for CrossLinksDB) or excludes them (for ModificationsDB).

```cpp
class OBODataProvider : public ModificationDataProvider {
public:
  explicit OBODataProvider(const String& filename, bool cross_links_only = false);
  std::vector<std::unique_ptr<ResidueModification>> loadModifications() override;
private:
  String filename_;
  bool cross_links_only_;
};
```

### InMemoryDataProvider

For testing. Accepts pre-built modifications.

```cpp
class InMemoryDataProvider : public ModificationDataProvider {
public:
  explicit InMemoryDataProvider(std::vector<std::unique_ptr<ResidueModification>> mods);
  std::vector<std::unique_ptr<ResidueModification>> loadModifications() override;
private:
  std::vector<std::unique_ptr<ResidueModification>> mods_;
};
```

## ModificationsDB Changes

### Single constructor path: providers only

```cpp
class ModificationsDB {
public:
  // Construct from providers (no file I/O in this constructor)
  explicit ModificationsDB(std::vector<std::unique_ptr<ModificationDataProvider>> providers);

  // Existing singleton API unchanged
  static ModificationsDB* getInstance();
  static ModificationsDB* initializeModificationsDB(
    String unimod_file = "", String custommod_file = "",
    String psimod_file = "", String xlmod_file = "");
  static bool isInstantiated();

  // All query methods unchanged
  // ...

private:
  void loadFromProviders_(std::vector<std::unique_ptr<ModificationDataProvider>>& providers);
  // Indexing logic stays here (registering lookup keys in modification_names_)
};
```

### initializeModificationsDB creates providers internally

```cpp
ModificationsDB* ModificationsDB::initializeModificationsDB(...) {
  std::vector<std::unique_ptr<ModificationDataProvider>> providers;
  if (!unimod_file.empty())
    providers.push_back(make_unique<UnimodXMLDataProvider>(unimod_file));
  if (!custommod_file.empty())
    providers.push_back(make_unique<UnimodXMLDataProvider>(custommod_file));
  if (!psimod_file.empty())
    providers.push_back(make_unique<OBODataProvider>(psimod_file));
  if (!xlmod_file.empty())
    providers.push_back(make_unique<OBODataProvider>(xlmod_file));
  return new ModificationsDB(std::move(providers));
}
```

### Test usage

```cpp
std::vector<std::unique_ptr<ModificationDataProvider>> providers;
providers.push_back(make_unique<InMemoryDataProvider>(my_test_mods));
ModificationsDB db(std::move(providers));
// Fast, isolated, no files needed
```

## CrossLinksDB Adaptation

CrossLinksDB creates its own provider set via a helper, no more clear-then-reload:

```cpp
CrossLinksDB::CrossLinksDB() :
    ModificationsDB(makeCrossLinkProviders_())
{
}

std::vector<std::unique_ptr<ModificationDataProvider>>
CrossLinksDB::makeCrossLinkProviders_() {
    std::vector<std::unique_ptr<ModificationDataProvider>> providers;
    providers.push_back(
        make_unique<OBODataProvider>("CHEMISTRY/XLMOD.obo", true));
    return providers;
}
```

## File Organization

### New files
- `src/openms/include/OpenMS/CHEMISTRY/ModificationDataProvider.h` -- interface + InMemoryDataProvider
- `src/openms/include/OpenMS/CHEMISTRY/UnimodXMLDataProvider.h`
- `src/openms/source/CHEMISTRY/UnimodXMLDataProvider.cpp`
- `src/openms/include/OpenMS/CHEMISTRY/OBODataProvider.h`
- `src/openms/source/CHEMISTRY/OBODataProvider.cpp`

### Modified files
- `ModificationsDB.h` / `.cpp` -- replace file-path constructor with provider constructor
- `CrossLinksDB.h` / `.cpp` -- adapt to provider-based parent constructor
- `CMakeLists.txt` -- register new source files

### What moves where
- ~400 lines of OBO parsing from `ModificationsDB::readFromOBOFile()` -> `OBODataProvider::loadModifications()`
- UnimodXMLFile delegation from `readFromUnimodXMLFile()` -> `UnimodXMLDataProvider::loadModifications()`
- Indexing logic (registering lookup keys in `modification_names_`) stays in `ModificationsDB::loadFromProviders_()`

## Backwards Compatibility

- `getInstance()` works identically (creates file providers internally)
- `initializeModificationsDB()` works identically (creates file providers from paths)
- All 126 call sites using `getInstance()` require zero changes
- `readFromUnimodXMLFile()` and `readFromOBOFile()` are deprecated but kept as thin wrappers
