# Remove Qt6::Core from libOpenMS — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove Qt6::Core from the OpenMS core library (`libOpenMS`) so it can be built and distributed without Qt.

**Architecture:** 7 incremental PRs (PR 5 split into 5a/5b), each independently compilable and testable. Leaf dependencies first (JSON, DateTime), then infrastructure (filesystem, process), then string API cleanup, then CMake cut. PRs 1-4 are mostly independent; 5a depends on all of 1-4; 6 depends on 5a+5b.

**Tech Stack:** C++20, std::filesystem, std::chrono, boost::process (Boost 1.85), nlohmann/json (already in-tree), sscanf/snprintf for date parsing

**Spec:** `docs/superpowers/specs/2026-03-18-remove-qt-from-core-design.md`

---

## Task 1: PR 1 — QJson → nlohmann/json (OMSFileLoad)

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/OMSFileLoad.h:15,179,186,202`
- Modify: `src/openms/source/FORMAT/OMSFileLoad.cpp:15-18,31-42,1104-1173`
- Test: `src/tests/class_tests/openms/source/OMSFile_test.cpp`

- [ ] **Step 1: Run existing OMSFile tests to establish baseline**

Run: `ctest --test-dir OpenMS-build -R OMSFile -V`
Expected: All pass

- [ ] **Step 2: Replace QJson includes and types in OMSFileLoad.h**

In `src/openms/include/OpenMS/FORMAT/OMSFileLoad.h`:

Replace:
```cpp
#include <QtCore/QJsonArray>
```
With:
```cpp
#include <nlohmann/json.hpp>
```

Replace the method signature (line 179):
```cpp
QJsonArray exportTableToJSON_(const QString& table, const QString& order_by);
```
With:
```cpp
nlohmann::ordered_json exportTableToJSON_(const String& table, const String& order_by);
```

Replace member (line 186):
```cpp
QString subquery_score_;
```
With:
```cpp
String subquery_score_;
```

Replace member (line 202):
```cpp
std::map<QString, QString> export_order_by_;
```
With:
```cpp
std::map<String, String> export_order_by_;
```

Remove any remaining `#include <QtCore/QString>` or `class QString;` forward declarations in this header.

- [ ] **Step 3: Replace QJson usage in OMSFileLoad.cpp**

In `src/openms/source/FORMAT/OMSFileLoad.cpp`:

Remove Qt includes (lines 15-18):
```cpp
#include <QtCore/QString>
#include <QtCore/QJsonDocument>
#include <QtCore/QJsonObject>
```

Add (if not already present):
```cpp
#include <nlohmann/json.hpp>
```

Replace the static `export_order_by_` initialization (lines 31-42) — change `QString` literals to `String` or `std::string` literals.

Rewrite `exportTableToJSON_()` (lines 1104-1140): Replace `QJsonObject`/`QJsonArray` with `nlohmann::ordered_json`. Use `ordered_json` to preserve insertion order matching Qt's output.

Rewrite `exportToJSON()` (lines 1143-1173): Replace `QJsonDocument(obj).toJson()` with `obj.dump(2)` for indented output.

Replace all `.toQString()` calls in this file with direct `String`/`std::string` usage.

- [ ] **Step 4: Build and run tests**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc) && ctest --test-dir OpenMS-build -R OMSFile -V`
Expected: All pass

- [ ] **Step 5: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/OMSFileLoad.h src/openms/source/FORMAT/OMSFileLoad.cpp
git commit -m "refactor: replace QJson with nlohmann/json in OMSFileLoad

Removes Qt6 JSON dependency from core library. Uses ordered_json
to preserve insertion-order key output matching Qt's default."
```

---

## Task 2: PR 2 — DateTime/Date internal rewrite

**Files:**
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h:21,200`
- Modify: `src/openms/source/DATASTRUCTURES/DateTime.cpp` (full rewrite of internals)
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/Date.h:18,45,123`
- Modify: `src/openms/source/DATASTRUCTURES/Date.cpp` (full rewrite of internals)
- Test: `src/tests/class_tests/openms/source/DateTime_test.cpp`
- Test: `src/tests/class_tests/openms/source/Date_test.cpp`

### Sub-task 2a: DateTime rewrite

- [ ] **Step 1: Run existing DateTime tests to establish baseline**

Run: `ctest --test-dir OpenMS-build -R DateTime_test -V`
Expected: All pass. Note the exact test output for round-trip verification.

- [ ] **Step 2: Replace DateTime internal representation in header**

In `src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h`:

Remove forward declaration (line 21):
```cpp
class QDateTime;
```

Replace the private member (line 200):
```cpp
std::unique_ptr<QDateTime> dt_;
```
With:
```cpp
struct Fields {
    int year = 0, month = 0, day = 0;
    int hour = 0, minute = 0, second = 0, millisecond = 0;
    bool valid = false;
} fields_;
```

Remove `#include <memory>` if it was only needed for `unique_ptr<QDateTime>` (check other uses first).

The class no longer needs the Rule-of-5 implementations (copy/move ctor, assignment) since the struct is trivially copyable. Simplify or default them.

- [ ] **Step 3: Rewrite DateTime.cpp — remove Qt, implement with sscanf/snprintf**

In `src/openms/source/DATASTRUCTURES/DateTime.cpp`:

Remove:
```cpp
#include <QtCore/QDateTime> // very expensive to include!
```

Add:
```cpp
#include <cstdio>    // sscanf, snprintf
#include <cstring>   // for legacy fallback parsing
#include <chrono>
#include <ctime>
```

Rewrite all methods using the `fields_` struct:

**Constructors/assignment:** Trivial — default-initialize `fields_` or copy struct.

**`isValid()`:** Return `fields_.valid`.

**`isNull()`:** Return `!fields_.valid` (matches Qt: default-constructed QDateTime is null and invalid).

**`clear()`:** `fields_ = Fields{};`

**`operator==`:** Compare all fields (year through millisecond) and valid flag.

**`operator<`:** Lexicographic compare: year, month, day, hour, minute, second, millisecond.

**`set(const String& date)` — 8 parse branches:**

Use `sscanf` for each branch. Example for branch 1 (German format):
```cpp
if (date.has('.') && !date.has('T'))
{
    if (sscanf(date.c_str(), "%d.%d.%d %d:%d:%d",
        &fields_.day, &fields_.month, &fields_.year,
        &fields_.hour, &fields_.minute, &fields_.second) == 6)
    {
        fields_.valid = true;
        return;
    }
}
```

For branch 3 (ISO with ms + tz suffix stripped):
```cpp
// Strip timezone suffix at '+'
String stripped = date.has('+') ? date.prefix('+') : date;
if (sscanf(stripped.c_str(), "%d-%d-%dT%d:%d:%d.%d",
    &fields_.year, &fields_.month, &fields_.day,
    &fields_.hour, &fields_.minute, &fields_.second, &fields_.millisecond) == 7)
{
    fields_.valid = true;
    return;
}
```

For branch 5 (`"yyyy-MM-ddZ"` — Z is discarded, NOT UTC):
```cpp
if (sscanf(date.c_str(), "%d-%d-%dZ", &fields_.year, &fields_.month, &fields_.day) == 3)
{
    fields_.hour = fields_.minute = fields_.second = 0;
    fields_.valid = true;
    return;
}
```

For branch 6 (`"yyyy-MM-dd+hh:mm"` — `+` is separator, NOT timezone):
```cpp
if (sscanf(date.c_str(), "%d-%d-%d+%d:%d",
    &fields_.year, &fields_.month, &fields_.day,
    &fields_.hour, &fields_.minute) == 5)
{
    fields_.second = 0;
    fields_.valid = true;
    return;
}
```

For branch 8 (legacy `ddd MMM d YYYY` fallback):
```cpp
// Lookup table for month abbreviations
static const std::map<std::string, int> months = {
    {"Jan",1},{"Feb",2},{"Mar",3},{"Apr",4},{"May",5},{"Jun",6},
    {"Jul",7},{"Aug",8},{"Sep",9},{"Oct",10},{"Nov",11},{"Dec",12}
};
char day_name[4], month_name[4];
int d, y;
if (sscanf(date.c_str(), "%3s %3s %d %d", day_name, month_name, &d, &y) == 4)
{
    auto it = months.find(month_name);
    if (it != months.end())
    {
        fields_.year = y; fields_.month = it->second; fields_.day = d;
        fields_.valid = true;
        return;
    }
}
```

**`set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second)`:** Assign fields directly, validate ranges.

**`get()` — output formatting:**
```cpp
if (!fields_.valid) return "0000-00-00 00:00:00";
char buf[32];
std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d",
    fields_.year, fields_.month, fields_.day,
    fields_.hour, fields_.minute, fields_.second);
return String(buf);
```

**`getDate()`:** `snprintf` with `"%04d-%02d-%02d"`.

**`getTime()`:** `snprintf` with `"%02d:%02d:%02d"`.

**`toString(format)` / `fromString(date, format)` — format dispatcher:**
```cpp
String DateTime::toString(const std::string& format) const
{
    if (!fields_.valid) return "";
    char buf[64];
    if (format == "yyyy-MM-ddThh:mm:ss")
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02d", fields_.year, fields_.month, fields_.day, fields_.hour, fields_.minute, fields_.second);
    else if (format == "yyyy-MM-ddThh:mm:ss.zzz")
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02d.%03d", fields_.year, fields_.month, fields_.day, fields_.hour, fields_.minute, fields_.second, fields_.millisecond);
    else if (format == "yyyy-MM-dd hh:mm:ss")
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d", fields_.year, fields_.month, fields_.day, fields_.hour, fields_.minute, fields_.second);
    else if (format == "yyyy-MM-dd+hh:mm")
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d+%02d:%02d", fields_.year, fields_.month, fields_.day, fields_.hour, fields_.minute);
    else if (format == "yyyy-MM-ddThh:mm:ssZ")
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02dZ", fields_.year, fields_.month, fields_.day, fields_.hour, fields_.minute, fields_.second);
    else if (format == "yyyy-MM-dd")
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d", fields_.year, fields_.month, fields_.day);
    else if (format == "hh:mm:ss")
        std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d", fields_.hour, fields_.minute, fields_.second);
    else
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown DateTime format: " + format);
    return String(buf);
}
```

**`setDate(const String& date)` — 3 date-only parse branches (same pattern as DateTime::set):**
```cpp
void DateTime::setDate(const String& date)
{
    int y, m, d;
    if (date.has('-') && sscanf(date.c_str(), "%d-%d-%d", &y, &m, &d) == 3)
        { fields_.year = y; fields_.month = m; fields_.day = d; }
    else if (date.has('.') && sscanf(date.c_str(), "%d.%d.%d", &d, &m, &y) == 3)
        { fields_.year = y; fields_.month = m; fields_.day = d; }
    else if (date.has('/') && sscanf(date.c_str(), "%d/%d/%d", &m, &d, &y) == 3)
        { fields_.year = y; fields_.month = m; fields_.day = d; }
    else
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Could not set date");
    // validate ...
}
```

**`setTime(const String& time)`:** `sscanf(time.c_str(), "%d:%d:%d", &h, &m, &s)`.

**`setDate(UInt month, UInt day, UInt year)`:** Assign fields directly, validate ranges, throw on invalid.

**`setTime(UInt hour, UInt minute, UInt second)`:** Assign fields directly, validate ranges.

**`get(UInt& month, UInt& day, UInt& year, UInt& hour, UInt& minute, UInt& second)`:** Extract from `fields_`.

**`getDate(UInt& month, UInt& day, UInt& year)`:** Extract from `fields_`.

**`getTime(UInt& hour, UInt& minute, UInt& second)`:** Extract from `fields_`.

**`fromString(const std::string& date, const std::string& format)` — static parse with format:**
```cpp
DateTime DateTime::fromString(const std::string& date, const std::string& format)
{
    DateTime d;
    // Reuse the same format dispatcher as toString, but in reverse (parsing)
    if (format == "yyyy-MM-ddThh:mm:ss" || format.empty())
        sscanf(date.c_str(), "%d-%d-%dT%d:%d:%d", &d.fields_.year, &d.fields_.month, &d.fields_.day, &d.fields_.hour, &d.fields_.minute, &d.fields_.second);
    else if (format == "yyyy-MM-ddThh:mm:ss.zzz")
        sscanf(date.c_str(), "%d-%d-%dT%d:%d:%d.%d", &d.fields_.year, &d.fields_.month, &d.fields_.day, &d.fields_.hour, &d.fields_.minute, &d.fields_.second, &d.fields_.millisecond);
    else if (format == "yyyy-MM-dd hh:mm:ss")
        sscanf(date.c_str(), "%d-%d-%d %d:%d:%d", &d.fields_.year, &d.fields_.month, &d.fields_.day, &d.fields_.hour, &d.fields_.minute, &d.fields_.second);
    else if (format == "yyyy-MM-dd+hh:mm")
        sscanf(date.c_str(), "%d-%d-%d+%d:%d", &d.fields_.year, &d.fields_.month, &d.fields_.day, &d.fields_.hour, &d.fields_.minute);
    else if (format == "yyyy-MM-ddThh:mm:ssZ")
        sscanf(date.c_str(), "%d-%d-%dT%d:%d:%dZ", &d.fields_.year, &d.fields_.month, &d.fields_.day, &d.fields_.hour, &d.fields_.minute, &d.fields_.second);
    else if (format == "yyyy-MM-dd")
        sscanf(date.c_str(), "%d-%d-%d", &d.fields_.year, &d.fields_.month, &d.fields_.day);
    else if (format == "hh:mm:ss")
        sscanf(date.c_str(), "%d:%d:%d", &d.fields_.hour, &d.fields_.minute, &d.fields_.second);
    else
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown DateTime format: " + format);
    d.fields_.valid = true; // TODO: validate ranges
    return d;
}
```

**ABI note:** The class changes from pimpl (`unique_ptr<QDateTime>`) to a value struct. The Rule-of-5 members (copy/move ctors, assignment ops) declared with `OPENMS_DLLAPI` in the header must be updated — either `= default` them or remove the explicit declarations. The `std::hash<DateTime>` specialization (DateTime.h ~line 206) calls `toString("yyyy-MM-ddThh:mm:ss.zzz")` — this format is supported by the dispatcher above. Verify pyOpenMS bindings in `src/pyOpenMS/bindings/` don't expose DateTime internals that would break.

**`now()` — local time:**
```cpp
DateTime DateTime::now()
{
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf{};
#ifdef _WIN32
    localtime_s(&tm_buf, &t);
#else
    localtime_r(&t, &tm_buf);
#endif
    DateTime d;
    d.fields_ = {tm_buf.tm_year + 1900, tm_buf.tm_mon + 1, tm_buf.tm_mday,
                  tm_buf.tm_hour, tm_buf.tm_min, tm_buf.tm_sec, 0, true};
    return d;
}
```

**`nowUTC()` — UTC time:**
Same but use `gmtime_r` / `_gmtime64_s`.

**`addSecs(int s)` — use chrono for carry arithmetic:**
```cpp
DateTime& DateTime::addSecs(int s)
{
    using namespace std::chrono;
    auto ymd = year{fields_.year} / month{unsigned(fields_.month)} / day{unsigned(fields_.day)};
    auto tp = sys_days{ymd} + hours{fields_.hour} + minutes{fields_.minute} + seconds{fields_.second};
    tp += seconds{s};
    auto dp = floor<days>(tp);
    year_month_day new_ymd{dp};
    hh_mm_ss hms{tp - dp};
    fields_.year = int(new_ymd.year());
    fields_.month = unsigned(new_ymd.month());
    fields_.day = unsigned(new_ymd.day());
    fields_.hour = hms.hours().count();
    fields_.minute = hms.minutes().count();
    fields_.second = hms.seconds().count();
    return *this;
}
```

- [ ] **Step 4: Run DateTime tests**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc) && ctest --test-dir OpenMS-build -R DateTime_test -V`
Expected: All pass

### Sub-task 2b: Date rewrite

- [ ] **Step 5: Replace Date internal representation in header**

In `src/openms/include/OpenMS/DATASTRUCTURES/Date.h`:

Remove forward declaration `class QDate;` (line 18).
Remove `Date(const QDate& date);` constructor (line 45).
Replace private member (line 123):
```cpp
std::unique_ptr<QDate> date_;
```
With:
```cpp
struct Fields { int year = 0, month = 0, day = 0; bool valid = false; } fields_;
```

Simplify Rule-of-5 (trivially copyable struct).

- [ ] **Step 6: Rewrite Date.cpp**

Same pattern as DateTime — use `sscanf`/`snprintf` for the 3 date formats (German `dd.MM.yyyy`, US `MM/dd/yyyy`, ISO `yyyy-MM-dd`). Remove `#include <QtCore/QDate>`.

`Date::today()`: Use `std::chrono::system_clock::now()` → `localtime_r` → extract year/month/day.

- [ ] **Step 7: Run Date tests**

Run: `ctest --test-dir OpenMS-build -R Date_test -V`
Expected: All pass

### Sub-task 2c: BuildInfo.h

- [ ] **Step 8: Replace QSysInfo in BuildInfo.h**

In `src/openms/include/OpenMS/SYSTEM/BuildInfo.h`:

Remove includes (lines 14-15):
```cpp
#include <QtCore/QSysInfo>
#include <QtCore/QString>
```

Replace `QSysInfo::WordSize` check (line 96) with `sizeof(void*) * 8` or call `getBinaryArchitecture()`.

Replace `QSysInfo::productVersion()` (line 93) with platform-specific:
```cpp
#if defined(_WIN32)
    // Use RtlGetVersion via ntdll
    info.os_version_ = Internal::getWindowsVersion();
#elif defined(__APPLE__)
    // sysctlbyname("kern.osproductversion", ...)
    info.os_version_ = Internal::getMacOSVersion();
#elif defined(__linux__)
    // Parse /etc/os-release VERSION_ID=
    info.os_version_ = Internal::getLinuxVersion();
#endif
```

Implement the three helper functions in `BuildInfo.cpp` (or a new `BuildInfo.cpp` if it only has SIMD code currently).

- [ ] **Step 9: Build and verify**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc)`
Expected: Compiles without Qt includes in DateTime.h, Date.h, BuildInfo.h

### Sub-task 2d: DimMapper QLocale and TOPPBase QDateTime

- [ ] **Step 10: Replace QLocale in DimMapper.cpp**

In `src/openms/source/KERNEL/DimMapper.cpp`:

Remove `#include <QtCore/QLocale>` (line 13).

Replace `QLocale::c().toString(value, 'f', valuePrecision())` (line 32) with a helper:
```cpp
// In DimMapper.cpp (local helper or in a shared utility)
static String formatWithGroupSeparators(double value, int precision)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", precision, value);
    String s(buf);
    // Insert thousand separators in integer part
    auto dot = s.find('.');
    String intpart = (dot != String::npos) ? s.substr(0, dot) : s;
    String fracpart = (dot != String::npos) ? s.substr(dot) : "";
    String result;
    int count = 0;
    for (int i = (int)intpart.size() - 1; i >= 0; --i)
    {
        if (count > 0 && count % 3 == 0 && intpart[i] != '-') result = "," + result;
        result = String(1, intpart[i]) + result;
        ++count;
    }
    return result + fracpart;
}
```

- [ ] **Step 11: Replace direct QDateTime in TOPPBase.cpp**

In `src/openms/source/APPLICATIONS/TOPPBase.cpp`:

Remove `#include <QtCore/QDateTime>` (or the line that pulls it in).

Replace all `QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString()` calls (lines ~1641, 1648, 1655, 1664, 1673, 1678, 1873) with:
```cpp
DateTime::now().get()
```

This requires `#include <OpenMS/DATASTRUCTURES/DateTime.h>` (may already be included).

- [ ] **Step 12: Build full and run all DateTime/Date/DimMapper tests**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc) && ctest --test-dir OpenMS-build -R "DateTime_test|Date_test" -V`
Expected: All pass

- [ ] **Step 13: Commit PR 2**

```bash
git add src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h \
        src/openms/source/DATASTRUCTURES/DateTime.cpp \
        src/openms/include/OpenMS/DATASTRUCTURES/Date.h \
        src/openms/source/DATASTRUCTURES/Date.cpp \
        src/openms/include/OpenMS/SYSTEM/BuildInfo.h \
        src/openms/source/SYSTEM/BuildInfo.cpp \
        src/openms/source/KERNEL/DimMapper.cpp \
        src/openms/source/APPLICATIONS/TOPPBase.cpp
git commit -m "refactor: replace QDateTime/QDate/QSysInfo/QLocale with std::chrono and platform APIs

DateTime and Date now use plain structs instead of unique_ptr<QDateTime/QDate>.
Parsing uses sscanf, formatting uses snprintf — locale-independent.
BuildInfo uses platform-specific APIs for OS version detection.
All 8 DateTime::set() parse formats and 7 toString formats preserved."
```

---

## Task 3: PR 3 — QDir/QFile/QFileInfo → std::filesystem

**Files:**
- Create: `src/openms/include/OpenMS/SYSTEM/PathUtils.h`
- Modify: `src/openms/source/SYSTEM/File.cpp:21-23` (remove Qt includes, full rewrite of Qt calls)
- Modify: `src/openms/include/OpenMS/SYSTEM/File.h` (signature changes for `copyDirRecursively`, `removeDir`)
- Modify: ~15 other files with QDir/QFileInfo includes (see list below)
- Test: `src/tests/class_tests/openms/source/File_test.cpp`

**Other files to modify (remove QDir/QFile/QFileInfo includes, replace usage):**
- `src/openms/source/SYSTEM/UpdateCheck.cpp`
- `src/openms/source/APPLICATIONS/ToolHandler.cpp:329-391`
- `src/openms/source/APPLICATIONS/TOPPBase.cpp:43`
- `src/openms/source/SYSTEM/PythonInfo.cpp:15` (QDir only — `QDir::isRelativePath`)
- `src/openms/source/SYSTEM/JavaInfo.cpp:15` (QDir only)
- `src/openms/source/CONCEPT/ClassTest.cpp:32`
- `src/openms/source/CONCEPT/FuzzyStringComparator.cpp:14`
- `src/openms/source/METADATA/ExperimentalDesign.cpp:22`
- `src/openms/source/METADATA/DocumentIdentifier.cpp:14`
- `src/openms/source/FORMAT/QcMLFile.cpp:21`
- `src/openms/source/FORMAT/ExperimentalDesignFile.cpp:17`
- `src/openms/source/FORMAT/MascotGenericFile.cpp:20`
- `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:18`
- `src/openms/source/FORMAT/DATAACCESS/SiriusFragmentAnnotation.cpp:13`
- `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:25`
- `src/openms/source/QC/MQEvidenceExporter.cpp:18`
- `src/openms/source/QC/MQMsmsExporter.cpp:18`
- `src/openms/source/MATH/STATISTICS/PosteriorErrorProbabilityModel.cpp:23`
- `src/openms/source/ANALYSIS/ID/IDRipper.cpp:13`
- `src/openms/source/ANALYSIS/ID/FIAMSScheduler.cpp:16`

- [ ] **Step 1: Create PathUtils.h with to_path() helper**

Create `src/openms/include/OpenMS/SYSTEM/PathUtils.h`:
```cpp
#pragma once
#include <filesystem>
#include <string>

namespace OpenMS
{
    /// Convert a UTF-8 std::string to std::filesystem::path safely on all platforms.
    inline std::filesystem::path to_path(const std::string& s)
    {
        return std::filesystem::path(
            std::u8string(reinterpret_cast<const char8_t*>(s.data()), s.size()));
    }
} // namespace OpenMS
```

- [ ] **Step 2: Run existing File tests to establish baseline**

Run: `ctest --test-dir OpenMS-build -R File_test -V`
Expected: All pass

- [ ] **Step 3: Rewrite File.cpp — replace Qt with std::filesystem**

In `src/openms/source/SYSTEM/File.cpp`:

Remove Qt includes (lines 21-23):
```cpp
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
```

Add:
```cpp
#include <OpenMS/SYSTEM/PathUtils.h>
```

Replace each Qt call per the mapping in the spec. Key replacements:

- `QDir::mkpath(path)` → `std::filesystem::create_directories(to_path(path))`
- `QFileInfo(f).exists()` → `std::filesystem::exists(to_path(f))`
- `QFile(f).size()` → `std::filesystem::file_size(to_path(f))`
- `QFileInfo(f).canonicalFilePath()` → `std::filesystem::canonical(to_path(f)).string()`
- `QFile::rename(a, b)` → `std::filesystem::rename(to_path(a), to_path(b))`
- `QFile::copy(a, b)` → `std::filesystem::copy_file(to_path(a), to_path(b))`
- `QDir::cleanPath(p)` → `to_path(p).lexically_normal().string()` (lexical only!)
- `QDir::tempPath()` → `std::filesystem::temp_directory_path().string()`
- `QDir::homePath()` → platform-specific `getenv("HOME")` / `getenv("USERPROFILE")`
- `QDir::separator()` → `std::filesystem::path::preferred_separator`

For `fileList()` with glob pattern and name-sorting:
```cpp
bool File::fileList(const String& dir, const String& file_pattern, StringList& output, bool full_path)
{
    namespace fs = std::filesystem;
    output.clear();
    auto dir_path = to_path(dir);
    if (!fs::is_directory(dir_path)) return false;
    for (const auto& entry : fs::directory_iterator(dir_path))
    {
        if (!entry.is_regular_file()) continue;
        String fname = entry.path().filename().string();
#ifdef _WIN32
        if (!PathMatchSpecA(fname.c_str(), file_pattern.c_str())) continue;
#else
        if (fnmatch(file_pattern.c_str(), fname.c_str(), 0) != 0) continue;
#endif
        output.push_back(full_path ? entry.path().string() : fname);
    }
    std::sort(output.begin(), output.end());
    return !output.empty();
}
```

For permission checks on Windows, use `_access_s()`:
```cpp
bool File::readable(const String& file)
{
    auto p = to_path(file);
    if (!std::filesystem::exists(p)) return false;
#ifdef _WIN32
    return _access_s(file.c_str(), 4) == 0; // 4 = read
#else
    auto perms = std::filesystem::status(p).permissions();
    return (perms & std::filesystem::perms::owner_read) != std::filesystem::perms::none;
#endif
}
```

Change `copyDirRecursively` and `removeDir` signatures from `const QString&` to `const String&`.

`removeDirRecursively()` simplifies to:
```cpp
bool File::removeDirRecursively(const String& dir_name)
{
    std::error_code ec;
    std::filesystem::remove_all(to_path(dir_name), ec);
    if (ec) { std::cerr << "Could not remove directory " << dir_name << ": " << ec.message() << std::endl; return false; }
    return true;
}
```

- [ ] **Step 4: Update File.h signatures**

Change `copyDirRecursively` and `removeDir` parameter types from `const QString&` to `const String&`.

- [ ] **Step 5: Fix secondary files — remove QDir/QFileInfo includes**

For each file in the list above, remove the Qt include and replace the specific usage:
- Files that only used `QDir` for `isRelativePath()` → `std::filesystem::path(p).is_relative()`
- Files that used `QFileInfo` for `exists()` → `std::filesystem::exists(to_path(p))` or `File::exists(p)`
- `ToolHandler.cpp` lines 329-391: Replace `QDir(path, "*.ttd")` + `entryList()` with `File::fileList()` or `std::filesystem::directory_iterator` + sort
- `OpenSwathOSWParquetWriter.cpp:400` — uses `File::copyDirRecursively()` with `.toQString()`. After signature change, pass `String` directly
- **`UpdateCheck.cpp` — cross-cutting file:** Has QDateTime/QSysInfo (PR 2 scope) AND QDir/QFileInfo. Handle ALL its Qt removal in this PR to avoid cross-PR conflicts. Replace QDateTime calls with `DateTime::now()` equivalents, QSysInfo with BuildInfo replacements from PR 2.

Add `#include <fnmatch.h>` on Unix / `#include <Shlwapi.h>` on Windows for glob matching in `fileList()`. On Windows, link `Shlwapi.lib` if not already linked.

- [ ] **Step 6: Build and run File tests + full test suite**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc) && ctest --test-dir OpenMS-build -R File_test -V`
Then: `ctest --test-dir OpenMS-build -j$(nproc)` (full suite to catch regressions)
Expected: All pass

- [ ] **Step 7: Commit PR 3**

```bash
git add src/openms/include/OpenMS/SYSTEM/PathUtils.h \
        src/openms/include/OpenMS/SYSTEM/File.h \
        src/openms/source/SYSTEM/File.cpp \
        # ... all other modified files
git commit -m "refactor: replace QDir/QFile/QFileInfo with std::filesystem

Introduces PathUtils.h with to_path() for UTF-8 safe path construction.
File::fileList() results are explicitly sorted (directory_iterator is unordered).
Uses lexically_normal() for cleanPath (lexical, not filesystem-touching).
Windows: uses _access_s() for ACL-aware permission checks."
```

---

## Task 4: PR 4 — QProcess/QObject → boost::process

**Files:**
- Modify: `src/openms/include/OpenMS/SYSTEM/ExternalProcess.h` (remove QObject, change API)
- Modify: `src/openms/source/SYSTEM/ExternalProcess.cpp` (full rewrite)
- Modify: `src/openms/include/OpenMS/SYSTEM/RWrapper.h` (change API)
- Modify: `src/openms/source/SYSTEM/RWrapper.cpp`
- Modify: `src/openms/source/SYSTEM/PythonInfo.cpp`
- Modify: `src/openms/source/SYSTEM/JavaInfo.cpp`
- Modify: `src/openms/include/OpenMS/APPLICATIONS/TOPPBase.h:882,886` (change `runExternalProcess_` signatures)
- Modify: `src/openms/source/APPLICATIONS/TOPPBase.cpp` (implementation)
- Modify: ~11 TOPP adapter files in `src/topp/` (update `runExternalProcess_` calls)
- Modify: `src/openms_gui/source/VISUAL/MISC/ExternalProcessMBox.cpp` (add idle callback)
- Modify: `src/openms_gui/include/OpenMS/VISUAL/MISC/ExternalProcessMBox.h`
- Modify: `cmake/cmake_findExternalLibs.cmake:26` (add boost process component)
- Modify: `src/openms/CMakeLists.txt` (add Boost::process link)
- Modify: `vcpkg.json` (add boost-process)
- Test: `src/tests/class_tests/openms/source/ExternalProcess_test.cpp`

- [ ] **Step 1: Add boost::process to build system**

**NOTE:** `boost::process` (v1) is **header-only** — it has no compiled library component. Do NOT add `process` to the `find_boost()` COMPONENTS list (that would cause cmake to fail looking for a non-existent compiled library).

In `vcpkg.json`: Add `"boost-process"` to the dependencies array (this ensures the headers are available; it will pull in transitive deps like `boost-asio`, `boost-system`, `boost-filesystem`).

No changes needed to `cmake/cmake_findExternalLibs.cmake` — the headers are found via `Boost::boost` (already linked at `src/openms/CMakeLists.txt:72`). No new link target needed.

Verify: `#include <boost/process.hpp>` compiles in a test source file.

- [ ] **Step 2: Rewrite ExternalProcess.h — remove QObject, add idle callback**

```cpp
#pragma once
#include <OpenMS/DATASTRUCTURES/String.h>
#include <functional>
#include <map>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI ExternalProcess
  {
  public:
    enum class RETURNSTATE { SUCCESS, NONZERO_EXIT, CRASH, FAILED_TO_START };
    enum class IO_MODE { NO_IO, READ_ONLY, WRITE_ONLY, READ_WRITE };

    ExternalProcess();
    ExternalProcess(std::function<void(const String&)> callbackStdOut,
                    std::function<void(const String&)> callbackStdErr);
    ~ExternalProcess();

    void setCallbacks(std::function<void(const String&)> callbackStdOut,
                      std::function<void(const String&)> callbackStdErr);

    RETURNSTATE run(const String& exe, const std::vector<String>& args,
                    const String& working_dir, bool verbose, String& error_msg,
                    IO_MODE io_mode = IO_MODE::READ_WRITE,
                    const std::map<String, String>& env = {},
                    std::function<void()> idle_callback = nullptr);

    RETURNSTATE run(const String& exe, const std::vector<String>& args,
                    const String& working_dir, bool verbose,
                    IO_MODE io_mode = IO_MODE::READ_WRITE,
                    const std::map<String, String>& env = {},
                    std::function<void()> idle_callback = nullptr);
  private:
    std::function<void(const String&)> callbackStdOut_;
    std::function<void(const String&)> callbackStdErr_;
  };
} // namespace OpenMS
```

- [ ] **Step 3: Rewrite ExternalProcess.cpp with boost::process**

Use `boost::process` v1 API (available in Boost 1.85):
```cpp
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <boost/process.hpp>
#include <boost/asio.hpp>

namespace bp = boost::process;
namespace OpenMS
{
    // ... constructors/setCallbacks same as before minus QProcess ...

    ExternalProcess::RETURNSTATE ExternalProcess::run(
        const String& exe, const std::vector<String>& args,
        const String& working_dir, bool verbose, String& error_msg,
        IO_MODE io_mode, const std::map<String, String>& env,
        std::function<void()> idle_callback)
    {
        error_msg.clear();
        boost::asio::io_context io;
        bp::async_pipe out_pipe(io), err_pipe(io);
        std::string out_buf, err_buf;

        // Build environment
        auto sys_env = boost::this_process::environment();
        for (const auto& [k, v] : env) sys_env[k] = v;

        // Build args
        std::vector<std::string> bp_args(args.begin(), args.end());

        // Start process
        std::error_code ec;
        bp::child child(
            exe, bp::args(bp_args),
            bp::std_out > out_pipe, bp::std_err > err_pipe,
            bp::env(sys_env),
            working_dir.empty() ? bp::start_dir(".") : bp::start_dir(working_dir),
            ec
        );

        if (ec || !child.valid()) {
            error_msg = "Process '" + exe + "' failed to start.";
            return RETURNSTATE::FAILED_TO_START;
        }

        // Async read with callbacks
        // ... (set up async reads on out_pipe/err_pipe that call callbackStdOut_/callbackStdErr_)

        // Poll loop with idle callback for GUI responsiveness
        while (child.running()) {
            io.poll();
            if (idle_callback) idle_callback();
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        io.run(); // drain remaining

        child.wait();
        int exit_code = child.exit_code();

        if (exit_code < 0) { /* CRASH */ return RETURNSTATE::CRASH; }
        if (exit_code != 0) { /* NONZERO */ return RETURNSTATE::NONZERO_EXIT; }
        return RETURNSTATE::SUCCESS;
    }
}
```

- [ ] **Step 4: Rewrite PythonInfo.cpp, JavaInfo.cpp, RWrapper.cpp**

Replace `QProcess` with `boost::process::child` (sync, no asio needed):
```cpp
// PythonInfo::canRun example
bp::child qp(python_executable, "--version", bp::std_out > bp::null, bp::std_err > bp::null);
if (!qp.wait_for(std::chrono::seconds(30))) {
    qp.terminate();
    // report timeout
}
```

Replace `QDir::isRelativePath()` with `std::filesystem::path(p).is_relative()`.

Update `RWrapper.h` public API: `QStringList` → `std::vector<String>`, `QString` → `String`.

- [ ] **Step 5: Update TOPPBase.h runExternalProcess_ signatures**

Change `const QString&` / `QStringList` params to `const String&` / `std::vector<String>`.
Change `std::map<QString, QString>` to `std::map<String, String>`.

- [ ] **Step 6: Update TOPP adapter files and core callers**

In each of the ~11 TOPP adapter files (MaRaClusterAdapter, CometAdapter, SageAdapter, etc.), replace `QStringList` argument construction with `std::vector<String>`.

**Also update these core files that directly call ExternalProcess::run() or RWrapper::runScript() with Qt types:**
- `src/openms/source/QC/DBSuitability.cpp:345` — calls `ep.run(adapter_name.toQString(), QStringList() << ...)`. Change to `ep.run(adapter_name, std::vector<String>{...}, ...)`
- `src/openms/source/PROCESSING/CALIBRATION/InternalCalibration.cpp:291,447,510` — calls `RWrapper::runScript()` with `QStringList`. Change to `std::vector<String>`

- [ ] **Step 7: Update ExternalProcessMBox in GUI**

In `src/openms_gui/source/VISUAL/MISC/ExternalProcessMBox.cpp`, pass a lambda as `idle_callback`:
```cpp
auto result = ep_.run(exe, args, working_dir, verbose, error_msg,
    ExternalProcess::IO_MODE::READ_WRITE, {},
    []() { QCoreApplication::processEvents(); });
```

Update `ExternalProcessMBox.h` signatures: `QString`/`QStringList` → `String`/`std::vector<String>`.

- [ ] **Step 8: Build and run tests**

Run: `cmake --build OpenMS-build -j$(nproc) && ctest --test-dir OpenMS-build -R "ExternalProcess|PythonInfo|JavaInfo" -V`
Expected: All pass

- [ ] **Step 9: Commit PR 4**

```bash
git add vcpkg.json cmake/ src/openms/ src/openms_gui/ src/topp/
git commit -m "refactor: replace QProcess/QObject with boost::process

ExternalProcess no longer inherits QObject. No Q_OBJECT macro in core.
Uses boost::process for subprocess management with async_pipe for streaming.
Adds idle_callback parameter for GUI event loop pumping.
Preserves 30s timeout behavior for PythonInfo/JavaInfo.
Updates all TOPP adapter files and ExternalProcessMBox."
```

---

## Task 5: PR 5a — Remove Qt string API from core

**Files:**
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/String.h:18,73,335`
- Modify: `src/openms/source/DATASTRUCTURES/String.cpp:43-46,271-273`
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h:80,306,349`
- Modify: `src/openms/source/DATASTRUCTURES/DataValue.cpp:130-134,357-363,764-767`
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h:203-211`
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/StringListUtils.h:42`
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/StringView.h:20`
- Delete: `src/openms/include/OpenMS/CONCEPT/Qt5Port.h`
- Modify: ~20 remaining core files with `toQString()` calls
- Test: `src/tests/class_tests/openms/source/String_test.cpp`
- Test: `src/tests/class_tests/openms/source/DataValue_test.cpp`

- [ ] **Step 1: Run baseline tests**

Run: `ctest --test-dir OpenMS-build -R "String_test|DataValue_test" -V`
Expected: All pass

- [ ] **Step 2: Remove QString from String class**

In `String.h`: Remove `class QString;` forward declaration (line 18), the constructor declaration (line 73), and `toQString()` declaration (line 335).

In `String.cpp`: Remove `String(const QString& s)` implementation and `toQString()` implementation. Remove Qt includes.

- [ ] **Step 3: Remove QString from DataValue class**

In `DataValue.h`: Remove `DataValue(const QString&)`, `operator=(const QString&)`, `toQString()`.
In `DataValue.cpp`: Remove implementations and Qt includes.

- [ ] **Step 4: Fix StringUtils::number() and remove toQString()**

In `StringUtils.h`, replace:
```cpp
static String number(double d, UInt n)
{
    return QString::number(d, 'f', n);
}
static QString toQString(const String & this_s)
{
    return QString(this_s.c_str());
}
```
With:
```cpp
static String number(double d, UInt n)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", static_cast<int>(n), d);
    return String(buf);
}
```

Remove Qt includes from this header.

- [ ] **Step 5: Remove fromQStringList from StringListUtils**

Remove the method declaration and implementation. Remove Qt includes.

- [ ] **Step 6: Remove forward declaration from StringView.h**

Remove `class QString;` from line 20.

- [ ] **Step 7: Delete Qt5Port.h**

Remove the file and its entry from `src/openms/include/OpenMS/CONCEPT/sources.cmake`.

- [ ] **Step 8: Fix all remaining toQString() calls in core files**

**Complete list of remaining core files with Qt usage after PRs 1-4:**

| File | Qt usage | Fix |
|------|----------|-----|
| `source/APPLICATIONS/INIUpdater.cpp:32` | `name.toQString().count(':')` | `std::count(name.begin(), name.end(), ':')` |
| `source/APPLICATIONS/ToolHandler.cpp:334-376` | `toQString()` for path ops | Remove — paths already String after PR 3 |
| `source/CONCEPT/FuzzyStringComparator.cpp:194-213` | 8x `toQString()` calls | Replace with direct `String` operations |
| `source/CONCEPT/ProgressLogger.cpp:60` | `QString::number(float, 'f', 2)` | `StringUtils::number(value, 2)` (now snprintf-based) |
| `source/FORMAT/MzTabFile.cpp:1049,1417` | `cells[i].toQString()` | Remove `.toQString()` — was round-tripping through `String(QString)` |
| `source/FORMAT/QcMLFile.cpp:63,68,73,141,146,151` | `toQString()` for path/name ops | Remove — use `String` directly |
| `source/FORMAT/ExperimentalDesignFile.cpp:28,36,38` | `toQString()` | Remove |
| `source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:933` | `toQString()` | Remove |
| `source/FORMAT/DATAACCESS/SiriusFragmentAnnotation.cpp:304` | `toQString()` | Remove |
| `source/METADATA/DocumentIdentifier.cpp:42` | `toQString()` | Remove |
| `source/ANALYSIS/ID/IDRipper.cpp:77` | `toQString()` | Remove |
| `source/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp:400` | `toQString()` | Remove (should be partially done in PR 3) |
| `source/QC/DBSuitability.cpp:345` | `toQString()`, `QStringList` | Should be done in PR 4; verify clean |
| `source/PROCESSING/CALIBRATION/InternalCalibration.cpp:291,447,510` | `QStringList` | Should be done in PR 4; verify clean |
| `source/MATH/STATISTICS/PosteriorErrorProbabilityModel.cpp:99,289` | `toQString()` for QDir path | Remove — QDir gone after PR 3 |
| `include/OpenMS/APPLICATIONS/TOPPBase.h:25,27` | `#include <QtCore/QString>`, `qcontainerfwd.h` | Remove includes, remove forward decls |
| `include/OpenMS/APPLICATIONS/ToolHandler.h:16` | `#include <QtCore/QStringList>` | Remove, change return types to `StringList` |
| `include/OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h:15` | `#include <QtCore/QStringList>` | Remove, replace with `#include <vector>` + `std::vector<String>` |
| `source/DATASTRUCTURES/StringListUtils.cpp` | `#include <QStringList>` for `fromQStringList` | Already deleted in step 5 |
| `source/CONCEPT/ClassTest.cpp:32` | `#include <QtCore/QFileInfo>` | Remove (QFileInfo usage gone after PR 3) |

- [ ] **Step 9: Verify no QtCore includes or Qt references remain in src/openms/**

Run these greps — all must return zero matches:
```bash
grep -r "QtCore\|Q_OBJECT\|qcontainerfwd" src/openms/include/ src/openms/source/
grep -r "class QString" src/openms/include/ src/openms/source/
grep -r "toQString\|fromQString\|QStringList\|QString" src/openms/include/ src/openms/source/
```
Expected: Zero matches for all three

- [ ] **Step 10: Build and run full test suite**

Run: `cmake --build OpenMS-build -j$(nproc) && ctest --test-dir OpenMS-build -j$(nproc)`
Expected: All pass

- [ ] **Step 11: Commit PR 5a**

```bash
git add src/openms/
git commit -m "refactor: remove all Qt string API from core library

Breaking change: String(QString), String::toQString(), DataValue(QString),
DataValue::toQString() removed. Use QString::fromStdString()/toStdString()
in GUI code instead. StringUtils::number() now uses snprintf."
```

---

## Task 6: PR 5b — Mechanical replacement in GUI/TOPP/tests

**Files:**
- Modify: ~59 files in `src/openms_gui/`
- Modify: ~20 files in `src/topp/`
- Modify: ~6 files in `src/tests/`

- [ ] **Step 1: Automated find-and-replace for toQString()**

```bash
# In src/openms_gui/, src/topp/, and src/tests/:
# Replace .toQString() with QString::fromStdString() wrapping
# This needs careful regex since the receiver varies:
#   foo.toQString() → QString::fromStdString(foo)
#   expr.toQString() → QString::fromStdString(expr)
```

This is best done per-file with editor support. The pattern is:
- `X.toQString()` → `QString::fromStdString(X)` where X is a String/std::string expression
- `String(qstr)` → `qstr.toStdString()` (or just `std::string(qstr.toStdString())`)

- [ ] **Step 2: Move Qt5Port.h toQSet helper into openms_gui**

The 3 GUI files that used `Qt5Port.h` (TableView.cpp, TreeView.cpp, FilterableList.cpp) need a local replacement. Either inline the helper or add it to a GUI utility header.

- [ ] **Step 3: Build and run full test suite**

Run: `cmake --build OpenMS-build -j$(nproc) && ctest --test-dir OpenMS-build -j$(nproc)`
Expected: All pass

- [ ] **Step 4: Commit PR 5b**

```bash
git add src/openms_gui/ src/topp/ src/tests/
git commit -m "refactor: replace toQString()/String(QString) in GUI, TOPP, and tests

Mechanical replacement: .toQString() → QString::fromStdString()
and String(qstr) → qstr.toStdString(). ~591 occurrences across 85 files."
```

---

## Task 7: PR 6 — CMake cleanup, cut Qt6::Core from core

**Files:**
- Modify: `src/openms/CMakeLists.txt:61`
- Modify: `cmake/cmake_findExternalLibs.cmake:285-299`

- [ ] **Step 1: Remove Qt6::Core from core link dependencies**

In `src/openms/CMakeLists.txt` line 61, remove `Qt6::Core` from the `OPENMS_DEP_LIBRARIES` list.

- [ ] **Step 2: Conditionalize Qt6 Core find_package**

In `cmake/cmake_findExternalLibs.cmake`, wrap the Qt6::Core section (lines 285-299) so it only runs when `WITH_GUI` is ON:

```cmake
if(WITH_GUI)
  find_package(Qt6 6.1.0 COMPONENTS Core REQUIRED)
  # ... existing GUI Qt6 discovery ...
endif()
```

- [ ] **Step 3: Verify zero Qt references in core**

Run: `grep -r "QtCore\|Qt6::Core" src/openms/`
Expected: Zero matches

Run: `grep -r "Q_OBJECT\|QObject\|QString\|QProcess\|QDateTime\|QDir\|QFile\|QJson\|QLocale\|QSysInfo" src/openms/include/ src/openms/source/`
Expected: Zero matches

- [ ] **Step 4: Build with WITH_GUI=OFF**

```bash
cmake -B OpenMS-build-nogui -DWITH_GUI=OFF <other flags>
cmake --build OpenMS-build-nogui -j$(nproc)
ctest --test-dir OpenMS-build-nogui -j$(nproc)
```
Expected: Builds and tests pass without Qt installed or found.

- [ ] **Step 5: Build with WITH_GUI=ON**

```bash
cmake --build OpenMS-build -j$(nproc)
ctest --test-dir OpenMS-build -j$(nproc)
```
Expected: Full build including GUI still works.

- [ ] **Step 5b: Verify pyOpenMS builds and tests pass**

```bash
cmake --build OpenMS-build --target pyopenms -j$(nproc)
cd /tmp && PYTHONPATH=/path/to/OpenMS-build/pyOpenMS python3 -m pytest /path/to/src/pyOpenMS/tests/ -v
```
Expected: pyOpenMS builds and tests pass. The DateTime/Date ABI change (pimpl → value struct) may require binding updates in `src/pyOpenMS/bindings/` if these classes are exposed. Check and fix if needed.

- [ ] **Step 6: Commit PR 6**

```bash
git add src/openms/CMakeLists.txt cmake/cmake_findExternalLibs.cmake
git commit -m "build: remove Qt6::Core from libOpenMS link dependencies

libOpenMS can now be built without Qt (WITH_GUI=OFF).
Qt6 is only required when building the GUI library."
```
