# Qt to Standard C++ Replacement Patterns

This document provides quick reference patterns for replacing Qt classes with standard C++ equivalents.

## Table of Contents
- [QString → std::string](#qstring--stdstring)
- [QStringList → std::vector<std::string>](#qstringlist--vectorstdstring)
- [QFile → std::fstream](#qfile--stdfstream)
- [QFileInfo → std::filesystem](#qfileinfo--stdfilesystem)
- [QDir → std::filesystem](#qdir--stdfilesystem)
- [QDateTime → std::chrono](#qdatetime--stdchrono)
- [QByteArray → std::vector/std::string](#qbytearray--stdvectorstdstring)
- [QProcess → boost::process](#qprocess--boostprocess)
- [QList/QVector → std::vector](#qlistvector--stdvector)
- [QMap → std::map/unordered_map](#qmap--stdmapunordered_map)

---

## QString → std::string

### Before (Qt)
```cpp
QString str = "Hello";
QString str2 = QString::fromStdString("World");
str.append(" World");
bool empty = str.isEmpty();
int len = str.length();
QString upper = str.toUpper();
QString trimmed = str.trimmed();
QStringList parts = str.split(" ");
```

### After (std::string)
```cpp
std::string str = "Hello";
std::string str2 = "World";
str.append(" World");
bool empty = str.empty();
size_t len = str.length();

// For upper case (need custom or use boost)
#include <algorithm>
#include <cctype>
std::string upper = str;
std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

// For trimming
#include <algorithm>
std::string trimmed = str;
trimmed.erase(trimmed.begin(), std::find_if(trimmed.begin(), trimmed.end(), 
    [](unsigned char ch) { return !std::isspace(ch); }));
trimmed.erase(std::find_if(trimmed.rbegin(), trimmed.rend(),
    [](unsigned char ch) { return !std::isspace(ch); }).base(), trimmed.end());

// For splitting (need custom function)
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
std::vector<std::string> parts = split(str, ' ');
```

### Helper Utilities (Recommended)
Create in `StringUtils.h`:
```cpp
namespace OpenMS {
    inline std::string toUpper(const std::string& s) {
        std::string result = s;
        std::transform(result.begin(), result.end(), result.begin(), ::toupper);
        return result;
    }
    
    inline std::string toLower(const std::string& s) {
        std::string result = s;
        std::transform(result.begin(), result.end(), result.begin(), ::tolower);
        return result;
    }
    
    inline std::string trim(const std::string& s) {
        auto start = std::find_if(s.begin(), s.end(),
            [](unsigned char ch) { return !std::isspace(ch); });
        auto end = std::find_if(s.rbegin(), s.rend(),
            [](unsigned char ch) { return !std::isspace(ch); }).base();
        return (start < end) ? std::string(start, end) : std::string();
    }
}
```

---

## QStringList → std::vector<std::string>

### Before (Qt)
```cpp
QStringList list;
list << "one" << "two" << "three";
QString joined = list.join(", ");
int count = list.count();
```

### After (std::vector<std::string>)
```cpp
std::vector<std::string> list = {"one", "two", "three"};

// For joining
#include <sstream>
std::string join(const std::vector<std::string>& v, const std::string& delimiter) {
    std::ostringstream result;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i > 0) result << delimiter;
        result << v[i];
    }
    return result.str();
}
std::string joined = join(list, ", ");
size_t count = list.size();
```

---

## QFile → std::fstream

### Before (Qt)
```cpp
QFile file("data.txt");
if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    QByteArray data = file.readAll();
    file.close();
}

if (file.exists("data.txt")) {
    qint64 size = file.size();
}
```

### After (std::fstream)
```cpp
#include <fstream>
#include <filesystem>

std::ifstream file("data.txt");
if (file.is_open()) {
    std::string data((std::istreambuf_iterator<char>(file)),
                     std::istreambuf_iterator<char>());
    file.close();
}

if (std::filesystem::exists("data.txt")) {
    auto size = std::filesystem::file_size("data.txt");
}
```

---

## QFileInfo → std::filesystem

### Before (Qt)
```cpp
QFileInfo fi("path/to/file.txt");
QString absolutePath = fi.absoluteFilePath();
QString fileName = fi.fileName();
QString baseName = fi.baseName();
QString extension = fi.suffix();
qint64 size = fi.size();
QDateTime modified = fi.lastModified();
bool exists = fi.exists();
bool isFile = fi.isFile();
bool isDir = fi.isDir();
bool readable = fi.isReadable();
bool writable = fi.isWritable();
```

### After (std::filesystem)
```cpp
#include <filesystem>
namespace fs = std::filesystem;

fs::path p("path/to/file.txt");
std::string absolutePath = fs::absolute(p).string();
std::string fileName = p.filename().string();
std::string baseName = p.stem().string();
std::string extension = p.extension().string();
auto size = fs::file_size(p);
auto modified = fs::last_write_time(p);
bool exists = fs::exists(p);
bool isFile = fs::is_regular_file(p);
bool isDir = fs::is_directory(p);

// For permissions
auto perms = fs::status(p).permissions();
bool readable = (perms & fs::perms::owner_read) != fs::perms::none;
bool writable = (perms & fs::perms::owner_write) != fs::perms::none;
```

---

## QDir → std::filesystem

### Before (Qt)
```cpp
QDir dir("/path/to/directory");
if (dir.exists()) {
    QStringList files = dir.entryList(QDir::Files);
    for (const QString& file : files) {
        // process file
    }
}

QDir::mkpath("path/to/new/directory");
bool isRelative = QDir::isRelativePath("relative/path");
QString absolute = QDir::toNativeSeparators("/unix/path");
```

### After (std::filesystem)
```cpp
#include <filesystem>
namespace fs = std::filesystem;

fs::path dir("/path/to/directory");
if (fs::exists(dir) && fs::is_directory(dir)) {
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.is_regular_file()) {
            std::string file = entry.path().filename().string();
            // process file
        }
    }
}

fs::create_directories("path/to/new/directory");
bool isRelative = fs::path("relative/path").is_relative();
std::string absolute = fs::path("/unix/path").make_preferred().string();
```

---

## QDateTime → std::chrono

### Before (Qt)
```cpp
QDateTime now = QDateTime::currentDateTime();
QString formatted = now.toString("yyyy-MM-dd HH:mm:ss");
QDateTime parsed = QDateTime::fromString("2024-01-15 10:30:00", "yyyy-MM-dd HH:mm:ss");
qint64 timestamp = now.toSecsSinceEpoch();
```

### After (std::chrono)
```cpp
#include <chrono>
#include <iomanip>
#include <sstream>

auto now = std::chrono::system_clock::now();

// For formatting (C++20)
#include <format> // C++20
std::string formatted = std::format("{:%Y-%m-%d %H:%M:%S}", now);

// For C++17 (without format):
#include <ctime>
auto now_c = std::chrono::system_clock::to_time_t(now);
std::stringstream ss;
ss << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S");
std::string formatted = ss.str();

// Timestamp
auto timestamp = std::chrono::duration_cast<std::chrono::seconds>(
    now.time_since_epoch()).count();

// Parsing (more complex, consider boost::date_time if needed)
```

### Alternative: boost::date_time
```cpp
#include <boost/date_time/posix_time/posix_time.hpp>

boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
std::string formatted = boost::posix_time::to_simple_string(now);
```

---

## QByteArray → std::vector/std::string

### Before (Qt)
```cpp
QByteArray data;
data.append("Hello");
data.append(someCharArray, length);
int size = data.size();
const char* rawData = data.constData();
```

### After (std::vector<char> or std::string)
```cpp
// Option 1: std::vector<char>
std::vector<char> data;
std::string str = "Hello";
data.insert(data.end(), str.begin(), str.end());
data.insert(data.end(), someCharArray, someCharArray + length);
size_t size = data.size();
const char* rawData = data.data();

// Option 2: std::string (if text data)
std::string data;
data.append("Hello");
data.append(someCharArray, length);
size_t size = data.size();
const char* rawData = data.c_str();
```

---

## QProcess → boost::process

### Before (Qt)
```cpp
QProcess process;
process.start("python", QStringList() << "script.py" << "--arg");
process.waitForFinished();
QByteArray output = process.readAllStandardOutput();
QByteArray error = process.readAllStandardError();
int exitCode = process.exitCode();
```

### After (boost::process)
```cpp
#include <boost/process.hpp>
namespace bp = boost::process;

std::string output;
std::string error;

bp::ipstream out_stream;
bp::ipstream err_stream;

bp::child process(
    "python script.py --arg",
    bp::std_out > out_stream,
    bp::std_err > err_stream
);

std::string line;
while (std::getline(out_stream, line)) {
    output += line + "\n";
}
while (std::getline(err_stream, line)) {
    error += line + "\n";
}

process.wait();
int exitCode = process.exit_code();
```

---

## QList/QVector → std::vector

### Before (Qt)
```cpp
QList<int> numbers;
numbers << 1 << 2 << 3;
numbers.append(4);
int first = numbers.first();
numbers.removeAt(0);
```

### After (std::vector)
```cpp
std::vector<int> numbers = {1, 2, 3};
numbers.push_back(4);
int first = numbers.front();
numbers.erase(numbers.begin());
```

---

## QMap → std::map/unordered_map

### Before (Qt)
```cpp
QMap<QString, int> map;
map["key1"] = 10;
map.insert("key2", 20);
bool contains = map.contains("key1");
int value = map.value("key1", -1);
```

### After (std::map or std::unordered_map)
```cpp
std::map<std::string, int> map;
map["key1"] = 10;
map.insert({"key2", 20});
bool contains = map.find("key1") != map.end();
int value = (map.find("key1") != map.end()) ? map["key1"] : -1;

// Or use unordered_map for hash-based lookup
std::unordered_map<std::string, int> fast_map;
```

---

## Common Patterns in OpenMS

### Pattern 1: String Conversion at Boundaries
```cpp
// When OpenMS String interacts with Qt (in GUI code only)
#ifdef WITH_GUI
QString toQt(const OpenMS::String& s) {
    return QString::fromStdString(s);
}

OpenMS::String fromQt(const QString& s) {
    return OpenMS::String(s.toStdString());
}
#endif
```

### Pattern 2: File Path Handling
```cpp
// Before
QFileInfo fi(filename.toQString());
if (fi.exists() && fi.isReadable()) {
    // process
}

// After
namespace fs = std::filesystem;
fs::path p(filename);
if (fs::exists(p) && (fs::status(p).permissions() & fs::perms::owner_read) != fs::perms::none) {
    // process
}
```

### Pattern 3: Process Execution
```cpp
// Before
QProcess p;
p.start(executable.toQString(), args);
p.waitForFinished();

// After
#include <boost/process.hpp>
namespace bp = boost::process;
bp::child p(executable, bp::args = args);
p.wait();
```

---

## Build System Changes

### CMakeLists.txt

Before:
```cmake
set(OPENMS_DEP_LIBRARIES Evergreen LibSVM::LibSVM XercesC::XercesC Qt6::Core Qt6::Network)
```

After:
```cmake
set(OPENMS_DEP_LIBRARIES Evergreen LibSVM::LibSVM XercesC::XercesC)
# Qt is still needed for openms_gui
if(WITH_GUI)
  find_package(Qt6 REQUIRED COMPONENTS Core Widgets Network)
endif()
```

---

## Testing After Changes

```bash
# Build core library
cmake -B build -S . -DWITH_GUI=OFF
cmake --build build --target OpenMS -j4

# Run specific tests
ctest --test-dir build -R String_test -V

# Verify no Qt linkage
ldd build/lib/libOpenMS.so | grep -i qt  # Should return nothing
```

---

## Reference Links

- [std::filesystem](https://en.cppreference.com/w/cpp/filesystem)
- [std::chrono](https://en.cppreference.com/w/cpp/chrono)
- [Boost.Process](https://www.boost.org/doc/libs/1_83_0/doc/html/process.html)
- [Boost.DateTime](https://www.boost.org/doc/libs/1_83_0/doc/html/date_time.html)
