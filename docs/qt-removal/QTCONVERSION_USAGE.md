# QtConversion.h - Bridge Between Core and GUI

## Purpose

The `QtConversion.h` header provides lightweight conversion utilities between OpenMS core types (C++ standard library) and Qt types used in the GUI.

## Key Design Principles

### 1. **Conditional Compilation**
```cpp
#ifdef WITH_GUI
// Only compiled when building with GUI support
#endif
```

This ensures the core library remains Qt-free while GUI components can use these converters.

### 2. **Zero Runtime Overhead**
All functions are `static` and will typically be inlined by the compiler, resulting in zero function call overhead.

### 3. **Explicit Conversions**
Conversions are explicit (must be called), not implicit (no automatic type conversion). This prevents accidental mixing of Qt and standard types.

### 4. **UTF-8 Handling**
All conversions use UTF-8 encoding via `QString::fromStdString()` and `QString::toStdString()`, ensuring proper Unicode support.

## Usage Guidelines

### ✅ **CORRECT Usage (GUI Code)**

```cpp
// In src/openms_gui/ files only:

#include <OpenMS/DATASTRUCTURES/QtConversion.h>

// Converting core result to display in GUI
void SomeGUIClass::displayResult()
{
  std::string core_result = core_algorithm.getResult();
  QString display_text = QtConversion::toQt(core_result);
  label->setText(display_text);
}

// Converting user input to pass to core
void SomeGUIClass::processInput()
{
  QString user_text = lineEdit->text();
  std::string core_input = QtConversion::fromQt(user_text);
  core_algorithm.setParameter(core_input);
}
```

### ❌ **INCORRECT Usage (Core Code)**

```cpp
// NEVER in src/openms/ files:

#include <OpenMS/DATASTRUCTURES/QtConversion.h>  // ❌ NO!

void CoreClass::someMethod()
{
  QString qstr = ...;  // ❌ NO Qt types in core!
  std::string str = QtConversion::fromQt(qstr);  // ❌ NO!
}
```

## API Reference

### QString ↔ std::string

```cpp
// std::string → QString
QString toQt(const std::string& str);

// QString → std::string
std::string fromQt(const QString& qstr);
```

**Example:**
```cpp
std::string filename = "data.mzML";
QString qt_filename = QtConversion::toQt(filename);
// Use qt_filename with Qt widgets

QString user_path = QFileDialog::getOpenFileName(...);
std::string core_path = QtConversion::fromQt(user_path);
// Pass core_path to core library functions
```

### QStringList ↔ std::vector<std::string>

```cpp
// std::vector<std::string> → QStringList
QStringList toQt(const std::vector<std::string>& vec);

// QStringList → std::vector<std::string>
std::vector<std::string> fromQt(const QStringList& qlist);
```

**Example:**
```cpp
std::vector<std::string> file_list = core_function.getFileList();
QStringList qt_list = QtConversion::toQt(file_list);
listWidget->addItems(qt_list);
```

### QByteArray ↔ std::vector<char> / std::string

```cpp
// Binary data conversions
QByteArray toQt(const std::vector<char>& vec);
std::vector<char> fromQtToVector(const QByteArray& qarray);

// String-based binary data
QByteArray toQtBytes(const std::string& str);
std::string fromQtBytes(const QByteArray& qarray);
```

**Example:**
```cpp
// For binary data
std::vector<char> binary_data = read_binary_file();
QByteArray qt_data = QtConversion::toQt(binary_data);
network->post(url, qt_data);
```

### Helper Utilities

```cpp
// Join vector to QString with delimiter
QString joinToQt(const std::vector<std::string>& vec, const QString& delimiter);

// Split QString and convert to vector
std::vector<std::string> splitFromQt(const QString& qstr, const QString& delimiter);
```

**Example:**
```cpp
std::vector<std::string> params = {"param1", "param2", "param3"};
QString display = QtConversion::joinToQt(params, ", ");
// Result: "param1, param2, param3"

QString csv_line = "apple,banana,cherry";
std::vector<std::string> items = QtConversion::splitFromQt(csv_line, ",");
// Result: {"apple", "banana", "cherry"}
```

## Integration Pattern

### Before Qt Removal (Old Code)
```cpp
// In core library (OLD - uses Qt directly)
QString CoreClass::getName() const
{
  return name_;  // name_ is QString
}

// In GUI (OLD - direct Qt usage)
void GUIClass::display()
{
  QString name = core_obj.getName();
  label->setText(name);
}
```

### After Qt Removal (New Code)
```cpp
// In core library (NEW - Qt-free)
std::string CoreClass::getName() const
{
  return name_;  // name_ is now std::string
}

// In GUI (NEW - with conversion)
#include <OpenMS/DATASTRUCTURES/QtConversion.h>

void GUIClass::display()
{
  std::string name = core_obj.getName();
  QString qt_name = QtConversion::toQt(name);
  label->setText(qt_name);
}
```

## Build System Integration

### CMakeLists.txt for Core Library
```cmake
# src/openms/CMakeLists.txt
# Qt is NO LONGER linked here
set(OPENMS_DEP_LIBRARIES Evergreen LibSVM::LibSVM XercesC::XercesC)
# No Qt6::Core, no Qt6::Network
```

### CMakeLists.txt for GUI Library
```cmake
# src/openms_gui/CMakeLists.txt
# Qt IS linked here
find_package(Qt6 REQUIRED COMPONENTS Core Widgets Network)
target_link_libraries(OpenMS_GUI PUBLIC OpenMS Qt6::Core Qt6::Widgets Qt6::Network)
```

## Performance Considerations

### Zero-Cost Abstraction
The conversion functions are extremely lightweight:

```cpp
static QString toQt(const std::string& str)
{
  return QString::fromStdString(str);  // Single function call
}
```

This typically inlines to the same machine code as calling `QString::fromStdString()` directly.

### When to Cache Conversions
If you're converting the same string repeatedly in a tight loop:

```cpp
// ❌ BAD - converts in every iteration
for (int i = 0; i < 1000; ++i)
{
  QString qt_str = QtConversion::toQt(my_string);
  someFunction(qt_str);
}

// ✅ GOOD - convert once
QString qt_str = QtConversion::toQt(my_string);
for (int i = 0; i < 1000; ++i)
{
  someFunction(qt_str);
}
```

## Testing

### Unit Test Example
```cpp
// In tests/class_tests/openms/source/QtConversion_test.cpp

#ifdef WITH_GUI

START_TEST(QtConversion, toQt_string)
{
  std::string std_str = "Hello, World!";
  QString qt_str = QtConversion::toQt(std_str);
  TEST_EQUAL(qt_str.toStdString(), std_str);
}

START_TEST(QtConversion, fromQt_string)
{
  QString qt_str = "Hello, Qt!";
  std::string std_str = QtConversion::fromQt(qt_str);
  TEST_EQUAL(std_str, "Hello, Qt!");
}

START_TEST(QtConversion, toQt_vector)
{
  std::vector<std::string> vec = {"one", "two", "three"};
  QStringList qlist = QtConversion::toQt(vec);
  TEST_EQUAL(qlist.size(), 3);
  TEST_EQUAL(qlist[0].toStdString(), "one");
}

START_TEST(QtConversion, unicode_handling)
{
  std::string utf8_str = u8"Hello 世界 🌍";
  QString qt_str = QtConversion::toQt(utf8_str);
  std::string back = QtConversion::fromQt(qt_str);
  TEST_EQUAL(back, utf8_str);
}

#endif // WITH_GUI
```

## Migration Checklist

When updating a GUI file that interacts with the core:

- [ ] Include `QtConversion.h` at the top
- [ ] Identify all calls to core library functions
- [ ] For each function that now returns `std::string` instead of `QString`:
  - [ ] Add conversion: `QtConversion::toQt(core_result)`
- [ ] For each function parameter that now takes `std::string` instead of `QString`:
  - [ ] Add conversion: `QtConversion::fromQt(qt_value)`
- [ ] Test the GUI functionality
- [ ] Verify no performance regression

## Common Patterns

### Pattern 1: Display Core Data in QLabel
```cpp
std::string core_message = algorithm.getMessage();
label->setText(QtConversion::toQt(core_message));
```

### Pattern 2: Get User Input for Core
```cpp
QString user_input = lineEdit->text();
core_algorithm.setParameter(QtConversion::fromQt(user_input));
```

### Pattern 3: Populate QComboBox from Core List
```cpp
std::vector<std::string> options = core_function.getOptions();
comboBox->addItems(QtConversion::toQt(options));
```

### Pattern 4: Get Selected Item from QListWidget
```cpp
QString selected_qt = listWidget->currentItem()->text();
std::string selected = QtConversion::fromQt(selected_qt);
core_function.processSelection(selected);
```

## Rationale

### Why Not Just Use Qt Everywhere?
1. **Build Time**: Qt adds 2-4 hours to build time
2. **Package Size**: pyOpenMS packages are 3x larger with Qt
3. **Licensing**: LGPL concerns for some users
4. **Dependencies**: Headless servers don't need GUI libraries
5. **Conflicts**: Qt can clash with other Python GUI frameworks

### Why Not Just Use std::string Everywhere (Including GUI)?
Qt's `QString` is specifically designed for:
- Efficient UTF-16 handling
- Integration with Qt widgets
- Internationalization (i18n)
- Rich text formatting

For GUI code, `QString` is the right choice. The conversion layer lets us use the right tool in each layer.

## Future Improvements

Potential enhancements for future consideration:

1. **String View Support** (C++17)
   ```cpp
   static QString toQt(std::string_view sv);
   ```

2. **Move Semantics**
   ```cpp
   static QString toQt(std::string&& str);  // Take ownership
   ```

3. **Locale-Aware Conversions**
   ```cpp
   static QString toQtLocale(const std::string& str, const std::locale& loc);
   ```

## Questions?

If you encounter issues or edge cases during the Qt removal process, please:
1. Document the issue in the PR
2. Add a test case to `QtConversion_test.cpp`
3. Propose an enhancement to this utility class
