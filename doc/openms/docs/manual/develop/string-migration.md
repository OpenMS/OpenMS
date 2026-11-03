Migrating from OpenMS::String
=============================

OpenMS 3.6 removed its own string class. C++ code written against OpenMS 3.5, such as an external project or a TOPP
tool, needs these changes:

- `OpenMS::String` is replaced by `std::string`, and `OpenMS::StringView` by `std::string_view`.
- `StringList` still exists. It is now `std::vector<std::string>` (`OpenMS/DATASTRUCTURES/TypeAliases.h`).
- The member functions of `String` are free functions in the namespace `OpenMS::StringUtils`, declared in
  `OpenMS/DATASTRUCTURES/StringUtils.h`. They take the string as their first argument.
- The headers `String.h`, `StringConversions.h` and `StringUtilsSimple.h` are gone. Include `StringUtils.h` instead.
  The namespace `StringConversions` is still there, as a thin layer over `StringUtils`, but new code should call
  `StringUtils` directly.

In most places, replacing `String` with `std::string` is all it takes. The table lists the member functions and their
replacements; the section after it lists the changes that the compiler does not catch.

## Replacements

| OpenMS 3.5 | OpenMS 3.6 |
| --- | --- |
| `String(42)`, `String(3.14)` | `StringUtils::toStr(42)`, `StringUtils::toStr(3.14)` |
| `String(x, false)` (at most three decimals) | `StringUtils::toStr(x, false)` |
| `String(data_value)` | `StringUtils::toStr(data_value)`; see [DataValue](#datavalue-to-string) below |
| `s.toInt32()`, `s.toInt64()`, `s.toFloat()`, `s.toDouble()` | `StringUtils::toInt32(s)`, `StringUtils::toInt64(s)`, `StringUtils::toFloat(s)`, `StringUtils::toDouble(s)` |
| `s.toInt()` | `StringUtils::toInt32(s)`; there is no `toInt` |
| `s.hasPrefix(p)`, `s.hasSuffix(p)`, `s.hasSubstring(p)` | `s.starts_with(p)`, `s.ends_with(p)`, `s.contains(p)`, or `StringUtils::hasPrefix(s, p)`, `hasSuffix(s, p)`, `hasSubstring(s, p)` |
| `s.has(c)` | `StringUtils::has(s, c)` |
| `s.prefix(n)`, `s.suffix(n)`, `s.substr(pos, n)`, `s.chop(n)` | `StringUtils::prefix(s, n)`, `suffix(s, n)`, `substr(s, pos, n)`, `chop(s, n)` |
| `s.prefix(c)`, `s.suffix(c)` with a character `c` | `StringUtils::prefix(s, c)`, `suffix(s, c)`; see [prefix and suffix](#prefix-and-suffix-no-longer-throw) below |
| `s.trim()`, `s.toUpper()`, `s.toLower()`, `s.simplify()`, `s.firstToUpper()`, `s.reverse()`, `s.fillLeft(c, n)`, `s.fillRight(c, n)`, `s.substitute(a, b)`, `s.remove(c)`, `s.ensureLastChar(c)`, `s.removeWhitespaces()` | `StringUtils::trim(s)`, `toUpper(s)`, ... with the same arguments after `s`. Like before, they change `s` and return a reference to it. |
| a changed copy of `s` | `StringUtils::trimmed(s)`, `toUppered(s)`, `toLowered(s)`, `substituted(s, a, b)` return one and leave `s` unchanged |
| `s.quote(q, String::ESCAPE)`, `s.unquote(q, String::ESCAPE)`, `s.isQuoted(q)` | `StringUtils::quote(s, q, QuotingMethod::ESCAPE)`, `unquote(s, q, QuotingMethod::ESCAPE)`, `isQuoted(s, q)` |
| `String::QuotingMethod` (`NONE`, `ESCAPE`, `DOUBLE`) | `OpenMS::QuotingMethod`, an `enum class` with the same values |
| `s.split(c, parts)`, `s.split(separator, parts)`, `s.split_quoted(separator, parts)` | `StringUtils::split(s, c, parts)`, `split(s, separator, parts)`, `split_quoted(s, separator, parts)`, with `parts` a `std::vector<std::string>` |
| `s.concatenate(first, last, glue)` | `s = StringUtils::concatenate(first, last, glue)`, or `StringUtils::concatenate(container, glue)` for a whole container |
| `String::number(d, n)`, `String::numberLength(d, n)`, `String::random(n)` | `StringUtils::number(d, n)`, `numberLength(d, n)`, `random(n)` |
| `s.toQString()`, `String(q_string)` | `QString::fromStdString(s)`, `q_string.toStdString()` |
| `ListUtils::create<String>("a,b")` | `ListUtils::create<std::string>("a,b")` |

`std::string::contains` is part of C++23, which OpenMS requires.

## Changes the compiler does not catch

### DataValue to string

The removed constructor `String(const DataValue&)` turned any value into text: numbers, lists and empty values alike.
Its replacement is `StringUtils::toStr(data_value)`, which never throws.

The conversion operator of `DataValue` to `std::string` is strict: it throws `Exception::ConversionError` unless the
value holds a string. It is also used implicitly, so code that compiled and worked with `String` can compile with
`std::string` and throw at run time:

```cpp
// throws Exception::ConversionError if "score" holds a number
std::string score = hit.getMetaValue("score");

// works for every value type, like String(hit.getMetaValue("score")) did
std::string score = StringUtils::toStr(hit.getMetaValue("score"));
```

The same holds for `ParamValue`; use `StringUtils::toStr(param_value)`.

### Appending numbers

`String` appended numbers as text. For `std::string`, `StringUtils.h` declares the operators `+` and `+=` for numbers,
so `s + 5` and `s += 1.5` append `"5"` and `"1.5"` as before. Include `StringUtils.h` wherever you do this: without it,
`s + 5` does not compile, but `s += 5` still does, and appends the character with code 5 instead of the text `"5"`.

`std::to_string()` is no replacement for numbers you write to files. In C++23 it writes floating-point numbers with
six decimals, so `std::to_string(1e-7)` is `"0.000000"`. `StringUtils::toStr()` writes floating-point numbers as `String` did: with full
precision, in scientific notation for very small and very large values, and `NaN` for not-a-number values.

### prefix and suffix no longer throw

`s.prefix(c)` and `s.suffix(c)` threw `Exception::ElementNotFound` when `s` did not contain the character `c`.
`StringUtils::prefix(s, c)` and `StringUtils::suffix(s, c)` return the whole of `s` instead. Code that relied on the
exception has to check for the character itself, for example with `s.find(c) == std::string::npos`.
