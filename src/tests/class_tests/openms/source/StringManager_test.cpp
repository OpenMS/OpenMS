// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Init.h>
#include <OpenMS/FORMAT/HANDLERS/StringManager.h>
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>
#include <future>
#include <limits>
#include <type_traits>

using namespace OpenMS;
using namespace OpenMS::Internal;

static_assert(!std::is_constructible_v<XMLAttributes, XMLAttributes::Values&&>);
static_assert(!std::is_constructible_v<XMLAttributes, const XMLAttributes::Values&&>);
static_assert(std::is_copy_constructible_v<XMLAttributes>);
static_assert(std::is_copy_assignable_v<XMLAttributes>);
static_assert(!std::is_constructible_v<XMLAttributes, XMLAttributes::NativeStore&&>);

START_TEST(StringManager, "$Id$")
START_SECTION(UTF8 and UTF16 conversions own their storage and preserve Unicode)
{
  const std::string utf8 = "ASCII Grüße β \U0001F600 \U0010FFFF";
  const std::u16string utf16 = u"ASCII Grüße β \U0001F600 \U0010FFFF";
  TEST_TRUE(StringManager::convert(utf8) == utf16)
  TEST_EQUAL(StringManager::convert(utf16.c_str()), utf8)
  auto owned = StringManager::convertPtr(utf8);
  TEST_TRUE(std::u16string(owned.get()) == utf16)
  auto moved = std::move(owned);
  TEST_TRUE(owned == nullptr)
  TEST_EQUAL(StringManager::convert(moved.get()), utf8)
  TEST_TRUE(StringManager::convert("").empty())
  TEST_TRUE(StringManager::convert(u"").empty())
  TEST_TRUE(StringManager::convert(static_cast<const char*>(nullptr)).empty())
  TEST_TRUE(StringManager::convert(static_cast<const char16_t*>(nullptr)).empty())
  const std::string ascii = "abcdefghijklmnopqrstuvwxyz0123456789";
  const auto wide = StringManager::convert(ascii);
  TEST_TRUE(StringManager::isASCII(wide.c_str(), wide.size()))
  TEST_EQUAL(StringManager::convert(wide.c_str()), ascii)
}
END_SECTION

START_SECTION(Malformed Unicode is rejected without silently losing data)
{
  for (const char* invalid : {"\x80", "\xc0\xaf", "\xe2\x82", "\xed\xa0\x80", "\xf4\x90\x80\x80"})
  {
    TEST_EXCEPTION(Exception::ConversionError, StringManager::convert(invalid))
  }
  const char16_t high[]{0xd800, 0};
  const char16_t low[]{0xdc00, 0};
  const char16_t bad_pair[]{0xd800, 'a', 0};
  TEST_EXCEPTION(Exception::ConversionError, StringManager::convert(high))
  TEST_EXCEPTION(Exception::ConversionError, StringManager::convert(low))
  TEST_EXCEPTION(Exception::ConversionError, StringManager::convert(bad_pair))
}
END_SECTION

START_SECTION(UTF16 integers preserve signs whitespace and overflow checks)
{
  TEST_EQUAL(StringManager::parseInt(u" \t+42\n"), 42)
  TEST_EQUAL(StringManager::parseInt(u"-2147483648"), std::numeric_limits<Int>::min())
  TEST_EQUAL(StringManager::parseInt(u"2147483647"), std::numeric_limits<Int>::max())
  TEST_EXCEPTION(Exception::ConversionError, StringManager::parseInt(u"2147483648"))
  TEST_EXCEPTION(Exception::ConversionError, StringManager::parseInt(u"12x"))
  TEST_EXCEPTION(Exception::ConversionError, StringManager::parseInt(u""))
}
END_SECTION

START_SECTION(Native attributes retain ASCII fast lookup and UTF8 fallback)
{
  const XMLAttributes::Values values{{u"id", u"one"}, {u"idRef", u"two"}, {u"β", u"\U0001F600"}};
  const XMLAttributes attributes(values);
  TEST_TRUE(attributes.value("id") == attributes.value(u"id"))
  TEST_TRUE(attributes.value("idRef") == attributes.value(u"idRef"))
  TEST_TRUE(attributes.value("idR") == nullptr)
  TEST_TRUE(attributes.value("idReference") == nullptr)
  TEST_TRUE(attributes.value("missing") == nullptr)
  TEST_TRUE(attributes.value("β") != nullptr)
  TEST_TRUE(std::u16string(attributes.value("β")) == u"\U0001F600")
  TEST_EQUAL(attributes.index(u"idRef"), 1)
  TEST_TRUE(attributes.valueByIndex(1) == attributes.value("idRef"))
  TEST_TRUE(attributes.valueByIndex(3) == nullptr)
  TEST_EXCEPTION(Exception::ConversionError, attributes.value("\xc0\xaf"))
}
END_SECTION

START_SECTION(Attribute views copy and assign without copying their backing storage)
{
  XMLAttributes::Values first{{u"id", u"first"}}, second{{u"id", u"second"}};
  XMLAttributes original(first);
  XMLAttributes copied(original);
  XMLAttributes assigned(second);
  assigned = original;
  TEST_TRUE(copied.value(u"id") == first[0].second.c_str())
  TEST_TRUE(assigned.value(u"id") == first[0].second.c_str())
  first[0].second = u"changed";
  TEST_TRUE(std::u16string(copied.value(u"id")) == u"changed")
  TEST_TRUE(std::u16string(assigned.value(u"id")) == u"changed")
  TEST_TRUE(std::u16string(XMLAttributes(second).value(u"id")) == u"second")
}
END_SECTION

START_SECTION(UTF8 slices convert into reused storage)
{
  std::u16string target = u"previous content that is longer";
  const std::string text = "Grüße \U0001F600 tail";
  StringManager::convert(text.data(), 6, target); // the first six bytes are "Grüß"
  TEST_TRUE(target == u"Grüß")
  StringManager::convert(text.data(), text.size(), target);
  TEST_TRUE(target == u"Grüße \U0001F600 tail")
  StringManager::convert(text.data(), 0, target);
  TEST_TRUE(target.empty())
  TEST_EXCEPTION(Exception::ConversionError, StringManager::convert("a\xc3", 2, target))
}
END_SECTION

START_SECTION(Parser-owned UTF8 attributes convert lazily and are reusable)
{
  const std::string id = "first", accent = "é \U0001F600";
  XMLAttributes::NativeStore store;
  store.add({}, "id", id);
  store.add("p", "attr", accent);
  store.add({}, "β", "beta");
  const XMLAttributes attributes(store);
  TEST_EQUAL(attributes.size(), 3)
  TEST_EQUAL(attributes.index(u"id"), 0)
  TEST_EQUAL(attributes.index(u"p:attr"), 1)
  TEST_EQUAL(attributes.index(u"β"), 2)
  TEST_EQUAL(attributes.index(u"attr"), -1)
  TEST_EQUAL(attributes.index(u"i"), -1)
  TEST_EQUAL(attributes.index(u"idx"), -1)
  TEST_TRUE(std::u16string(attributes.value(u"p:attr")) == u"é \U0001F600")
  // conversion happens once, the pointer stays valid for the callback
  TEST_TRUE(attributes.value(u"p:attr") == attributes.valueByIndex(1))
  TEST_TRUE(attributes.value("p:attr") == attributes.valueByIndex(1))
  TEST_TRUE(std::u16string(attributes.value("β")) == u"beta")
  TEST_TRUE(attributes.value("missing") == nullptr)
  TEST_TRUE(attributes.valueByIndex(3) == nullptr)
  TEST_EXCEPTION(Exception::ConversionError, attributes.value("\xc0\xaf"))
  std::string utf8 = "untouched";
  TEST_FALSE(attributes.valueUTF8(u"missing", utf8))
  TEST_EQUAL(utf8, "untouched")
  TEST_TRUE(attributes.valueUTF8(u"p:attr", utf8))
  TEST_EQUAL(utf8, accent)
  TEST_TRUE(attributes.valueUTF8("id", utf8))
  TEST_EQUAL(utf8, "first")

  // the next element reuses the entries, including their cached UTF-16 values
  store.clear();
  store.add("p", "attr", "second");
  TEST_EQUAL(attributes.size(), 1)
  TEST_TRUE(std::u16string(attributes.value(u"p:attr")) == u"second")
  TEST_TRUE(attributes.value(u"id") == nullptr)

  // legacy UTF-16 attributes give the same UTF-8 values
  const XMLAttributes::Values values{{u"p:attr", u"é \U0001F600"}};
  const XMLAttributes legacy(values);
  TEST_TRUE(legacy.valueUTF8("p:attr", utf8))
  TEST_EQUAL(utf8, accent)
  TEST_EQUAL(legacy.size(), 1)
}
END_SECTION

START_SECTION(Concurrent initialization and independent conversions)
{
  std::vector<std::future<bool>> workers;
  for (Size i = 0; i < 8; ++i)
  {
    workers.emplace_back(std::async(std::launch::async, []
    {
      for (Size iteration = 0; iteration < 100; ++iteration)
      {
        libxmlInitialize();
        const auto text = StringManager::convert("parallel β \U0001F600");
        if (StringManager::convert(text.c_str()) != "parallel β \U0001F600") return false;
      }
      return true;
    }));
  }
  for (auto& worker : workers) { TEST_TRUE(worker.get()) }
}
END_SECTION
END_TEST
