// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/test_config.h>
#include <fstream>
#include <future>
#include <map>
#include <stdexcept>
#include <sstream>

using namespace OpenMS;
namespace
{
  struct Reader : Internal::XMLFile
  {
    using XMLFile::parseBuffer_;
    using XMLFile::parse_;
    using XMLFile::enforceEncoding_;
  };
  struct NamespaceReader : Reader
  {
    void parseBuffer_(const std::string& xml, Internal::XMLHandler* handler)
    {
      Internal::XMLFile::parseBuffer_(xml, handler, true);
    }
  };
  struct Handler : Internal::XMLHandler
  {
    Handler() : XMLHandler("memory", "1.0") {}
    std::u16string text, attribute, unicode_attribute;
    Size starts = 0;
    bool stop = false, fail = false;
    void onStartElement(const char16_t*, const Internal::XMLAttributes& attrs) override
    {
      ++starts;
      if (stop) { throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
      if (fail) { throw std::runtime_error("callback exception"); }
      if (attrs.value("clé")) { unicode_attribute = attrs.value("clé"); }
      if (attrs.value("value"))
      {
        attribute = attrs.value("value");
        if (attrs.index(u"value") < 0 || !attrs.valueByIndex(attrs.index(u"value")))
        {
          throw std::runtime_error("attribute lookup mismatch");
        }
      }
    }
    void onCharacters(const char16_t* value, Size length) override { text.append(value, length); }
    void onEndElement(const char16_t*) override {}
  };
}
START_TEST(XMLParserBackend, "$Id$")
START_SECTION((UTF8 Unicode, entities, CDATA, attributes, and chunk boundaries))
  Reader reader;
  Handler handler;
  const std::string repeated(70000, 'x');
  reader.parseBuffer_("<root clé='é🧬' value='Grüße &amp; &#x1F9EC;'>&lt;<![CDATA[&]]>" + repeated + "🧬</root>", &handler);
  TEST_EQUAL(handler.attribute == u"Grüße & 🧬", true)
  TEST_EQUAL(handler.unicode_attribute == u"é🧬", true)
  TEST_EQUAL(handler.text == std::u16string(u"<&") + std::u16string(70000, u'x') + u"🧬", true)
  TEST_EQUAL(handler.starts, 1)
END_SECTION
START_SECTION(Malformed XML and exception propagation)
  Reader reader;
  Handler malformed;
  TEST_EXCEPTION(Exception::ParseError, reader.parseBuffer_("<root><broken></root>", &malformed))
  Handler failure;
  failure.fail = true;
  bool callback_failed = false;
  try { reader.parseBuffer_("<root/>", &failure); }
  catch (const std::runtime_error& error)
  {
    callback_failed = true;
    TEST_STRING_EQUAL(error.what(), "callback exception")
  }
  TEST_TRUE(callback_failed)
  Handler stopped;
  stopped.stop = true;
  reader.parseBuffer_("<root><ignored/></root>", &stopped);
  TEST_EQUAL(stopped.starts, 1)
END_SECTION
START_SECTION(Parser errors retain the actual diagnostic and location)
  Reader reader;
  Handler handler;
  bool caught = false;
  try
  {
    reader.parseBuffer_("<root>\n<broken></root>", &handler);
  }
  catch (const Exception::ParseError& error)
  {
    caught = true;
    const std::string message = error.getMessage();
    TEST_TRUE(message.find("Opening and ending tag mismatch") != std::string::npos)
    TEST_TRUE(message.find("broken") != std::string::npos)
    TEST_TRUE(message.find("root") != std::string::npos)
    TEST_TRUE(message.find("( in line 2 column ") != std::string::npos)
  }
  TEST_TRUE(caught)
END_SECTION
START_SECTION(Generic readers preserve legacy namespace-disabled qualified names)
  struct Legacy : Internal::XMLHandler
  {
    Legacy() : XMLHandler("legacy", "1.0") {}
    std::u16string name, attribute;
    void onStartElement(const char16_t* qname, const Internal::XMLAttributes& attrs) override
    {
      name = qname;
      if (const auto* value = attrs.value("OpenXQuest:log_occupancy")) attribute = value;
    }
    void onEndElement(const char16_t*) override {}
    void onCharacters(const char16_t*, Size) override {}
  } qualified, namespace_uri;
  Reader reader;
  reader.parseBuffer_("<OpenXQuest:root OpenXQuest:log_occupancy='legacy'/>", &qualified);
  TEST_TRUE(qualified.name == u"OpenXQuest:root")
  TEST_TRUE(qualified.attribute == u"legacy")
  reader.parseBuffer_("<root xmlns:xsi='http://www.w3.org/2001/XMLSche ma-instance'/>", &namespace_uri);
  TEST_TRUE(namespace_uri.name == u"root")
END_SECTION
START_SECTION(Comments and processing instructions are ignored across streaming chunks)
  std::string filename;
  NEW_TMP_FILE(filename)
  {
    std::ofstream output(filename, std::ios::binary);
    output << "<!-- before --><?before ignored?><root value='kept'>begin";
    for (Size i = 0; i < 10000; ++i)
    {
      output << "<!-- ignored comment --><?ignored processing instruction?>";
    }
    output << "end</root><?after ignored?><!-- after -->";
  }
  for (const bool namespace_aware : {false, true})
  {
    Reader reader;
    Handler handler;
    reader.parse_(filename, &handler, namespace_aware);
    TEST_EQUAL(handler.starts, 1)
    TEST_TRUE(handler.attribute == u"kept")
    TEST_TRUE(handler.text == u"beginend")
  }
END_SECTION
START_SECTION(SAX2 warning-only namespace diagnostics do not abort parsing)
  NamespaceReader reader;
  Handler handler;
  std::ostringstream warnings;
  getThreadLocalLogWarn().insert(warnings);
  getThreadLocalLogDebug().insert(warnings);
  try
  {
    reader.parseBuffer_("<root xmlns='relative'>valid text</root>", &handler);
  }
  catch (...)
  {
    getThreadLocalLogWarn().remove(warnings);
    getThreadLocalLogDebug().remove(warnings);
    throw;
  }
  getThreadLocalLogWarn().remove(warnings);
  getThreadLocalLogDebug().remove(warnings);
  TEST_EQUAL(handler.starts, 1)
  TEST_TRUE(handler.text == u"valid text")
  TEST_TRUE(warnings.str().find("URI relative is not absolute") != std::string::npos)
END_SECTION
START_SECTION(SAX2 rejects namespace-invalid XML)
  for (const std::string xml : {
      "<root ghost:flag='1'/>",
      "<ghost:root/>",
      "<root xmlns:a='urn:same' xmlns:b='urn:same' a:value='1' b:value='2'/>",
      "<root xmlns:xml='urn:wrong'/>",
      "<root xmlns:xmlns='urn:forbidden'/>",
      "<root xmlns:p='http://www.w3.org/XML/1998/namespace'/>"})
  {
    NamespaceReader reader;
    Handler handler;
    TEST_EXCEPTION(Exception::ParseError, reader.parseBuffer_(xml, &handler))
  }
END_SECTION
START_SECTION(SAX2 preserves qualified names namespace declarations and bounded attribute values)
  struct Namespaces : Internal::XMLHandler
  {
    Namespaces() : XMLHandler("namespaces", "1.0") {}
    std::vector<std::u16string> starts, ends;
    std::map<std::u16string, std::u16string> values;
    void onStartElement(const char16_t* name, const Internal::XMLAttributes& attrs) override
    {
      starts.emplace_back(name);
      for (const char16_t* key : {u"xmlns", u"xmlns:p", u"xmlns:xsi", u"value", u"xsi:schemaLocation",
                                  u"mzML:version", u"p:attr", u"tail"})
      {
        if (const auto* value = attrs.value(key)) values[key] = value;
      }
    }
    void onEndElement(const char16_t* name) override { ends.emplace_back(name); }
    void onCharacters(const char16_t*, Size) override {}
  } handler;
  NamespaceReader reader;
  reader.parseBuffer_("<!DOCTYPE root [<!ENTITY named 'entity'>]>"
      "<root xmlns='urn:default' xmlns:p='urn:prefix' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' "
      "xmlns:mzML='urn:mzml' value='&named; &amp; &#x1F9EC;' xsi:schemaLocation='urn:default schema.xsd' mzML:version='1.1'>"
      "<p:child p:attr='é &lt; &#946;' tail='last'/></root>", &handler);
  TEST_TRUE(handler.starts == std::vector<std::u16string>({u"root", u"p:child"}))
  TEST_TRUE(handler.ends == std::vector<std::u16string>({u"p:child", u"root"}))
  TEST_TRUE(handler.values[u"xmlns"] == u"urn:default")
  TEST_TRUE(handler.values[u"xmlns:p"] == u"urn:prefix")
  TEST_TRUE(handler.values[u"xmlns:xsi"] == u"http://www.w3.org/2001/XMLSchema-instance")
  TEST_TRUE(handler.values[u"value"] == u"entity & 🧬")
  TEST_TRUE(handler.values[u"xsi:schemaLocation"] == u"urn:default schema.xsd")
  TEST_TRUE(handler.values[u"mzML:version"] == u"1.1")
  TEST_TRUE(handler.values[u"p:attr"] == u"é < β")
  TEST_TRUE(handler.values[u"tail"] == u"last")
  Handler stopped;
  stopped.stop = true;
  reader.parseBuffer_("<p:root xmlns:p='urn:prefix'><p:ignored/></p:root>", &stopped);
  TEST_EQUAL(stopped.starts, 1)
  Handler failed;
  failed.fail = true;
  TEST_EXCEPTION(std::runtime_error, reader.parseBuffer_("<p:root xmlns:p='urn:prefix'/>", &failed))
END_SECTION
START_SECTION(External entity loading is disabled)
  std::string filename;
  NEW_TMP_FILE(filename)
  { std::ofstream out(filename); out << "external secret"; }
  Reader reader;
  Handler handler;
  const std::string xml = "<!DOCTYPE root [<!ENTITY external SYSTEM 'file://" + filename + "'>]><root>&external;</root>";
  reader.parseBuffer_(xml, &handler);
  TEST_TRUE(handler.text.empty())
END_SECTION
START_SECTION(Validator stream exceptions are rethrown after C callbacks return)
  struct FailingBuffer : std::streambuf
  {
    int_type overflow(int_type) override { return traits_type::eof(); }
  } buffer;
  std::ostream output(&buffer);
  output.exceptions(std::ios::badbit | std::ios::failbit);
  XMLValidator validator;
  TEST_EXCEPTION(std::ios_base::failure,
    validator.isValid(OPENMS_GET_TEST_DATA_PATH("XMLValidator_missing_attribute.xml"),
                      OPENMS_GET_TEST_DATA_PATH("XMLValidator.xsd"), output))
END_SECTION
START_SECTION(Encoding overrides are per reader)
  Reader latin, utf8;
  latin.enforceEncoding_("ISO-8859-1");
  utf8.enforceEncoding_("UTF-8");
  Handler latin_handler, utf8_handler;
  latin.parseBuffer_("<root>Gr\xfc\xdf" "e</root>", &latin_handler);
  utf8.parseBuffer_("<root>Grüße</root>", &utf8_handler);
  TEST_EQUAL(latin_handler.text == u"Grüße", true)
  TEST_EQUAL(utf8_handler.text == u"Grüße", true)
END_SECTION
START_SECTION(ISO-8859-1 declarations decode text and attributes across read boundaries)
  for (const std::string declaration : {std::string("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"),
                                        std::string("<?xml version='1.0' encoding = 'iso-8859-1' standalone='yes'?>\n")})
  {
    Reader reader;
    Handler handler;
    std::string xml = declaration + "<root value='Gr\xfc\xdf" "e'>";
    // Put non-ASCII bytes on both sides of the 64 KiB read boundary.
    const Size padding = 65535 - xml.size();
    xml += std::string(padding, 'x') + "\xe9\xe9\xff</root>";
    reader.parseBuffer_(xml, &handler);
    TEST_TRUE(handler.attribute == u"Grüße")
    TEST_TRUE(handler.text == std::u16string(padding, u'x') + u"ééÿ")
  }
  // Other declarations keep libxml2's decoding: these bytes are invalid UTF-8.
  Reader reader;
  Handler handler;
  TEST_EXCEPTION(Exception::ParseError, reader.parseBuffer_("<?xml version='1.0' encoding='UTF-8'?><root>\xfc</root>", &handler))
END_SECTION
START_SECTION(Unsupported encoding diagnostics retain the input filename)
  std::string filename;
  NEW_TMP_FILE(filename)
  { std::ofstream output(filename); output << "<root/>"; }
  Reader reader;
  reader.enforceEncoding_("not-an-encoding");
  Handler handler;
  bool caught = false;
  try { reader.parse_(filename, &handler); }
  catch (const Exception::ParseError& error)
  {
    caught = true;
    const std::string message = error.getMessage();
    TEST_TRUE(message.find(filename) != std::string::npos)
    TEST_TRUE(message.find("not-an-encoding") != std::string::npos)
  }
  TEST_TRUE(caught)
END_SECTION
START_SECTION(Indexed decoder supports a binary text node larger than 10 MB)
  const std::string large_payload = std::string(13 * 1024 * 1024, ' ') + "AAAAAAAAAAA=";
  const auto array = [](const std::string& kind, const std::string& payload)
  {
    return "<binaryDataArray><cvParam accession='MS:1000523' name='64-bit float'/>"
           "<cvParam accession='" + kind + "' name='" + (kind == "MS:1000514" ? "m/z array" : "intensity array") +
           "'/><cvParam accession='MS:1000576' name='no compression'/><binary>" + payload + "</binary></binaryDataArray>";
  };
  const std::string xml = "<spectrum id='large' defaultArrayLength='1'><!-- </spectrum> --><binaryDataArrayList count='2'>" +
    array("MS:1000514", large_payload) + array("MS:1000515", "AAAAAAAAAAA=") + "</binaryDataArrayList></spectrum></spectrumList><chromatogramList>";
  MzMLSpectrumDecoder decoder;
  MSSpectrum spectrum;
  decoder.domParseSpectrum(xml, spectrum);
  TEST_EQUAL(spectrum.size(), 1)
  if (!spectrum.empty())
  {
    TEST_EQUAL(spectrum[0].getMZ(), 0.0)
    TEST_EQUAL(spectrum[0].getIntensity(), 0.0)
  }
END_SECTION
START_SECTION(Indexed decoder rejects malformed XML inside the root)
  MzMLSpectrumDecoder decoder;
  MSSpectrum spectrum;
  TEST_EXCEPTION(Exception::ParseError, decoder.domParseSpectrum("<spectrum defaultArrayLength='0'><broken></spectrum>", spectrum))
  TEST_EXCEPTION(Exception::ParseError, decoder.domParseSpectrum("<spectrum defaultArrayLength='0'>raw &-</spectrum>", spectrum))
END_SECTION
START_SECTION(Indexed decoder retains the entity amplification guard under HUGE)
  std::string dtd = "<!DOCTYPE spectrum [<!ENTITY e0 '0123456789'>";
  for (int i = 1; i <= 7; ++i)
  {
    dtd += "<!ENTITY e" + std::to_string(i) + " '";
    for (int j = 0; j < 10; ++j) { dtd += "&e" + std::to_string(i - 1) + ";"; }
    dtd += "'>";
  }
  const std::string xml = dtd + "]><spectrum defaultArrayLength='0'>&e7;</spectrum>";
  MzMLSpectrumDecoder decoder;
  MSSpectrum spectrum;
  // The same complete DTD must parse below the budget, including expansion.
  // This prevents a missing/stripped DTD from making the negative test pass.
  decoder.domParseSpectrum(dtd + "]><spectrum defaultArrayLength='0' id='&e2;'/>", spectrum);
  std::string expanded;
  for (int i = 0; i < 100; ++i) expanded += "0123456789";
  TEST_STRING_EQUAL(spectrum.getNativeID(), expanded)
  TEST_EXCEPTION(Exception::ParseError, decoder.domParseSpectrum(xml, spectrum))
END_SECTION
START_SECTION(Independent parser contexts on eight threads)
  std::vector<std::future<bool>> workers;
  for (int i = 0; i < 8; ++i)
  {
    workers.emplace_back(std::async(std::launch::async, []
    {
      for (int j = 0; j < 100; ++j)
      {
        Reader reader;
        Handler handler;
        reader.parseBuffer_("<root value='parallel'>text</root>", &handler);
        if (handler.text != u"text" || handler.attribute != u"parallel") { return false; }
      }
      return true;
    }));
  }
  for (auto& worker : workers) { TEST_TRUE(worker.get()) }
END_SECTION
START_SECTION(Native ASCII attribute lookup and UTF8 fallback)
  const std::vector<std::pair<std::u16string, std::u16string>> values{
      {u"id", u"one"}, {u"idRef", u"two"}, {u"β", u"unicode"}, {u"naïve", u"accent"}};
  const Internal::XMLAttributes attributes(values);
  TEST_TRUE(attributes.value(u"id") != nullptr)
  TEST_TRUE(attributes.value("id") == attributes.value(u"id"))
  TEST_TRUE(attributes.value("idRef") == attributes.value(u"idRef"))
  TEST_TRUE(attributes.value("idR") == nullptr)
  TEST_TRUE(attributes.value("idReference") == nullptr)
  TEST_TRUE(attributes.value("") == nullptr)
  TEST_TRUE(attributes.value("missing") == nullptr)
  TEST_TRUE(attributes.value(u"β") != nullptr)
  TEST_TRUE(attributes.value("β") == attributes.value(u"β"))
  TEST_TRUE(attributes.value("naïve") == attributes.value(u"naïve"))
  TEST_TRUE(attributes.value("βmissing") == nullptr)
END_SECTION
START_SECTION(Large attributes and deep nesting remain readable as with Xerces)
  Reader reader;
  Handler handler;
  const std::string large(10000001, 'x');
  reader.parseBuffer_("<root value='" + large + "'>" + large + "</root>", &handler);
  TEST_EQUAL(handler.attribute.size(), large.size())
  TEST_EQUAL(handler.text.size(), large.size())
  std::string deep;
  for (Size i = 0; i < 300; ++i) { deep += "<e>"; }
  for (Size i = 0; i < 300; ++i) { deep += "</e>"; }
  Handler nested;
  reader.parseBuffer_(deep, &nested);
  TEST_EQUAL(nested.starts, 300)
  // The entity amplification budget still applies.
  std::string bomb = "<!DOCTYPE r [<!ENTITY a '0123456789'>";
  std::string previous = "a";
  for (char level = 'b'; level <= 'j'; ++level)
  {
    bomb += std::string("<!ENTITY ") + level + " '" ;
    for (int i = 0; i < 10; ++i) { bomb += "&" + previous + ";"; }
    bomb += "'>";
    previous = std::string(1, level);
  }
  bomb += "]><r>&j;</r>";
  Handler exploded;
  TEST_EXCEPTION(Exception::ParseError, reader.parseBuffer_(bomb, &exploded))
END_SECTION
START_SECTION(Compressed input is read ahead without changing results or blocking early stops)
  // Several MiB of decompressed XML: more than the read-ahead queue holds.
  MSExperiment experiment;
  for (Size i = 0; i < 300; ++i)
  {
    MSSpectrum spectrum;
    spectrum.setRT(static_cast<double>(i));
    spectrum.setNativeID("scan=" + std::to_string(i));
    for (Size j = 0; j < 2000; ++j) { spectrum.emplace_back(100.0 + static_cast<double>(j) / 7.0, static_cast<float>(i + j)); }
    experiment.addSpectrum(spectrum);
  }
  std::string gz;
  NEW_TMP_FILE_EXT(gz, ".mzML.gz")
  MzMLFile().store(gz, experiment);
  MSExperiment loaded;
  MzMLFile().load(gz, loaded);
  TEST_EQUAL(loaded.size(), experiment.size())
  bool same_peaks = loaded.size() == experiment.size();
  for (Size i = 0; same_peaks && i < loaded.size(); ++i)
  {
    same_peaks = loaded[i].getNativeID() == experiment[i].getNativeID() && loaded[i].size() == experiment[i].size();
    for (Size j = 0; same_peaks && j < loaded[i].size(); ++j)
    {
      same_peaks = loaded[i][j].getMZ() == experiment[i][j].getMZ() && loaded[i][j].getIntensity() == experiment[i][j].getIntensity();
    }
  }
  TEST_TRUE(same_peaks)

  // Stopping at the first element must not wait for the rest of the stream.
  Reader reader;
  Handler stopped;
  stopped.stop = true;
  reader.parse_(gz, &stopped);
  TEST_EQUAL(stopped.starts, 1)

  // Decompression errors on the read-ahead thread reach the caller.
  Handler handler;
  TEST_EXCEPTION(Exception::ParseError, reader.parse_(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1_corrupt.bz2"), &handler))
END_SECTION
START_SECTION(File input handles empty and one-byte files without reading uninitialized magic)
  std::string filename;
  NEW_TMP_FILE(filename)
  Reader reader;
  Handler handler;
  { std::ofstream out(filename, std::ios::binary); }
  TEST_EXCEPTION(Exception::ParseError, reader.parse_(filename, &handler))
  { std::ofstream out(filename, std::ios::binary); out << "<"; }
  TEST_EXCEPTION(Exception::ParseError, reader.parse_(filename, &handler))
  { std::ofstream out(filename, std::ios::binary); out << "<root value='file'>text</root>"; }
  reader.parse_(filename, &handler);
  TEST_EQUAL(handler.text == u"text", true)
  TEST_EQUAL(handler.attribute == u"file", true)
END_SECTION
END_TEST
