// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ImzMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLWriter.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/IMAGING/MSImagingGeometry.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLException.hpp>
#include <xercesc/sax/SAXException.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <vector>

namespace OpenMS
{

namespace
{
  // Parse an imzML UUID string (with or without dashes) into 16 raw bytes.
  // Returns false if the string is not exactly 32 hex digits. Mirrors the
  // writer's ImzMLWriter::uuidStringToBytes_ so reader and writer agree on the
  // byte order written to / read from the .ibd header.
  bool uuidStringToBytes_(const std::string& uuid, unsigned char out[16])
  {
    std::string hex = uuid;
    // Strip the separators imzML UUIDs may carry: dashes and (some writers) { } braces.
    hex.erase(std::remove_if(hex.begin(), hex.end(),
                             [](char c) { return c == '-' || c == '{' || c == '}'; }),
              hex.end());
    if (hex.size() != 32) { return false; }
    for (int i = 0; i < 16; ++i)
    {
      const std::string byte_str = hex.substr(static_cast<std::string::size_type>(i) * 2, 2);
      char* end = nullptr;
      const unsigned long value = std::strtoul(byte_str.c_str(), &end, 16);
      if (end != byte_str.c_str() + byte_str.size() || value > 255) { return false; }
      out[i] = static_cast<unsigned char>(value);
    }
    return true;
  }

  std::string bytesToHex_(const unsigned char* bytes, size_t n)
  {
    static const char* digits = "0123456789abcdef";
    std::string s;
    s.reserve(n * 2);
    for (size_t i = 0; i < n; ++i)
    {
      s.push_back(digits[bytes[i] >> 4]);
      s.push_back(digits[bytes[i] & 0x0F]);
    }
    return s;
  }

  // imzML integrity check: the .ibd must begin with the 16-byte UUID declared in
  // the .imzML XML (IMS:1000080). A mismatch means the .imzML and .ibd do not
  // belong together (or the .ibd is truncated/corrupt) — reject loudly instead of
  // silently decoding garbage offsets. If the XML carries no (or an unparsable)
  // UUID we can only warn, since some non-conformant writers omit it.
  // Advisory imzML integrity check: the .ibd should begin with the 16-byte UUID declared
  // in the .imzML XML (IMS:1000080). A mismatch usually means the .imzML and .ibd do not
  // belong together. This is reported as a loud WARNING rather than a hard error so that
  // legacy / non-conformant datasets (e.g. an .ibd that stores array data from offset 0
  // with no UUID prefix) still load; callers who require strict conformance can inspect
  // the warning log.
  void verifyIbdUuid_(const std::string& ibd_path, const std::string& xml_uuid, const std::string& imzml_path)
  {
    unsigned char expected[16];
    if (xml_uuid.empty() || !uuidStringToBytes_(xml_uuid, expected))
    {
      OPENMS_LOG_WARN << "imzML: missing or invalid universally unique identifier "
                         "(IMS:1000080) in '" << imzml_path
                      << "'; cannot verify the companion .ibd UUID header." << std::endl;
      return;
    }
    std::FILE* f = std::fopen(ibd_path.c_str(), "rb");
    if (f == nullptr)
    {
      OPENMS_LOG_WARN << "imzML: cannot open companion .ibd '" << ibd_path
                      << "' to verify its UUID header." << std::endl;
      return;
    }
    unsigned char actual[16];
    const size_t n = std::fread(actual, 1, sizeof(actual), f);
    std::fclose(f);
    if (n != sizeof(actual))
    {
      OPENMS_LOG_WARN << "imzML: companion .ibd '" << ibd_path
                      << "' is shorter than the 16-byte UUID header; skipping UUID verification."
                      << std::endl;
      return;
    }
    if (std::memcmp(expected, actual, sizeof(actual)) != 0)
    {
      OPENMS_LOG_WARN << "imzML: .ibd UUID header (" << bytesToHex_(actual, sizeof(actual))
                      << ") does not match the universally unique identifier declared in '"
                      << imzml_path << "' (" << bytesToHex_(expected, sizeof(expected))
                      << "); the .imzML and .ibd may not belong together." << std::endl;
      return;
    }
  }

  void attachImzMLMeta_(MSExperiment& exp, const ImzMLMeta& meta)
  {
    exp.setMetaValue("imzml:imaging_mode", meta.imaging_mode);
    exp.setMetaValue("imzml:ibd_path", meta.ibd_file_path);
    exp.setMetaValue("imzml:max_count_x", static_cast<UInt>(meta.max_count_x));
    exp.setMetaValue("imzml:max_count_y", static_cast<UInt>(meta.max_count_y));
    exp.setMetaValue("imzml:max_count_z", static_cast<UInt>(meta.max_count_z));
    if (meta.pixel_size_x > 0)
    {
      exp.setMetaValue("imzml:pixel_size_x", meta.pixel_size_x);
    }
    if (meta.pixel_size_y > 0)
    {
      exp.setMetaValue("imzml:pixel_size_y", meta.pixel_size_y);
    }
    if (meta.max_dim_x > 0)
    {
      exp.setMetaValue("imzml:max_dim_x", meta.max_dim_x);
    }
    if (meta.max_dim_y > 0)
    {
      exp.setMetaValue("imzml:max_dim_y", meta.max_dim_y);
    }
    if (!meta.uuid.empty())
    {
      exp.setMetaValue("imzml:uuid", meta.uuid);
    }
    if (!meta.ibd_sha1.empty())
    {
      exp.setMetaValue("imzml:ibd_sha1", meta.ibd_sha1);
    }
    if (!meta.ibd_md5.empty())
    {
      exp.setMetaValue("imzml:ibd_md5", meta.ibd_md5);
    }
    if (!meta.mz_data_type.empty())
    {
      exp.setMetaValue("imzml:mz_data_type", meta.mz_data_type);
    }
    if (!meta.int_data_type.empty())
    {
      exp.setMetaValue("imzml:int_data_type", meta.int_data_type);
    }
    if (!meta.scan_pattern.empty())
    {
      exp.setMetaValue("imzml:scan_pattern", meta.scan_pattern);
    }
    if (!meta.scan_direction.empty())
    {
      exp.setMetaValue("imzml:scan_direction", meta.scan_direction);
    }
    if (!meta.line_scan_direction.empty())
    {
      exp.setMetaValue("imzml:line_scan_direction", meta.line_scan_direction);
    }
    if (!meta.polarity.empty())
    {
      exp.setMetaValue("imzml:polarity", meta.polarity);
    }
  }

  /// Balance Xerces Initialize()/Terminate() for each load invocation.
  struct XercesPlatformGuard
  {
    XercesPlatformGuard()
    {
      xercesc::XMLPlatformUtils::Initialize();
    }

    ~XercesPlatformGuard()
    {
      xercesc::XMLPlatformUtils::Terminate();
    }
  };
} // namespace

ImzMLFile::ImzMLFile() :
  Internal::XMLFile("/SCHEMAS/mzML_1_10.xsd", "1.1.0")
{
}

PeakFileOptions& ImzMLFile::getOptions()             { return options_; }
const PeakFileOptions& ImzMLFile::getOptions() const { return options_; }
void ImzMLFile::setOptions(const PeakFileOptions& o) { options_ = o; }

void ImzMLFile::load(const std::string& filename, MSImagingExperiment& exp)
{
  // imzML is a mass spectrometry imaging format, so MSImagingExperiment is the
  // canonical in-memory target (cf. BrukerTimsImagingFile). Decode the spectra into a
  // plain MSExperiment, derive the pixel geometry, and hand both to the imaging model.
  // Callers that only need the flat spectra can use exp.getMSExperiment() afterwards.
  MSExperiment ms_exp;
  ms_exp.reset();
  ImzMLMeta meta;
  std::vector<ImzMLSpectrumIndex> index;
  loadImpl_(filename, nullptr, ms_exp, &meta, &index, false);

  // Build the geometry directly from the parsed per-spectrum index (the source of truth),
  // not by reading imzml:x/y MetaValues back off the spectra. index[i] corresponds to
  // spectrum i in ms_exp (both driven by the same delivery order).
  MSImagingGeometry geom;
  buildImagingGeometry(index, meta, geom);

  exp.setMSExperiment(std::move(ms_exp));
  exp.setGeometry(std::move(geom));
  exp.validate();
}

void ImzMLFile::buildImagingGeometry(const MSExperiment& exp, MSImagingGeometry& geom)
{
  geom.clear();

  UInt width = 0;
  UInt height = 0;
  if (exp.metaValueExists("imzml:max_count_x"))
  {
    width = static_cast<UInt>(exp.getMetaValue("imzml:max_count_x"));
  }
  if (exp.metaValueExists("imzml:max_count_y"))
  {
    height = static_cast<UInt>(exp.getMetaValue("imzml:max_count_y"));
  }
  if (width > 0 && height > 0)
  {
    geom.setDimensions(width, height);
  }

  UInt max_x = 0;
  UInt max_y = 0;
  const Size n_spectra = exp.getNrSpectra();
  for (Size i = 0; i < n_spectra; ++i)
  {
    const MSSpectrum& spec = exp[i];
    if (!spec.metaValueExists("imzml:x") || !spec.metaValueExists("imzml:y"))
    {
      continue;
    }

    const Int x_imz = spec.getMetaValue("imzml:x");
    const Int y_imz = spec.getMetaValue("imzml:y");
    if (x_imz < 1 || y_imz < 1)
    {
      // imzML coordinates are 1-based; a <1 value is non-conformant. Warn and skip the
      // spectrum (it stays accessible by index) rather than aborting the whole load.
      OPENMS_LOG_WARN << "imzML: pixel coordinates must be >= 1; skipping spectrum " << i
                      << " at (" << OpenMS::StringConversions::toString(x_imz) << ","
                      << OpenMS::StringConversions::toString(y_imz) << ")." << std::endl;
      continue;
    }

    Int z_imz = 1;
    if (spec.metaValueExists("imzml:z"))
    {
      z_imz = spec.getMetaValue("imzml:z");
    }
    if (z_imz != 1)
    {
      continue;
    }

    const UInt x = static_cast<UInt>(x_imz - 1);
    const UInt y = static_cast<UInt>(y_imz - 1);
    // Soften the geometry's strict in-bounds rule to a warning: a pixel beyond the
    // declared grid (IMS:1000042/043) is dropped from the geometry (the spectrum stays
    // accessible by index) instead of aborting the whole load.
    if (width > 0 && height > 0 && (x >= width || y >= height))
    {
      OPENMS_LOG_WARN << "imzML: pixel (" << OpenMS::StringConversions::toString(x_imz) << ","
                      << OpenMS::StringConversions::toString(y_imz) << ") at spectrum " << i
                      << " lies outside the declared " << width << "x" << height
                      << " grid; excluding it from the imaging geometry." << std::endl;
      continue;
    }
    max_x = std::max(max_x, x);
    max_y = std::max(max_y, y);
    geom.addPixel(x, y, i); // still throws Exception::InvalidValue on duplicate coordinates
  }

  if (width == 0 && (max_x > 0 || geom.getNumberOfPixels() > 0))
  {
    width = max_x + 1;
  }
  if (height == 0 && (max_y > 0 || geom.getNumberOfPixels() > 0))
  {
    height = max_y + 1;
  }
  if (width > 0 && height > 0
      && (geom.getWidth() != width || geom.getHeight() != height))
  {
    geom.setDimensions(width, height);
  }

  if (exp.metaValueExists("imzml:pixel_size_x") && exp.metaValueExists("imzml:pixel_size_y"))
  {
    const double px = static_cast<double>(exp.getMetaValue("imzml:pixel_size_x"));
    const double py = static_cast<double>(exp.getMetaValue("imzml:pixel_size_y"));
    geom.setPixelSize(px, py, "micrometer");
  }
}

void ImzMLFile::buildImagingGeometry(const std::vector<ImzMLSpectrumIndex>& index,
                                     const ImzMLMeta& meta,
                                     MSImagingGeometry& geom)
{
  // Source-of-truth geometry builder: coordinates come straight from the parsed index
  // (1-based imzML), not from imzml:x/y MetaValues. Shared by ImzMLFile::load(MSImagingExperiment&)
  // and OnDiscImzMLExperiment so the in-memory and on-disc paths cannot diverge.
  geom.clear();

  UInt width = meta.max_count_x;
  UInt height = meta.max_count_y;
  if (width > 0 && height > 0)
  {
    geom.setDimensions(width, height);
  }

  UInt max_x = 0;
  UInt max_y = 0;
  for (Size i = 0; i < index.size(); ++i)
  {
    const ImzMLSpectrumIndex& e = index[i];
    if (e.z != 1) // MSImagingGeometry is 2D: only the first plane is mapped
    {
      continue;
    }
    if (e.x < 1 || e.y < 1)
    {
      // imzML coordinates are 1-based; a <1 value is non-conformant. Warn and skip the
      // spectrum (it stays reachable by linear index) rather than aborting the load.
      OPENMS_LOG_WARN << "imzML: pixel coordinates must be >= 1; skipping spectrum " << i
                      << " at (" << OpenMS::StringConversions::toString(e.x) << ","
                      << OpenMS::StringConversions::toString(e.y) << ")." << std::endl;
      continue;
    }
    const UInt x = static_cast<UInt>(e.x - 1);
    const UInt y = static_cast<UInt>(e.y - 1);
    // Soften the geometry's strict in-bounds rule to a warning: a pixel beyond the declared
    // grid is dropped from the geometry (the spectrum stays reachable by index).
    if (width > 0 && height > 0 && (x >= width || y >= height))
    {
      OPENMS_LOG_WARN << "imzML: pixel (" << OpenMS::StringConversions::toString(e.x) << ","
                      << OpenMS::StringConversions::toString(e.y) << ") at spectrum " << i
                      << " lies outside the declared " << width << "x" << height
                      << " grid; excluding it from the imaging geometry." << std::endl;
      continue;
    }
    max_x = std::max(max_x, x);
    max_y = std::max(max_y, y);
    geom.addPixel(x, y, i); // throws Exception::InvalidValue on duplicate coordinates
  }

  if (width == 0 && (max_x > 0 || geom.getNumberOfPixels() > 0))
  {
    width = max_x + 1;
  }
  if (height == 0 && (max_y > 0 || geom.getNumberOfPixels() > 0))
  {
    height = max_y + 1;
  }
  if (width > 0 && height > 0 && (geom.getWidth() != width || geom.getHeight() != height))
  {
    geom.setDimensions(width, height);
  }

  if (meta.pixel_size_x > 0 && meta.pixel_size_y > 0)
  {
    geom.setPixelSize(meta.pixel_size_x, meta.pixel_size_y, "micrometer");
  }
}

void ImzMLFile::load(const std::string& filename, Interfaces::IMSDataConsumer& consumer)
{
  MSExperiment meta_only;
  loadImpl_(filename, &consumer, meta_only, nullptr, nullptr, false);
}

void ImzMLFile::loadSpectraIndex(const std::string& filename, ImzMLMeta& meta, std::vector<ImzMLSpectrumIndex>& index, const std::string& ibd_path)
{
  MSExperiment dummy;
  loadImpl_(filename, nullptr, dummy, &meta, &index, true, ibd_path);
}

bool ImzMLFile::isValid(const std::string& filename, std::ostream& os)
{
  return Internal::XMLFile::isValid(filename, os);
}

void ImzMLFile::store(const std::string& filename, const MSExperiment& exp) const
{
  Internal::ImzMLWriter::store(filename, exp, options_, const_cast<ImzMLFile&>(*this));
}

void ImzMLFile::store(const std::string& filename, const MSImagingExperiment& exp) const
{
  // Geometry-driven store: the MSImagingGeometry is the source of truth for pixel
  // coordinates and grid dimensions, so this works for any MSImagingExperiment — including
  // ones built without imzml:x/y MetaValues (e.g. from BrukerTimsImagingFile). We synthesize
  // those MetaValues from the geometry onto a copy of the experiment and reuse the writer.
  const MSImagingGeometry& geom = exp.getGeometry();
  MSExperiment out = exp.getMSExperiment();

  for (const MSImagingGeometry::Pixel& px : geom.getPixels())
  {
    if (px.spectrum_index >= out.getNrSpectra())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "imaging geometry references a spectrum index outside the experiment",
        OpenMS::StringConversions::toString(px.spectrum_index));
    }
    MSSpectrum& s = out[px.spectrum_index];
    s.setMetaValue("imzml:x", static_cast<Int>(px.x) + 1); // 0-based geometry -> 1-based imzML
    s.setMetaValue("imzml:y", static_cast<Int>(px.y) + 1);
    if (!s.metaValueExists("imzml:z"))
    {
      s.setMetaValue("imzml:z", 1);
    }
  }

  if (geom.getWidth() > 0)  { out.setMetaValue("imzml:max_count_x", static_cast<UInt>(geom.getWidth())); }
  if (geom.getHeight() > 0) { out.setMetaValue("imzml:max_count_y", static_cast<UInt>(geom.getHeight())); }
  if (geom.getPixelSizeX() > 0) { out.setMetaValue("imzml:pixel_size_x", geom.getPixelSizeX()); }
  if (geom.getPixelSizeY() > 0) { out.setMetaValue("imzml:pixel_size_y", geom.getPixelSizeY()); }

  store(filename, out);
}

void ImzMLFile::loadImpl_(const std::string& filename,
                          Interfaces::IMSDataConsumer* consumer,
                          MSExperiment& meta_exp,
                          ImzMLMeta* out_meta,
                          std::vector<ImzMLSpectrumIndex>* out_index,
                          bool index_only,
                          const std::string& ibd_path_override)
{
  XercesPlatformGuard xerces_guard;

  PeakMap& peak_map = static_cast<PeakMap&>(meta_exp);
  ProgressLogger logger;
  logger.setLogType(getLogType());

  Internal::ImzMLHandler handler(peak_map, filename, logger);

  PeakFileOptions handler_opts = options_;
  if (index_only)
  {
    handler_opts.setFillData(false);
    handler_opts.setAlwaysAppendData(false);
  }
  else if (consumer == nullptr)
  {
    handler_opts.setAlwaysAppendData(true);
  }
  handler.setOptions(handler_opts);

  const bool append_to_map = !index_only && (consumer == nullptr);
  handler.connectDecodeConsumer(consumer, append_to_map, !index_only);

  // Honor a caller-supplied .ibd path (e.g. OnDiscImzMLExperiment::open(imzml, ibd_override))
  // so the index load AND the UUID integrity check target the file we will actually read from,
  // not an inferred sibling that may be missing or stale.
  const std::string ibd_path = ibd_path_override.empty() ? inferIbdPath_(filename) : ibd_path_override;
  handler.openIBD(ibd_path);
  handler.getImzMLMeta().ibd_file_path = ibd_path;

  Internal::StringManager sm;
  xercesc::LocalFileInputSource src(sm.convert(filename).c_str());

  std::unique_ptr<xercesc::SAX2XMLReader> reader{xercesc::XMLReaderFactory::createXMLReader()};
  reader->setFeature(xercesc::XMLUni::fgSAX2CoreValidation, false);
  reader->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces, true);
  reader->setFeature(xercesc::XMLUni::fgXercesLoadExternalDTD, false);
  reader->setFeature(xercesc::XMLUni::fgXercesDisableDefaultEntityResolution, true);
  reader->setFeature(xercesc::XMLUni::fgXercesSchema, false);
  reader->setFeature(xercesc::XMLUni::fgXercesSchemaFullChecking, false);
  reader->setContentHandler(&handler);
  reader->setErrorHandler(&handler);

  try
  {
    reader->parse(src);
  }
  catch (const Internal::XMLHandler::EndParsingSoftly&)
  {
  }
  catch (const xercesc::XMLException& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
      std::string("XML error: ") + Internal::StringManager::convert(e.getMessage()));
  }
  catch (const xercesc::SAXException& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
      std::string("SAX error: ") + Internal::StringManager::convert(e.getMessage()));
  }
  catch (const Exception::BaseException&)
  {
    throw;
  }
  catch (const std::exception& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, e.what());
  }

  // The XML UUID (IMS:1000080) is only known after parsing; verify it against the
  // 16-byte header of the companion .ibd now. Covers every read path: the in-memory
  // loads and the on-disc reader (OnDiscImzMLExperiment::open() routes through
  // loadSpectraIndex() -> loadImpl_()).
  verifyIbdUuid_(ibd_path, handler.getImzMLMeta().uuid, filename);

  if (out_meta)
  {
    *out_meta = handler.getImzMLMeta();
  }
  if (out_index)
  {
    *out_index = handler.getIndex();
  }

  if (!index_only)
  {
    attachImzMLMeta_(meta_exp, handler.getImzMLMeta());
  }
}

std::string ImzMLFile::inferIbdPath_(const std::string& imzml_path)
{
  std::string p = imzml_path;
  std::string lower = p;
  StringUtils::toLower(lower);

  if (StringUtils::hasSuffix(lower, ".imzml"))
    return p.substr(0, p.size() - 6) + ".ibd";

  return p + ".ibd";
}

} // namespace OpenMS
