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
#include <memory>
#include <vector>

namespace OpenMS
{

namespace
{
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

  /// Restores global WARN/ERROR log levels and stderr sinks on scope exit.
  struct GlobalLogSilencer
  {
    std::string warn_level_;
    std::string err_level_;

    GlobalLogSilencer()
    {
      warn_level_ = getGlobalLogWarn().getLevel();
      getGlobalLogWarn().setLevel("FATAL");
      getGlobalLogWarn().removeAllStreams();
      err_level_ = getGlobalLogError().getLevel();
      getGlobalLogError().setLevel("FATAL");
      getGlobalLogError().removeAllStreams();
    }

    ~GlobalLogSilencer()
    {
      getGlobalLogWarn().setLevel(warn_level_);
      getGlobalLogWarn().insert(std::cerr);
      getGlobalLogError().setLevel(err_level_);
      getGlobalLogError().insert(std::cerr);
    }
  };

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

void ImzMLFile::load(const String& filename, MSExperiment& exp)
{
  exp.reset();
  loadImpl_(filename, nullptr, exp, nullptr, nullptr, false);
}

void ImzMLFile::load(const String& filename, MSImagingExperiment& exp)
{
  MSExperiment ms_exp;
  load(filename, ms_exp);

  MSImagingGeometry geom;
  buildImagingGeometry(ms_exp, geom);

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
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "imzML pixel coordinates must be >= 1",
                                      String("(") + x_imz + "," + y_imz + ") at spectrum " + i);
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
    max_x = std::max(max_x, x);
    max_y = std::max(max_y, y);
    geom.addPixel(x, y, i);
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

void ImzMLFile::load(const String& filename, Interfaces::IMSDataConsumer& consumer)
{
  MSExperiment meta_only;
  loadImpl_(filename, &consumer, meta_only, nullptr, nullptr, false);
}

void ImzMLFile::loadSpectraIndex(const String& filename, ImzMLMeta& meta, std::vector<ImzMLSpectrumIndex>& index)
{
  MSExperiment dummy;
  loadImpl_(filename, nullptr, dummy, &meta, &index, true);
}

bool ImzMLFile::isValid(const String& filename, std::ostream& os)
{
  return Internal::XMLFile::isValid(filename, os);
}

void ImzMLFile::store(const String& filename, const MSExperiment& exp) const
{
  Internal::ImzMLWriter::store(filename, exp, options_, const_cast<ImzMLFile&>(*this));
}

void ImzMLFile::loadImpl_(const String& filename,
                          Interfaces::IMSDataConsumer* consumer,
                          MSExperiment& meta_exp,
                          ImzMLMeta* out_meta,
                          std::vector<ImzMLSpectrumIndex>* out_index,
                          bool index_only)
{
  XercesPlatformGuard xerces_guard;

  const GlobalLogSilencer log_silencer;

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

  const String ibd_path = inferIbdPath_(filename);
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
      String("XML error: ") + xercesc::XMLString::transcode(e.getMessage()));
  }
  catch (const xercesc::SAXException& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
      String("SAX error: ") + xercesc::XMLString::transcode(e.getMessage()));
  }
  catch (const Exception::BaseException&)
  {
    throw;
  }
  catch (const std::exception& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, e.what());
  }

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

String ImzMLFile::inferIbdPath_(const String& imzml_path)
{
  String p = imzml_path;
  String lower = p;
  lower.toLower();

  if (lower.hasSuffix(".imzml"))
    return p.substr(0, p.size() - 6) + ".ibd";

  return p + ".ibd";
}

} // namespace OpenMS
