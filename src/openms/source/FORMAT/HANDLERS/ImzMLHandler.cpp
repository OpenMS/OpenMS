// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ImzMLHandler.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <xercesc/util/XMLString.hpp>

#include <cstring>
#include <limits>
#include <stdexcept>

namespace OpenMS
{
namespace Internal
{

namespace
{
  [[noreturn]] void throwImsParseError_(const std::string& file, const std::string& acc, const std::string& val, const std::string& reason)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file,
                                "Invalid IMS CV value for " + acc + " ('" + val + "'): " + reason);
  }

  uint32_t parseImsUInt32_(const std::string& file, const std::string& acc, const std::string& val)
  {
    if (val.empty())
    {
      throwImsParseError_(file, acc, val, "empty value");
    }
    try
    {
      size_t pos = 0;
      const unsigned long long v = std::stoull(val, &pos);
      if (pos != val.size())
      {
        throwImsParseError_(file, acc, val, "trailing characters");
      }
      if (v > std::numeric_limits<uint32_t>::max())
      {
        throwImsParseError_(file, acc, val, "out of range for uint32");
      }
      return static_cast<uint32_t>(v);
    }
    catch (const Exception::ParseError&)
    {
      throw;
    }
    catch (...)
    {
      throwImsParseError_(file, acc, val, "not a valid unsigned integer");
    }
  }

  uint64_t parseImsUInt64_(const std::string& file, const std::string& acc, const std::string& val)
  {
    if (val.empty())
    {
      throwImsParseError_(file, acc, val, "empty value");
    }
    if (val.find('-') != std::string::npos)
    {
      // std::stoull silently wraps a negative value to a huge unsigned number (no range
      // ceiling exists for uint64, unlike parseImsUInt32_), so reject it explicitly.
      throwImsParseError_(file, acc, val, "negative value not allowed");
    }
    try
    {
      size_t pos = 0;
      const unsigned long long v = std::stoull(val, &pos);
      if (pos != val.size())
      {
        throwImsParseError_(file, acc, val, "trailing characters");
      }
      return v;
    }
    catch (const Exception::ParseError&)
    {
      throw;
    }
    catch (...)
    {
      throwImsParseError_(file, acc, val, "not a valid unsigned integer");
    }
  }

  double parseImsDouble_(const std::string& file, const std::string& acc, const std::string& val)
  {
    if (val.empty())
    {
      throwImsParseError_(file, acc, val, "empty value");
    }
    try
    {
      size_t pos = 0;
      const double v = std::stod(val, &pos);
      if (pos != val.size())
      {
        throwImsParseError_(file, acc, val, "trailing characters");
      }
      return v;
    }
    catch (const Exception::ParseError&)
    {
      throw;
    }
    catch (...)
    {
      throwImsParseError_(file, acc, val, "not a valid floating-point number");
    }
  }
} // namespace

// ===========================================================================
// ImzMLInterceptConsumer — internal consumer wired between MzMLHandler and
// the user's IMSDataConsumer (or PeakMap collector).
//
// MzMLHandler fires consumeSpectrum() for every finalised spectrum.  This
// consumer attaches the IMS pixel coordinates and decodes the binary arrays
// from the companion .ibd file using the offsets captured by ImzMLHandler.
// ===========================================================================
class ImzMLInterceptConsumer final : public Interfaces::IMSDataConsumer
{
public:
  ImzMLInterceptConsumer(ImzMLHandler&                handler,
                         Interfaces::IMSDataConsumer* downstream)
    : handler_(handler), downstream_(downstream)
  {}

  void setExpectedSize(size_t n_spectra, size_t n_chromatograms) override
  {
    if (downstream_) downstream_->setExpectedSize(n_spectra, n_chromatograms);
  }

  void setExperimentalSettings(const ExperimentalSettings& es) override
  {
    if (downstream_) downstream_->setExperimentalSettings(es);
  }

  void consumeSpectrum(MSSpectrum& s) override
  {
    // MzMLHandler batches delivery — consumeSpectrum fires after the entire
    // <spectrumList> is parsed, so live cur_x_/y_/z_ in handler_ are 0.
    // Use the per-spectrum snapshot stored at endElement("spectrum").
    const int idx = spec_idx_++;
    const ImzMLHandler::SpecIMS* ims = nullptr;
    if (idx < static_cast<int>(handler_.spec_ims_.size()))
      ims = &handler_.spec_ims_[idx];

    // Annotate pixel coordinates as MetaValues
    if (ims)
    {
      s.setMetaValue("imzml:x", static_cast<Int>(ims->x));
      s.setMetaValue("imzml:y", static_cast<Int>(ims->y));
      s.setMetaValue("imzml:z", static_cast<Int>(ims->z));

      // Update dataset-level bounding box
      if (ims->x > handler_.meta_.max_count_x) handler_.meta_.max_count_x = ims->x;
      if (ims->y > handler_.meta_.max_count_y) handler_.meta_.max_count_y = ims->y;
      if (ims->z > handler_.meta_.max_count_z) handler_.meta_.max_count_z = ims->z;

      // Decode binary arrays from .ibd when either array is external; keep inline arrays from MzMLHandler.
      if (handler_.ibd_ && handler_.decode_ibd_ && (ims->mz_meta.is_ext || ims->int_meta.is_ext))
      {
        std::vector<double> mz_vec;
        std::vector<float>  int_vec;
        const std::string& ibd_path = handler_.meta_.ibd_file_path;
        if (ims->mz_meta.is_ext)
        {
          ImzMLBinaryIO::readMzArray(handler_.ibd_, ims->mz_meta.offset, ims->mz_meta.count,
                                     ims->mz_meta.dt, mz_vec, ibd_path);
        }
        else
        {
          mz_vec.resize(s.size());
          for (Size i = 0; i < s.size(); ++i)
          {
            mz_vec[i] = s[i].getMZ();
          }
        }
        if (ims->int_meta.is_ext)
        {
          ImzMLBinaryIO::readIntArray(handler_.ibd_, ims->int_meta.offset, ims->int_meta.count,
                                      ims->int_meta.dt, int_vec, ibd_path);
        }
        else
        {
          int_vec.resize(s.size());
          for (Size i = 0; i < s.size(); ++i)
          {
            int_vec[i] = s[i].getIntensity();
          }
        }

        if (mz_vec.size() != int_vec.size())
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path,
            "m/z and intensity array length mismatch at pixel (" + StringConversions::toString(ims->x) + "," + StringConversions::toString(ims->y) + "," + StringConversions::toString(ims->z)
            + "): mz=" + StringConversions::toString(mz_vec.size()) + " intensity=" + StringConversions::toString(int_vec.size()));
        }

        s.resize(mz_vec.size());
        for (Size i = 0; i < s.size(); ++i)
        {
          s[i].setMZ(mz_vec[i]);
          s[i].setIntensity(int_vec[i]);
        }

        // Record first-seen data types at dataset level
        if (handler_.meta_.mz_data_type.empty())
          handler_.meta_.mz_data_type  = dtStr_(ims->mz_meta.dt);
        if (handler_.meta_.int_data_type.empty())
          handler_.meta_.int_data_type = dtStr_(ims->int_meta.dt);
      }

      // Build full spectrum index entry
      ImzMLSpectrumIndex entry;
      entry.index      = idx;
      entry.x          = ims->x;
      entry.y          = ims->y;
      entry.z          = ims->z;
      entry.mz_offset  = ims->mz_meta.offset;
      entry.mz_length  = ims->mz_meta.count;
      entry.mz_type    = ims->mz_meta.dt;
      entry.int_offset = ims->int_meta.offset;
      entry.int_length = ims->int_meta.count;
      entry.int_type   = ims->int_meta.dt;
      handler_.index_.push_back(entry);
    }

    if (downstream_) downstream_->consumeSpectrum(s);
  }

  void consumeChromatogram(MSChromatogram& c) override
  {
    if (downstream_) downstream_->consumeChromatogram(c);
  }

private:
  ImzMLHandler&                handler_;
  Interfaces::IMSDataConsumer* downstream_;
  int                          spec_idx_  {0};

  static std::string dtStr_(ImzMLSpectrumIndex::DataType dt)
  {
    switch (dt)
    {
      case ImzMLSpectrumIndex::DataType::FLOAT32: return "float32";
      case ImzMLSpectrumIndex::DataType::FLOAT64: return "float64";
      case ImzMLSpectrumIndex::DataType::INT32:   return "int32";
      case ImzMLSpectrumIndex::DataType::INT64:   return "int64";
      default:                                    return "unknown";
    }
  }
};

// ===========================================================================
// ImzMLHandler — constructor / destructor
// ===========================================================================

ImzMLHandler::ImzMLHandler(PeakMap&              exp,
                            const std::string&    filename,
                            const ProgressLogger& logger)
  : MzMLHandler(exp, filename, "1.1.0", logger)
{}

ImzMLHandler::~ImzMLHandler()
{
  decode_bridge_.reset();
  if (ibd_) { fclose(ibd_); ibd_ = nullptr; }
}

void ImzMLHandler::connectDecodeConsumer(Interfaces::IMSDataConsumer* downstream,
                                         bool append_spectra_to_map,
                                         bool decode_ibd)
{
  decode_ibd_ = decode_ibd;
  decode_bridge_ = std::make_unique<ImzMLInterceptConsumer>(*this, downstream);
  setMSDataConsumer(decode_bridge_.get());
  PeakFileOptions opt = getOptions();
  opt.setAlwaysAppendData(append_spectra_to_map);
  setOptions(opt);
}

void ImzMLHandler::openIBD(const std::string& ibd_path)
{
  FILE* new_ibd = fopen(ibd_path.c_str(), "rb");
  if (!new_ibd)
  {
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path);
  }
  if (ibd_)
  {
    fclose(ibd_);
  }
  ibd_ = new_ibd;
}

// ===========================================================================
// SAX2 overrides
// ===========================================================================

void ImzMLHandler::startElement(const XMLCh*               uri,
                                 const XMLCh*               localname,
                                 const XMLCh*               qname,
                                 const xercesc::Attributes& attrs)
{
  const std::string tag = sm_.convert(qname);

  // --- Track IMS parsing state BEFORE handing off to MzMLHandler base ---

  if (tag == "spectrum")
  {
    in_spectrum_ = true;
    cur_x_ = cur_y_ = 0;
    cur_z_ = 1;
    cur_x_seen_ = false;
    cur_y_seen_ = false;
    cur_mz_meta_.reset();
    cur_int_meta_.reset();
  }
  if (tag == "scan"            && in_spectrum_) in_scan_ = true;
  if (tag == "binaryDataArray" && in_spectrum_)
  {
    in_bda_ = true;
    cur_array_.reset();
  }

  if (tag == "referenceableParamGroup")
  {
    in_ref_group_ = true;
    std::string id_om;
    optionalAttributeAsString_(id_om, attrs, "id");
    cur_ref_id_ = id_om;
    ref_groups_[cur_ref_id_]; // ensure key exists
  }
  if (tag == "referenceableParamGroupRef")
  {
    std::string ref_id_om;
    optionalAttributeAsString_(ref_id_om, attrs, "ref");
    applyRefGroup_(ref_id_om);
  }

  // Intercept <cvParam> before MzMLHandler: capture IMS:* terms
  if (tag == "cvParam")
  {
    std::string acc_om, val_om;
    optionalAttributeAsString_(acc_om, attrs, "accession");
    optionalAttributeAsString_(val_om, attrs, "value");

    if (in_ref_group_)
      ref_groups_[cur_ref_id_].emplace_back(acc_om, val_om);
    else
      handleIMSCvParam_(acc_om, val_om);
    // Fall through: MzMLHandler still processes MS:* terms below
  }

  MzMLHandler::startElement(uri, localname, qname, attrs);
}

void ImzMLHandler::endElement(const XMLCh* uri,
                               const XMLCh* localname,
                               const XMLCh* qname)
{
  const std::string tag = sm_.convert(qname);

  if (tag == "binaryDataArray" && in_spectrum_)
  {
    if      (cur_array_.is_mz)  cur_mz_meta_  = cur_array_;
    else if (cur_array_.is_int) cur_int_meta_ = cur_array_;

    // For external (.ibd) arrays the inline <binary> is an empty placeholder, but the
    // base MzMLHandler captured this array's expected length from defaultArrayLength
    // (BinaryData::size). MzMLHandlerHelper::decodeBase64Arrays() — run inside
    // MzMLHandler::populateSpectraWithData_()'s `#pragma omp parallel for` — compares
    // that length against the 0 decoded values and logs a per-array warning. Those
    // concurrent writes to the non-thread-safe LogStream corrupt the heap on Windows
    // (non-deterministic 0xC0000374). Zero the captured length so the empty inline
    // array decodes cleanly and silently; the real peaks come from the .ibd later.
    if (cur_array_.is_ext && !bin_data_.empty())
    {
      bin_data_.back().size = 0;
    }

    cur_array_.reset();
    in_bda_ = false;
  }
  if (tag == "scan")     in_scan_ = false;
  if (tag == "spectrum")
  {
    if (!cur_x_seen_ || !cur_y_seen_)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_,
                                  "imzML spectrum missing required IMS pixel coordinate (x and/or y)");
    }
    // Snapshot per-spectrum IMS state NOW — before cur_* are reset.
    // MzMLHandler delivers spectra to consumeSpectrum() after the entire
    // spectrumList is parsed, so live cur_* fields would all be 0 by then.
    //
    // Only snapshot spectra that MzMLHandler will actually deliver. skip_spectrum_
    // (protected base member, set during this spectrum's content parsing when it
    // is filtered out via PeakFileOptions) holds this spectrum's final decision
    // here — it is the same value MzMLHandler::endElement() uses below to decide
    // whether to collect/deliver the spectrum, and is only reset for the *next*
    // spectrum inside that base call. ImzMLInterceptConsumer::consumeSpectrum()
    // re-associates delivered spectra with spec_ims_ by a running counter, so an
    // unconditional push would desync coordinates and .ibd offsets as soon as any
    // spectrum is skipped.
    if (!skip_spectrum_)
    {
      SpecIMS snap;
      snap.x        = cur_x_;
      snap.y        = cur_y_;
      snap.z        = cur_z_;
      snap.mz_meta  = cur_mz_meta_;
      snap.int_meta = cur_int_meta_;
      spec_ims_.push_back(snap);
    }

    // imzML peak data lives in the companion .ibd, so the inline <binary> arrays are
    // empty placeholders while the spectrum still declares defaultArrayLength = N (the
    // external element count). Left untouched, MzMLHandler::populateSpectraWithData_()
    // would see the 0-vs-N size mismatch and, for *every* spectrum, emit "array has
    // size 0 but should have size N" / "Fixing faulty defaultArrayLength" warnings —
    // from inside its `#pragma omp parallel for`. The shared LogStreamBuf is not
    // thread-safe (one streambuf, a static line buffer, an incomplete_line_ string and
    // a message cache), so those concurrent log writes corrupt the heap (observed as a
    // non-deterministic 0xC0000374 abort in the Windows pyOpenMS wheel tests). Zeroing
    // the length here — before MzMLHandler::endElement() snapshots it for the batch —
    // makes the base treat the spectrum as having no inline peaks: no mismatch, no
    // warning, no wasted populate work. ImzMLInterceptConsumer fills the peaks from the
    // .ibd afterwards using the external offsets/counts, which are independent of this.
    if (cur_mz_meta_.is_ext || cur_int_meta_.is_ext)
    {
      default_array_length_ = 0;
    }
    in_spectrum_ = false;
  }
  if (tag == "referenceableParamGroup")
  {
    in_ref_group_ = false;
    cur_ref_id_.clear();
  }

  MzMLHandler::endElement(uri, localname, qname);

  if (tag == "spectrum")
  {
    cur_x_ = cur_y_ = 0;
    cur_z_ = 1;
    cur_mz_meta_.reset();
    cur_int_meta_.reset();
  }
}

// ===========================================================================
// IMS CV dispatch
// ===========================================================================

void ImzMLHandler::handleIMSCvParam_(const std::string& acc, const std::string& val)
{
  // Pixel coordinates (inside <scan> inside <spectrum>)
  if (in_scan_ && in_spectrum_)
  {
    if (acc == "IMS:1000050") { cur_x_ = parseImsUInt32_(file_, acc, val); cur_x_seen_ = true; return; }
    if (acc == "IMS:1000051") { cur_y_ = parseImsUInt32_(file_, acc, val); cur_y_seen_ = true; return; }
    if (acc == "IMS:1000052") { cur_z_ = parseImsUInt32_(file_, acc, val); return; }
  }

  // Inside a binaryDataArray
  if (in_bda_)
  {
    using DT = ImzMLSpectrumIndex::DataType;
    if (acc == "MS:1000521") { cur_array_.dt = DT::FLOAT32; return; }
    if (acc == "MS:1000523") { cur_array_.dt = DT::FLOAT64; return; }
    if (acc == "MS:1000519") { cur_array_.dt = DT::INT32;   return; }
    if (acc == "MS:1000522") { cur_array_.dt = DT::INT64;   return; }
    if (acc == "MS:1000514") { cur_array_.is_mz  = true;    return; }
    if (acc == "MS:1000515") { cur_array_.is_int = true;    return; }
    if (acc == "IMS:1000101"){ cur_array_.is_ext = true;    return; }
    if (acc == "IMS:1000102") { cur_array_.offset = parseImsUInt64_(file_, acc, val); return; }
    if (acc == "IMS:1000103") { cur_array_.count  = parseImsUInt64_(file_, acc, val); return; }
  }

  // Dataset-level imaging metadata
  if (acc == "IMS:1000042") { meta_.max_count_x = parseImsUInt32_(file_, acc, val); return; }
  if (acc == "IMS:1000043") { meta_.max_count_y = parseImsUInt32_(file_, acc, val); return; }
  if (acc == "IMS:1000030") { meta_.imaging_mode        = "continuous";  return; }
  if (acc == "IMS:1000031") { meta_.imaging_mode        = "processed";   return; }
  if (acc == "IMS:1000091") { meta_.ibd_sha1            = val;           return; }
  if (acc == "IMS:1000090") { meta_.ibd_md5             = val;           return; }
  if (acc == "IMS:1000046") { meta_.pixel_size_x  = parseImsDouble_(file_, acc, val); return; }
  if (acc == "IMS:1000047") { meta_.pixel_size_y  = parseImsDouble_(file_, acc, val); return; }
  if (acc == "IMS:1000080") { if (!val.empty()) meta_.uuid = val;        return; }
  if (acc == "IMS:1000044") { meta_.max_dim_x     = parseImsDouble_(file_, acc, val); return; }
  if (acc == "IMS:1000045") { meta_.max_dim_y     = parseImsDouble_(file_, acc, val); return; }
  if (acc == "MS:1000129")  { meta_.polarity            = "negative";    return; }
  if (acc == "MS:1000130")  { meta_.polarity            = "positive";    return; }
  if (acc == "IMS:1000401") { meta_.scan_pattern        = "top down";    return; }
  if (acc == "IMS:1000402") { meta_.scan_pattern        = "bottom up";   return; }
  if (acc == "IMS:1000413") { meta_.scan_direction      = "flyback";     return; }
  if (acc == "IMS:1000412") { meta_.scan_direction      = "meander";     return; }
  if (acc == "IMS:1000480") { meta_.scan_direction      = "horizontal";  return; }
  if (acc == "IMS:1000481") { meta_.scan_direction      = "vertical";    return; }
  if (acc == "IMS:1000491") { meta_.line_scan_direction = "left-right";  return; }
  if (acc == "IMS:1000492") { meta_.line_scan_direction = "right-left";  return; }
}

void ImzMLHandler::applyRefGroup_(const std::string& id)
{
  auto it = ref_groups_.find(id);
  if (it == ref_groups_.end()) return;
  for (const auto& kv : it->second)
    handleIMSCvParam_(kv.first, kv.second);
}

} // namespace Internal
} // namespace OpenMS
