// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <cstdio>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{
namespace Internal
{

  class ImzMLInterceptConsumer; // defined in ImzMLHandler.cpp

  /**
    @brief SAX2 handler for imzML 1.1.0 files, extending OpenMS MzMLHandler.

    @ingroup FileIO

    imzML (HUPO-PSI v1.1.0) stores MSI data across two files:
    - @c .imzML  — mzML 1.1.0 XML with IMS ontology CV extensions
    - @c .ibd    — raw binary companion (m/z and intensity arrays)

    This handler extends @p MzMLHandler so that all standard mzML metadata
    (instrument, data processing, scan settings, RT, MS level) is parsed by the
    base class at no extra cost.  Only the ~15 IMS:* CV terms that MzMLHandler
    does not know are intercepted here.

    @b Delivery timing.  MzMLHandler batches spectrum delivery: it calls
    @p consumeSpectrum on the registered @p IMSDataConsumer only after the
    entire @c \<spectrumList\> has been parsed.  ImzMLHandler therefore records
    per-spectrum IMS state in @p spec_ims_ at @p endElement("spectrum") and
    uses a spectrum counter inside @p ImzMLInterceptConsumer to re-associate
    each delivered @p MSSpectrum with the correct IMS entry.

    @b Usage.  Do not use this class directly.  Use @p ImzMLFile::load() or
    @p OnDiscImzMLExperiment::open() instead.

    @b Inheritance diagram:
    @code
    xercesc::DefaultHandler
      └─ OpenMS::Internal::XMLHandler
            └─ OpenMS::Internal::MzMLHandler
                  └─ OpenMS::Internal::ImzMLHandler  ← this class
                        └─ ImzMLInterceptConsumer : IMSDataConsumer
    @endcode
  */
  class OPENMS_DLLAPI ImzMLHandler final : public MzMLHandler
  {
    friend class ImzMLInterceptConsumer;

  public:

    /**
      @brief Construct an ImzMLHandler.

      @param exp       PeakMap that receives @p ExperimentalSettings from MzMLHandler.
      @param filename  Path to the .imzML file (for error messages and Xerces URI).
      @param logger    ProgressLogger inherited from the owning @p ImzMLFile.
    */
    ImzMLHandler(PeakMap&              exp,
                 const std::string&    filename,
                 const ProgressLogger& logger);

    ~ImzMLHandler() override;  ///< Closes the .ibd FILE* if open.

    ImzMLHandler(const ImzMLHandler&)            = delete;
    ImzMLHandler& operator=(const ImzMLHandler&) = delete;

    /**
      @brief Open the companion .ibd file for random-access binary reads.

      Must be called before @p parse_() when peak data is needed.
      When loading metadata only (no IBD), skip this call.

      @param ibd_path  Absolute path to the .ibd file.

      @throws Exception::FileNotFound if the file cannot be opened.
    */
    void openIBD(const std::string& ibd_path);

    /// Imaging metadata accumulated during SAX parse (read after parse_() returns).
    const ImzMLMeta& getImzMLMeta() const noexcept { return meta_; }

    /// Non-const access for ImzMLFile to set the .ibd path before parse.
    ImzMLMeta&       getImzMLMeta()       noexcept { return meta_; }

    /// Per-spectrum binary index built during parse (one entry per spectrum).
    const std::vector<ImzMLSpectrumIndex>& getIndex() const noexcept { return index_; }

    /**
      @brief Wire @p ImzMLInterceptConsumer between MzMLHandler and an optional downstream consumer.

      Decodes .ibd peak data when spectra are delivered. When @p append_spectra_to_map
      is @c true, decoded spectra are also appended to the bound @p PeakMap.
    */
    void connectDecodeConsumer(Interfaces::IMSDataConsumer* downstream,
                               bool append_spectra_to_map,
                               bool decode_ibd = true);

    // ------------------------------------------------------------------
    // SAX2 overrides (MzMLHandler is called through for MS:* terms)
    // ------------------------------------------------------------------

    void startElement(const XMLCh*               uri,
                      const XMLCh*               localname,
                      const XMLCh*               qname,
                      const xercesc::Attributes& attrs) override;

    void endElement(const XMLCh* uri,
                    const XMLCh* localname,
                    const XMLCh* qname) override;

  private:

    // ------------------------------------------------------------------
    // IMS CV dispatch
    // ------------------------------------------------------------------
    void handleIMSCvParam_(const std::string& acc, const std::string& val);
    void applyRefGroup_   (const std::string& id);

    // ------------------------------------------------------------------
    // Per-binaryDataArray IMS metadata
    // ------------------------------------------------------------------
    struct ArrayMeta
    {
      ImzMLSpectrumIndex::DataType dt      { ImzMLSpectrumIndex::DataType::UNKNOWN };
      bool     is_mz  { false };
      bool     is_int { false };
      bool     is_ext { false };  ///< IMS:1000101 — external (in .ibd)
      uint64_t offset { 0 };     ///< IMS:1000102 — byte offset
      uint64_t count  { 0 };     ///< IMS:1000103 — element count

      void reset() noexcept { *this = ArrayMeta{}; }
    };

    // ------------------------------------------------------------------
    // Per-spectrum IMS state snapshot (stored at endElement("spectrum"))
    // ------------------------------------------------------------------
    struct SpecIMS
    {
      uint32_t  x {0}, y {0}, z {1};
      ArrayMeta mz_meta;
      ArrayMeta int_meta;
    };

    // ------------------------------------------------------------------
    // State
    // ------------------------------------------------------------------
    FILE*      ibd_          { nullptr };
    ImzMLMeta  meta_;

    bool       in_spectrum_  { false };
    bool       in_scan_      { false };
    bool       in_bda_       { false };
    uint32_t   cur_x_        { 0 };
    uint32_t   cur_y_        { 0 };
    uint32_t   cur_z_        { 1 };
    ArrayMeta  cur_array_;
    ArrayMeta  cur_mz_meta_;
    ArrayMeta  cur_int_meta_;

    std::vector<SpecIMS>             spec_ims_;  ///< IMS state per spectrum (document order)
    std::vector<ImzMLSpectrumIndex>  index_;     ///< Full index for OnDiscImzMLExperiment

    using CvPair = std::pair<std::string, std::string>;
    std::unordered_map<std::string, std::vector<CvPair>> ref_groups_;
    std::string cur_ref_id_;
    bool   in_ref_group_ { false };
    bool   decode_ibd_ { true };
    bool   cur_x_seen_ { false };
    bool   cur_y_seen_ { false };

    std::unique_ptr<ImzMLInterceptConsumer> decode_bridge_;
  };

} // namespace Internal
} // namespace OpenMS
