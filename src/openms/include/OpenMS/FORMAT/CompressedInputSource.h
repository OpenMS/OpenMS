// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <xercesc/sax/InputSource.hpp>

namespace OpenMS
{
  /**
      @brief This class is based on xercesc::LocalFileInputSource
  */
  class OPENMS_DLLAPI CompressedInputSource :
    public xercesc::InputSource
  {
public:
    ///Constructor
    CompressedInputSource(const   String & file_path, const String & header, xercesc::MemoryManager * const manager = xercesc::XMLPlatformUtils::fgMemoryManager);
    ///Constructor
    CompressedInputSource(const   XMLCh * const file_path, const String & header, xercesc::MemoryManager * const manager = xercesc::XMLPlatformUtils::fgMemoryManager);
    ///Constructor
    ~CompressedInputSource() override;

    /**
       @brief Depending on the header in the Constructor a Bzip2InputStream or a GzipInputStream object is returned
       @note InputSource interface implementation
    */
    xercesc::BinInputStream * makeStream() const override;

private:
    String head_;
    /// private CTor - not implemented
    CompressedInputSource();
    CompressedInputSource(const CompressedInputSource & source);
    CompressedInputSource & operator=(const CompressedInputSource & source);
  };

} // namespace OpenMS

