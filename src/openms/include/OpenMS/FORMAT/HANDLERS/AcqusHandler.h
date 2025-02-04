// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Guillaume Belz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <map>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Read-only acqus File handler for XMass Analysis.

      acqus File contains meta data about calibration (conversion for time to mz ratio),
      instrument specification and acquisition method.

      @note Do not use this class directly. It is only needed for XMassFile.
    */
    class OPENMS_DLLAPI AcqusHandler
    {
public:
      /**
        @brief Constructor with filename.

        Open acqus File as stream and import params.

        @param filename to acqus File.

        @exception Exception::FileNotFound is thrown if the file could not be opened.
        @exception Exception::ConversionError is thrown if error conversion from String to calibration param.
      */
      explicit AcqusHandler(const String & filename);

      /// Destructor
      virtual ~AcqusHandler();

      /// Conversion from index to MZ ratio using internal calibration params
      double getPosition(Size index) const;

      /// Read param as string
      String getParam(const String & param);

      /// Get size of spectrum
      Size getSize() const;

private:
      /// Private default constructor
      AcqusHandler();

      /// Map for params saving
      std::map<String, String> params_;

      /**@name Internal params for calibration */
      //@{
      double dw_;
      Size delay_;
      double ml1_;
      double ml2_;
      double ml3_;
      Size td_;
      //@}
    };
  }   // namespace Internal
} // namespace OpenMS

