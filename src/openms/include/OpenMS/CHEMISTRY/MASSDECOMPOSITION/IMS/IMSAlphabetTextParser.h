// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetParser.h>

namespace OpenMS
{

  namespace ims
  {

    /**
      @brief Implements abstract @c AlphabetParser to read data from the plain text format.

      @c AlphabetTextParser parses the data source using overridden @c parse(std::istream&)
      and stores the parsed data permanently. That can be retrieved by @c getElements() function.
    */
    class OPENMS_DLLAPI IMSAlphabetTextParser :
      public IMSAlphabetParser<>
    {
private:
      /**
        The parsed data.
      */
      ContainerType elements_;
public:
      /**
        Gets the parsed data.

        @return The parsed data.
      */
      ContainerType & getElements() override { return elements_; }

      /**
        Parses the input stream @c is.

        @param is The input stream to be parsed
      */
      void parse(std::istream & is) override;
    };

  } // namespace ims
} // namespace OpenMS

