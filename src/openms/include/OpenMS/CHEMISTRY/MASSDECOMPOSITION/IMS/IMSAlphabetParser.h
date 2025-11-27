// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <fstream>
#include <istream>
#include <map>
#include <string>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  namespace ims
  {

    /**
      @brief An abstract templatized parser to load the data that is used to initialize @c Alphabet objects.

      @c AlphabetParser reads the input source, which is given as a template parameter @c InputSource , by
      @c load (const std::string& fname) function where @c fname is the source name.
      Loaded data can be retrieved by calling @c getElements().

      @see Alphabet
    */
    template <typename AlphabetElementType = double,
              typename Container = std::map<std::string, AlphabetElementType>,
              typename InputSource = std::istream>
    class IMSAlphabetParser
    {
public:
      /**
        Type of data to be loaded.
       */
      typedef Container ContainerType;

      /**
        Loads the data from the InputSource with the name @c fname.
        If there is an error occurred while reading data from InputSource,
        @c IOException is thrown.

        @param fname The name of the input source.
       */
      void load(const std::string & fname);

      /**
        Gets the data that was loaded.

        @return The data.
       */
      virtual ContainerType & getElements() = 0;

      /**
        Parses the the given input source @c is .

        @param is The InputSource
       */
      virtual void parse(InputSource & is) = 0;

      /// Destructor
      virtual ~IMSAlphabetParser() {}
    };

    template <typename AlphabetElementType, typename Container, typename InputSource>
    void IMSAlphabetParser<AlphabetElementType, Container, InputSource>::load(const std::string & fname)
    {
      std::ifstream ifs(fname.c_str());
      if (!ifs)
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, fname);
      }
      this->parse(ifs);
    }

  } // namespace ims
} // namespace OpenMS

