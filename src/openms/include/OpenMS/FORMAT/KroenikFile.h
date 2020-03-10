// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  class String;

  /**
    @brief File adapter for Kroenik (HardKloer sibling) files.

    The first line is the header and contains the column names:<br>
    File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications

    Every subsequent line is a feature.

    All properties in the file are converted to Feature properties, whereas "First Scan", "Last Scan", "Num of Scans" and "Modifications" are stored as
    metavalues with the following names "FirstScan", "LastScan", "NumOfScans" and "AveragineModifications".

    The width in m/z of the overall convex hull of each feature is set to 3 Th in lack of a value provided by the Kroenik file.

    @note Kroenik files are Tab (\\t) separated files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI KroenikFile
  {
public:
    /// Default constructor
    KroenikFile();
    /// Destructor
    virtual ~KroenikFile();

    /**
      @brief Loads a Kroenik file into a featureXML.

      The content of the file is stored in @p features.

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, FeatureMap& feature_map);

    /**
      @brief Stores a featureXML as a Kroenik file.

      NOT IMPLEMENTED

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    template <typename SpectrumType>
    void store(const String& filename, const SpectrumType& spectrum) const
    {
      std::cerr << "Store() for KroenikFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

  };
} // namespace OpenMS

