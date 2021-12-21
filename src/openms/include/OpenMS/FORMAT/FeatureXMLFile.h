// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <iosfwd>

namespace OpenMS
{
  class Feature;
  class FeatureMap;

  /**
    @brief This class provides Input/Output functionality for feature maps

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    @todo Take care that unique ids are assigned properly by TOPP tools before
    calling FeatureXMLFile::store().  There will be a message on OPENMS_LOG_INFO but
    we will make no attempt to fix the problem in this class.  (all developers)

    @note This format will eventually be replaced by the HUPO-PSI AnalysisXML
    (mzIdentML and mzQuantML) formats!

    @ingroup FileIO
  */
  class OPENMS_DLLAPI FeatureXMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {

public:

    /** @name Constructors and Destructor */
    //@{
    ///Default constructor
    FeatureXMLFile();
    ///Destructor
    ~FeatureXMLFile() override;
    //@}

    /**
        @brief loads the file with name @p filename into @p map and calls updateRanges().

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, FeatureMap& feature_map);

    Size loadSize(const String& filename);

    /**
        @brief stores the map @p feature_map in file with name @p filename.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const FeatureMap& feature_map);

    /// Mutable access to the options for loading/storing
    FeatureFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const FeatureFileOptions& getOptions() const;

    /// setter for options for loading/storing
    void setOptions(const FeatureFileOptions&);

protected:

    /// Options that can be set
    FeatureFileOptions options_;

  };

} // namespace OpenMS

