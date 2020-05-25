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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
      @brief Description of the meta data arrays of MSSpectrum.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI MetaInfoDescription :
    public MetaInfoInterface
  {
public:
    /// Constructor
    MetaInfoDescription() = default;
    /// Copy constructor
    MetaInfoDescription(const MetaInfoDescription &) = default;
    /// Move constructor
    MetaInfoDescription(MetaInfoDescription&&) = default;
    /// Destructor
    ~MetaInfoDescription();

    /// Assignment operator
    MetaInfoDescription & operator=(const MetaInfoDescription &) = default;
    /// Move assignment operator
    MetaInfoDescription& operator=(MetaInfoDescription&&) & = default;

    /// Equality operator
    bool operator==(const MetaInfoDescription & rhs) const;

    /// returns the name of the peak annotations
    const String & getName() const;
    /// sets the name of the peak annotations
    void setName(const String & name);

    /// returns a const reference to the description of the applied processing
    const std::vector<ConstDataProcessingPtr> & getDataProcessing() const;
    /// returns a mutable reference to the description of the applied processing
    std::vector<DataProcessingPtr> & getDataProcessing();
    /// sets the description of the applied processing
    void setDataProcessing(const std::vector<DataProcessingPtr> & data_processing);

protected:
    String comment_;
    String name_;
    std::vector<DataProcessingPtr> data_processing_;
  };
} // namespace OpenMS

