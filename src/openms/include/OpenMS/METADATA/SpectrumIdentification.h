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
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/IdentificationHit.h>

namespace OpenMS
{
  /**
    @brief Represents a object which can store the information of an analysisXML instance

        //@todo docu (Andreas)

        @ingroup Metadata
  */
  class OPENMS_DLLAPI SpectrumIdentification :
    public MetaInfoInterface
  {
public:

    /// @name constructors,destructors,assignment operator
    //@{
    /// Default constructor
    SpectrumIdentification() = default;
    /// Destructor
    virtual ~SpectrumIdentification();
    /// Copy constructor
    SpectrumIdentification(const SpectrumIdentification &) = default;
    /// Move constructor
    SpectrumIdentification(SpectrumIdentification&&) = default;
    /// Assignment operator
    SpectrumIdentification & operator=(const SpectrumIdentification &) = default;
    /// Move assignment operator
    SpectrumIdentification& operator=(SpectrumIdentification&&) & = default;
    /// Equality operator
    bool operator==(const SpectrumIdentification & rhs) const;
    /// Inequality operator
    bool operator!=(const SpectrumIdentification & rhs) const;
    //@}

    // @name Accessors
    //@{
    /// sets the identification hits of this spectrum identification (corresponds to single peptide hit in the list)
    void setHits(const std::vector<IdentificationHit> & hits);

    /// adds a single identification hit to the hits
    void addHit(const IdentificationHit & hit);

    /// returns the identification hits of this spectrum identification
    const std::vector<IdentificationHit> & getHits() const;
    //@}

protected:

    String id_; ///< Identifier
    std::vector<IdentificationHit> hits_; ///< Single peptide hits
  };

} //namespace OpenMS
