// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IDENTIFICATION_H
#define OPENMS_METADATA_IDENTIFICATION_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/SpectrumIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Represents a object which can store the information of an analysisXML instance

        //@todo docu (Andreas)

        @ingroup Metadata
  */
  class OPENMS_DLLAPI Identification :
    public MetaInfoInterface
  {
public:

    /// @name constructors,destructors,assignment operator
    //@{
    /// default constructor
    Identification();
    /// destructor
    virtual ~Identification();
    /// copy constructor
    Identification(const Identification & source);
    /// assignment operator
    Identification & operator=(const Identification & source);
    /// Equality operator
    bool operator==(const Identification & rhs) const;
    /// Inequality operator
    bool operator!=(const Identification & rhs) const;
    //@}

    /// @name Accessors
    //@{
    /// sets the date and time the file was written
    void setCreationDate(const DateTime & date);

    /// returns the date and time the file was created
    const DateTime & getCreationDate() const;

    /// sets the spectrum identifications
    void setSpectrumIdentifications(const std::vector<SpectrumIdentification> & ids);

    /// adds a spectrum identification
    void addSpectrumIdentification(const SpectrumIdentification & id);

    /// returns the spectrum identifications stored
    const std::vector<SpectrumIdentification> & getSpectrumIdentifications() const;
    //@}
protected:

    String id_;                                     ///< Identifier
    DateTime creation_date_;            ///< Date and time the search was performed
    std::vector<SpectrumIdentification> spectrum_identifications_;

  };

} //namespace OpenMS
#endif // OPENMS_METADATA_IDENTIFICATION_H
