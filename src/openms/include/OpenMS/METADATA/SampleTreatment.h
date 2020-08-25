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
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
      @brief Base class for sample treatments (Digestion, Modification, Tagging, ...)

      Virtual base class for all sample treatments.

      The type of the treatment can be determined through the getType() method.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI SampleTreatment :
    public MetaInfoInterface
  {
public:

    /**
        @brief Constructor.

        Use a unique type string for each treatment type
    */
    explicit SampleTreatment(const String & type);

    /**
        @brief Copy constructor

        @note Do not forget to call it when you derive a class from SampleTreatment!
    */
    SampleTreatment(const SampleTreatment &) = default;

    /// Move constructor
    SampleTreatment(SampleTreatment&&) = default;

    /// Destructor
    virtual ~SampleTreatment();

    /**
        @brief Assignment operator

        @note Do not forget to call it when you derive a class from SampleTreatment!
    */
    SampleTreatment & operator=(const SampleTreatment &) = default;

    /// Move assignment operator
    SampleTreatment& operator=(SampleTreatment&&) & = default;

    /**
        @brief Equality operator

        The equality operators of derived classes also take a SampleTreatment reference as argument.
        They check the type and cast the reference to the right type if the type matches.

    @note Do not forget to call it when you derive a class from SampleTreatment!
    */
    virtual bool operator==(const SampleTreatment & rhs) const;

    /**
        @brief return the treatment type

        The type_ has to be set in the default constructor.
        It is used to determine the kind of sample treatment, when only a pointer to this base class is available.
    */
    const String & getType() const;

    /// returns the description of the sample treatment
    const String & getComment() const;

    /// sets the description of the sample treatment
    void setComment(const String & comment);

    /**
        @brief A clone method

        clone method that creates a copy and returns a pointer (base class pointer).
        Used to copy sample treatments when only a pointer to this base class is available.
    */
    virtual SampleTreatment * clone() const = 0;

protected:
    String type_;
    String comment_;

private:
    /// Default constructor hidden to force setting of a type
    SampleTreatment();

  };
} // namespace OpenMS

