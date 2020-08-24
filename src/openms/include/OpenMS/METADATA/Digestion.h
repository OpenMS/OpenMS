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
#include <OpenMS/METADATA/SampleTreatment.h>

namespace OpenMS
{
  /**
      @brief Meta information about digestion of a sample

      Representation of a digestion.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Digestion :
    public SampleTreatment
  {
public:
    /// Default constructor
    Digestion();
    /// Copy constructor
    Digestion(const Digestion &) = default;
    /// Move constructor
    Digestion(Digestion&&) = default;
    /// Destructor
    ~Digestion() override;

    /// Assignment operator
    Digestion & operator=(const Digestion &) = default;
    /// Move assignment operator
    Digestion& operator=(Digestion&&) & = default;

    /**
      @brief Equality operator

      Although this operator takes a reference to a SampleTreatment as argument
      it tests for the equality of Tagging instances!
    */
    bool operator==(const SampleTreatment & rhs) const override;

    /// clone method. See SampleTreatment
    SampleTreatment * clone() const override;

    /// returns the enzyme name (default is "")
    const String & getEnzyme() const;
    /// sets the enzyme name
    void setEnzyme(const String & enzyme);

    /// returns the digestion time in minutes (default is 0.0)
    double getDigestionTime() const;
    /// sets the digestion time in minutes
    void setDigestionTime(double digestion_time);

    /// return the temperature during digestion in degree C (default is 0.0)
    double getTemperature() const;
    /// sets the temperature during digestion in degree C
    void setTemperature(double temperature);

    /// returns the pH value (default is 0.0)
    double getPh() const;
    /// sets the pH value
    void setPh(double ph);

protected:
    String enzyme_;
    double digestion_time_;
    double temperature_;
    double ph_;
  };
} // namespace OpenMS

