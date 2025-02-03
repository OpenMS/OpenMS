// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

