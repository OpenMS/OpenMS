// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <list>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  class SampleTreatment;

  /**
      @brief Meta information about the sample

      It contains basic descriptions like name, number (i.e. order number), mass,
      volume, concentration, state and a comment.

      Additionally sample treatments like Digestion, Modification or Tagging can be added.

      A Sample can be composed of other samples.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Sample :
    public MetaInfoInterface
  {
public:
    ///state of aggregation of the sample
    enum SampleState {SAMPLENULL, SOLID, LIQUID, GAS, SOLUTION, EMULSION, SUSPENSION, SIZE_OF_SAMPLESTATE};
    /// Names of sample states
    static const std::string NamesOfSampleState[SIZE_OF_SAMPLESTATE];

    /// Default constructor
    Sample();
    /// Copy constructor
    Sample(const Sample & source);
    /// Move constructor
    Sample(Sample&&) = default;
    /// Destructor
    ~Sample();

    /// Assignment operator
    Sample & operator=(const Sample & source);
    /// Move assignment operator
    Sample& operator=(Sample&&) & = default;

    /// Equality operator
    bool operator==(const Sample & rhs) const;

    /// returns the sample name (default: "")
    const String & getName() const;
    /// sets the sample name
    void setName(const String & name);

    /// returns the sample name (default: "")
    const String & getOrganism() const;
    /// sets the sample name
    void setOrganism(const String & organism);

    /// returns the sample number (default: "")
    const String & getNumber() const;
    /// sets the sample number (e.g. sample ID)
    void setNumber(const String & number);

    /// returns the comment (default: "")
    const String & getComment() const;
    /// sets the comment (may contain newline characters)
    void setComment(const String & comment);

    /// returns the state of aggregation (default: SAMPLENULL)
    SampleState getState() const;
    /// sets the state of aggregation
    void setState(SampleState state);

    /// returns the mass (in gram) (default: 0.0)
    double getMass() const;
    /// sets the mass (in gram)
    void setMass(double mass);

    /// returns the volume (in ml) (default: 0.0)
    double getVolume() const;
    /// sets the volume (in ml)
    void setVolume(double volume);

    /// returns the concentration (in g/l) (default: 0.0)
    double getConcentration() const;
    /// sets the concentration (in g/l)
    void setConcentration(double concentration);

    /// returns a mutable reference to the vector of subsamples that were combined to create this sample
    std::vector<Sample> & getSubsamples();
    /// returns a const reference to the vector of subsamples that were combined to create this sample
    const std::vector<Sample> & getSubsamples() const;
    /// sets the vector of subsamples that were combined to create this sample
    void setSubsamples(const std::vector<Sample> & subsamples);

    /**
        @brief adds a sample treatment before the given position (default is
        the end of the list). Sample treatments are ordered in the order of
        application to the sample. If before_position is smaller than 0, the
        sample treatment is appended to the list.

        @exception Exception::IndexOverflow is thrown if the position is invalid.
    */
    void addTreatment(const SampleTreatment & treatment, Int before_position = -1);

    /**
        @brief returns a mutable reference to the sample treatment at the given position

        @exception Exception::IndexOverflow is thrown if the position is invalid.
    */
    SampleTreatment & getTreatment(UInt position);

    /**
        @brief returns a const reference to the sample treatment at the given position

        @exception Exception::IndexOverflow is thrown if the position is invalid.
    */
    const SampleTreatment & getTreatment(UInt position) const;

    /**
        @brief removes the sample treatment at the given position

        @exception Exception::IndexOverflow is thrown if the position is invalid.
    */
    void removeTreatment(UInt position);

    /// returns the number of sample treatments
    Int countTreatments() const;

protected:
    String name_;
    String number_;
    String comment_;
    String organism_;
    SampleState state_;
    double mass_;
    double volume_;
    double concentration_;
    std::vector<Sample> subsamples_;

    // note: default move constructor / assignment operator will work on this,
    // since it will move the whole vector over to the new object who will then
    // own the memory.
    std::list<SampleTreatment *> treatments_;

  };
} // namespace OpenMS

