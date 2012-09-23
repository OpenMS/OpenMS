// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SAMPLE_H
#define OPENMS_METADATA_SAMPLE_H

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

    ///default constructor
    Sample();
    ///copy constructor
    Sample(const Sample & source);
    ///desctuctor
    ~Sample();

    ///assignment operator
    Sample & operator=(const Sample & source);

    /// Equality operator
    bool operator==(const Sample & rhs) const;

    ///retuns the sample name (default: "")
    const String & getName() const;
    ///sets the sample name
    void setName(const String & name);

    ///retuns the sample name (default: "")
    const String & getOrganism() const;
    ///sets the sample name
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
    DoubleReal getMass() const;
    /// sets the mass (in gram)
    void setMass(DoubleReal mass);

    /// returns the volume (in ml) (default: 0.0)
    DoubleReal getVolume() const;
    /// sets the volume (in ml)
    void setVolume(DoubleReal volume);

    /// returns the concentration (in g/l) (default: 0.0)
    DoubleReal getConcentration() const;
    /// sets the concentration (in g/l)
    void setConcentration(DoubleReal concentration);

    /// returns a mutable reference to the vector of subsamples that were combined to create this sample
    std::vector<Sample> & getSubsamples();
    /// returns a const referenct to the vector of subsamples that were combined to create this sample
    const std::vector<Sample> & getSubsamples() const;
    /// sets the vector of subsamples that were combined to create this sample
    void setSubsamples(const std::vector<Sample> & subsamples);

    /**
        @brief adds a sample treatment before the given postion (default is the end of the list). Sample treatments are ordered in the order of application to the sample. If before_position is smaller than 0, the sample treatment is appended to the list.

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
    DoubleReal mass_;
    DoubleReal volume_;
    DoubleReal concentration_;
    std::vector<Sample> subsamples_;
    std::list<SampleTreatment *> treatments_;

  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SAMPLE_H
