// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_DIGESTIONENZYME_H
#define OPENMS_CHEMISTRY_DIGESTIONENZYME_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <iosfwd>
#include <set>
#include <vector>

namespace OpenMS
{
  /**
     @ingroup Chemistry

     @brief Abstract base class for digestion enzymes
  */
  class OPENMS_DLLAPI DigestionEnzyme
  {
  public:

    /** @name Constructors
    */
    //@{
    /// copy constructor
    DigestionEnzyme(const DigestionEnzyme& enzyme);

    /// detailed constructor
    explicit DigestionEnzyme(const String& name,
                             const String& cleavage_regex,
                             const std::set<String>& synonyms = std::set<String>(),
                             String regex_description = "");

    /// destructor
    virtual ~DigestionEnzyme();
    //@}

    /** @name Assignment
     */
    //@{
    /// assignment operator
    DigestionEnzyme& operator=(const DigestionEnzyme& enzyme);
    //@}

    /** Accessors
    */
    //@{
    /// sets the name of the enzyme
    void setName(const String& name);

    /// returns the name of the enzyme
    String getName() const;

    /// sets the synonyms
    void setSynonyms(const std::set<String>& synonyms);

    /// adds a synonym
    void addSynonym(const String& synonym);

    /// returns the synonyms
    const std::set<String>& getSynonyms() const;

    /// sets the cleavage regex
    void setRegEx(const String& cleavage_regex);

    /// returns the cleavage regex
    String getRegEx() const;

    /// sets the regex description
    void setRegExDescription(const String& value);

    /// returns the regex description
    String getRegExDescription() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const DigestionEnzyme& enzyme) const;

    /// inequality operator
    bool operator!=(const DigestionEnzyme& enzyme) const;

    /// equality operator for regex
    bool operator==(const String& cleavage_regex) const;

    /// equality operator for regex
    bool operator!=(const String& cleavage_regex) const;

    /// order operator
    bool operator<(const DigestionEnzyme& enzyme) const;
    //@}

    /**
       @brief Set the value of a member variable based on an entry from an input file

       Returns whether the key was recognized and the value set successfully.
    */
    virtual bool setValueFromFile(const String& key, const String& value);

    /// ostream iterator to write the enzyme to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzyme& enzyme);

  protected:
    /// default constructor
    DigestionEnzyme();

    // basic
    String name_;

    String cleavage_regex_;

    std::set<String> synonyms_;

    String regex_description_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzyme& enzyme);
}

#endif
