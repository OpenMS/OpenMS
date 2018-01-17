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
// $Maintainer: Timo Sachsenberg $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_TRYPTICITERATOR_H
#define OPENMS_CHEMISTRY_TRYPTICITERATOR_H

#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>

namespace OpenMS
{
/**
@brief finds all tryptic Peptides with every missed cleavage

*/
  class OPENMS_DLLAPI TrypticIterator :
    public PepIterator
  {

public:

    typedef std::pair<String, String> FASTAEntry;
    /**
    @brief Constructor
    */
    TrypticIterator();
    /**
    @brief Copy Constructor
    */
    TrypticIterator(const TrypticIterator &);
    /**
    @brief Destructor
    */
    ~TrypticIterator() override;

    /**
    @brief * operator for getting the value of the iterator
    @return FASTAEntry with specific candidate
    @throw InvalidIterator if iterator has not been initialized
    */
    FASTAEntry operator*() override;

    /**
    @brief operator ++ for post-increment
    @return Reference to PepIterator
    @throw InvalidIterator if iterator has not been initialized
    */
    PepIterator & operator++() override;

    /**
    @brief operator ++ for pre-increment
    @return pointer to PepIterator
    @throw InvalidIterator if iterator has not been initialized
    */
    PepIterator * operator++(int i) override;

    /**
    @brief setter for fasta file
    @param f String with fasta file location
    @throw FileNotFound if file could not be found
    */
    void setFastaFile(const String & f) override;

    /**
    @brief getter for FASTA file
    @return String with file location
    */
    String getFastaFile() override;

    /**
    @brief setter for tolerance
    @throw NotImplemented because its not available for tryptic iterator
    */
    void setTolerance(double) override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
    @brief getter for tolerance
    @return tolerance
    @throw NotImplemented because its not available for tryptic iterator
    */
    double getTolerance() override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
    @brief setter for spectrum
    @throw NotImplemented because its not available for tryptic iterator
    */
    void setSpectrum(const std::vector<double> &) override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
    @brief getter for spectrum
    @return the used spectrum
    @throw NotImplemented because its not available for tryptic iterator
    */
    const std::vector<double> & getSpectrum() override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
    @brief initializing iterator
    @return true if everything was ok
    @throw throws InvalidIterator if begin iterator is not valid
    */
    bool begin() override;

    /**
    @brief indicates whether iterator is at end
    @return true if iterator is at end
    @see hasNext
    */
    bool isAtEnd() override;

    /**
    @brief indicated if a digesting enzyme will cut at this position
    @return true if digesting enzyme cuts the sequence
    */
    virtual bool isDigestingEnd(char aa1, char aa2);

    /**
    @brief needed by Factory
    @return const string name of class
    */
    static const String getProductName()
    {
      return "TrypticIterator";
    }

    /**
    @brief needed by Factory
    @return pointer to new object
    */
    static PepIterator * create()
    {
      return new TrypticIterator;
    }

protected:
    /**
    @brief getting the next candidate
    @return string with next sequence
    */
    virtual std::string next_();

    /**
    @brief indicates if there will be a next element
    @return true if iterator has more elements
    */
    bool hasNext_();

    /**
    @brief finds the next starting position where a digesting enzyme will cut the sequence
    */
    void goToNextAA_();

    String f_file_;     ///< fasta file location

    std::string actual_pep_;     ///< actual peptide

    bool is_at_end_;     ///< indicates if iterator is at end

    PepIterator * f_iterator_;    ///< FastaIterator

    FASTAEntry f_entry_;     ///< actual fasta entry

    unsigned int b_, e_;    ///< to ints representing a position within the actual string (b = begin, e = end)



  };

}
#endif //OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATOR_H
