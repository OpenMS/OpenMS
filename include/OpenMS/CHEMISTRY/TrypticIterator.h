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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
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
    virtual ~TrypticIterator();

    /**
    @brief * operator for getting the value of the iterator
    @return FASTAEntry with specific candidate
    @throw InvalidIterator if iterator has not been initialized
    */
    virtual FASTAEntry operator*();

    /**
    @brief opperator ++ for postincrement
    @return Reference to PepIterator
    @throw InvalidIterator if iterator has not been initialized
    */
    virtual PepIterator & operator++();

    /**
    @brief opperator ++ for preincrement
    @return pointer to PepIterator
    @throw InvalidIterator if iterator has not been initialized
    */
    virtual PepIterator * operator++(int i);

    /**
    @brief setter for fasta file
    @param f String with fasta file location
    @throw FileNotFound if file could not be found
    */
    virtual void setFastaFile(const String & f);

    /**
    @brief getter for FASTA file
    @return String with file location
    */
    virtual String getFastaFile();

    /**
    @brief setter for tolerance
    @throw NotImplemented because its not available for tryptic iterator
    */
    virtual void setTolerance(DoubleReal)
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

    /**
    @brief getter for tolerance
    @return tolerance
    @throw NotImplemented because its not available for tryptic iterator
    */
    virtual DoubleReal getTolerance()
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

    /**
    @brief setter for spectrum
    @throw NotImplemented because its not available for tryptic iterator
    */
    virtual void setSpectrum(const std::vector<DoubleReal> &)
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

    /**
    @brief getter for spectrum
    @return the used spectrum
    @throw NotImplemented because its not available for tryptic iterator
    */
    virtual const std::vector<DoubleReal> & getSpectrum()
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

    /**
    @brief initializing iterator
    @return true if everything was ok
    @throw throws InvalidIterator if begin iterator is not valid
    */
    virtual bool begin();

    /**
    @brief indicates whether iterator is at end
    @return true if iterator is at end
    @see hasNext
    */
    virtual bool isAtEnd();

    /**
    @brief indicated if a digesting enzyme will cut at this position
    @return true if digenting enzym cuts the sequence
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
    @return poiter to new object
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
