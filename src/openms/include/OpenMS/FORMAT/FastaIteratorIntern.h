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

#ifndef OPENMS_FORMAT_FASTAITERATORINTERN_H
#define OPENMS_FORMAT_FASTAITERATORINTERN_H

#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <vector>

namespace OpenMS
{

  /**
      @brief Iterator for a FASTA file

      In comparison to FastaIterator the FASTA file will be loaded first and stored to RAM, while the FastaIterator just iterates over the FASTA file without loading it completely to memory.

      @see FastaIterator
  */
  class OPENMS_DLLAPI FastaIteratorIntern :
    public PepIterator
  {

public:

    typedef std::pair<String, String> FASTAEntry;

    /**
    @brief constructor
    */
    FastaIteratorIntern();

    /**
    @brief copy constructor
    */
    FastaIteratorIntern(const FastaIteratorIntern &);

    /**
    @brief constructor
    */
    ~FastaIteratorIntern() override;

    /**
    @brief * Operator for dereferencing of iterator
    @return FASTEEntry iterator at actual position
    @throw Exception::InvalidIterator if iterator has not been initialized
    */
    FASTAEntry operator*() override;

    /**
    @brief ++ Operator for the iterator
    @return reference to PepIterator
    @throw Exception::InvalidIterator if iterator has not been initialized
    */
    PepIterator & operator++() override;

    /**
    @brief ++ Operator for the iterator
    @return pointer to PepIterator
    @throw Exception::InvalidIterator if iterator has not been initialized
    */
    PepIterator * operator++(int i) override;

    /**
    @brief setter for FASTA file
    @param f const String reference representing file location
    @throw Exception::FileNotFound
    @throw Exception::ParseError
    */
    void setFastaFile(const String & f) override;

    /**
    @brief getter for FASTA file
    @return String with file location
    */
    String getFastaFile() override;

    /**
    @brief setter for spectrum
    @note note available for FastaIterator
    @throw Exception::NotImplemented
    */
    void setSpectrum(const std::vector<double> & /*spec*/) override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
    @brief getter for spectrum
    @note note available for FastaIterator
    @throw Exception::NotImplemented
    */
    const std::vector<double> & getSpectrum() override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
    @brief setter for tolerance
    @note note available for FastaIterator
    @throw Exception::NotImplemented
    */
    void setTolerance(double /* t */) override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
    @brief getter for tolerance
    @note note available for FastaIterator
    @return tolerance
    @throw Exception::NotImplemented
    */
    double getTolerance() override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
    @brief initializing of iterator
    @return true if everything went right
    @throw Exception::InvalidIterator if fastaFile was not set
    */
    bool begin() override;

    /**
    @brief indicates whether iterator is at end
    @return bool true if iterator is at end
    */
    bool isAtEnd() override;

    /**
    @brief needed by Factory
    @return const string name of class
    */
    static const std::string getProductName()
    {
      return "FastaIteratorIntern";
    }

    /**
    @brief needed by Factory
    @return pointer to new object
    */
    static PepIterator * create()
    {
      return new FastaIteratorIntern;
    }

protected:

    String fasta_file_;         ///< location of the fasta file

    std::vector<FASTAEntry> entrys_;          ///< content of fasta file

    std::vector<FASTAEntry>::iterator it_;          ///< iterator over fasta file content

  };

}

#endif //OpenMS/FORMAT/FastaIteratorIntern
