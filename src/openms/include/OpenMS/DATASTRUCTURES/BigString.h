// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_BIGSTRING_H
#define OPENMS_DATASTRUCTURES_BIGSTRING_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <vector>

namespace OpenMS
{

/**
@brief concatenates Proteins given as FASTAEntry to one big string separated by a unique character (by default $)

Concatenates the strings given as FASTAEntry separating them with a unique character and storing the headers of FASTAEntry as well as the position of separator characters. So a substring can be accessed easily and the corresponding header can be found fast by using binary search.
*/
  class OPENMS_DLLAPI BigString
  {

public:

    typedef std::pair<String, String> FASTAEntry;

    /**
    @brief constructor
    */
    BigString();

    /**
    @brief copy constructor
    */
    BigString(const BigString & bs);

    /**
    @brief destructor
    */
    virtual ~BigString();

    /**
    @brief add new string to bigString
    @param new_entry FASTAEntry to be added to big_string
    */
    void add(FASTAEntry const & new_entry);

    /**
    @brief setter for separator character by default $
    @param sep separator character
    */
    void setSeparator(const char sep);

    /**
    @brief getter for separator character
    @return separator character
    */
    char getSeparator();

    /**
    @brief returns the number of strings
    @return int with number of strings
    */
    Size size();

    /**
    @brief length of bigString
    @return int with length of the created bigString
    */
    Size length();

    /**
    @brief getPeptide from start position with given length this includes FASTAHeader
    @param entry contains the entry of the given range after calling
    @param start start index
    @param length length of desired substring
    @return FASTAEntry describing the protein
    @throw InvalidValue if a peptide is part of two different fasta entrys
    */
    void getPeptide(FASTAEntry & entry, Size start, Size length);

    /**
    @brief returns bigString
    @return const reference to bigString
    */
    const String & getBigString() const;

protected:

    /**
    @brief private function to implement binary search
    @param index
    @param start start index
    @param end end index
    @return int with index
    */
    Size getIndex_(Size index, Size start, Size end);

    /**
    @brief retrieves index of inserted protein by bigStringPosition
    @param index
    @return int with index
    */
    Size getIndex_(Size index);

    String big_string_; ///< concatenated String

    char separator_; ///< separator sign

    Size count_; ///< number of Strings added to big_string

    Size len_; ///< length of the big_string

    std::vector<Size> sep_indices_; ///< indices of separators

    std::vector<String> FASTA_header_; ///< vector with headers of FASTAEntry

  };
}
#endif // OPENMS_DATASTRUCTURES_BIGSTRING_H
