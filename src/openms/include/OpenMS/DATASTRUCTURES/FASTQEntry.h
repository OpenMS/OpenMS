// --------------------------------------------------------------------------
//                       OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Authors: Marie Hoffmann $
// --------------------------------------------------------------------------

/**
  @brief FASTQ entry type (identifier, description, sequence, quality)

  The first string corresponds to the identifier that is written after 
  the @ in the FASTQ file. The part after the first whitespace is stored 
  as a description and the text from the next line until the next linebreak is stored
  as a sequence string. A newline starting with + (and optionally an identifier) is then 
  followed by a quality score of the same length as the sequence.
*/


#ifndef OPENMS_DATASTRUCTURES_FASTQENTRY_H
#define OPENMS_DATASTRUCTURES_FASTQENTRY_H

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <algorithm>

struct FASTQEntry
{
  seqan::CharString identifier;
  seqan::CharString description;
  seqan::CharString sequence;
  seqan::CharString quality;

  FASTQEntry() :
    identifier(""),
    description(""),
    sequence(""),
    quality("")
  {
  }

  FASTQEntry(seqan::CharString id, seqan::CharString desc, seqan::CharString seq, seqan::CharString qual) :
    identifier(id),
    description(desc),
    sequence(seq),
    quality(qual)
  {
  }

  bool operator==(const FASTQEntry& rhs) const
  {
    return identifier == rhs.identifier
           && description == rhs.description
           && sequence == rhs.sequence
           && quality == rhs.quality;
  }
  
  // Illumina 1.8+ Phred+33 with score ranges ['!', 'J'] corresponding to [0, 41] 
  std::vector<int> qual2phred()
  {
    std::vector<int> ps;
    std::transform(begin(this->quality), end(this->quality), std::back_inserter(ps), std::bind2nd(std::minus<int>(), '!'));
    return ps;
  }
  
};

#endif  // OPENMS_DATASTRUCTURES_FASTQENTRY_H