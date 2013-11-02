// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifndef OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATORTRYPTIC_H
#define OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATORTRYPTIC_H

#include <OpenMS/CHEMISTRY/EdwardsLippertIterator.h>

namespace OpenMS
{
/**
@brief EdwardsLippertIterator that only retrieves tryptic sequences
*/
  class OPENMS_DLLAPI EdwardsLippertIteratorTryptic :
    public EdwardsLippertIterator
  {
public:

    /// default constructor
    EdwardsLippertIteratorTryptic();

    /// copy constructor
    EdwardsLippertIteratorTryptic(const EdwardsLippertIteratorTryptic & rhs);

    /// destructor
    virtual ~EdwardsLippertIteratorTryptic();

    /// assignment operator
    EdwardsLippertIteratorTryptic & operator=(const EdwardsLippertIteratorTryptic & rhs);

    /**
    @brief indicates if trypsin will cat between the two amino acids
    @param aa1 first amino acid
    @param aa2 second amino acid
    */
    virtual bool isDigestingEnd(char aa1, char aa2);

    /**
    @brief needed by Factory
    @return const string name of class
    */
    static const String getProductName()
    {
      return "EdwardsLippertIteratorTryptic";
    }

    /**
    @brief needed by Factory
    @return pointer to new object
    */
    static PepIterator * create()
    {
      return new EdwardsLippertIteratorTryptic;
    }

  };

}

#endif //OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATORTRYPTIC_H
