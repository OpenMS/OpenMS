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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_HELPERS_H
#define OPENMS_CONCEPT_HELPERS_H

#include <vector>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  namespace Helpers 
  {

    /**
        @brief Helper function to add constness to a vector of shared pointers
    */
    template <class T>
    const std::vector<boost::shared_ptr<const T> >&
    constifyPointerVector(const std::vector<boost::shared_ptr<T> >& vec) 
    {
      return reinterpret_cast<const std::vector<boost::shared_ptr<const T> >&>(vec);
    }


    /**
      * @brief Helper comparing two pointers for equality (taking NULL into account)
    */
    template <class PtrType>
    inline bool cmpPtrSafe(const PtrType& a, const PtrType& b)
    {
       // We are not interested whether the pointers are equal but whether the
       // contents are equal
      if (a == nullptr && b == nullptr)
      {
        return true;
      }
      else if (a == nullptr || b == nullptr)
      {
        return false; // one is null the other is not
      }
      else
      {
        // compare the internal object
        return (*a == *b);
      }
    }

    /**
      * @brief Helper function to compare two pointer-containers for equality of all elements
    */
    template <class ContainerType>
    inline bool cmpPtrContainer(const ContainerType& a, const ContainerType& b)
    {
      if (a.size() != b.size()) return false;

      // check that all elements of a and b are equal using safe comparison
      // (taking NULL into account)
      for (Size i = 0; i < a.size(); i++)
      {
        if (!cmpPtrSafe(a[i], b[i]))
        {
          return false;
        }
      }
      return true;
    }

  }
}


#endif //OPENMS_CONCEPT_HELPERS_H
