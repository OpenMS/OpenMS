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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DoubleList.h>
using namespace std;

namespace OpenMS
{

  DoubleList::DoubleList()
  {
  }

  DoubleList::DoubleList(const DoubleList & rhs) :
    vector<DoubleReal>(rhs)
  {

  }

  DoubleList::DoubleList(const vector<DoubleReal> & rhs) :
    vector<DoubleReal>(rhs)
  {

  }

  DoubleList::DoubleList(const vector<Real> & rhs)
  {
    this->resize(rhs.size());
    for (Size i = 0; i < rhs.size(); ++i)
    {
      (*this)[i] = (DoubleReal)rhs[i];
    }
  }

  DoubleList & DoubleList::operator=(const DoubleList & rhs)
  {
    vector<DoubleReal>::operator=(rhs);
    return *this;
  }

  DoubleList & DoubleList::operator=(const vector<Real> & rhs)
  {
    this->resize(rhs.size());
    for (Size i = 0; i < rhs.size(); ++i)
    {
      (*this)[i] = (DoubleReal)rhs[i];
    }
    return *this;
  }

  DoubleList & DoubleList::operator=(const vector<DoubleReal> & rhs)
  {
    vector<DoubleReal>::operator=(rhs);
    return *this;
  }

  DoubleList DoubleList::create(const String & list)
  {
    DoubleList ret;
    vector<String> out;
    list.split(',', out);
    ret.resize(out.size());
    for (Size i = 0; i < out.size(); ++i)
    {
      ret[i] = out[i].toDouble();
    }
    return ret;
  }

  DoubleList DoubleList::create(const StringList & list)
  {
    DoubleList ret;
    for (UInt i = 0; i < list.size(); ++i)
    {
      ret.push_back(list[i].toDouble());
    }
    return ret;
  }

  bool DoubleList::contains(DoubleReal s, DoubleReal tolerance) const
  {
    for (Size i = 0; i < this->size(); ++i)
    {
      if (std::fabs(this->operator[](i) - s) < tolerance)
        return true;
    }
    return false;
  }

  // ----------------- Output operator ----------------------

  ostream & operator<<(std::ostream & os, const DoubleList & p)
  {
    os << "[";
    if (p.size() > 0)
    {
      os << precisionWrapper(p[0]);
    }

    for (Size i = 1; i < p.size(); ++i)
    {
      os << ", " << precisionWrapper(p[i]);
    }
    os << "]";
    return os;
  }

} // namespace OpenMS
