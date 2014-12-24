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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/CVReference.h>

#include <OpenMS/DATASTRUCTURES/String.h>

using namespace std;

namespace OpenMS
{
  // CV reference implementation
  CVReference::CVReference()
  {
  }

  CVReference::~CVReference()
  {
  }

  CVReference::CVReference(const CVReference & rhs) :
    name_(rhs.name_),
    identifier_(rhs.identifier_)
  {
  }

  CVReference & CVReference::operator=(const CVReference & rhs)
  {
    if (this != &rhs)
    {
      name_ = rhs.name_;
      identifier_ = rhs.identifier_;
    }
    return *this;
  }

  bool CVReference::operator==(const CVReference & rhs) const
  {
    return name_ == rhs.name_ && identifier_ == rhs.identifier_;
  }

  bool CVReference::operator!=(const CVReference & rhs) const
  {
    return !(*this == rhs);
  }

  void CVReference::setName(const String & name)
  {
    name_ = name;
  }

  const String & CVReference::getName() const
  {
    return name_;
  }

  void CVReference::setIdentifier(const String & identifier)
  {
    identifier_ = identifier;
  }

  const String & CVReference::getIdentifier() const
  {
    return identifier_;
  }

} // namespace OpenMS
