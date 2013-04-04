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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Identification.h>

using namespace std;

namespace OpenMS
{

  Identification::Identification() :
    MetaInfoInterface(),
    id_(),
    creation_date_()
  {
  }

  Identification::Identification(const Identification & rhs) :
    MetaInfoInterface(rhs),
    id_(rhs.id_),
    creation_date_(rhs.creation_date_),
    spectrum_identifications_(rhs.spectrum_identifications_)
  {
  }

  Identification::~Identification()
  {
  }

  Identification & Identification::operator=(const Identification & rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }

    MetaInfoInterface::operator=(rhs);
    id_ = rhs.id_;
    creation_date_ = rhs.creation_date_;
    spectrum_identifications_ = rhs.spectrum_identifications_;

    return *this;
  }

  // Equality operator
  bool Identification::operator==(const Identification & rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && creation_date_ == rhs.creation_date_
           && spectrum_identifications_ == rhs.spectrum_identifications_;
  }

  // Inequality operator
  bool Identification::operator!=(const Identification & rhs) const
  {
    return !(*this == rhs);
  }

  void Identification::setCreationDate(const DateTime & date)
  {
    creation_date_ = date;
  }

  const DateTime & Identification::getCreationDate() const
  {
    return creation_date_;
  }

  void Identification::setSpectrumIdentifications(const vector<SpectrumIdentification> & ids)
  {
    spectrum_identifications_ = ids;
  }

  void Identification::addSpectrumIdentification(const SpectrumIdentification & id)
  {
    spectrum_identifications_.push_back(id);
  }

  const vector<SpectrumIdentification> & Identification::getSpectrumIdentifications() const
  {
    return spectrum_identifications_;
  }

} // namespace OpenMS
