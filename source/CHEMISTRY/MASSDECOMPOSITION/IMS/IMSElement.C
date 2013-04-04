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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>

namespace OpenMS
{

  namespace ims
  {

    /**
      @note Value for electron mass is taken from
      @link www.mcelwee.net/html/table_of_physical_constants.html
    */
    const IMSElement::mass_type IMSElement::ELECTRON_MASS_IN_U = 0.00054858;

    IMSElement & IMSElement::operator=(const IMSElement & element)
    {
      // if one doesn't assign object to itself,
      // assign all object elements to the elements of the given object
      if (this != &element)
      {
        name_ = element.name_;
        sequence_ = element.sequence_;
        isotopes_ = element.isotopes_;
      }
      return *this;
    }

    bool IMSElement::operator==(const IMSElement & element) const
    {
      return this == &element ||
             (name_ == element.name_ &&
              sequence_ == element.sequence_ &&
              isotopes_ == element.isotopes_);
    }

    bool IMSElement::operator!=(const IMSElement & element) const
    {
      return !this->operator==(element);
    }

    std::ostream & operator<<(std::ostream & os, const IMSElement & element)
    {
      os << "name:\t" << element.getName() << "\nsequence:\t" << element.getSequence()
      << "\nisotope distribution:\n" << element.getIsotopeDistribution() << '\n';
      return os;
    }

  } // namespace ims
} // namespace OpenMS
