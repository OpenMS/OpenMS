// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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


#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FORMAT/FastaIterator.h>
#include <OpenMS/CHEMISTRY/EdwardsLippertIterator.h>
#include <OpenMS/FORMAT/FastaIteratorIntern.h>
#include <OpenMS/CHEMISTRY/EdwardsLippertIteratorTryptic.h>
#include <OpenMS/CHEMISTRY/TrypticIterator.h>

namespace OpenMS
{

  typedef std::pair<String, String> FASTAEntry;

  PepIterator::PepIterator()
  {
  }

  PepIterator::PepIterator(const PepIterator & /*source*/)
  {
  }

  PepIterator::~PepIterator()
  {
  }

  void PepIterator::registerChildren()
  {
    //register new products here
    Factory<PepIterator>::registerProduct(EdwardsLippertIterator::getProductName(), &EdwardsLippertIterator::create);
    Factory<PepIterator>::registerProduct(EdwardsLippertIteratorTryptic::getProductName(), &EdwardsLippertIteratorTryptic::create);
    Factory<PepIterator>::registerProduct(FastaIterator::getProductName(), &FastaIterator::create);
    Factory<PepIterator>::registerProduct(FastaIteratorIntern::getProductName(), &FastaIteratorIntern::create);
    Factory<PepIterator>::registerProduct(TrypticIterator::getProductName(), &TrypticIterator::create);
  }

}
