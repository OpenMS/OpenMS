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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

namespace OpenMS
{

  const int PeptideEvidence::UNKNOWN_POSITION = -1;
  const int PeptideEvidence::N_TERMINAL_POSITION = 0;
  const char PeptideEvidence::UNKNOWN_AA = 'X';
  const char PeptideEvidence::N_TERMINAL_AA = '[';
  const char PeptideEvidence::C_TERMINAL_AA = ']';

  PeptideEvidence::PeptideEvidence()
   : accession_(),
     start_(UNKNOWN_POSITION),
     end_(UNKNOWN_POSITION),
     aa_before_(UNKNOWN_AA),
     aa_after_(UNKNOWN_AA)
  {
  }

  PeptideEvidence::PeptideEvidence(const String& accession, Int start, Int end, char aa_before, char aa_after)
    : accession_(accession),
      start_(start),
      end_(end),
      aa_before_(aa_before),
      aa_after_(aa_after)
  {
  }

  PeptideEvidence::PeptideEvidence(const PeptideEvidence& rhs)
  {
    accession_ = rhs.accession_;
    start_ = rhs.start_;
    end_ = rhs.end_;
    aa_before_ = rhs.aa_before_;
    aa_after_ = rhs.aa_after_;
  }

  PeptideEvidence& PeptideEvidence::operator=(const PeptideEvidence& rhs)
  {
    accession_ = rhs.accession_;
    start_ = rhs.start_;
    end_ = rhs.end_;
    aa_before_ = rhs.aa_before_;
    aa_after_ = rhs.aa_after_;
    return *this;
  }

  bool PeptideEvidence::operator==(const PeptideEvidence& rhs) const
  {
    return accession_ == rhs.accession_ &&
           start_ == rhs.start_ &&
           end_ == rhs.end_ &&
           aa_before_ == rhs.aa_before_ &&
           aa_after_ == rhs.aa_after_;
  }

  bool PeptideEvidence::operator<(const PeptideEvidence& rhs) const
  {
    if (accession_ != rhs.accession_) return accession_ < rhs.accession_;
    if (start_ != rhs.start_) return start_ < rhs.start_;
    if (end_ != rhs.end_) return end_ < rhs.end_;
    if (aa_before_ != rhs.aa_before_) return aa_before_ < rhs.aa_before_;
    if (aa_after_ != rhs.aa_after_) return aa_after_ < rhs.aa_after_;
    return false;
  }

  
  bool PeptideEvidence::operator!=(const PeptideEvidence& rhs) const
  {
    return !operator==(rhs);
  }
  
  bool PeptideEvidence::hasValidLimits() const
  {
    return !(
      getStart() == UNKNOWN_POSITION ||
      getEnd() == UNKNOWN_POSITION ||
      getEnd() == N_TERMINAL_POSITION);
  }

  void PeptideEvidence::setProteinAccession(const String& s)
  {
    accession_ = s;
  }

  const String& PeptideEvidence::getProteinAccession() const
  {
    return accession_;
  }

  void PeptideEvidence::setStart(const Int a)
  {
    start_ = a;
  }

  Int PeptideEvidence::getStart() const
  {
    return start_;
  }

  void PeptideEvidence::setEnd(const Int a)
  {
    end_ = a;
  }

  Int PeptideEvidence::getEnd() const
  {
    return end_;
  }

  void PeptideEvidence::setAABefore(const char acid)
  {
    aa_before_ = acid;
  }

  char PeptideEvidence::getAABefore() const
  {
    return aa_before_;
  }

  void PeptideEvidence::setAAAfter(const char acid)
  {
    aa_after_ = acid;
  }

  char PeptideEvidence::getAAAfter() const
  {
    return aa_after_;
  }

} // namespace OpenMS
