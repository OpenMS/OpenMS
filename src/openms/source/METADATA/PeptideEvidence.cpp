// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

  PeptideEvidence::PeptideEvidence(const String& accession, Int start, Int end, char aa_before, char aa_after) :
      accession_(accession),
      start_(start),
      end_(end),
      aa_before_(aa_before),
      aa_after_(aa_after)
  {
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
    if (accession_ != rhs.accession_)
    {
      return accession_ < rhs.accession_;
    }
    if (start_ != rhs.start_)
    {
      return start_ < rhs.start_;
    }
    if (end_ != rhs.end_)
    {
      return end_ < rhs.end_;
    }
    if (aa_before_ != rhs.aa_before_)
    {
      return aa_before_ < rhs.aa_before_;
    }
    if (aa_after_ != rhs.aa_after_)
    {
      return aa_after_ < rhs.aa_after_;
    }
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

