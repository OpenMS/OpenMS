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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace std;

namespace OpenMS
{
  ModificationDefinition::ModificationDefinition() :
    mod_(nullptr),
    fixed_modification_(true),
    max_occurrences_(0)
  {
  }

  ModificationDefinition::ModificationDefinition(const ModificationDefinition& rhs) :
    mod_(rhs.mod_),
    fixed_modification_(rhs.fixed_modification_),
    max_occurrences_(rhs.max_occurrences_)
  {
  }

  ModificationDefinition::ModificationDefinition(const String& mod, bool fixed, UInt max_occur) :
    mod_(nullptr),
    fixed_modification_(fixed),
    max_occurrences_(max_occur)
  {
    setModification(mod);
  }

  ModificationDefinition::ModificationDefinition(const ResidueModification& mod, bool fixed, UInt max_occur) :
    mod_(&mod),
    fixed_modification_(fixed),
    max_occurrences_(max_occur)
  {
  }

  ModificationDefinition& ModificationDefinition::operator=(const ModificationDefinition& rhs)
  {
    if (this != &rhs)
    {
      mod_ = rhs.mod_;
      fixed_modification_ = rhs.fixed_modification_;
      max_occurrences_ = rhs.max_occurrences_;
    }
    return *this;
  }

  bool ModificationDefinition::operator==(const ModificationDefinition& rhs) const
  {
    return mod_ == rhs.mod_ &&
           fixed_modification_ == rhs.fixed_modification_ &&
           max_occurrences_ == rhs.max_occurrences_;
  }

  bool ModificationDefinition::operator!=(const ModificationDefinition& rhs) const
  {
    return !(*this == rhs);
  }

  ModificationDefinition::~ModificationDefinition()
  {
  }

  bool ModificationDefinition::operator<(const ModificationDefinition& rhs) const
  {
    return this->getModificationName() < rhs.getModificationName();
  }

  void ModificationDefinition::setFixedModification(bool fixed_mod)
  {
    fixed_modification_ = fixed_mod;
  }

  bool ModificationDefinition::isFixedModification() const
  {
    return fixed_modification_;
  }

  void ModificationDefinition::setModification(const String& modification)
  {
    //cerr << "setModification(" << modification << ")" << endl;
    mod_ = &ModificationsDB::getInstance()->getModification(modification);
    //cerr << "setModification: id=" << mod_->getId() << ", full_id=" << mod_->getFullId() << ", UniMod=" << mod_->getUniModAccession() << ", origin=" << mod_->getOrigin() << ", PSI-MOD=" << mod_->getPSIMODAccession() << endl;
  }

  const ResidueModification& ModificationDefinition::getModification() const
  {
    if (!mod_)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "No modification defined", nullptr);
    }
    return *mod_;
  }

  String ModificationDefinition::getModificationName() const
  {
    if (mod_ != nullptr)
    {
      return mod_->getFullId();
    }
    return "";
  }

  void ModificationDefinition::setMaxOccurrences(UInt max_occurrences)
  {
    max_occurrences_ = max_occurrences;
  }

  UInt ModificationDefinition::getMaxOccurrences() const
  {
    return max_occurrences_;
  }

} // namespace OpenMS
