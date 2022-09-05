// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/VarMod.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
    VarMod::VarMod (const ResidueModification *mod , Int binary_group, Int max_variable_mods_in_peptide)
    {

      String residues = mod->getOrigin();

      //TODO support mod-specific limit (default for now is the overall max per peptide)
      int max_current_mod_per_peptide = max_variable_mods_in_peptide;
      //TODO support term-distances?
      int term_distance = -1;
      int nc_term = 0;

      //TODO support agglomeration of Modifications to same AA. Watch out for nc_term value then.
      if (mod->getTermSpecificity() == ResidueModification::C_TERM)
      {
        if (mod->getOrigin() == 'X')
        {
          residues = "c";
        } // else stays mod.getOrigin()
        term_distance = 0;
        // Since users need to specify mods that apply to multiple residues/terms separately
        // 3 and -1 should be equal for now.
        nc_term = 3;
      }
      else if (mod->getTermSpecificity() == ResidueModification::N_TERM)
      {
        if (mod->getOrigin() == 'X')
        {
          residues = "n";
        } // else stays mod.getOrigin()
        term_distance = 0;
        // Since users need to specify mods that apply to multiple residues/terms separately
        // 2 and -1 should be equal for now.
        nc_term = 2;
      }
      else if (mod->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
      {
        if (mod->getOrigin() == 'X')
        {
          residues = "n";
        } // else stays mod.getOrigin()
        term_distance = 0;
        nc_term = 0;
      }
      else if (mod->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
      {
        if (mod->getOrigin() == 'X')
        {
          residues = "c";
        } // else stays mod.getOrigin()
        term_distance = 0;
        nc_term = 1;
      }

      this->mod_mass_ = mod->getDiffMonoMass();
      this->mod_residue_ = residues;
      this->mod_term_specificity_ = mod->getTermSpecificity();
      this->mod_binary_group_ = binary_group;
      this->mod_max_current_mod_per_peptide_ = max_current_mod_per_peptide;
      this->mod_term_distance_ = term_distance;
      this->mod_nc_term_ = nc_term;
      this->mod_required_ = false;       //TODO: support required variable mods
      this->mod_neutral_loss_ = 0.0;    // TODO: add neutral losses (from Residue or user defined?)
    }

    VarMod::~VarMod()
    {
    }

    bool VarMod::isMergeableWith(const VarMod& otherMod)
    {
        return   this->mod_mass_ == otherMod.mod_mass_
              && this->mod_binary_group_ == otherMod.mod_binary_group_
              && this->mod_max_current_mod_per_peptide_ == otherMod.mod_max_current_mod_per_peptide_
              && this->mod_term_distance_ == otherMod.mod_term_distance_;
    }

    void VarMod::merge(VarMod& modification)
    {
        this->mod_residue_.append(modification.mod_residue_);
        modification.mod_mass_ = 0.0;
    }

    String VarMod::toCometString() const
    {
        String cometParam = String(this->mod_mass_) + " " + this->mod_residue_ + " "
        + String(this->mod_binary_group_) + " "
        + String(this->mod_max_current_mod_per_peptide_) + " "
        + String(this->mod_term_distance_) + " "
        + String(this->mod_nc_term_) + " "
        + String(this->mod_required_) + " "
        + String(this->mod_neutral_loss_);

        return cometParam;
    }
}
