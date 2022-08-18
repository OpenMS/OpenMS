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

#include "OpenMS/CHEMISTRY/EmpiricalFormula.h"
#include <OpenMS/CHEMISTRY/VarMod.h>

namespace OpenMS
{

    bool VarMod::isMergeableWith(const VarMod& otherMod)
    {
        if(    this->mod_mass_ == otherMod.mod_mass_
            && this->mod_binary_group_ == otherMod.mod_binary_group_
            && this->mod_max_current_mod_per_peptide_ == otherMod.mod_max_current_mod_per_peptide_
            && this->mod_term_distance_ == otherMod.mod_term_distance_
            && this->mod_nc_term_ == otherMod.mod_nc_term_
            && this->mod_term_specificity_ == otherMod.mod_term_specificity_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    void VarMod::merge(VarMod& modification)
    {
        this->mod_residue_.append(modification.mod_residue_);
        modification.mod_mass_ = 0.0;
    }

    String VarMod::toCometString() const
    {
        String cometParam = String(this->mod_mass_) + " " + this->mod_residue_ + " "
        + this->mod_binary_group_ + " "
        + this->mod_max_current_mod_per_peptide_ + " "
        + this->mod_term_distance_ + " "
        + this->mod_nc_term_ + " "
        + this->mod_required_ + " "
        + String(this->mod_neutral_loss_);

        return cometParam;
    }
}
