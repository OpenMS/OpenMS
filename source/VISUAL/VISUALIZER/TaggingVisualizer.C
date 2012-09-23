// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/TaggingVisualizer.h>

#include <QtGui/QComboBox>
#include <QtGui/QTextEdit>
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>

#include <iostream>

using namespace std;

namespace OpenMS
{

  TaggingVisualizer::TaggingVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Tagging>()
  {
    addLabel_("Modify Tagging information");
    addSeparator_();
    addLineEdit_(treatmenttype_, "Treatment type");
    addTextEdit_(treatmentcomment_, "Comment");
    addLineEdit_(modificationname_, "Reagent name");
    addDoubleLineEdit_(modificationmass_, "Mass");

    addComboBox_(modificationspecificity_, "Specificity Type");
    addLineEdit_(modificationAA_, "Affected Amino Acids");

    addDoubleLineEdit_(taggingmass_shift_, "Mass_Shift");
    addComboBox_(taggingvariant_, "Variant");

    finishAdding_();
  }

  void TaggingVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(modificationspecificity_, &temp_.NamesOfSpecificityType[temp_.getSpecificityType()], 1);
      fillComboBox_(taggingvariant_, &temp_.NamesOfIsotopeVariant[temp_.getVariant()], 1);
    }
    else
    {
      fillComboBox_(modificationspecificity_, temp_.NamesOfSpecificityType, Tagging::SIZE_OF_SPECIFICITYTYPE);
      fillComboBox_(taggingvariant_, temp_.NamesOfIsotopeVariant, Tagging::SIZE_OF_ISOTOPEVARIANT);
      modificationspecificity_->setCurrentIndex(temp_.getSpecificityType());
      taggingvariant_->setCurrentIndex(temp_.getVariant());
    }
    treatmenttype_->setText(temp_.getType().c_str());
    treatmenttype_->setReadOnly(true);
    treatmentcomment_->setText(temp_.getComment().c_str());
    modificationname_->setText(temp_.getReagentName().c_str());
    modificationmass_->setText(String(temp_.getMass()).c_str());

    modificationAA_->setText(temp_.getAffectedAminoAcids().c_str());
    taggingmass_shift_->setText(String(temp_.getMassShift()).c_str());
  }

  void TaggingVisualizer::store()
  {
    ptr_->setComment(treatmentcomment_->toPlainText());
    ptr_->setReagentName(modificationname_->text());
    ptr_->setMass(modificationmass_->text().toDouble());
    ptr_->setSpecificityType((Modification::SpecificityType)modificationspecificity_->currentIndex());
    ptr_->setAffectedAminoAcids(modificationAA_->text());
    ptr_->setMassShift(taggingmass_shift_->text().toFloat());
    ptr_->setVariant((Tagging::IsotopeVariant)taggingvariant_->currentIndex());

    temp_ = (*ptr_);
  }

  void TaggingVisualizer::undo_()
  {
    update_();
  }

}
