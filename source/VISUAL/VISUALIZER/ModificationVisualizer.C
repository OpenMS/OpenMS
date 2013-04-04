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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/ModificationVisualizer.h>

//QT
#include <QtGui/QTextEdit>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>

#include <iostream>
#include <vector>

using namespace std;

namespace OpenMS
{

  ModificationVisualizer::ModificationVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Modification>()
  {
    addLabel_("Modify Modification information");
    addSeparator_();
    addLineEdit_(treatmenttype_, "Treatment type");
    addTextEdit_(treatmentcomment_, "Comment");
    addLineEdit_(modificationname_, "Reagent name");
    addDoubleLineEdit_(modificationmass_, "Mass change");

    addComboBox_(modificationspecificity_, "Specificity Type");
    addLineEdit_(modificationAA_, "Affected Amino Acids");

    finishAdding_();
  }

  void ModificationVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(modificationspecificity_, &temp_.NamesOfSpecificityType[temp_.getSpecificityType()], 1);
    }
    else
    {
      fillComboBox_(modificationspecificity_, temp_.NamesOfSpecificityType, Modification::SIZE_OF_SPECIFICITYTYPE);
      modificationspecificity_->setCurrentIndex(temp_.getSpecificityType());
    }
    treatmenttype_->setText(temp_.getType().c_str());
    treatmenttype_->setReadOnly(true);
    treatmentcomment_->setText(temp_.getComment().c_str());
    modificationname_->setText(temp_.getReagentName().c_str());
    modificationmass_->setText(String(temp_.getMass()).c_str());
    modificationAA_->setText(temp_.getAffectedAminoAcids().c_str());
  }

  void ModificationVisualizer::store()
  {
    try
    {
      ptr_->setComment(treatmentcomment_->toPlainText());
      ptr_->setReagentName(modificationname_->text());
      ptr_->setMass(modificationmass_->text().toFloat());
      ptr_->setSpecificityType((Modification::SpecificityType)modificationspecificity_->currentIndex());
      ptr_->setAffectedAminoAcids(modificationAA_->text());
      temp_ = (*ptr_);
    }
    catch (exception & e)
    {
      std::cout << "Error while trying to store the new modification data. " << e.what() << endl;
    }
  }

  void ModificationVisualizer::undo_()
  {
    update_();
  }

}
