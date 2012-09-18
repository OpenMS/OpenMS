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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/PeptideHitVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  PeptideHitVisualizer::PeptideHitVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<PeptideHit>()
  {
    addLineEdit_(peptidehit_score_, "Score");
    addLineEdit_(peptidehit_charge_, "Charge");
    addLineEdit_(peptidehit_rank_, "Rank");
    addTextEdit_(peptidehit_sequence_, "Sequence");

    finishAdding_();
  }

  void PeptideHitVisualizer::update_()
  {
    peptidehit_score_->setText(String(temp_.getScore()).c_str());
    peptidehit_score_->setReadOnly(true);
    peptidehit_charge_->setText(String(temp_.getCharge()).c_str());
    peptidehit_charge_->setReadOnly(true);
    peptidehit_rank_->setText(String(temp_.getRank()).c_str());
    peptidehit_rank_->setReadOnly(true);
    peptidehit_sequence_->setText(temp_.getSequence().toString().c_str());
    peptidehit_sequence_->setReadOnly(true);
  }

  void PeptideHitVisualizer::store()
  {
    //TODO?
    (*ptr_) = temp_;
  }

  void PeptideHitVisualizer::undo_()
  {
    update_();
  }

}
