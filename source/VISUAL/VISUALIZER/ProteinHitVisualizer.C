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
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/ProteinHitVisualizer.h>

//QT
#include <QtGui/QTextEdit>
#include <QtGui/QLineEdit>

#include <iostream>

using namespace std;

namespace OpenMS
{

  ProteinHitVisualizer::ProteinHitVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<ProteinHit>()
  {
    addLineEdit_(proteinhit_score_, "Score");
    addLineEdit_(proteinhit_rank_, "Rank");
    addLineEdit_(proteinhit_accession_, "Accession");
    addTextEdit_(proteinhit_sequence_, "Sequence");

    finishAdding_();
  }

  void ProteinHitVisualizer::update_()
  {
    proteinhit_score_->setText(String(temp_.getScore()).c_str());
    proteinhit_score_->setReadOnly(true);
    proteinhit_rank_->setText(String(temp_.getRank()).c_str());
    proteinhit_rank_->setReadOnly(true);
    proteinhit_accession_->setText(temp_.getAccession().c_str());
    proteinhit_accession_->setReadOnly(true);
    proteinhit_sequence_->setText(temp_.getSequence().c_str());
    proteinhit_sequence_->setReadOnly(true);
  }

  void ProteinHitVisualizer::store()
  {
    //TODO?
    (*ptr_) = temp_;
  }

  void ProteinHitVisualizer::undo_()
  {
    update_();
  }

}
