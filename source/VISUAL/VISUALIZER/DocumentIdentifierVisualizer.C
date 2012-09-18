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

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/DocumentIdentifierVisualizer.h>
#include <OpenMS/FORMAT/FileHandler.h>

//QT
#include <QtGui/QLineEdit>

using namespace std;

namespace OpenMS
{

  DocumentIdentifierVisualizer::DocumentIdentifierVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<DocumentIdentifier>()
  {
    addLabel_("Modify DocumentIdentifier information");
    addSeparator_();
    addLineEdit_(identifier_, "Identifier");
    addSeparator_();
    addLineEdit_(file_path_, "Loaded from file");
    addLineEdit_(file_type_, "File type");
    finishAdding_();
  }

  void DocumentIdentifierVisualizer::update_()
  {
    identifier_->setText(temp_.getIdentifier().c_str());
    file_path_->setText(temp_.getLoadedFilePath().c_str());
    file_type_->setText(FileHandler::typeToName(temp_.getLoadedFileType()).c_str());
    file_path_->setReadOnly(true);
    file_type_->setReadOnly(true);
  }

  void DocumentIdentifierVisualizer::store()
  {
    ptr_->setIdentifier(identifier_->text());

    temp_ = (*ptr_);
  }

  void DocumentIdentifierVisualizer::undo_()
  {
    update_();
  }

}
