// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/SwathTabWidget.h>
#include <ui_SwathTabWidget.h>

#include <OpenMS/VISUAL/DIALOGS/PythonModuleRequirement.h>


using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    SwathTabWidget::SwathTabWidget(QWidget *parent) :
        QTabWidget(parent),
        ui(new Ui::SwathTabWidget)
    {
        ui->setupUi(this);
        auto py_selector = ui->config->findChild<PythonSelector*>("py_selector", Qt::FindDirectChildrenOnly);

        auto py_pyprophet = ui->config->findChild<PythonModuleRequirement*>("py_pyprophet", Qt::FindDirectChildrenOnly);
        py_pyprophet->setRequiredModules({"pyprophet", "stats"});
        py_pyprophet->setFreeText("In order to run PyProphet after OpenSWATH, the above modules need to be installed\n" \
                                  "Once they are available, the 'pyProphet' tab will become active and configurable.");
        py_pyprophet->setTitle("External: PyProphet tool");
        connect(py_selector, &PythonSelector::valueChanged, py_pyprophet, &PythonModuleRequirement::validate);
        
        // call once to update py_pyprophet canvas 
        // alternative: load latest data from .ini and set py_selector (will update py_pyprophet via above signal/slot)
        py_pyprophet->validate(py_selector->getLastPython().toQString());
    }

    SwathTabWidget::~SwathTabWidget()
    {
        delete ui;
    }

  }   //namespace Internal
} //namspace OpenMS
