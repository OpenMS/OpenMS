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
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/InstrumentSettingsVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  InstrumentSettingsVisualizer::InstrumentSettingsVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<InstrumentSettings>()
  {
    addLabel_("Modify the settings of the instrument.");
    addSeparator_();
    addComboBox_(instrumentsettings_scan_mode_, "Scan mode");
    addBooleanComboBox_(zoom_scan_, "Zoom scan");
    addComboBox_(instrumentsettings_polarity_, "Polarity");

    finishAdding_();
  }

  void InstrumentSettingsVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(instrumentsettings_scan_mode_, &temp_.NamesOfScanMode[temp_.getScanMode()], 1);
      fillComboBox_(instrumentsettings_polarity_, &IonSource::NamesOfPolarity[temp_.getPolarity()], 1);

    }
    else
    {
      fillComboBox_(instrumentsettings_scan_mode_, InstrumentSettings::NamesOfScanMode, InstrumentSettings::SIZE_OF_SCANMODE);
      fillComboBox_(instrumentsettings_polarity_, IonSource::NamesOfPolarity, IonSource::SIZE_OF_POLARITY);


      instrumentsettings_scan_mode_->setCurrentIndex(temp_.getScanMode());
      zoom_scan_->setCurrentIndex(temp_.getZoomScan());
      instrumentsettings_polarity_->setCurrentIndex(temp_.getPolarity());
    }
  }

  void InstrumentSettingsVisualizer::store()
  {
    ptr_->setScanMode((InstrumentSettings::ScanMode)instrumentsettings_scan_mode_->currentIndex());
    ptr_->setZoomScan(zoom_scan_->currentIndex());
    ptr_->setPolarity((IonSource::Polarity)instrumentsettings_polarity_->currentIndex());

    temp_ = (*ptr_);
  }

  void InstrumentSettingsVisualizer::undo_()
  {
    update_();
  }

}
