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
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/MassAnalyzerVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  MassAnalyzerVisualizer::MassAnalyzerVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<MassAnalyzer>()
  {
    addLabel_("Modify massanalyzer information.");
    addSeparator_();

    addIntLineEdit_(order_, "Order");
    addComboBox_(type_, "Type");
    addComboBox_(res_method_, "Resolution method");
    addComboBox_(res_type_, "Resolution type");
    addComboBox_(scan_dir_, "Scan direction");
    addComboBox_(scan_law_, "Scan law");
    addComboBox_(reflectron_state_, "Reflectron state");

    addDoubleLineEdit_(res_, "Resolution");
    addDoubleLineEdit_(acc_, "Accuracy");
    addDoubleLineEdit_(scan_rate_, "Scan rate (in s)");
    addDoubleLineEdit_(scan_time_, "Scan time (in s)");
    addDoubleLineEdit_(TOF_, "TOF Total path length (in meter)");
    addDoubleLineEdit_(iso_, "Isolation width (in m/z)");
    addDoubleLineEdit_(final_MS_, "Final MS exponent");
    addDoubleLineEdit_(magnetic_fs_, "Magnetic field strength (in T)");

    finishAdding_();
  }

  void MassAnalyzerVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(type_, &temp_.NamesOfAnalyzerType[temp_.getType()], 1);
      fillComboBox_(res_method_, &temp_.NamesOfResolutionMethod[temp_.getResolutionMethod()], 1);
      fillComboBox_(res_type_, &temp_.NamesOfResolutionType[temp_.getResolutionType()], 1);
      fillComboBox_(scan_dir_, &temp_.NamesOfScanDirection[temp_.getScanDirection()], 1);
      fillComboBox_(scan_law_, &temp_.NamesOfScanLaw[temp_.getScanLaw()], 1);
      fillComboBox_(reflectron_state_, &temp_.NamesOfReflectronState[temp_.getReflectronState()], 1);
    }
    else
    {
      fillComboBox_(type_, temp_.NamesOfAnalyzerType, MassAnalyzer::SIZE_OF_ANALYZERTYPE);
      fillComboBox_(res_method_, temp_.NamesOfResolutionMethod, MassAnalyzer::SIZE_OF_RESOLUTIONMETHOD);
      fillComboBox_(res_type_, temp_.NamesOfResolutionType, MassAnalyzer::SIZE_OF_RESOLUTIONTYPE);
      fillComboBox_(scan_dir_, temp_.NamesOfScanDirection, MassAnalyzer::SIZE_OF_SCANDIRECTION);
      fillComboBox_(scan_law_, temp_.NamesOfScanLaw, MassAnalyzer::SIZE_OF_SCANLAW);
      fillComboBox_(reflectron_state_, temp_.NamesOfReflectronState, MassAnalyzer::SIZE_OF_REFLECTRONSTATE);

      type_->setCurrentIndex(temp_.getType());
      res_method_->setCurrentIndex(temp_.getResolutionMethod());
      res_type_->setCurrentIndex(temp_.getResolutionType());
      scan_dir_->setCurrentIndex(temp_.getScanDirection());
      scan_law_->setCurrentIndex(temp_.getScanLaw());
      reflectron_state_->setCurrentIndex(temp_.getReflectronState());
    }

    order_->setText(String(temp_.getOrder()).c_str());
    res_->setText(String(temp_.getResolution()).c_str());
    acc_->setText(String(temp_.getAccuracy()).c_str());
    scan_rate_->setText(String(temp_.getScanRate()).c_str());
    scan_time_->setText(String(temp_.getScanTime()).c_str());
    TOF_->setText(String(temp_.getTOFTotalPathLength()).c_str());
    iso_->setText(String(temp_.getIsolationWidth()).c_str());
    final_MS_->setText(String(temp_.getFinalMSExponent()).c_str());
    magnetic_fs_->setText(String(temp_.getMagneticFieldStrength()).c_str());
  }

  void MassAnalyzerVisualizer::store()
  {
    ptr_->setOrder(order_->text().toInt());
    ptr_->setType((MassAnalyzer::AnalyzerType)type_->currentIndex());
    ptr_->setResolutionMethod((MassAnalyzer::ResolutionMethod)res_method_->currentIndex());
    ptr_->setResolutionType((MassAnalyzer::ResolutionType)res_type_->currentIndex());
    ptr_->setScanDirection((MassAnalyzer::ScanDirection)scan_dir_->currentIndex());
    ptr_->setScanLaw((MassAnalyzer::ScanLaw)scan_law_->currentIndex());
    ptr_->setReflectronState((MassAnalyzer::ReflectronState)reflectron_state_->currentIndex());

    ptr_->setResolution(res_->text().toDouble());
    ptr_->setAccuracy(acc_->text().toDouble());
    ptr_->setScanRate(scan_rate_->text().toDouble());
    ptr_->setScanTime(scan_time_->text().toDouble());
    ptr_->setTOFTotalPathLength(TOF_->text().toDouble());
    ptr_->setIsolationWidth(iso_->text().toDouble());
    ptr_->setFinalMSExponent(final_MS_->text().toInt());
    ptr_->setMagneticFieldStrength(magnetic_fs_->text().toDouble());

    temp_ = (*ptr_);
  }

  void MassAnalyzerVisualizer::undo_()
  {
    update_();
  }

}
