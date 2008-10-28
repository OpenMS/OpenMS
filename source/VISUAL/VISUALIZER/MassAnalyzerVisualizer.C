// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free MassAnalyzer Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  this library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer:  Marc Sturm $
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
	
	MassAnalyzerVisualizer::MassAnalyzerVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<MassAnalyzer>()
	{
		addLabel("Modify massanalyzer information.");	
		addSeparator();  
		
		addComboBox(massanalyzer_type_, "Type");
		addComboBox(massanalyzer_res_method_, "Resolution method");
		addComboBox(massanalyzer_res_type_, "Resolution type");
		addComboBox(massanalyzer_scan_func_, "Scan function");
		addComboBox(massanalyzer_scan_dir_, "Scan direction");
		addComboBox(massanalyzer_scan_law_, "Scan law");
		addComboBox(massanalyzer_tandem_scan_method_, "Tandem scan maethod");
		addComboBox(massanalyzer_reflectron_state_, "Reflectron state");
			
		addDoubleLineEdit(massanalyzer_res_, "Resolution" );
		addDoubleLineEdit(massanalyzer_acc_, "Accuracy" );
		addDoubleLineEdit(massanalyzer_scan_rate_, "Scan rate (in s)" );
		addDoubleLineEdit(massanalyzer_scan_time_, "Scan time (in s)" );
		addDoubleLineEdit(massanalyzer_TOF_, "TOF Total path length (in mm)" );
		addDoubleLineEdit(massanalyzer_iso_, "Isolation width (in m/z)" );
		addDoubleLineEdit(massanalyzer_final_MS_, "Final MS exponent" );
		addDoubleLineEdit(massanalyzer_magnetic_fs_, "Magnetic field strength (in T)" );
		
		finishAdding_();
	}
	
	void MassAnalyzerVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox(massanalyzer_type_,& temp_.NamesOfAnalyzerType[temp_.getType()] , 1);
			fillComboBox(massanalyzer_res_method_,& temp_.NamesOfResolutionMethod[temp_.getResolutionMethod()] , 1);
			fillComboBox(massanalyzer_res_type_,& temp_.NamesOfResolutionType[temp_.getResolutionType()] , 1);
			fillComboBox(massanalyzer_scan_func_,& temp_.NamesOfScanFunction[temp_.getScanFunction()] , 1);
			fillComboBox(massanalyzer_scan_dir_,& temp_.NamesOfScanDirection[temp_.getScanDirection()] , 1);
			fillComboBox(massanalyzer_scan_law_,& temp_.NamesOfScanLaw[temp_.getScanLaw()] , 1);
			fillComboBox(massanalyzer_tandem_scan_method_,& temp_.NamesOfTandemScanningMethod[temp_.getTandemScanMethod()] , 1);
			fillComboBox(massanalyzer_reflectron_state_,& temp_.NamesOfReflectronState[temp_.getReflectronState()] , 1);
		}
		else
		{
			fillComboBox(massanalyzer_type_, temp_.NamesOfAnalyzerType , MassAnalyzer::SIZE_OF_ANALYZERTYPE);
			fillComboBox(massanalyzer_res_method_, temp_.NamesOfResolutionMethod , MassAnalyzer::SIZE_OF_RESOLUTIONMETHOD);
			fillComboBox(massanalyzer_res_type_, temp_.NamesOfResolutionType , MassAnalyzer::SIZE_OF_RESOLUTIONTYPE);
			fillComboBox(massanalyzer_scan_func_, temp_.NamesOfScanFunction , MassAnalyzer::SIZE_OF_SCANFUNCTION);
			fillComboBox(massanalyzer_scan_dir_, temp_.NamesOfScanDirection , MassAnalyzer::SIZE_OF_SCANDIRECTION);
			fillComboBox(massanalyzer_scan_law_, temp_.NamesOfScanLaw , MassAnalyzer::SIZE_OF_SCANLAW);
			fillComboBox(massanalyzer_tandem_scan_method_, temp_.NamesOfTandemScanningMethod , MassAnalyzer::SIZE_OF_TANDEMSCANNINGMETHOD);
			fillComboBox(massanalyzer_reflectron_state_, temp_.NamesOfReflectronState , MassAnalyzer::SIZE_OF_REFLECTRONSTATE);
			
			massanalyzer_type_->setCurrentIndex(temp_.getType()); 
			massanalyzer_res_method_->setCurrentIndex(temp_.getResolutionMethod()); 
			massanalyzer_res_type_->setCurrentIndex(temp_.getResolutionType()); 
			massanalyzer_scan_func_->setCurrentIndex(temp_.getScanFunction()); 
			massanalyzer_scan_dir_->setCurrentIndex(temp_.getScanDirection()); 
			massanalyzer_scan_law_->setCurrentIndex(temp_.getScanLaw()); 
			massanalyzer_tandem_scan_method_->setCurrentIndex(temp_.getTandemScanMethod()); 
			massanalyzer_reflectron_state_->setCurrentIndex(temp_.getReflectronState()); 
		}
		
		massanalyzer_res_->setText(String( temp_.getResolution() ).c_str());
		massanalyzer_acc_->setText(String( temp_.getAccuracy() ).c_str());
		massanalyzer_scan_rate_->setText(String( temp_.getScanRate() ).c_str());
		massanalyzer_scan_time_->setText(String( temp_.getScanTime() ).c_str());
		massanalyzer_TOF_->setText(String( temp_.getTOFTotalPathLength() ).c_str());
		massanalyzer_iso_->setText(String( temp_.getIsolationWidth() ).c_str());
		massanalyzer_final_MS_->setText(String( temp_.getFinalMSExponent() ).c_str());
		massanalyzer_magnetic_fs_->setText(String( temp_.getMagneticFieldStrength() ).c_str());
	}
	
	void MassAnalyzerVisualizer::store()
	{
		ptr_->setType((MassAnalyzer::AnalyzerType)massanalyzer_type_->currentIndex());		
		ptr_->setResolutionMethod((MassAnalyzer::ResolutionMethod)massanalyzer_res_method_->currentIndex());		
		ptr_->setResolutionType((MassAnalyzer::ResolutionType)massanalyzer_res_type_->currentIndex());		
		ptr_->setScanFunction((MassAnalyzer::ScanFunction)massanalyzer_scan_func_->currentIndex());		
		ptr_->setScanDirection((MassAnalyzer::ScanDirection)massanalyzer_scan_dir_->currentIndex());		
		ptr_->setScanLaw((MassAnalyzer::ScanLaw)massanalyzer_scan_law_->currentIndex());		
		ptr_->setTandemScanMethod((MassAnalyzer::TandemScanningMethod)	massanalyzer_tandem_scan_method_->currentIndex());		
		ptr_->setReflectronState((MassAnalyzer::ReflectronState)massanalyzer_reflectron_state_->currentIndex());		
		
		String m(massanalyzer_res_->text().toStdString());
		ptr_->setResolution(m.toFloat());
		String n(massanalyzer_acc_->text().toStdString());
		ptr_->setAccuracy(n.toFloat());
		String o(massanalyzer_scan_rate_->text().toStdString());
		ptr_->setScanRate(o.toFloat());
		String p(massanalyzer_scan_time_->text().toStdString());
		ptr_->setScanTime(p.toFloat());
		String q(massanalyzer_TOF_->text().toStdString());
		ptr_->setTOFTotalPathLength(q.toFloat());
		String r(massanalyzer_iso_->text().toStdString());
		ptr_->setIsolationWidth(r.toFloat());
		
		String s(massanalyzer_final_MS_->text().toStdString());
		ptr_->setFinalMSExponent(s.toInt());
		
		String t(massanalyzer_magnetic_fs_->text().toStdString());
		ptr_->setMagneticFieldStrength(t.toFloat() );
		
		temp_=(*ptr_);
	}
	
	void MassAnalyzerVisualizer::undo_()
	{
		update_();
	}

}
