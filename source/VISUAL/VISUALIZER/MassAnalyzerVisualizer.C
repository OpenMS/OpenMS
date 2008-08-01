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

//Constructor
MassAnalyzerVisualizer::MassAnalyzerVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
	type_="MassAnalyzer";
  
	addLabel("Modify massanalyzer information.");	
	addSeperator();  
	
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



void MassAnalyzerVisualizer::load(MassAnalyzer &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempmassanalyzer_=s;
			
  
	
	update_();
}

void MassAnalyzerVisualizer::update_()
{
		if(! isEditable())
		{
			
			fillComboBox(massanalyzer_type_, &tempmassanalyzer_.NamesOfAnalyzerType[tempmassanalyzer_.getType()] , 1);
			fillComboBox(massanalyzer_res_method_, &tempmassanalyzer_.NamesOfResolutionMethod[tempmassanalyzer_.getResolutionMethod()] , 1);
			fillComboBox(massanalyzer_res_type_, &tempmassanalyzer_.NamesOfResolutionType[tempmassanalyzer_.getResolutionType()] , 1);
			fillComboBox(massanalyzer_scan_func_, &tempmassanalyzer_.NamesOfScanFunction[tempmassanalyzer_.getScanFunction()] , 1);
			fillComboBox(massanalyzer_scan_dir_, &tempmassanalyzer_.NamesOfScanDirection[tempmassanalyzer_.getScanDirection()] , 1);
			fillComboBox(massanalyzer_scan_law_, &tempmassanalyzer_.NamesOfScanLaw[tempmassanalyzer_.getScanLaw()] , 1);
			fillComboBox(massanalyzer_tandem_scan_method_, &tempmassanalyzer_.NamesOfTandemScanningMethod[tempmassanalyzer_.getTandemScanMethod()] , 1);
			fillComboBox(massanalyzer_reflectron_state_, &tempmassanalyzer_.NamesOfReflectronState[tempmassanalyzer_.getReflectronState()] , 1);
			
		}
		else
		{
			fillComboBox(massanalyzer_type_, tempmassanalyzer_.NamesOfAnalyzerType , MassAnalyzer::SIZE_OF_ANALYZERTYPE);
			fillComboBox(massanalyzer_res_method_, tempmassanalyzer_.NamesOfResolutionMethod , MassAnalyzer::SIZE_OF_RESOLUTIONMETHOD);
			fillComboBox(massanalyzer_res_type_, tempmassanalyzer_.NamesOfResolutionType , MassAnalyzer::SIZE_OF_RESOLUTIONTYPE);
			fillComboBox(massanalyzer_scan_func_, tempmassanalyzer_.NamesOfScanFunction , MassAnalyzer::SIZE_OF_SCANFUNCTION);
			fillComboBox(massanalyzer_scan_dir_, tempmassanalyzer_.NamesOfScanDirection , MassAnalyzer::SIZE_OF_SCANDIRECTION);
			fillComboBox(massanalyzer_scan_law_, tempmassanalyzer_.NamesOfScanLaw , MassAnalyzer::SIZE_OF_SCANLAW);
			fillComboBox(massanalyzer_tandem_scan_method_, tempmassanalyzer_.NamesOfTandemScanningMethod , MassAnalyzer::SIZE_OF_TANDEMSCANNINGMETHOD);
			fillComboBox(massanalyzer_reflectron_state_, tempmassanalyzer_.NamesOfReflectronState , MassAnalyzer::SIZE_OF_REFLECTRONSTATE);
			
			massanalyzer_type_->setCurrentIndex(tempmassanalyzer_.getType()); 
			massanalyzer_res_method_->setCurrentIndex(tempmassanalyzer_.getResolutionMethod()); 
			massanalyzer_res_type_->setCurrentIndex(tempmassanalyzer_.getResolutionType()); 
			massanalyzer_scan_func_->setCurrentIndex(tempmassanalyzer_.getScanFunction()); 
			massanalyzer_scan_dir_->setCurrentIndex(tempmassanalyzer_.getScanDirection()); 
			massanalyzer_scan_law_->setCurrentIndex(tempmassanalyzer_.getScanLaw()); 
			massanalyzer_tandem_scan_method_->setCurrentIndex(tempmassanalyzer_.getTandemScanMethod()); 
			massanalyzer_reflectron_state_->setCurrentIndex(tempmassanalyzer_.getReflectronState()); 
		
		}
		
		massanalyzer_res_->setText(String( tempmassanalyzer_.getResolution() ).c_str());
		massanalyzer_acc_->setText(String( tempmassanalyzer_.getAccuracy() ).c_str());
		massanalyzer_scan_rate_->setText(String( tempmassanalyzer_.getScanRate() ).c_str());
		massanalyzer_scan_time_->setText(String( tempmassanalyzer_.getScanTime() ).c_str());
		massanalyzer_TOF_->setText(String( tempmassanalyzer_.getTOFTotalPathLength() ).c_str());
		massanalyzer_iso_->setText(String( tempmassanalyzer_.getIsolationWidth() ).c_str());
		massanalyzer_final_MS_->setText(String( tempmassanalyzer_.getFinalMSExponent() ).c_str());
		massanalyzer_magnetic_fs_->setText(String( tempmassanalyzer_.getMagneticFieldStrength() ).c_str());
}

void MassAnalyzerVisualizer::store_()
{
	try
	{
		
		(*ptr_).setType((MassAnalyzer::AnalyzerType)massanalyzer_type_->currentIndex());		
		(*ptr_).setResolutionMethod((MassAnalyzer::ResolutionMethod)massanalyzer_res_method_->currentIndex());		
		(*ptr_).setResolutionType((MassAnalyzer::ResolutionType)massanalyzer_res_type_->currentIndex());		
		(*ptr_).setScanFunction((MassAnalyzer::ScanFunction)massanalyzer_scan_func_->currentIndex());		
		(*ptr_).setScanDirection((MassAnalyzer::ScanDirection)massanalyzer_scan_dir_->currentIndex());		
		(*ptr_).setScanLaw((MassAnalyzer::ScanLaw)massanalyzer_scan_law_->currentIndex());		
		(*ptr_).setTandemScanMethod((MassAnalyzer::TandemScanningMethod)	massanalyzer_tandem_scan_method_->currentIndex());		
		(*ptr_).setReflectronState((MassAnalyzer::ReflectronState)massanalyzer_reflectron_state_->currentIndex());		
		
		String m(massanalyzer_res_->text().toStdString());
		(*ptr_).setResolution(m.toFloat());
		String n(massanalyzer_acc_->text().toStdString());
		(*ptr_).setAccuracy(n.toFloat());
		String o(massanalyzer_scan_rate_->text().toStdString());
		(*ptr_).setScanRate(o.toFloat());
		String p(massanalyzer_scan_time_->text().toStdString());
		(*ptr_).setScanTime(p.toFloat());
		String q(massanalyzer_TOF_->text().toStdString());
		(*ptr_).setTOFTotalPathLength(q.toFloat());
		String r(massanalyzer_iso_->text().toStdString());
		(*ptr_).setIsolationWidth(r.toFloat());
		
		String s(massanalyzer_final_MS_->text().toStdString());
		(*ptr_).setFinalMSExponent(s.toInt());
		
		String t(massanalyzer_magnetic_fs_->text().toStdString());
		(*ptr_).setMagneticFieldStrength(t.toFloat() );
		
		
		
		tempmassanalyzer_=(*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new mass analyzer data. "<<e.what()<<endl;
	}
	
}

void MassAnalyzerVisualizer::reject_()
{
	
	try
	{

		update_();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original mass analyzer data. "<<e.what()<<endl;
	}
	
}

}
