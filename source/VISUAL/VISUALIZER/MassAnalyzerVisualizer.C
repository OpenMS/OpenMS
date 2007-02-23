// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free massanalyzer; you can redistribute it and/or
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
// $Maintainer:  stefan_heess $
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
		
	addLineEdit(massanalyzer_res_, "Resolution" );
	addLineEdit(massanalyzer_acc_, "Accuracy" );
	addLineEdit(massanalyzer_scan_rate_, "Scan rate (in s)" );
	addLineEdit(massanalyzer_scan_time_, "Scan time (in s)" );
	addLineEdit(massanalyzer_TOF_, "TOF Total path length (in mm)" );
	addLineEdit(massanalyzer_iso_, "Isolation width (in m/z)" );
	addLineEdit(massanalyzer_final_MS_, "Final MS exponent" );
	addLineEdit(massanalyzer_magnetic_fs_, "Magnetic field strength (in T)" );
	
	finishAdding_();
	
	
	// A validator to check the input for the resolution
	QDoubleValidator *massanalyzer_res_vali_= new QDoubleValidator(massanalyzer_res_);
	massanalyzer_res_->setValidator(massanalyzer_res_vali_);
	// A validator to check the input for the accuracy
	QDoubleValidator *massanalyzer_acc_vali_ = new QDoubleValidator(massanalyzer_acc_);
	massanalyzer_acc_->setValidator(massanalyzer_acc_vali_);
	// A validator to check the input for the scan rate
	QDoubleValidator *massanalyzer_sr_vali_= new QDoubleValidator(massanalyzer_scan_rate_);
	massanalyzer_scan_rate_->setValidator(massanalyzer_sr_vali_);
	// A validator to check the input for the scan time
	QDoubleValidator *massanalyzer_st_vali_ = new QDoubleValidator(massanalyzer_scan_time_);
	massanalyzer_scan_time_->setValidator(massanalyzer_st_vali_);
	// A validator to check the input for the TOF total path length
	QDoubleValidator *massanalyzer_TOF_vali_ = new QDoubleValidator(massanalyzer_TOF_);
	massanalyzer_TOF_->setValidator(massanalyzer_TOF_vali_);
	// A validator to check the input for the isolation width
	QDoubleValidator *massanalyzer_iso_vali_ = new QDoubleValidator(massanalyzer_iso_);
	massanalyzer_iso_->setValidator(massanalyzer_iso_vali_);
	// A validator to check the input for the final MS exponent
	QIntValidator *massanalyzer_final_vali_ = new QIntValidator(massanalyzer_final_MS_);
	massanalyzer_final_MS_->setValidator(massanalyzer_final_vali_);
	// A validator to check the input for the magnetic field strngth
	QDoubleValidator *massanalyzer_fs_vali_ = new QDoubleValidator(massanalyzer_magnetic_fs_);
	massanalyzer_magnetic_fs_->setValidator(massanalyzer_fs_vali_);
			
}


void MassAnalyzerVisualizer::load(MassAnalyzer &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempmassanalyzer_=s;
			
  fillComboBox(massanalyzer_type_, s.NamesOfAnalyzerType , MassAnalyzer::SIZE_OF_ANALYZERTYPE);
	fillComboBox(massanalyzer_res_method_, s.NamesOfResolutionMethod , MassAnalyzer::SIZE_OF_RESOLUTIONMETHOD);
	fillComboBox(massanalyzer_res_type_, s.NamesOfResolutionType , MassAnalyzer::SIZE_OF_RESOLUTIONTYPE);
	fillComboBox(massanalyzer_scan_func_, s.NamesOfScanFunction , MassAnalyzer::SIZE_OF_SCANFUNCTION);
	fillComboBox(massanalyzer_scan_dir_, s.NamesOfScanDirection , MassAnalyzer::SIZE_OF_SCANDIRECTION);
	fillComboBox(massanalyzer_scan_law_, s.NamesOfScanLaw , MassAnalyzer::SIZE_OF_SCANLAW);
	fillComboBox(massanalyzer_tandem_scan_method_, s.NamesOfTandemScanningMethod , MassAnalyzer::SIZE_OF_TANDEMSCANNINGMETHOD);
	fillComboBox(massanalyzer_reflectron_state_, s.NamesOfReflectronState , MassAnalyzer::SIZE_OF_REFLECTRONSTATE);
	
	update_();
}

void MassAnalyzerVisualizer::update_()
{
		massanalyzer_type_->setCurrentIndex(tempmassanalyzer_.getType()); 
		massanalyzer_res_method_->setCurrentIndex(tempmassanalyzer_.getResolutionMethod()); 
		massanalyzer_res_type_->setCurrentIndex(tempmassanalyzer_.getResolutionType()); 
		massanalyzer_scan_func_->setCurrentIndex(tempmassanalyzer_.getScanFunction()); 
		massanalyzer_scan_dir_->setCurrentIndex(tempmassanalyzer_.getScanDirection()); 
		massanalyzer_scan_law_->setCurrentIndex(tempmassanalyzer_.getScanLaw()); 
		massanalyzer_tandem_scan_method_->setCurrentIndex(tempmassanalyzer_.getTandemScanMethod()); 
		massanalyzer_reflectron_state_->setCurrentIndex(tempmassanalyzer_.getReflectronState()); 
		
		massanalyzer_res_->setText(String( tempmassanalyzer_.getResolution() ).c_str());
		massanalyzer_acc_->setText(String( tempmassanalyzer_.getAccuracy() ).c_str());
		massanalyzer_scan_rate_->setText(String( tempmassanalyzer_.getScanRate() ).c_str());
		massanalyzer_scan_time_->setText(String( tempmassanalyzer_.getScanTime() ).c_str());
		massanalyzer_TOF_->setText(String( tempmassanalyzer_.getTOFTotalPathLength() ).c_str());
		massanalyzer_iso_->setText(String( tempmassanalyzer_.getIsolationWidth() ).c_str());
		massanalyzer_final_MS_->setText(String( tempmassanalyzer_.getFinalMSExponent() ).c_str());
		massanalyzer_magnetic_fs_->setText(String( tempmassanalyzer_.getMagneticFieldStrength() ).c_str());
}

void MassAnalyzerVisualizer::store()
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

void MassAnalyzerVisualizer::reject()
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
