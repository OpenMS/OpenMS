// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_ANDIHANDLER_H
#define OPENMS_FORMAT_HANDLERS_ANDIHANDLER_H

#include <netcdf.h>
#include <ms10.h>
#include <fstream>
#include <sstream>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
	namespace Internal
	{
	  /**
	  	@brief Read-only File handler for ANDI/MS format (Version 1.0)
	
			MapType has to be a MSExperiment or have the same interface.
	  	
	  	@note Do not use this class directly. It is only needed for ANDIFile.
	  */
		template <typename MapType>
	  class ANDIHandler
	  {
	    public:
	      /**@name Constructors and destructor */
	      //@{
	      ///
	      ANDIHandler(MapType& exp, const ProgressLogger& logger)
				:	exp_(exp),
					pol_(),
					logger_(logger)
	  		{
					MetaInfoRegistry& registry =	MetaInfo().registry();
					for (int i=0; i<NUM_PARAM; i++) 
					{
						registry.registerName(user_params_[i], description_[i]);
					}
				}
		    /// Destructor
		    virtual ~ANDIHandler()
		    {
		    }
		    //@}

				/// read the ANDIFile using the ANDI/MS-NETCDF library
				void parse(std::string file_name);

    	protected:
		    /**@name Handlers for each NETCDF struct */
		    //@{
		    /// fill administration data with @p admin_data
				void getAdminData_(MS_Admin_Data* admin_data);
		
				/// fill sample data with @p sample_data
				void getSampleData_(MS_Sample_Data* sample_data);
		
				/// fill test data with @p test_data
				void getTestData_(MS_Test_Data* test_data);
		
				/// fill instrument data with first instrument of @p inst_data
				void getInstrumentData_(MS_Instrument_Data* inst_data);
		
				/// fill scan data with @p scan_data and @p global_data
				void getRawPerScan_(UInt scan_number, MS_Raw_Per_Scan* scan_data, MS_Raw_Data_Global* global_data);
		    //@}
		
				/// map pointer for reading
				MapType& exp_;		
				///Scan polarity
				IonSource::Polarity pol_;
		
				/**@name meta value handling */
				//@{
				/// index of meta value string
				enum userParamsID {CONTACT=0, PROC, ANDI_ERROR, CALHIST, CALTIMES,
													INSTSERIAL, INSTCOMMENTS, INSTSOFTWARE, INSTFIRMWARE,
													INSTOS, INSTID, INLETTEMP, IONMODEADD, SRCTEMP, ACCPOT,
													INSTPARAMS, DETPOT, DETENTRPOT, NUM_PARAM};
				/// strings used as meta values
				static const std::string user_params_[];
				/// description of the meta values
				static const std::string description_[];
				//@}
		
				/// convert all char* struct members to string in case member is NULL
				inline String string_(char* input)
				{
					return (input == NULL)? "" : String(input);
				}
		
				/** 
					@brief check all float struct members in case member is negative
		
					unset member usually indicated by value -9999
				*/
				inline float float_(float input, float def=0.0f)
				{
					return (input < -1000)? def : input;
				}
		
				/** 
					@brief check all int struct members in case member is negative
		
					unset member usually indicated by value -9999
				*/
				inline int int_(int input, int def=0)
				{
					return (input < -1000)? def : input;
				}
				
				///Progress logging helper class
				const ProgressLogger& logger_;
				
				///DataProcessing auxilary variable
				DataProcessing data_processing_;
  };

	//---------------------------------------------------------------------------

	template <typename MapType>
	const std::string ANDIHandler<MapType>::user_params_[] = {"ContactPosition",
			"ProcessingNumer", "ErrorLog", "CalibrationHistory", "NumOfCalibrations",
			"InstSerial", "InstComments", "InstSoftware", "InstFirmware", "InstOS", "InstID",
			"InletTemp", "IonModeAdd", "SrcTemp", "AccPot", "InstParams", "DetPot", "DetEntrPot"};

	template <typename MapType>
	const std::string ANDIHandler<MapType>::description_[] = {"Position of the contact person",
			"number of times processed", "Processing Method error log", "history of calibration",
			"number of times calibrated", "Instrument serial number", "Instrument id comments",
			"Instrument software revision", "Instrument firmware revision",
			"Operating system revision", "Instrument ProteinIdentification code",
			"Spectrometer inlet temperature", "Additional ionization mode information",
			"Ionization source temperature", "Accelerating Potential",
			"Instrument parameter comments", "Detector potential", "Detector entrance potential"};

	template <typename MapType>
	void ANDIHandler<MapType>::parse(std::string file)
	{
		if (file=="") return;
	 
		char* file_name = (char*) file.c_str();
		
		std::ifstream infile(file_name);
		if (!infile)
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_name);
		}
		infile.close();

		ncopts = 0;  // Disable automatic error-abort of netcdf
		int file_id = ms_open_read(file_name);
		
		MS_Admin_Data      ms_admin;
		MS_Sample_Data     ms_sample;
		MS_Test_Data       ms_test;
		MS_Raw_Data_Global ms_raw_global;
		MS_Instrument_Data ms_inst;
		
		if (file_id == MS_ERROR)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_name,"Error: could not open ANDI/MS file. Invalid file!");
		}
		
		ms_init_global( 0, &ms_admin, &ms_sample, &ms_test, &ms_raw_global);
		
		if (ms_read_global( file_id, &ms_admin, &ms_sample, &ms_test, &ms_raw_global) == MS_ERROR) 
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_name,"Error: could not read global data of ANDI/MS file. Invalid file!");
		}
		
		// Global Data
		getAdminData_(&ms_admin);
		getSampleData_(&ms_sample);
		getTestData_(&ms_test);

		long index;
		long num_scans = ms_raw_global.nscans;
		logger_.startProgress(0,num_scans,"reading ANDI/MS file");
		exp_.reserve(num_scans);
		long num_inst = ms_admin.number_instrument_components;
		ms_admin_expt_t expt_type = ms_admin.experiment_type;
		bool is_library = (expt_library == expt_type );
		
		ms_init_instrument(0,&ms_inst);
		if (num_inst > 1) num_inst = 1; // Read only the first instrument
		for (index = 0; index < num_inst; index++)
		{
			logger_.setProgress(index);
			ms_inst.inst_no = index;

			if (ms_read_instrument( file_id, &ms_inst) == MS_ERROR)
			{
				std::cerr << "Warning: could not read instrument data of ANDI/MS file '" << file <<"'." << "\n";
			}
			else
			{
				//read data
				getInstrumentData_(&ms_inst);
			}
			//clean up
			ms_init_instrument(1,&ms_inst);
		}
		MS_Raw_Per_Scan		 ms_raw;
		MS_Raw_Library		 ms_lib;

		int err_code;
		
		ms_init_per_scan(0, &ms_raw, &ms_lib);
		
		for (index=0; index<num_scans; index++) 
		{
			ms_raw.scan_no = index;

			if (is_library)
			{
				ms_lib.scan_no = index;
				err_code = ms_read_per_scan(file_id, &ms_raw, &ms_lib);
			}
			else	
			{
				err_code = ms_read_per_scan(file_id, &ms_raw, NULL);
			}
			
			if (err_code == MS_ERROR)
			{
				std::cerr << "Warning: could not read scan " << index << " of ANDI/MS file '" << file << "'." << "\n";
			}
			else
			{
				//parse data
				getRawPerScan_(index, &ms_raw,&ms_raw_global);
			}	
			//clean up
			ms_init_per_scan( 1, &ms_raw, &ms_lib);					
		}

		ms_init_global( 1, &ms_admin, &ms_sample, &ms_test, &ms_raw_global);
		ms_close( file_id);
		logger_.endProgress();
	}

	template <typename MapType>
	void ANDIHandler<MapType>::getAdminData_(MS_Admin_Data* admin_data)
	{
		// Partition file reference into name and path
 		std::string file = string_(admin_data->source_file_reference);
		int last_slash = file.rfind("/",file.size());
		int last_backslash = file.rfind("\\",file.size());
		int cut = (last_slash > last_backslash)? last_slash : last_backslash;
		SourceFile sf;
		sf.setNameOfFile( file.substr(cut+1) ); 
		sf.setPathToFile( file.substr(0,cut+1) ); 
		sf.setFileType( string_(admin_data->source_file_format) );
		exp_.getSourceFiles().push_back(sf);
		
		ContactPerson contact;
		contact.setLastName( string_(admin_data->operator_name) );
		contact.setMetaValue(user_params_[CONTACT], String("Operator"));
		exp_.getContacts().push_back(contact);
		
		contact = ContactPerson();
		contact.setLastName( string_(admin_data->dataset_owner) );
		contact.setContactInfo( string_(admin_data->dataset_origin) );
		contact.setMetaValue(user_params_[CONTACT], String("Dataset owner"));
		exp_.getContacts().push_back(contact);
 		
 		data_processing_ = DataProcessing();
		data_processing_.getSoftware().setName( string_(admin_data->post_expt_program_name) );
		data_processing_.setMetaValue(user_params_[ANDI_ERROR], string_(admin_data->error_log));
		data_processing_.setMetaValue(user_params_[PROC], int_(admin_data->number_times_processed));
		std::stringstream buffer;
		buffer << string_(admin_data->calibration_history_0) << string_(admin_data->calibration_history_1)
					 << string_(admin_data->calibration_history_2) << string_(admin_data->calibration_history_3);
		data_processing_.setMetaValue(user_params_[CALHIST], String(buffer.str()));	
		data_processing_.setMetaValue(user_params_[CALTIMES], int_(admin_data->number_times_calibrated));	

		// unused MS_Admin_Data fields: comments, experiment_title, ms_template_revision, netcdf_revision, languages
		// experiment_date_time, netcdf_date_time, source_file_date_time, pre_expt_program_name
		// experiment_x_ref_0,	experiment_x_ref_1, experiment_x_ref_2, experiment_x_ref_3, external_file_ref_0
		// external_file_ref_1, external_file_ref_2, external_file_ref_3
	}


	template <typename MapType>
	void ANDIHandler<MapType>::getSampleData_(MS_Sample_Data* sample_data)
	{
		std::stringstream buffer;
		if (sample_data->internal_id != NULL) buffer << sample_data->internal_id;
		if (sample_data->external_id != NULL) buffer << "(" << sample_data->external_id << ")";
		exp_.getSample().setNumber(buffer.str());
		exp_.getSample().setComment( string_(sample_data->comments) );

		// Solid, Liquid, Gas, Supercritical Fluid, Plasma, Other
		int sample_map[] = {Sample::SOLID, Sample::LIQUID, Sample::GAS, 0, 0, 0, 0};
		exp_.getSample().setState( (Sample::SampleState) sample_map[sample_data->state - state_solid]);

		// unused sample_data fields: owner, receipt_date_time, procedure_name, prep_procedure, matrix, storage
		// disposal, history, manual_handling, prep_comments
	}

	template <typename MapType>
	void ANDIHandler<MapType>::getTestData_(MS_Test_Data* test_data)
	{
		exp_.getInstrument().setMetaValue(user_params_[INSTPARAMS], string_(test_data->comments));


		// Membrane Separator, Capillary Direct, Open Split, Jet Separator, Direct Inlet Probe, Septum, Particle Beam,
		// Reservoir, Moving Belt, Atmospheric Pressure Chemical Ionization, Flow Injection Analysis, Electrospray,
		// Infusion, Thermospray, Other Probe, Other
		exp_.getInstrument().getIonSources().resize(1);
		IonSource& src = exp_.getInstrument().getIonSources()[0];
		int inlet_map[] = {IonSource::MEMBRANESEPARATOR, 0, IonSource::OPENSPLIT, IonSource::JETSEPARATOR, IonSource::DIRECT, IonSource::SEPTUM,
		IonSource::PARTICLEBEAM, IonSource::RESERVOIR, IonSource::MOVINGBELT, 0, IonSource::FLOWINJECTIONANALYSIS, IonSource::ELECTROSPRAYINLET,
		IonSource::INFUSION, IonSource::THERMOSPRAYINLET, 0, 0};
		src.setInletType( (IonSource::InletType) inlet_map[test_data->ms_inlet - inlet_membrane]);
		src.setMetaValue(user_params_[INLETTEMP], float_(test_data->ms_inlet_temperature));

		// Electron Impact, Chemical Ionization, Fast Atom Bombardment, Field Desorption, Field Ionization,
		// Electrospray, Thermospray, Atmospheric Pressure Chemical Ionization, Plasma Desorption,
		// Laser Desorption, Spark Ionization, Thermal Ionization, Other
		int ion_map[] = {IonSource::EI, IonSource::CI, IonSource::FAB, IonSource::FD, IonSource::FI, IonSource::ESI, 
										 IonSource::TSP, IonSource::APCI, IonSource::PD, IonSource::LD, IonSource::SI, IonSource::TI, 0};
		src.setIonizationMethod( (IonSource::IonizationMethod) ion_map[test_data->ionization_mode - ionization_ei]);

		std::stringstream buffer;
		if (test_data->fab_type != NULL) buffer << "FABType=" << test_data->fab_type << " ";
		if (test_data->fab_matrix != NULL) buffer << "FABMatrix=" << test_data->fab_matrix << " ";
		if (test_data->reagent_gas != NULL) buffer << "ReagentGas=" << test_data->reagent_gas << " ";
		buffer << "ReagentGasPressure=" << test_data->reagent_gas_pressure << " ";
		buffer << "ElectronEnergy=" << test_data->electron_energy << " ";
		buffer << "LaserWaveLength=" << test_data->laser_wavelength << " ";
		buffer << "FilamentCurrent=" << test_data->filament_current << " ";
		buffer << "EmissionCurrent=" << test_data->emission_current << " ";
		src.setMetaValue(user_params_[IONMODEADD], String(buffer.str()));
		src.setMetaValue(user_params_[SRCTEMP], test_data->source_temperature);
		src.setMetaValue(user_params_[ACCPOT], test_data->accelerating_potential);

		pol_ = (IonSource::Polarity) (test_data->ionization_polarity - polarity_plus + 1);

		// Electron Multiplier, Photomultplier, Focal Plane Array, Faraday Cup, Conversion Dynode Electron Multiplier,
		// Conversion dynode Photomultiplier, Multicollector, Other
		exp_.getInstrument().getIonDetectors().resize(1);
		IonDetector& det = exp_.getInstrument().getIonDetectors()[0];
		int detector_map[] = {IonDetector::ELECTRONMULTIPLIER, IonDetector::PHOTOMULTIPLIER, IonDetector::FOCALPLANEARRAY, IonDetector::FARADAYCUP,
													IonDetector::CONVERSIONDYNODEELECTRONMULTIPLIER, IonDetector::CONVERSIONDYNODEPHOTOMULTIPLIER,
													IonDetector::MULTICOLLECTOR, 0};
		det.setType( (IonDetector::Type) detector_map[test_data->detector_type - detector_em]);
		det.setMetaValue(user_params_[DETPOT], float_(test_data->detector_potential));
		det.setMetaValue(user_params_[DETENTRPOT], float_(test_data->detector_entrance_potential));



		MassAnalyzer analyzer;
		int dir_map[] = {MassAnalyzer::UP, MassAnalyzer::DOWN, 0};
		analyzer.setScanDirection( (MassAnalyzer::ScanDirection) dir_map[test_data->scan_direction - direction_up]);
		
		// Linear, Exponential, Quadratic,  Other
		int law_map[] = {MassAnalyzer::LINEAR, MassAnalyzer::EXPONENTIAL, MassAnalyzer::QUADRATIC, 0};
		analyzer.setScanLaw( (MassAnalyzer::ScanLaw) law_map[test_data->scan_law - law_linear]);

		analyzer.setResolutionType( (MassAnalyzer::ResolutionType) (test_data->resolution_type - resolution_constant));
		analyzer.setScanTime(test_data->scan_time);

		if (test_data->resolution_method != NULL)
		{
			if (String(test_data->resolution_method) == "50% peak height")
			{
				analyzer.setResolutionMethod(MassAnalyzer::FWHM);
			}
			else
			{
				if (String(test_data->resolution_method) == "10% peak valley") 
				{
					analyzer.setResolutionMethod(MassAnalyzer::TENPERCENTVALLEY);
				}
			}
		}
		exp_.getInstrument().getMassAnalyzers().push_back(analyzer);

		// unused MS_Test_Data fields: mass_calibration_file, external_reference_file, internal_reference_file
		// test_data->separation_type
	}

	template <typename MapType>
	void ANDIHandler<MapType>::getInstrumentData_(MS_Instrument_Data* inst_data)
	{
 		exp_.getInstrument().setName( string_(inst_data->name) );
		exp_.getInstrument().setVendor( string_(inst_data->manufacturer) );
		exp_.getInstrument().setModel( string_(inst_data->model_number) );

		exp_.getInstrument().setMetaValue(user_params_[INSTSERIAL], string_(inst_data->serial_number));
		exp_.getInstrument().setMetaValue(user_params_[INSTCOMMENTS], string_(inst_data->comments));
		exp_.getInstrument().setMetaValue(user_params_[INSTSOFTWARE], string_(inst_data->software_version));
		exp_.getInstrument().setMetaValue(user_params_[INSTFIRMWARE], string_(inst_data->firmware_version));
		exp_.getInstrument().setMetaValue(user_params_[INSTOS], string_(inst_data->operating_system));
		exp_.getInstrument().setMetaValue(user_params_[INSTID], string_(inst_data->id));
	}

	template <typename MapType>
	void ANDIHandler<MapType>::getRawPerScan_(UInt scan_number, MS_Raw_Per_Scan* scan_data, MS_Raw_Data_Global* global_data)
	{
		float mass_factor = float_(global_data->mass_factor, 1.0f);
		float intens_factor = float_(global_data->intensity_factor, 1.0f);
		float intens_offset = float_(global_data->intensity_offset);

		// in case anyone set the factor accidentally to zero -> avoid all zero values
		if (mass_factor==0.0f) mass_factor = 1.0f;
		if (intens_factor==0.0f) intens_factor = 1.0f;
		
		//Abort if no masses are present or if times are stored in this scan
		if (global_data->has_masses!=1 || global_data->has_times==1)
		{
			return;
		}
		
		//resize experiment
		exp_.resize(exp_.size()+1);
		
		//set scan data (besides peaks)
		typename MapType::SpectrumType& spectrum = exp_.back();
		spectrum.resize(scan_data->points);
		spectrum.setRT( float_(scan_data->scan_acq_time));
		spectrum.setNativeID(String("index=")+ scan_number);
		spectrum.setMSLevel(1);
		spectrum.getDataProcessing().push_back(data_processing_); //assign general data processing info to all spectra
		ScanWindow window;
		window.begin = float_(scan_data->mass_range_min);
		window.end = float_(scan_data->mass_range_max);
		spectrum.getInstrumentSettings().getScanWindows().push_back(window);
		spectrum.getInstrumentSettings().setPolarity(pol_);

		//check intensity/mass format
		std::string intensity_format = ms_enum_to_string(global_data->intensity_format);
	  if (intensity_format!= "Short" && intensity_format!= "Long" && intensity_format!= "Float" && intensity_format!= "Double")
	  {
	  	throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"","ANDI/MS file parse error. Unknown intensity format '"+intensity_format+"'.");
		}
		std::string mass_format = ms_enum_to_string(global_data->mass_format);
	  if (mass_format!= "Short" && mass_format!= "Long" && mass_format!= "Float" && mass_format!= "Double")
	  {
	  	throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"", "ANDI/MS file parse error. Unknown mass format '"+mass_format+"'.");
		}
		
		//load peak data
		for (int i=0; i < scan_data->points; i++)
		{
			//parse intensity
			double intensity(0);
			if (intensity_format == "Short") intensity =  (double) ((short*)  scan_data->intensities)[i];
			else if (intensity_format == "Long") intensity =  (double) ((long*)   scan_data->intensities)[i];
			else if (intensity_format == "Float") intensity =  (double) ((float*)  scan_data->intensities)[i];
			else if (intensity_format == "Double") intensity =  ((double*) scan_data->intensities)[i];
			spectrum[i].setIntensity(intensity * intens_factor + intens_offset);
			
			//parse mass
			double mass(0);
			if (mass_format == "Short") mass = (double) ((short*)  scan_data->masses)[i];
			else if (mass_format == "Long") mass = (double) ((long*)   scan_data->masses)[i];
			else if (mass_format == "Float") mass = (double) ((float*)  scan_data->masses)[i];
			else if (mass_format == "Double") mass = ((double*) scan_data->masses)[i];
			spectrum[i].setPosition(mass * mass_factor);
		}

		// unused MS_Raw_Data_Global fields:	mass_axis_global_min, mass_axis_global_max, time_axis_global_min,
		// time_axis_global_max, intensity_axis_global_min, intensity_axis_global_max, calibrated_mass_min);
		// calibrated_mass_max, uniform_flag, mass_label, time_label, intensity_label,
		// mass_units, time_units, intensity_units, total_intensity_units

		// unused MS_Raw_Per_Scan fields:	flags, a_d_rate, a_d_coadditions, scan_duration, time_range_min,
		// time_range_max, inter_scan_time, resolution, actual_scan_no
	}

	} // namespace Internal
} // namespace OpenMS

#endif
