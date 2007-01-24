// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_ANDIHANDLER_H
#define OPENMS_FORMAT_HANDLERS_ANDIHANDLER_H

#include <netcdf.h>
#include <ms10.h>
#include <fstream>
#include <sstream>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
	namespace Internal
	{

  /**
  	@brief Write-only File handler for ANDIFile (Version 1.0)

		MapType has to be a MSExperiment or have the same interface.
  	Do not use this class. It is only needed in ANDIFile.
  	  	
  */
	template <typename MapType>
  class ANDIHandler
  {
    public:
      /**@name Constructors and destructor */
      //@{
      ///
      ANDIHandler(MapType& exp)
			: exp_(exp), peak_count_(0), peak_(), pol_()
  		{
				MetaInfoRegistry& registry =	MetaInfo().registry();
				for (int i=0; i<NUM_PARAM; i++) {
					registry.registerName(user_params_[i], description_[i]);
				}
			}
      ///
      virtual ~ANDIHandler(){}
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
		void getRawPerScan_(long index, MS_Raw_Per_Scan* scan_data, MS_Raw_Data_Global* global_data);
    //@}

		/// map pointer for reading
		MapType& exp_;

		/**@name temporary datastructures to hold parsed data */
		//@{
		Size peak_count_;
		typename MapType::SpectrumType::PeakType peak_;
		typename MapType::SpectrumType* spec_;
		IonSource::Polarity pol_;
		//@}

		/**@name meta value handling */
		//@{
		/// index of meta value string
		enum userParamsID {CONTACT=0, PROC, ERROR, CALHIST, CALTIMES,
											INSTSERIAL, INSTCOMMENTS, INSTSOFTWARE, INSTFIRMWARE,
											INSTOS, INSTID, INLETTEMP, IONMODEADD, SRCTEMP, ACCPOT,
											INSTPARAMS, DETPOT, DETENTRPOT, NUM_PARAM};
		/// strings used as meta values
		static const std::string user_params_[];
		/// description of the meta values
		static const std::string description_[];
		//@}

		/// convert all char* struct members to string in case member is NULL
		inline std::string string_(char* input)
		{
			return (input == NULL)? "" : std::string(input);
		}

		/** @brief check all float struct members in case member is negative

				unset member usually indicated by value -9999
		*/
		inline float float_(float input, float def=0.0f)
		{
			return (input < -1000)? def : input;
		}

		/** @brief check all int struct members in case member is negative

				unset member usually indicated by value -9999
		*/
		inline int int_(int input, int def=0)
		{
			return (input < -1000)? def : input;
		}
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
			"Operating system revision", "Instrument identification code",
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
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_name);
		infile.close();

		ncopts = 0;  // Disable automatic Errorabort of netcdf
		int file_id = ms_open_read(file_name);
		
		MS_Admin_Data      ms_admin;
		MS_Sample_Data     ms_sample;
		MS_Test_Data       ms_test;
		MS_Raw_Data_Global ms_raw_global;
		MS_Instrument_Data ms_inst;
		
		if ( MS_ERROR == file_id ) 
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_name,"invalid ANDIFile");

		ms_init_global( 0, &ms_admin, &ms_sample, &ms_test, &ms_raw_global);
		if ( MS_ERROR == ms_read_global( file_id, &ms_admin, &ms_sample, &ms_test, &ms_raw_global) ) 
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_name,"ms_read_global() failed");

		// Global Data
		getAdminData_(&ms_admin);
		getSampleData_(&ms_sample);
		getTestData_(&ms_test);

		long index;
		long num_scans = ms_raw_global.nscans;
		exp_.resize(num_scans);
		long num_inst = ms_admin.number_instrument_components;
		ms_admin_expt_t expt_type = ms_admin.experiment_type;
		bool is_library = (expt_library == expt_type );
		
		ms_init_instrument(0,&ms_inst);
		if (num_inst > 1) num_inst = 1; // Read only the first instrument
		for (index = 0; index < num_inst; index++)
		{
			ms_inst.inst_no = index;

			if ( MS_ERROR == ms_read_instrument( file_id, &ms_inst) ) 
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_name,"ms_read_instrument() failed");

			// Instrument Data
			getInstrumentData_(&ms_inst);
			ms_init_instrument(1,&ms_inst);
		}
		MS_Raw_Per_Scan		 ms_raw;
		MS_Raw_Library		 ms_lib;

		int err_code;
		
		ms_init_per_scan(0, &ms_raw, &ms_lib);
		
		for (index=0; index<num_scans; index++) {
			ms_raw.scan_no = index;

			if (is_library) {
				ms_lib.scan_no = index;
				err_code = ms_read_per_scan(file_id, &ms_raw, &ms_lib);
			}
			else	err_code = ms_read_per_scan(file_id, &ms_raw, NULL);
			
			if (MS_ERROR == err_code)
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_name,"ms_read_per_scan() failed");

			// Raw Data
			getRawPerScan_(index,&ms_raw,&ms_raw_global);
			ms_init_per_scan( 1, &ms_raw, &ms_lib);
		}

		ms_init_global( 1, &ms_admin, &ms_sample, &ms_test, &ms_raw_global);
		ms_close( file_id);
	}

	template <typename MapType>
	void ANDIHandler<MapType>::getAdminData_(MS_Admin_Data* admin_data)
	{
		// Partition file reference into name and path
 		std::string file = string_(admin_data->source_file_reference);
		int last_slash = file.rfind("/",file.size());
		int last_backslash = file.rfind("\\",file.size());
		int cut = (last_slash > last_backslash)? last_slash : last_backslash;
		exp_.getSourceFile().setNameOfFile( file.substr(cut+1) ); 
		exp_.getSourceFile().setPathToFile( file.substr(0,cut+1) ); 

		exp_.getSourceFile().setFileType( string_(admin_data->source_file_format) );

		ContactPerson contact;
		contact.setLastName( string_(admin_data->operator_name) );
		contact.setMetaValue(user_params_[CONTACT], std::string("Operator"));
		exp_.getContacts().push_back(contact);

		
		contact = ContactPerson();
		contact.setLastName( string_(admin_data->dataset_owner) );
		contact.setContactInfo( string_(admin_data->dataset_origin) );
		contact.setMetaValue(user_params_[CONTACT], std::string("Dataset owner"));
		exp_.getContacts().push_back(contact);
 		
		typedef ProcessingMethod pm;
		// Centroided Mass Spectrum, Continuum Mass Spectrum, Library Mass Spectrum
		int exp_map[] = {SpectrumSettings::PEAKS, SpectrumSettings::RAWDATA, 0};
		exp_.getProcessingMethod().setSpectrumType( (SpectrumSettings::SpectrumType) exp_map[admin_data->experiment_type - expt_centroid]);

		exp_.getSoftware().setName( string_(admin_data->post_expt_program_name) );
		exp_.getProcessingMethod().setMetaValue(user_params_[ERROR], string_(admin_data->error_log));
		exp_.getProcessingMethod().setMetaValue(user_params_[PROC], int_(admin_data->number_times_processed));
	 
		std::stringstream buffer;
		buffer << string_(admin_data->calibration_history_0) << string_(admin_data->calibration_history_1)
					 << string_(admin_data->calibration_history_2) << string_(admin_data->calibration_history_3);
		exp_.getProcessingMethod().setMetaValue(user_params_[CALHIST], buffer.str());	
		exp_.getProcessingMethod().setMetaValue(user_params_[CALTIMES], int_(admin_data->number_times_calibrated));	

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
		typedef IonSource is;
		IonSource& src = exp_.getInstrument().getIonSource();
		int inlet_map[] = {is::MEMBRANESEPARATOR, 0, is::OPENSPLIT, is::JETSEPARATOR, is::DIRECT, is::SEPTUM,
		is::PARTICLEBEAM, is::RESERVOIR, is::MOVINGBELT, 0, is::FLOWINJECTIONANALYSIS, is::ELECTROSPRAYINLET,
		is::INFUSION, is::THERMOSPRAYINLET, 0, 0};
		src.setInletType( (is::InletType) inlet_map[test_data->ms_inlet - inlet_membrane]);
		src.setMetaValue(user_params_[INLETTEMP], float_(test_data->ms_inlet_temperature));

		// Electron Impact, Chemical Ionization, Fast Atom Bombardment, Field Desorption, Field Ionization,
		// Electrospray, Thermospray, Atmospheric Pressure Chemical Ionization, Plasma Desorption,
		// Laser Desorption, Spark Ionization, Thermal Ionization, Other
		int ion_map[] = {is::EI, is::CI, is::FAB, is::FD, is::FI, is::ESI, 
										 is::TSP, is::APCI, is::PD, is::LD, is::SI, is::TI, 0};
		src.setIonizationMethod( (is::IonizationMethod) ion_map[test_data->ionization_mode - ionization_ei]);

		std::stringstream buffer;
		if (test_data->fab_type != NULL) buffer << "FABType=" << test_data->fab_type << " ";
		if (test_data->fab_matrix != NULL) buffer << "FABMatrix=" << test_data->fab_matrix << " ";
		if (test_data->reagent_gas != NULL) buffer << "ReagentGas=" << test_data->reagent_gas << " ";
		buffer << "ReagentGasPressure=" << test_data->reagent_gas_pressure << " ";
		buffer << "ElectronEnergy=" << test_data->electron_energy << " ";
		buffer << "LaserWaveLength=" << test_data->laser_wavelength << " ";
		buffer << "FilamentCurrent=" << test_data->filament_current << " ";
		buffer << "EmissionCurrent=" << test_data->emission_current << " ";
		src.setMetaValue(user_params_[IONMODEADD], buffer.str());
		src.setMetaValue(user_params_[SRCTEMP], test_data->source_temperature);
		src.setMetaValue(user_params_[ACCPOT], test_data->accelerating_potential);

		pol_ = (IonSource::Polarity) (test_data->ionization_polarity - polarity_plus + 1);

		// Electron Multiplier, Photomultplier, Focal Plane Array, Faraday Cup, Conversion Dynode Electron Multiplier,
		// Conversion dynode Photomultiplier, Multicollector, Other
		typedef IonDetector id;
		IonDetector& det = exp_.getInstrument().getIonDetector();
		int detector_map[] = {id::ELECTRONMULTIPLIER, id::PHOTOMULTIPLIER, id::FOCALPLANEARRAY, id::FARADAYCUP,
													id::CONVERSIONDYNODEELECTRONMULTIPLIER, id::CONVERSIONDYNODEPHOTOMULTIPLIER,
													id::MULTICOLLECTOR, 0};
		det.setType( (id::Type) detector_map[test_data->detector_type - detector_em]);
		det.setMetaValue(user_params_[DETPOT], float_(test_data->detector_potential));
		det.setMetaValue(user_params_[DETENTRPOT], float_(test_data->detector_entrance_potential));



		typedef MassAnalyzer ma;
		ma analyzer;

		int dir_map[] = {ma::UP, ma::DOWN, 0};
		analyzer.setScanDirection( (ma::ScanDirection) dir_map[test_data->scan_direction - direction_up]);
		
		// Linear, Exponential, Quadratic,  Other
		int law_map[] = {ma::LINEAR, ma::EXPONENTIAL, ma::QUADRATIC, 0};
		analyzer.setScanLaw( (ma::ScanLaw) law_map[test_data->scan_law - law_linear]);

		//Mass Scan, Selected Ion Detection, Other
		int function_map[] = {ma::MASSSCAN, ma::SELECTEDIONDETECTION, 0};
		analyzer.setScanFunction( (ma::ScanFunction) function_map[test_data->scan_function - function_scan]);

		analyzer.setResolutionType( (ma::ResolutionType) (test_data->resolution_type - resolution_constant));
		analyzer.setScanTime(test_data->scan_time);

		if (test_data->resolution_method != NULL)
		{
			if (test_data->resolution_method == "50% peak height") analyzer.setResolutionMethod(ma::FWHM);
			else if (test_data->resolution_method == "10% peak valley") analyzer.setResolutionMethod(ma::TENPERCENTVALLEY);
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
	void ANDIHandler<MapType>::getRawPerScan_(long index, MS_Raw_Per_Scan* scan_data, MS_Raw_Data_Global* global_data)
	{
		float mass_factor = float_(global_data->mass_factor, 1.0f);
		float intens_factor = float_(global_data->intensity_factor, 1.0f);
		float intens_offset = float_(global_data->intensity_offset);

		// in case anyone set the factor accidentally to zero -> avoid all zero values
		if (mass_factor==0.0f) mass_factor = 1.0f;
		if (intens_factor==0.0f) intens_factor = 1.0f;

		// Length of raw data array
		const long n = scan_data->points;

		double intensity;
		bool has_masses	=	(global_data->has_masses == 1);
		bool has_times = 	(global_data->has_times == 1);

		if (!has_masses || has_times) return;

		exp_[index].resize(n);
		spec_ = &exp_[index];

		spec_->setRetentionTime( float_(scan_data->scan_acq_time),
																			 float_(global_data->delay_time),
																			 float_(global_data->run_time));
		spec_->setMSLevel(1);
		spec_->getInstrumentSettings().setMzRangeStart(float_(scan_data->mass_range_min));
		spec_->getInstrumentSettings().setMzRangeStop(float_(scan_data->mass_range_max));
		spec_->getInstrumentSettings().setPolarity(pol_);


		for (int i=0; i < n; i++)
		{
			std::string format = ms_enum_to_string(global_data->intensity_format);

			if (format == "Short")       intensity =  (double) ((short*)  scan_data->intensities)[i];
			else if (format == "Long")   intensity =  (double) ((long*)   scan_data->intensities)[i];
			else if (format == "Float")  intensity =  (double) ((float*)  scan_data->intensities)[i];
			else if (format == "Double") intensity =  ((double*) scan_data->intensities)[i];
			else throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"",
											"ANDIFile parse error. Unknown intensity format.");

			intensity = intensity*intens_factor + intens_offset;

			double masses;
			format = ms_enum_to_string(global_data->mass_format);
			if (format == "Short")       masses = (double) ((short*)  scan_data->masses)[i];
			else if (format == "Long")   masses = (double) ((long*)   scan_data->masses)[i];
			else if (format == "Float")  masses = (double) ((float*)  scan_data->masses)[i];
			else if (format == "Double") masses = ((double*) scan_data->masses)[i];
			else throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"",
											"ANDIFile parse error. Unknown mass format.");

			masses *= mass_factor;

			// Build 1D peak
			spec_->getContainer()[i].getIntensity() = intensity;
			spec_->getContainer()[i].getPosition()[0] = masses;
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
