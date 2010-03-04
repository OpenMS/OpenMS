// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Flex series file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz$
// $Authors: Guillaume Belz$
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XMASSFILE_H
#define OPENMS_FORMAT_XMASSFILE_H

#include <OpenMS/FORMAT/HANDLERS/AcqusHandler.h>
#include <OpenMS/FORMAT/HANDLERS/FidHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
 	/**
 		@brief File adapter for 'XMass Analysis (fid)' files.
 		
 		XMass Analysis files is native format for Bruker spectrometer Flex Series.<br />
 		Each spectrum are saved in one directory. Each directory contains several files. 
 		We use 2 files for import in OpenMS :<br />
 		<b>acqus</b> : contains meta data about calibration (conversion for time to mz ratio), 
 		instrument specification and acquisition method.<br />
 		<b>fid</b> : contains intensity array. Intensity for each point are coded in 4 bytes integer.
 		
 		@note MZ ratio are calculed with forumla based on article :<br />
    <i>A database application for pre-processing, storage and comparison of mass spectra 
    derived from patients and controls</i><br />
    <b>Mark K Titulaer, Ivar Siccama, Lennard J Dekker, Angelique LCT van Rijswijk, 
    Ron MA Heeren, Peter A Sillevis Smitt, and Theo M Luider</b><br />
    BMC Bioinformatics. 2006; 7: 403<br />
    http://www.pubmedcentral.nih.gov/picrender.fcgi?artid=1594579&blobtype=pdf<br />
  	
  	@ingroup FileIO
  */
  
  class OPENMS_DLLAPI XMassFile
 	  :	public ProgressLogger
  {
    public:
      /// Default constructor
      XMassFile();
      /// Destructor
      virtual ~XMassFile();
      
			/**
				@brief Loads a spectrum from a XMass file.

				@param spectrum has to be a MSSpectrum or have the same interface.

				@exception Exception::FileNotFound is thrown if the file could not be read
			*/      
      template <class PeakType>
      void load(const String& filename, MSSpectrum<PeakType>& spectrum)
      {
        Internal::AcqusHandler acqus(filename.prefix(filename.length()-3) + String("acqus"));
        
        Internal::FidHandler fid(filename);
        if (!fid)
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}      
				
        //  Delete old spectrum
				spectrum.clear(true);
				
				//temporary variables
				PeakType p;
							
        while (spectrum.size() < acqus.getSize())
        {
				  //fill peak
				  p.setPosition( (typename PeakType::PositionType) acqus.getPosition(fid.getIndex()) );
				  p.setIntensity( (typename PeakType::IntensityType) fid.getIntensity() );
				  spectrum.push_back(p);
        }
        fid.close();
        
        // import metadata
        spectrum.setRT(0.0);
        spectrum.setMSLevel(1);
        spectrum.setName("Xmass analysis file " + acqus.getParam("$ID_raw"));
        spectrum.setType(SpectrumSettings::RAWDATA);
        spectrum.setNativeID("spectrum=xsd:" + acqus.getParam("$ID_raw").remove('<').remove('>'));
        spectrum.setComment("no comment");
        
        InstrumentSettings instrument_settings;
        instrument_settings.setScanMode(InstrumentSettings::MASSSPECTRUM);
        instrument_settings.setZoomScan(false);

        if(acqus.getParam(".IONIZATION MODE") == "LD+")
				{
          instrument_settings.setPolarity(IonSource::POSITIVE);
				}
        else if (acqus.getParam(".IONIZATION MODE") == "LD-")
				{
        	instrument_settings.setPolarity(IonSource::NEGATIVE);
				}
        else 
				{
					instrument_settings.setPolarity(IonSource::POLNULL);
				}
        spectrum.setInstrumentSettings(instrument_settings);
        
        AcquisitionInfo acquisition_info;
        acquisition_info.setMethodOfCombination("Sum of " + acqus.getParam("$NoSHOTS") + " raw spectrum");
        spectrum.setAcquisitionInfo(acquisition_info);
          
        SourceFile source_file;
        source_file.setNameOfFile("fid");
        source_file.setPathToFile(filename.prefix(filename.length()-3));
        source_file.setFileSize(4.0 * acqus.getSize() / 1024 / 1024); // 4 bytes / point
        source_file.setFileType("Xmass analysis file (fid)");
        spectrum.setSourceFile(source_file);
        
        DataProcessing data_processing;
        Software software;
        software.setName("FlexControl");
        String fc_ver = acqus.getParam("$FCVer"); // FlexControlVersion
        if(fc_ver.hasPrefix("<flexControl "))
				{
					fc_ver = fc_ver.suffix(' ');
				}
        if(fc_ver.hasSuffix(">"))
				{
					fc_ver = fc_ver.prefix('>');
				}
        software.setVersion(fc_ver);
        software.setMetaValue("Acquisition method", DataValue(acqus.getParam("$ACQMETH").remove('<').remove('>')));
       	data_processing.setSoftware(software);
        std::set<DataProcessing::ProcessingAction> actions;
        actions.insert(DataProcessing::SMOOTHING);
        actions.insert(DataProcessing::BASELINE_REDUCTION);
        actions.insert(DataProcessing::CALIBRATION);
        data_processing.setProcessingActions(actions);
        data_processing.setCompletionTime(DateTime::now());
          
        std::vector<DataProcessing> data_processing_vector;
          data_processing_vector.push_back(data_processing);
        spectrum.setDataProcessing(data_processing_vector);
      }

			/**
				@brief Import settings from a XMass file.

				@param exp has to be a MSExperiment or have the same interface.

				@exception Exception::FileNotFound is thrown if the file could not be opened
			*/              
      template <class PeakType>
      void importExperimentalSettings(const String& filename, MSExperiment<PeakType>& exp)
      {
        Internal::AcqusHandler acqus(filename.prefix(filename.length()-3) + String("acqus"));
  
        ExperimentalSettings& experimental_settings = exp.getExperimentalSettings();
        
        Instrument& instrument = experimental_settings.getInstrument();
        instrument.setName(acqus.getParam("SPECTROMETER/DATASYSTEM"));
        instrument.setVendor(acqus.getParam("ORIGIN"));
        instrument.setModel(acqus.getParam("$InstrID").remove('<').remove('>'));
          
        std::vector<IonSource>& ionSourceList = instrument.getIonSources();
        ionSourceList.clear();
        ionSourceList.resize(1);
        if (acqus.getParam(".INLET") == "DIRECT")
				{
        	ionSourceList[0].setInletType(IonSource::DIRECT);
				}
        else
				{
        	ionSourceList[0].setInletType(IonSource::INLETNULL);
          ionSourceList[0].setIonizationMethod(IonSource::MALDI);
				}
        if (acqus.getParam(".IONIZATION MODE") == "LD+")
				{
        	ionSourceList[0].setPolarity(IonSource::POSITIVE);
				}
        else if (acqus.getParam(".IONIZATION MODE") == "LD-")
				{
        	ionSourceList[0].setPolarity(IonSource::NEGATIVE);
				}
        else 
				{
        	ionSourceList[0].setPolarity(IonSource::POLNULL);
				}
        ionSourceList[0].setMetaValue("MALDI target reference", DataValue(acqus.getParam("$TgIDS").remove('<').remove('>')));
        ionSourceList[0].setOrder(0);
            
        std::vector<MassAnalyzer>& massAnalyzerList = instrument.getMassAnalyzers();
        massAnalyzerList.clear();
        massAnalyzerList.resize(1);  
        if(acqus.getParam(".SPECTROMETER TYPE") == "TOF")
				{
        	massAnalyzerList[0].setType(MassAnalyzer::TOF);
				}
        else
				{
        	massAnalyzerList[0].setType(MassAnalyzer::ANALYZERNULL);
				}
        
        DateTime date;
        date.set(acqus.getParam("$AQ_DATE").remove('<').remove('>') );
        experimental_settings.setDateTime(date);
      }

			/**
				@brief Stores a spectrum in a XMass file (not avaible)

				@exception Exception::FileNotWritable is thrown
			*/
      template <typename SpectrumType>
      void store(const String& /*filename*/, const SpectrumType& /*spectrum*/)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_XMASSFILE_H

