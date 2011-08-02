// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_DTA2DFILE_H
#define OPENMS_FORMAT_DTA2DFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <fstream>
#include <iostream>

namespace OpenMS
{
  class String;
  /**
			@brief DTA2D File adapter.

			File adapter for files with three tab/space-separated columns.

			The default format is: retention time (seconds) , m/z , intensity.

			If the first line starts with '#', a different order is defined by the
			the order of the keywords 'MIN' (retention time in minutes) or 'SEC' (retention time in seconds), 'MZ', and 'INT'.

			Example: '\#MZ MIN INT'

			The peaks of one retention time have to be in subsequent lines.

			@ingroup FileIO
  */
  class OPENMS_DLLAPI DTA2DFile
  	: public ProgressLogger
  {
	 private:
		PeakFileOptions options_;

	 public:

		/** @name Constructors and Destructor */
		//@{
		/// Default constructor
		DTA2DFile();
		/// Destructor
		~DTA2DFile();
		//@}

		/// Mutable access to the options for loading/storing
		PeakFileOptions& getOptions();

		/// Non-mutable access to the options for loading/storing
		const PeakFileOptions& getOptions() const;

		/**
			 @brief Loads a map from a DTA2D file.

			 @p map has to be a MSExperiment or have the same interface.

			@exception Exception::FileNotFound is thrown if the file could not be opened
			@exception Exception::ParseError is thrown if an error occurs during parsing
		*/
		template <typename MapType>
		void load(const String& filename, MapType& map)
		{
			startProgress(0,0,"loading DTA2D file");

			//try to open file
			std::ifstream is(filename.c_str());
			if (!is)
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
			}

			map.reset();

			//set DocumentIdentifier
			map.setLoadedFileType(filename);
			map.setLoadedFilePath(filename);

			// temporary variables to store the data in
			std::vector<String> strings(3);
			typename MapType::SpectrumType spec;
			spec.setRT(-1.0); //to make sure the first RT is different from the the initialized value
			typename MapType::SpectrumType::PeakType p;
			DoubleReal rt(0.0);
			char delimiter;

			// default dimension of the data
			Size rt_dim = 0;
			Size mz_dim = 1;
			Size int_dim = 2;

			//RT unit (default is seconds)
			bool time_in_minutes = false;

			// string to store the current line in
			String line;

			// native ID (numbers from 0)
			UInt native_id = 0;

			// line number counter
			Size line_number = 0;

			while (getline(is,line,'\n'))
			{
				++line_number;
				line.trim();

				if ( line.empty() ) continue;

				//test which delimiter is used in the line
				if (line.has('\t'))
				{
					delimiter = '\t';
				}
				else
				{
					delimiter = ' ';
				}

				//is header line
				if (line.hasPrefix("#"))
				{
					line = line.substr(1).trim();
					line.split(delimiter,strings);

					// flags to check if dimension is set correctly
					bool rt_set = false;
					bool mz_set = false;
					bool int_set = false;

					//assign new order
					for (Size i = 0 ; i<3;++i)
					{
						if ( strings[i]=="RT" || strings[i]=="RETENTION_TIME" || strings[i]=="MASS-TO-CHARGE" || strings[i]=="IT" || strings[i]=="INTENSITY")
						{
							std::cerr << "Warning: This file contains the deprecated keyword '" << strings[i] << "'." << "\n";
							std::cerr << "         Please use only the new keywords SEC/MIN, MZ, INT." << "\n";
						}
						if ( ( strings[i]=="SEC" || strings[i]=="RT" || strings[i]=="RETENTION_TIME" ) && rt_set==false)
						{
							rt_dim = i;
							rt_set = true;
						}
						else if ( ( strings[i]=="MIN") && rt_set==false)
						{
							rt_dim = i;
							rt_set = true;
							time_in_minutes = true;
						}
						else if ( ( strings[i]=="MZ" || strings[i]=="MASS-TO-CHARGE" ) && mz_set==false)
						{
							mz_dim = i;
							mz_set = true;
						}
						else if ( ( strings[i]=="INT" || strings[i]=="IT" || strings[i]=="INTENSITY" ) && int_set==false )
						{
							int_dim = i;
							int_set = true;
						}
						else
						{
							throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Misformatted header line!" ,filename);
						}
					}
					continue;
				}

				try
				{
					line.split(delimiter,strings);
					if (strings.size()!=3)
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line (" + String(line_number) + "): \"")+line+"\" (got  " + String(strings.size()) + ", expected 3 entries)" ,filename);
					}
					p.setIntensity(strings[int_dim].toFloat());
					p.setMZ(strings[mz_dim].toDouble());
					rt = (strings[rt_dim].toDouble()) * (time_in_minutes ? 60.0 : 1.0);
				}
				// conversion to double or something else could have gone wrong
				catch ( Exception::BaseException & /*e*/ )
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line (" + String(line_number) + "): \"")+line+"\"" ,filename);
				}

				// Retention time changed -> new Spectrum
				if (fabs(rt-spec.getRT())>0.0001)
				{
					if ( spec.size()!=0
						 	 &&
						 	 (!options_.hasRTRange() || options_.getRTRange().encloses(DPosition<1>(spec.getRT())))) // RT restriction fulfilled
					{
						map.push_back(spec);
					}
					setProgress(0);
					spec.clear(true);
					spec.setRT(rt);
					spec.setNativeID(String("index=")+native_id);
					++native_id;
				}

				//Skip peaks with invalid m/z or intensity value
				if (
						(!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(p.getMZ())))
						&&
						(!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(p.getIntensity())))
					 )
				{
					spec.push_back(p);
				}
			}

			// add last Spectrum
			if (
				  spec.size()!=0
				 	&&
				 	(!options_.hasRTRange() || options_.getRTRange().encloses(DPosition<1>(spec.getRT()))) // RT restriction fulfilled
				 )
			{
				map.push_back(spec);
			}

			is.close();
			endProgress();
		}
    
		/**
     @brief Stores a map in a DTA2D file.
     
     @p map has to be a MSExperiment or have the same interface.
     
     @exception Exception::UnableToCreateFile is thrown if the file could not be created
     */
		template <typename MapType>
		void store(const String& filename, const MapType& map) const
		{
			startProgress(0,map.size(),"storing DTA2D file");
      
			std::ofstream os(filename.c_str());
			if (!os)
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
			}
      
      // write header
      os	<< "RT" << "\t" << "MZ" << "\t" << "INT" << "\n";
      
			// Iterate over all peaks of each spectrum and
			// write one line for each peak of the spectrum.
			UInt count = 0;
			for (typename MapType::const_iterator spec=map.begin(); spec!=map.end(); ++spec)
			{
				setProgress(count++);
				for (typename MapType::SpectrumType::ConstIterator it = spec->begin(); it != spec->end(); ++it)
				{
					// Write rt, m/z and intensity.
					os	<< precisionWrapper(spec->getRT()) << "\t" << precisionWrapper(it->getPos()) << "\t" << precisionWrapper(it->getIntensity()) << "\n";
				}
        
			}
			os.close();
			endProgress();
		}
    
		/**
     @brief Stores the TIC of a map in a DTA2D file.
     
     @p map has to be a MSExperiment or have the same interface.
     
     @exception Exception::UnableToCreateFile is thrown if the file could not be created
     */
		template <typename MapType>
		void storeTIC(const String& filename, const MapType& map) const
		{
			startProgress(0,map.size(),"storing DTA2D file");
      
			std::ofstream os(filename.c_str());
			if (!os)
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
			}
      
      // write header (Always MZ=0 for chromatograms in DTA2D.)
      os	<< "RT" << "\t" << "MZ" << "\t" << "INT" << "\n";
      
      typename MapType::ChromatogramType TIC = map.getTIC();
      for (typename MapType::ChromatogramType::ConstIterator it = TIC.begin(); it != TIC.end(); ++it)
      {
        // write rt and intensity.
        os	<< precisionWrapper(it->getRT()) << "\t" << precisionWrapper(0) << "\t" << precisionWrapper(it->getIntensity()) << "\n";
      }
      
      
			// Iterate over spectra (retention times) and write summed intensities for each spectrum, i.e. the total ion chromatogram (TIC)
			/*UInt count = 0;
			for (typename MapType::const_iterator spec=map.begin(); spec!=map.end(); ++spec)
			{
				setProgress(count++);
        if (spec->getMSLevel() == 1)
        {
          DoubleReal totalIntensity = 0;
          for (typename MapType::SpectrumType::ConstIterator it = spec->begin(); it != spec->end(); ++it)
          {
            // sum intensities of each spectrum
            totalIntensity = totalIntensity + it->getIntensity();
          }
          // Write rt and intensity.
          os	<< precisionWrapper(spec->getRT()) << "\t" << precisionWrapper(0) << "\t" << precisionWrapper(totalIntensity) << "\n";
        }
			}*/
			os.close();
			endProgress();
		}
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_DTA2DFILE_H
