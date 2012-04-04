// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTGENERICFILE_H
#define OPENMS_FORMAT_MASCOTGENERICFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <fstream>

namespace OpenMS
{
	/**
		@brief Mascot input file adapter.
		
		Creates a file that can be used for Mascot search from a peak list or a whole experiment.

		@htmlinclude OpenMS_MascotGenericFile.parameters
	
  	@ingroup FileIO
	*/
  class OPENMS_DLLAPI MascotGenericFile
		: public ProgressLogger,
			public DefaultParamHandler
  {
    public:

			/// constructor
			MascotGenericFile();

			/// constructor
			virtual ~MascotGenericFile();

			/// stores the experiment data in a MascotGenericFile that can be used as input for MASCOT shell execution
			void store(const String& filename, const PeakMap& experiment);

			/// store the experiment data in a MascotGenericFile; the output is written to the given stream, the filename will be noted in the file
			void store(std::ostream& os, const String& filename, const PeakMap& experiment);

			/** loads a Mascot Generic File into a PeakMap
					
					@param filename file name which the map should be read from
					@param exp the map which is filled with the data from the given file
					@throw FileNotFound is thrown if the given file could not be found
			*/
			template <typename MapType> void load(const String& filename, MapType& exp)
      {
				if (!File::exists(filename))
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}

				exp.reset();

				std::ifstream is(filename.c_str());
				std::vector<std::pair<double, double> > spec;
				UInt charge(0);
				double pre_mz(0), pre_int(0), rt(-1);
				String title;
				UInt spectrum_number = 0;
				Size line_number=0;
				while (getNextSpectrum_(is, spec, charge, pre_mz, pre_int, rt, title, line_number))
				{
					typename MapType::SpectrumType spectrum;
					for (std::vector<std::pair<double, double> >::const_iterator it = spec.begin(); it != spec.end(); ++it)
					{
						typename MapType::PeakType p;
						p.setPosition(it->first);
						p.setIntensity(it->second);
						spectrum.push_back(p);
					}
					spectrum.setMSLevel(2);
					spectrum.getPrecursors().resize(1);
					spectrum.getPrecursors()[0].setMZ(pre_mz);
					spectrum.getPrecursors()[0].setIntensity(pre_int);
					spectrum.getPrecursors()[0].setCharge(charge);
					spectrum.setRT(rt);
					if (title != "")
					{
						spectrum.setMetaValue("TITLE", title);
						title = "";
					}

					spectrum.setNativeID(String("index=")+(spectrum_number++));
					exp.push_back(spectrum);
					
					// clean up
					spec.clear();
					charge = 0;
					pre_mz = 0;
					pre_int = 0;
				}
      }

      /**
        @brief enclosing Strings of the peak list body for HTTP submission
        
          Can be used to embed custom content into HTTP submission (when writing only the MGF header in HTTP format and then
          adding the peaks (in whatever format, e.g. mzXML) enclosed in this body.
          The @p filename can later be found in the Mascot response.
      */
      std::pair<String,String> getHTTPPeakListEnclosure(const String& filename) const;

    protected:

			/// writes a parameter header
			void writeParameterHeader_(const String& name, std::ostream& os);

			/// writes the full header
			void writeHeader_(std::ostream& os);
			
			/// writes the spectrum
			void writeSpectrum_(std::ostream& os, const PeakSpectrum& spec);
						
			/// writes the MSExperiment
			void writeMSExperiment_(std::ostream& os, const String& filename, const PeakMap& experiment);

			/// reads a spectrum block, the section between 'BEGIN IONS' and 'END IONS' of a mgf file
			bool getNextSpectrum_(std::istream& is, std::vector<std::pair<double, double> >& spectrum, UInt& charge, double& precursor_mz, double& precursor_int, double& rt, String& title, Size& line_number);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTGENERICFILE_H
