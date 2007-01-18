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

#ifndef OPENMS_FORMAT_DTAFILE_H
#define OPENMS_FORMAT_DTAFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <fstream>
#include <vector>

namespace OpenMS
{
 	/**
 		@brief File adapter for DTA files.
 		
  	The first line contains the singly protonated peptide mass (MH+) and the peptide charge state separated by a space.
  	Subsequent lines contain space separated pairs of fragment ion m/z and intensity values.
  	
  	From precusor mass and charge state the mass-charge-ratio is calculated and stored in the spectrum as precursor mass.
  	
  	@ingroup FileIO
  */
  class DTAFile
  {
    public:
      /// Default constructor
      DTAFile();
			/// Destructor
      ~DTAFile();
      
      /**
 				@brief Loads a DTA file to a spectrum.
 				
 				The content of the file is stored in @p spectrum.
 				@p spectrum has to be a DSpectrum<1>/MSSpectrum<1> or have the same interface.
      */
      template <typename SpectrumType>
      void load(const String& filename, SpectrumType& spectrum) throw (Exception::FileNotFound,Exception::ParseError)
      {
				std::ifstream is(filename.c_str());
				if (!is)
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}
		
				//  Delete old spectrum
				spectrum.getContainer().clear();
				
				//temporary variables
				String line;
				std::vector<String> strings(2);
				typename SpectrumType::PeakType p;
				char delimiter;
				
				//read first line and store precursor m/z and charge
				getline(is,line,'\n');
				line.trim();
		
				//test which delimiter is used in the line
				if (line.has('\t'))
				{
					delimiter = '\t';
				}
				else
				{
					delimiter = ' ';
				}
				
				try
				{
					line.split(delimiter,strings);
					if (strings.size()!=2)
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" ,filename);
					}
					double mh_mass = strings[0].toDouble();
					SignedInt charge = strings[1].toInt();
					if (charge != 0)
					{
						spectrum.getPrecursorPeak().getPosition()[0] = (mh_mass - 1.0) / charge + 1.0;
					}
					else
					{
						spectrum.getPrecursorPeak().getPosition()[0] = mh_mass;
					}
					spectrum.getPrecursorPeak().setCharge(charge);
				}
				catch(...)
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" ,filename);
				}
				
		    while (getline(is,line,'\n'))
		    {
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
		
					try 
					{
						line.split(delimiter,strings);
						if (strings.size()!=2)
						{
							throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" ,filename);
						}				
						
						//fill peak
						p.setPosition(strings[0].toDouble());
						p.setIntensity(strings[1].toDouble());
					} 
					catch ( Exception::Base & e )
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" ,filename);
					}
					spectrum.getContainer().push_back(p);
				}
				
				is.close();  	
      }

      /**
      	@brief Stores a spectrum in a DTA file.
      	
      	The content of @p spectrum is stored in a file.
      	@p spectrum has to be a DSpectrum<1>/MSSpectrum<1> or have the same interface.
      */
      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum) const throw (Exception::UnableToCreateFile)
      {
				std::ofstream os(filename.c_str());
				if (!os)
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}
		
				// Write mh+ mass
				if (spectrum.getPrecursorPeak().getCharge()==0)
				{
					//unknown charge
					os << spectrum.getPrecursorPeak().getPosition()[0];
				}
				else
				{
					//known charge
					os << ((spectrum.getPrecursorPeak().getPosition()[0] - 1.0) * spectrum.getPrecursorPeak().getCharge() +1.0);
				}
				 
				//charge
				os << " " << spectrum.getPrecursorPeak().getCharge() << std::endl;
		
				// Iterate over all peaks of the spectrum and
				// write one line for each peak of the spectrum.
				typename SpectrumType::ConstIterator it(spectrum.begin());
				for (; it != spectrum.end(); ++it)
				{
					// Write m/z and intensity.
					os << it->getPosition()[0] << " " << it->getIntensity() << std::endl;
				}
		
				// Done.
				os.close();		
      }
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_DTAFILE_H

