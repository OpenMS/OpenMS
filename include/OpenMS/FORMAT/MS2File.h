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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MS2FILE_H
#define OPENMS_FORMAT_MS2FILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>

#include <vector>
#include <fstream>

namespace OpenMS
{
	/**
		@brief MS2 input file adapter.

		For the format description take a look at:
		Rapid Commun Mass Spectrom. 2004;18(18):2162-8.

		MS1, MS2, and SQT-three unified, compact, and easily parsed file formats for the
		storage of shotgun proteomic spectra and identifications.

		McDonald WH, Tabb DL, Sadygov RG, MacCoss MJ, Venable J, Graumann J, Johnson JR,
		Cociorva D, Yates JR 3rd.

		PMID: 15317041

  	@ingroup FileIO
	*/
  class OPENMS_DLLAPI MS2File
		: public ProgressLogger
  {
    public:

			/// constructor
			MS2File();

			/// constructor
			virtual ~MS2File();

			template <typename MapType> void load(const String& filename, MapType& exp)
			{
	      //startProgress(0,0,"loading DTA2D file");

      	if (!File::exists(filename))
      	{
        	throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      	}
				if (!File::readable(filename))
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}

				exp.reset();

				//set DocumentIdentifier
				exp.setLoadedFileType(filename);
				exp.setLoadedFilePath(filename);

				std::ifstream in(filename.c_str());

				UInt spectrum_number = 0;
				typename MapType::SpectrumType spec;
      	typename MapType::SpectrumType::PeakType p;

				String line;
				bool first_spec(true);
				while (getline(in, line, '\n'))
				{
					line.trim();
					if (line.size() == 0) continue;

					// header
					if (line[0] == 'H')
					{
						continue;
					}

					// scan
					if (line[0] == 'S')
					{
						if (!first_spec)
						{
							spec.setMSLevel(2);
							spec.setNativeID(String("index=")+(spectrum_number++));
							exp.push_back(spec);
						}
						else
						{
							first_spec = false;
						}
						spec.clear(true);
						line.simplify();
						std::vector<String> split;
						line.split(' ', split);
						if (split.size() != 4)
						{
							throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "line '" + line  + "' should contain four values!", "");
						}
						spec.getPrecursors().resize(1);
						spec.getPrecursors()[0].setMZ(split[3].toDouble());
						continue;
					}

					// charge-independent analysis
					if (line[0] == 'I')
					{
						continue;
					}

					// charge specification
					if (line[0] == 'Z')
					{
						continue;
					}

					// charge-dependent analysis
					if (line[0] == 'D')
					{
						continue;
					}

					// yet another peak, hopefully
					line.simplify();
					std::vector<String> split;
					line.split(' ', split);
					if (split.size() != 2)
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "line '" + line  + "' should contain two values!", "");
					}

					// TODO catch exceptions
					p.setPosition(split[0].toDouble());
					p.setIntensity(split[1].toFloat());
					spec.push_back(p);
				}

				if (!first_spec)
				{
					spec.setMSLevel(2);
					spec.setNativeID(String("index=")+(spectrum_number++));
					exp.push_back(spec);
				}
			}

			/*
			template <typename MapType> void store(const String& filename, MapType& map)
			{

			}
			*/

    protected:

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MS2FILE_H
