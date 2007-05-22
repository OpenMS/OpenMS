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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MSPFILE_H
#define OPENMS_FORMAT_MSPFILE_H

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <cctype>
#include <fstream>

namespace OpenMS
{
	/**
		@brief File adapter for MzData files
	
		@todo add ProgressLogger to File? (Andreas)
		@todo is TextFile a good idea? (Andreas)
		@ingroup FileIO
	*/
	class MSPFile /*
			public ProgressLogger*/
	{
		public:

			///Default constructor
			MSPFile();

			///Destructor
			virtual ~MSPFile();
			
			/**
				@brief Loads a map from a MSPFile file.

				@p map has to be a MSExperiment or have the same interface.
			*/
			template <typename MapType>
			void load(const String& filename, std::vector<PeptideIdentification>& ids, MapType& map) throw (Exception::FileNotFound, Exception::ParseError)
			{
				if (!File::exists(filename))
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}
				if (!File::readable(filename))
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}
				String line;
				std::ifstream is(filename.c_str());

				typename MapType::SpectrumType spec;
				while (getline(is, line))
				{
					if (line.hasPrefix("Name:"))
					{
						std::vector<String> split, split2;
						line.split(' ', split);
						split[1].split('/', split2);
						PeptideIdentification id;
						id.insertHit(PeptideHit(0, "", 0, split2[1].toInt(), split2[0]));
						ids.push_back(id);
					}
					if (line.hasPrefix("MW:"))
					{
						// skip that as it is not necessary and might not be available at all
					}
					if (line.hasPrefix("Comment:"))
					{
						//spec.setMetaValue("MSPComment", line);
					}
					if (line.hasPrefix("Num peaks:"))
					{
						while (getline(is, line) && line.size() > 0 && std::isdigit(line[0]))
						{
							std::vector<String> split;
							line.split('\t', split);
							typename MapType::SpectrumType::PeakType peak;
							if (split.size() != 3)
							{
								throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, line, "not <mz> <intensity> <comment>");
							}
							peak.setMZ(split[0].toFloat());
							peak.setIntensity(split[1].toFloat());
							//peak.setMetaValue("MSPPeakInfo", split[2]);
							spec.push_back(peak);
						}
						
						map.push_back(spec);
						spec.clear();
					}
				}
			}

			/**
				@brief Stores a map in a MSPFile file.

				@p map has to be a MSExperiment or have the same interface.
			*/
			template <typename MapType>
			void store(const String& filename, const MapType& map)
			const throw (Exception::UnableToCreateFile)
			{
				// TODO
			}
		
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_MSPFILE_H
