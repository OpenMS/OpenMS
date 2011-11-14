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

#ifndef OPENMS_FORMAT_FILETYPES_H
#define OPENMS_FORMAT_FILETYPES_H

#include <OpenMS/config.h>

namespace OpenMS
{
	/**
		@brief Centralizes the file types recognized by FileHandler.

		FileType separate from FileHandler to avoid circular inclusions by DocumentIdentifier, ExperimentalSettings and FileHandler and respective fileclasses (e.g. DTA2DFile). See also: FileHandler::nameToType, FileHandler::typeToName and FileHandler::NameOfTypes .

		@ingroup FileIO
	*/
	struct FileTypes
	{
		//NOTE: if you change/add something here, do not forget to change FileHandler::NameOfTypes[]

		///Actual file types enum.
		enum Type
		{
			UNKNOWN,        		///< Unknown file extension
			DTA,            		///< DTA file (.dta)
			DTA2D,          		///< DTA2D file (.dta2d)
			MZDATA,         		///< MzData file (.MzData)
			MZXML,          		///< MzXML file (.MzXML)
			FEATUREXML,     		///< %OpenMS feature file (.featureXML)
			IDXML,  						///< %OpenMS identification format (.idXML)
			CONSENSUSXML,  			///< %OpenMS consensus map format (.consensusXML)
			MGF,								///< Mascot Generic Format (.mgf)
			INI,          			///< %OpenMS parameters file (.ini)
      TOPPAS,          		///< %OpenMS parameters file with workflow information (.toppas)
			TRANSFORMATIONXML,  ///< Tranformation description file (.trafoXML)
			MZML,								///< MzML file (.mzML)
			MS2,								///< MS2 file (.ms2)
			PEPXML,							///< TPP pepXML file (.pepXML)
			PROTXML,						///< TPP protXML file (.protXML)
			MZIDENTML,				  ///< mzIdentML (HUPO PSI AnalysisXML format) (.mzid)
			GELML,							///< GelML (HUPO PSI format) (.GelML)
			TRAML,							///< TraML (HUPO PSI format) for transitions (.TraML)
			MSP,								///< NIST spectra library file format (.msp)
			OMSSAXML,						///< OMSSA XML file format for peptide identifications (.xml)
			MASCOTXML,						///< Mascot XML file format for peptide identifications (.xml)
			PNG,                ///< Portable Network Graphics (.png)
			XMASS,              ///< XMass Analysis file (fid)
      TSV,                ///< msInspect file (.tsv)
      PEPLIST,            ///< specArray file (.pepList)
      HARDKLOER,          ///< hardkloer file (.hardkloer)
      KROENIK,            ///< kroenik file (.kroenik)
      FASTA,              ///< FASTA file (.fasta)
      EDTA,               ///< enhanced comma separated files (RT, m/z, Intensity, [meta])
			SIZE_OF_TYPE    		///< No file type. Simply stores the number of types
		};
	};

} //namespace OpenMS

#endif //OPENMS_FORMAT_FILETYPES_H
