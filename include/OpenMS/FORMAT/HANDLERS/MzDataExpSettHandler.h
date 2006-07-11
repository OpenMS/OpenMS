// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZDATAEXPSETTHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZDATAEXPSETTHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>


namespace OpenMS
{
	namespace Internal
	{

  /**
  	@brief XML handler for MzDataFile

		MapType has to be a MSExperiment or have the same interface.
  	Do not use this class. It is only needed in MzDataFile.
  */
  class MzDataExpSettHandler
		: public SchemaHandler
  {
    public:
      /**@name Constructors and destructor */
      //@{
      ///
      MzDataExpSettHandler(ExperimentalSettings& exp);

      ///
      MzDataExpSettHandler(const ExperimentalSettings& exp);

      ///
      virtual ~MzDataExpSettHandler();
      //@}

      ///
      virtual bool endElement( const QString & uri, const QString & local_name,
															 const QString & qname );
      ///
      virtual bool startElement(const QString & uri, const QString & local_name,
																const QString & qname, const QXmlAttributes & attributes );
		  ///
      virtual bool characters( const QString & chars );

  		///Writes the contents to a stream
			void writeTo(std::ostream& os);

    protected:
		/// map pointer for reading
		ExperimentalSettings* exp_;
		/// map pointer for writing
		const ExperimentalSettings* cexp_;

		/** @brief indices for tags used by mzData

			Used to access is_parser_in_tag_.
			If you add tags, also add them to XMLSchemes.h.
			Add no elements to the enum after TAG_NUM.
		*/
		enum Tags { TAGNULL, DESCRIPTION, ADMIN, SAMPLENAME_TAG,
								SAMPLEDESCRIPTION, SRCFILE, NAMEOFFILE,
								PATHTOFILE, FILETYPE,  CONTACT, NAME, CONTACTINST,
								CONTACTINFO, INSTRUMENT, INSTNAME,
								INSTSRC, DETECTOR, ANALYZERLIST, ANALYZER, INSTADDITIONAL,
								DATAPROCESSING, SOFTWARE, SWVERSION,
								PROCMETHOD, COMMENTS, CVPARAM, USERPARAM, TAG_NUM};

		/** @brief indices for ontology terms used by mzData

			If you add terms, also add them to XMLSchemes.h
			Caution: SAMPLENAME_TAG != SAMPLENAME_ONT
		*/
		enum Ontology { ONTNULL, SAMPLENUMBER, SAMPLENAME_ONT, SAMPLESTATE, SAMPLEMASS,
										SAMPLEVOLUME, SAMPLECONC,	INLETTYPE, IONTYPE, IONMODE,
										ANALYZTYPE, RESOLUTION, RESMETHOD, RESTYPE, ACCURACY,
										SCANRATE, SCANTIME, SCANFCT,SCANDIR, SCANLAW,
										TANDEM, REFLECTRON, TOFLENGTH, ISOWIDTH, FINALMSEXP, MAGSTRENGTH,
										DETECTTYPE, ACQMODE, DETECTRES, ADCFREQ,
										VENDOR, MODEL, CUSTOM, DEISOTOPED, DECONVOLVED, PEAKPROC};

		/** @brief indices for enum2str-maps used by mzData

			Used to access enum2str_().
			If you add maps, also add them to XMLSchemes.h.
			Add no elements to the enum after MAP_NUM.
			Each map corresponds to a string in XMLSchemes.h.
		*/
		enum MapTypes {	IONMODEMAP, RESMETHODMAP, RESTYPEMAP, SCANFUNCTIONMAP,
										SCANDIRECTIONMAP, SCANLAWMAP,
										SAMPLESTATEMAP,	PEAKPROCMAP, REFLECTRONMAP,
										ACQMODEMAP, IONTYPEMAP, INLETTYPEMAP, TANDEMMAP, TYPEMAP,
										ANALYZERTYPEMAP, ONTOLOGYMAP, TAGMAP, MAP_NUM};

		/**@name temporary datastructures to hold parsed data */
    //@{
		ContactPerson* contact_;
		MassAnalyzer* analyzer_;
    //@}

		/** @brief read attributes of MzData's cvParamType

			Example:
			&lt;cvParam cvLabel="psi" accession="PSI:1000001" name="@p name" value="@p value"/&gt;
			@p name and sometimes @p value are defined in the MzData ontology.
		*/
		void cvParam_(QString name, QString value);

		/** @brief read attributes of MzData's userParamType

			Example:
			&lt;userParam name="@p name" value="@p value"/&gt;
			@p name and @p value are stored as MetaValues
		*/
		void userParam_(QString name, QString value);
  };

	} // namespace Internal

} // namespace OpenMS

#endif //OPENMS_FORMAT_HANDLERS_MZDATAEXPSETTHANDLER_H
