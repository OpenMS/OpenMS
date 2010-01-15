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
// $Maintainer: Andreas Bertsch$
// $Authors: Andreas Bertsch, Daniel Jameson$
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTREMOTEQUERY_H
#define OPENMS_FORMAT_MASCOTREMOTEQUERY_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtCore/QFile>
#include <QtNetwork/QHttpRequestHeader>


namespace OpenMS
{
	/**
		@brief Class which handles the communication between OpenMS and the Mascot server

		This class provides a communication interface which is able to query the Mascot
		server and reports the identifications provided be the Mascot server

		@htmlinclude OpenMS_MascotRemoteQuery.parameters

	*/
	class OPENMS_DLLAPI MascotRemoteQuery 
		: public QObject,
			public DefaultParamHandler
	{
		Q_OBJECT
	
		public:
		
			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			MascotRemoteQuery(QObject *parent=0);
	
			/// destructor
			virtual ~MascotRemoteQuery();		
			//@}


			/// sets the query spectra, given in MGF file format
			void setQuerySpectra(const String& exp);
			
			/// returns the Mascot XML response which contains the identifications
			const QByteArray& getMascotXMLResponse() const;
		
			/// predicate which returns true if an error occurred during the query
			bool hasError() const;

			/// returns the error message, if hasError can be used to check whether an error has occurred
			const String& getErrorMessage() const;	
			
		protected:

			virtual void updateMembers_();
	
		public slots:

			void run();
		
		private slots:

			/** slot connected to signal requestFinished of QHttp: "This signal is emitted 
				  when processing the request identified by id has finished. error is true 
					if an error occurred during the processing; otherwise error is false"
			*/
			void httpRequestFinished(int request_id, bool error);

			/// slot connected to signal dataReadProgress of QHttp
			void httpDataReadProgress(int bytes_read, int bytes_total);

			/// slot connected to signal dataSendProgress of QHttp
			void httpDataSendProgress(int bytes_sent, int bytes_total);

			/// slot connected to signal requestStarted of QHttp, which indicates that the processing of request request_id has been started
			void httpRequestStarted(int request_id);

			/** slot connected to signal stateChanged of QHttp, which is emitted if 
		 			the http state changed. See 'enum QHttp::State' of Qt docu for more 
					info.
			*/
			void httpStateChanged(int state);

			/// slot connected to signal done of QHttp
			void httpDone(bool error);

			/// slot connect to responseHeaderRecieved, which indicates that a new response header is available
			void readResponseHeader(const QHttpResponseHeader& response_header);

			void login();
		
			void execQuery();

			void getResults();
		
		
		signals:
		
			void done();
		
			void loginDone();

			void queryDone();

		private:

			String query_spectra_;

			QByteArray mascot_xml_;

			QHttp* http_;

			QString results_path_;

			QString cookie_;

			String error_message_;

			/// assignment operator
      MascotRemoteQuery& operator = (const MascotRemoteQuery& rhs);
			
			/// copy constructor
      MascotRemoteQuery(const MascotRemoteQuery& rhs);

};

}
#endif /*OPENMS_FORMAT_MASCOTREMOTEQUERY_H*/
