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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTREMOTEQUERY_H
#define OPENMS_FORMAT_MASCOTREMOTEQUERY_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtCore/QFile>
#include <QtNetwork/QHttpRequestHeader>
#include <QTimer>


namespace OpenMS
{
	/**
		@brief Class which handles the communication between OpenMS and the Mascot server

		This class provides a communication interface which is able to query the Mascot
		server and reports the identifications provided be the Mascot server

		@htmlinclude OpenMS_MascotRemoteQuery.parameters

	*/
	class MascotRemoteQuery 
		: public QObject,
			public DefaultParamHandler
	{
		Q_OBJECT
	
		public:
		
			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			OPENMS_DLLAPI MascotRemoteQuery(QObject *parent=0);
	
			/// destructor
			OPENMS_DLLAPI virtual ~MascotRemoteQuery();		
			//@}


			/// sets the query spectra, given in MGF file format
			OPENMS_DLLAPI void setQuerySpectra(const String& exp);
			
			/// returns the Mascot XML response which contains the identifications
			OPENMS_DLLAPI const QByteArray& getMascotXMLResponse() const;
		
			/// predicate which returns true if an error occurred during the query
			OPENMS_DLLAPI bool hasError() const;

			/// returns the error message, if hasError can be used to check whether an error has occurred
			OPENMS_DLLAPI const String& getErrorMessage() const;	
			
		protected:

			OPENMS_DLLAPI virtual void updateMembers_();
	
		public slots:

			OPENMS_DLLAPI void run();

		private slots:

      OPENMS_DLLAPI void timedOut();
      
      OPENMS_DLLAPI void readyReadSlot ( const QHttpResponseHeader & resp );

			/** slot connected to signal requestFinished of QHttp: "This signal is emitted 
				  when processing the request identified by id has finished. error is true 
					if an error occurred during the processing; otherwise error is false"
			*/
			OPENMS_DLLAPI void httpRequestFinished(int request_id, bool error);

			/// slot connected to signal dataReadProgress of QHttp
			OPENMS_DLLAPI void httpDataReadProgress(int bytes_read, int bytes_total);

			/// slot connected to signal dataSendProgress of QHttp
			OPENMS_DLLAPI void httpDataSendProgress(int bytes_sent, int bytes_total);

			/// slot connected to signal requestStarted of QHttp, which indicates that the processing of request request_id has been started
			OPENMS_DLLAPI void httpRequestStarted(int request_id);

			/** slot connected to signal stateChanged of QHttp, which is emitted if 
		 			the http state changed. See 'enum QHttp::State' of Qt docu for more 
					info.
			*/
			OPENMS_DLLAPI void httpStateChanged(int state);

			/// slot connected to signal done of QHttp
			OPENMS_DLLAPI void httpDone(bool error);

			/// slot connect to responseHeaderRecieved, which indicates that a new response header is available
			OPENMS_DLLAPI void readResponseHeader(const QHttpResponseHeader& response_header);

			OPENMS_DLLAPI void login();
		
			OPENMS_DLLAPI void execQuery();

			OPENMS_DLLAPI void getResults();
		
      OPENMS_DLLAPI void loginSuccess();

		signals:
		
			OPENMS_DLLAPI void done();
		
			OPENMS_DLLAPI void loginDone();

			OPENMS_DLLAPI void queryDone();

		private:
			/// assignment operator
      OPENMS_DLLAPI MascotRemoteQuery& operator = (const MascotRemoteQuery& rhs);
			/// copy constructor
      OPENMS_DLLAPI MascotRemoteQuery(const MascotRemoteQuery& rhs);

      OPENMS_DLLAPI void endRun_();

			String query_spectra_;
			QByteArray mascot_xml_;
			QHttp* http_;
			QString results_path_;
			QString cookie_;
			String error_message_;
      QTimer timeout_;
      Int to_;
};

}
#endif /*OPENMS_FORMAT_MASCOTREMOTEQUERY_H*/
