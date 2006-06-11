// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: OpenMS.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
//
//
//

/////////////////////////
// OpenMS PYTHON SUPPORT
/////////////////////////
//
// module definition file
//

%Module OpenMS

%Include copying.sip


// std
%Include std__string.sip
%Include std__exception.sip

// CONCEPT
%Include Exception.sip
%Include MetaInfo.sip
%Include MetaInfoRegistry.sip
%Include MetaInfoInterface.sip
%Include TimeStamp.sip
%Include Types.sip
%Include VersionInfo.sip

// KERNEL
%Include DataValue.sip
%Include DPeak.sip         // template DPeak<D> expands to DPeak1, DPeak2, DPeak3
%Include DRawDataPoint.sip // template DRawDataPoint<D> expands to DRDP1...3
%Include DPosition.sip 		 // template DPosition<D> expands to DPosition1
%Include DDiscreteSpectrum.sip
%Include DSpectrumInfo.sip
%Include DPeakArray.sip

// FORMAT
%Include PersistenceManager.sip
%Include PersistenceManagerInterface.sip
%Include DTAFile.sip
