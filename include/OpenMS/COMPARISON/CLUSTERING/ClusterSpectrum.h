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
// $Maintainer: Andreas Bertsch  $
// --------------------------------------------------------------------------

#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERSPECTRUM_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERSPECTRUM_H

#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/PeptideHit.h>

namespace OpenMS
{
	class DBAdapter;
	
	/** @brief This class represents a spectrum used for clustering with additional information associated with

  this class allows the use of Spectra without worrying about the
  used Representation<br>
  usually a ClusterSpectrum should only contain 1 Representation
  if the other Representation is requested 2 things can happen:<br>
   - an Exception is thrown (wrong Representation)<br>
   - the appropriate Representation is created, this is done if a DBAdapter*<br>
     is given at Construction ( and needed ).
     (needs Database support)<br>
  note that ClusterSpectrum takes Possession of the given pointers (except DBAdapter)
  they are deleted on Destruction
  */
  class ClusterSpectrum
  {
  	
  public:

    
    
		/** @brief Exception which is thrown if the spectra a incompatible
		
				this Exception indicates that a the BinnedRep and the DSpectrum in the Constructor dont represent the same Spectrum
		*/
    class DifferentSpectra : public Exception::Base
    {
    public:
      DifferentSpectra(const char* file, int line, const char* function) throw();
      virtual ~DifferentSpectra() throw();
    };

    /** @brief Exception which is thrown if the representation is not available
	
		    the requested Representation is not available
    */
    class WrongRepresentation : public Exception::Base
    {
    public:
      WrongRepresentation(const char* file, int line, const char* function, const char* message
          = "ClusterSpectrum didnt contain what was requested and no DBAdapter was given at Construction") throw();
      virtual ~WrongRepresentation() throw();
    };

    
		/** @name Constructors and destructors
		*/
		//@{
		/** @brief constructor with database and binning parameters
		
    		allows a more convenient way of getting spectra from the db<br>
    		although it is more efficient to get the spectra in large chunks directly
        with DBAdapter (ca 5% in simple benchmarks(10k spectra)), getting them on demand allows using less memory<br>
    */
    ClusterSpectrum(long id, double binsize_ = 0 , uint binspread_ = 0);

    /// default constructor
    ClusterSpectrum();

    /// detailed constructor with PeakSpectrum, DBAdapter and binning parameters
    ClusterSpectrum(const PeakSpectrum& spec, double binsize = 0, uint binspread = 0);

    /// detailed constructor with PeakSpectrum pointer, DBAdapter and binning parameters
    ClusterSpectrum(PeakSpectrum* specp, double binsize = 0, uint binspread = 0);

    /// detailed constructor with BinnedRep pointer and DBAdapter
    ClusterSpectrum(BinnedRep* binrepp);

    /// detailed constructor with PeakSpectrum pointer and BinnedRep pointer
    ClusterSpectrum(PeakSpectrum* specp, BinnedRep* binrepp);

    /// copy constructor
    ClusterSpectrum(const ClusterSpectrum& source);

    /// destructor
    virtual ~ClusterSpectrum();
		//@}

		/** @name Accessors
		*/
		//@{
    /// assignment operator
    ClusterSpectrum& operator=(const ClusterSpectrum& source);

		/// returns the id of the spectrum
    int id() const;

		/// returns the retention time of the spectrum
    const double& getRetention() const;

		/// returns the parent mass of the parent ion
    const double& getParentMass() const;

		/// returns the charge of the parent ion
    const uint& getParentionCharge() const;

		/// returns the binning size of the binned spectrum
    const double& getBinSize() const {return binsize_;}

		/// returns the spreading of the binned spectrum
    const uint& getBinSpread() const { return binspread_;}

		/// return the top hit
    PeptideHit getTophit() const;

		/// returns the binned Representation of the spectrum
    const BinnedRep& getBinrep() const;

		/// returns the peak spectrum 
    const PeakSpectrum& getSpec() const;

    /// mutable access to the peak spectrum
    PeakSpectrum& spec();

    /// delete pointers to stick and bin spectrum
    void strip() const;

    /// access to peptide annotations
    const std::vector<Identification>& getIdentification() const;
		//@}

  private:

    void updatecache_() const;
    void clearcache_() const;

    //declared mutable so the appropriate representation can created inside
    //const accessor
    mutable PeakSpectrum* specp_;
    mutable BinnedRep* binrepp_;
    double binsize_;
    uint binspread_;
    long id_;

    // caches for global information
    mutable bool cached_;
    mutable double retention_;
    mutable double parent_mass_;
    mutable uint parentioncharge_;

  };

}
#endif //OPENMS_COMPARISON_CLUSTERING_CLUSTERSPECTRUM_H
