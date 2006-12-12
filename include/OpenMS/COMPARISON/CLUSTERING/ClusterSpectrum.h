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
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERSPECTRUM_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERSPECTRUM_H

#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/PeptideHit.h>

namespace OpenMS
{
	class DBAdapter;
	
/**
  this class allows the use of Spectra without worrying about the
  used Representation<br>
  usually a ClusterSpectrum should only contain 1 Representation
  if the other Representation is requested 2 things can happen:<Br>
   - an Exception is thrown (wrong Representation)<br>
   - the appropriate Representation is created, this is done if a DBAdapter*<br>
     is given at Construction ( and needed ).
     (needs Database support)<br>
  note that ClusterSpectrum takes Possession of the given pointers ( except DBAdapter )
  they are deleted on Destruction
  */
  class ClusterSpectrum
  {
  	
  public:

    /**
    this Exception indicates that a the BinnedRep and the DSpectrum in the Constructor dont represent the same Spectrum
    */
    class DifferentSpectra : public Exception::Base
    {
    public:
      DifferentSpectra(const char* file, int line, const char* function) throw();
      ~DifferentSpectra() throw();
    };

    /**
    the requested Representation is not available
    */
    class WrongRepresentation : public Exception::Base
    {
    public:
      WrongRepresentation(const char* file, int line, const char* function, const char* message
          = "ClusterSpectrum didnt contain what was requested and no DBAdapter was given at Construction") throw();
      ~WrongRepresentation() throw();
    };

    /**
    allows a more convenient way of getting spectra from the db<br>
    although it is more efficient to get the spectra in large chunks directly
        with DBAdapter (ca 5% in simple benchmarks(10k spectra)), getting them on demand allows using less memory<br>
    */
    ClusterSpectrum(long id, DBAdapter* adapterp , double binsize_ = 0 , uint binspread_ = 0);

    /** @brief standard constructor <br> */
    ClusterSpectrum();

    /** @brief copy spec <br> */
    ClusterSpectrum(const PeakSpectrum& spec, DBAdapter* adapterp = 0, double binsize = 0, uint binspread = 0);

    /** @brief use specp <br> */
    ClusterSpectrum(PeakSpectrum* specp, DBAdapter* adapterp = 0, double binsize = 0, uint binspread = 0);

    /** @brief use binrepp <br> */
    ClusterSpectrum(BinnedRep* binrepp, DBAdapter* adapterp = 0);

    /** @brief use specp and binrepp <br> */
    ClusterSpectrum(PeakSpectrum* specp, BinnedRep* binrepp);

    /** @brief copy constructor <br> */
    ClusterSpectrum(const ClusterSpectrum& source);

    /** @brief destructor <br> */
    virtual ~ClusterSpectrum();

    /** @brief assignment operator <br> */
    ClusterSpectrum& operator=(const ClusterSpectrum& source);

    /** @name read accessors <br> */
    //@{
    int id() const;
    const double& getRetention() const;
    const double& getParentMass() const;
    const uint& getParentionCharge() const;
    const double& getBinSize() const {return binsize_;}
    const uint& getBinSpread() const { return binspread_;}
    PeptideHit getTophit() const;
    const BinnedRep& getBinrep() const;
    const PeakSpectrum& getSpec() const;
    //@}

    /** @brief write accessor for stick spectrum <br> */
    PeakSpectrum& spec();

    /** @brief delete pointers to stick and bin spectrum <br> */
    void strip() const;

    /** @brief access to peptide annotations <br> */
    const std::vector<Identification>& getIdentification() const;
  private:

    void updatecache_() const;
    void clearcache_() const;

    //declared mutable so the appropriate representation can created inside
    //const accessor
    mutable PeakSpectrum* specp_;
    mutable BinnedRep* binrepp_;
    DBAdapter* adapterp_;
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
