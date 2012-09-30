/*
 * ISpectrumAccess.h
 *
 *  Created on: Jul 11, 2012
 *      Author: wolski
 */

#ifndef ISPECTRUMACCESS_H_
#define ISPECTRUMACCESS_H_

#include "DataStructures.h"
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

namespace OpenSwath
{

  class ISpectrumAccess
  {
  public:
    virtual ~ISpectrumAccess(){};
    virtual SpectrumPtr getSpectrumById(int id) const = 0;
    virtual std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const = 0;
    virtual size_t getNrSpectra() const = 0;
    virtual SpectrumMeta getSpectrumMetaById(int id) const = 0;

    virtual ChromatogramPtr getChromatogramById(int id) const = 0;
    //virtual std::vector<std::size_t> getChromatogramByPrecursorMZ(double mz,
    //    double deltaMZ) const = 0;
    virtual std::size_t getNrChromatograms() const = 0;
    virtual std::string getChromatogramNativeID(int id) const = 0;
  };

  typedef boost::shared_ptr<ISpectrumAccess> SpectrumAccessPtr;
}

#endif /* ISPECTRUMACCESS_H_ */

