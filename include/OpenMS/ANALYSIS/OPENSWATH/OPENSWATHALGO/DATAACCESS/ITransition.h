/*
 * ITransition.h
 *
 *  Created on: Jul 27, 2012
 *      Author: witold
 */

#ifndef ITRANSITION_H_
#define ITRANSITION_H_
#include <vector>
#include <boost/shared_ptr.hpp>
namespace OpenSwath{
  // Datastructures for Scoring
  class IFeature
  {
    public:
    virtual ~IFeature(){};
    virtual void getRT(std::vector<double> & rt) = 0;
    virtual void getIntensity(std::vector<double> & intens) = 0;
    virtual float getIntensity() = 0;
    virtual double getRT() = 0;
  };

  class IMRMFeature
  {
    public:
    virtual ~IMRMFeature(){};
    virtual boost::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID) = 0;
    virtual float getIntensity() = 0;
    virtual double getRT() = 0;
  };

  struct ITransitionGroup
  {
    virtual ~ITransitionGroup() {};
    virtual std::size_t size() = 0;
    virtual std::vector<std::string> getNativeIDs() = 0;
    virtual void getLibraryIntensities(std::vector<double> & intensities) = 0;
  };

  struct ISignalToNoise
  {
    virtual ~ISignalToNoise() {};
    virtual double getValueAtRT(double RT) = 0;
  };
  typedef boost::shared_ptr<ISignalToNoise> ISignalToNoisePtr;


} //end Namespace OpenSwath



#endif /* ITRANSITION_H_ */
