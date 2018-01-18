// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  class FeatureFinder;

  namespace Internal
  {
    //-------------------------------------------------------------------

    /**
        @brief Comparator that allows to compare the indices of two peaks by their intensity.
    */
    template <typename FeaFiModuleType>
    struct IntensityLess :
      std::binary_function<typename FeatureFinderDefs::IndexPair, typename FeatureFinderDefs::IndexPair, bool>
    {
      /// Constructor that takes a FeaFiModule reference
      IntensityLess(const FeaFiModuleType & module) :
        module_(module)
      {
      }

      /// Copy ctor
      IntensityLess(const IntensityLess & rhs) :
        module_(rhs.module_)
      {
      }

      /// Compare with respect to intensity
      inline bool operator()(const typename FeatureFinderDefs::IndexPair & left, const typename FeatureFinderDefs::IndexPair & right) const
      {
        return module_.getPeakIntensity(left) < module_.getPeakIntensity(right);
      }

private:
      /// Reference to the FeaFiModule
      const FeaFiModuleType & module_;
      /// Default ctor undefined since we cannot compare without a FeaFiModule.
      IntensityLess();
    };     // struct IntensityLess

    //-------------------------------------------------------------------

    ///Intensity iterator for a FeatureFinderDefs::IndexSet
    template <typename FeaFiModuleType>
    struct IntensityIterator :
      FeatureFinderDefs::IndexSet::const_iterator
    {
      IntensityIterator(const FeatureFinderDefs::IndexSet::const_iterator & iter, const FeaFiModuleType * module) :
        FeatureFinderDefs::IndexSet::const_iterator(iter),
        module_(module)
      {
      }

      typename FeaFiModuleType::IntensityType operator*() const
      {
        return module_->getPeakIntensity(FeatureFinderDefs::IndexSet::const_iterator::operator*());
      }

protected:
      const FeaFiModuleType * module_;
    };

    //-------------------------------------------------------------------

    ///m/z iterator for a FeatureFinderDefs::IndexSet
    template <typename FeaFiModuleType>
    struct MzIterator :
      FeatureFinderDefs::IndexSet::const_iterator
    {
      MzIterator(const FeatureFinderDefs::IndexSet::const_iterator & iter, const FeaFiModuleType * module) :
        FeatureFinderDefs::IndexSet::const_iterator(iter),
        module_(module)
      {
      }

      typename FeaFiModuleType::IntensityType operator*() const
      {
        return module_->getPeakMz(FeatureFinderDefs::IndexSet::const_iterator::operator*());
      }

protected:
      const FeaFiModuleType * module_;
    };

    //-------------------------------------------------------------------

    ///Retention time iterator for a FeatureFinderDefs::IndexSet
    template <typename FeaFiModuleType>
    struct RtIterator :
      FeatureFinderDefs::IndexSet::const_iterator
    {
      RtIterator(const FeatureFinderDefs::IndexSet::const_iterator & iter, const FeaFiModuleType * module) :
        FeatureFinderDefs::IndexSet::const_iterator(iter),
        module_(module)
      {
      }

      typename FeaFiModuleType::IntensityType operator*() const
      {
        return module_->getPeakRt(FeatureFinderDefs::IndexSet::const_iterator::operator*());
      }

protected:
      const FeaFiModuleType * module_;
    };

    //-------------------------------------------------------------------
  }   // namespace Internal

  /**
  @brief Implements a module of the FeatureFinder algorithm.
  */
  template <class PeakType>
  class FeaFiModule :
    public DefaultParamHandler
  {
public:
    ///Input map type
    typedef PeakMap MapType;
    ///Input spectrum type
    typedef typename MapType::SpectrumType SpectrumType;
    ///Input intensity type
    typedef typename PeakType::IntensityType IntensityType;
    ///Input coordinate type
    typedef typename PeakType::CoordinateType CoordinateType;

    ///Constructor
    FeaFiModule(const PeakMap * map, FeatureMap* features, FeatureFinder * ff) :
      DefaultParamHandler("FeaFiModule"),
      map_(nullptr),
      features_(nullptr),
      ff_(nullptr)
    {
      map_ = map;
      features_ = features;
      ff_ = ff;
    }

    /// destructor
    ~FeaFiModule() override
    {
    }

    /// Returns the intensity of a peak
    inline IntensityType getPeakIntensity(const FeatureFinderDefs::IndexPair & index) const
    {
      //Corrupt index
      OPENMS_PRECONDITION(index.first < map_->size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(index.second < (*map_)[index.first].size(), "Peak index outside of scan!");

      return (*map_)[index.first][index.second].getIntensity();
    }

    /// Returns the m/z of a peak
    inline CoordinateType getPeakMz(const FeatureFinderDefs::IndexPair & index) const
    {
      //Corrupt index
      OPENMS_PRECONDITION(index.first < map_->size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(index.second < (*map_)[index.first].size(), "Peak index outside of scan!");

      return (*map_)[index.first][index.second].getMZ();
    }

    /// Returns the retention time of a peak
    inline CoordinateType getPeakRt(const FeatureFinderDefs::IndexPair & index) const
    {
      //Corrupt index
      OPENMS_PRECONDITION(index.first < map_->size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(index.second < (*map_)[index.first].size(), "Peak index outside of scan!");

      return (*map_)[index.first].getRT();
    }

    /**
    @brief fills @p index with the index of next peak in m/z dimension

    @exception FeatureFinderDefs::NoSuccessor is thrown if there is no next peak
    @exception Exception::Precondition is thrown if an invalid index is given
    */
    inline void getNextMz(FeatureFinderDefs::IndexPair & index) const
    {
      //Corrupt index
      OPENMS_PRECONDITION(index.first < map_->size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(index.second < (*map_)[index.first].size(), "Peak index outside of scan!");

      //At the last peak of this spectrum
      if (index.second + 1 >= (*map_)[index.first].size())
      {
        throw FeatureFinderDefs::NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getNextMz", index);
      }

      ++index.second;
    }

    /**
    @brief fills @p index with the index of previous peak in m/z dimension

    @exception FeatureFinderDefs::NoSuccessor is thrown if there is no previous peak
    @exception Exception::Precondition is thrown if an invalid index is given
    */
    inline void getPrevMz(FeatureFinderDefs::IndexPair & index) const
    {
      //Corrupt index
      OPENMS_PRECONDITION(index.first < map_->size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(index.second < (*map_)[index.first].size(), "Peak index outside of scan!");

      //begin of scan
      if (index.second == 0)
      {
        throw FeatureFinderDefs::NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getPrevMz", index);
      }

      --index.second;
    }

    /**
    @brief fills @p index with the index of the nearest peak in the next scan

    @exception FeatureFinderDefs::NoSuccessor is thrown if there is no next spectrum
    @exception Exception::Precondition is thrown if an invalid index is given
    */
    void getNextRt(FeatureFinderDefs::IndexPair & index)
    {
      //Corrupt index
      OPENMS_PRECONDITION(index.first  < map_->size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(index.second < (*map_)[index.first].size(), "Peak index outside of scan!");

      CoordinateType mz_pos = (*map_)[index.first][index.second].getMZ();       // mz value we want to find
      Size index_first_tmp = index.first;

      ++index.first;
      while (index.first < map_->size() &&
             (*map_)[index.first].empty())
      {
        ++index.first;
      }
      //last scan
      if (index.first >= map_->size())
      {
        throw FeatureFinderDefs::NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getNextRt", index);
      }
      // now we have a spectrum with scans in it ...

      // perform binary search to find the neighbour in mz dimension
      typename SpectrumType::ConstIterator it = lower_bound((*map_)[index.first].begin(), (*map_)[index.first].end(), (*map_)[index_first_tmp][index.second], typename PeakType::PositionLess());

      // if the found peak is at the end of the spectrum, there is not much we can do...
      if (it == (*map_)[index.first].end())
      {
        index.second = (*map_)[index.first].size() - 1;
      }
      // if the found peak is at the beginning of the spectrum, there is also not much we can do !
      else if (it == (*map_)[index.first].begin())
      {
        index.second = 0;
      }
      // see if the next smaller one fits better
      else
      {
        // peak to the right is closer (in m/z dimension)
        if (it->getMZ() - mz_pos < mz_pos - (it - 1)->getMZ())
        {
          index.second = it - (*map_)[index.first].begin();
        }
        else            // left one is closer
        {
          index.second = --it - (*map_)[index.first].begin();
        }
      }
    }

    /**
    @brief fills @p index with the index of the nearest peak in the previous scan

    @exception FeatureFinderDefs::NoSuccessor is thrown if there is no previous spectrum
    @exception Exception::Precondition is thrown if an invalid index is given
    */
    void getPrevRt(FeatureFinderDefs::IndexPair & index)
    {
      //Corrupt index
      OPENMS_PRECONDITION(index.first < map_->size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(index.second < (*map_)[index.first].size(), "Peak index outside of scan!");

      // TODO: this seems useless (at least for debug mode) given preconditions above... (and why not in getNextRt()??)
      if (index.first >= map_->size())
      {
        std::cout << "Scan index outside of map!" << std::endl;
        std::cout << index.first << " " << index.second << std::endl;
        return;
      }
      if (index.second >= (*map_)[index.first].size())
      {
        std::cout << "Peak index outside of scan!" << std::endl;
        std::cout << index.first << " " << index.second << std::endl;
        return;
      }

      CoordinateType mz_pos = (*map_)[index.first][index.second].getMZ();
      Size index_first_tmp = index.first;

      // first scan
      if (index.first == 0)
      {
        throw FeatureFinderDefs::NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getPrevRt", index);
      }

      --index.first;
      while ((index.first > 0) && ((*map_)[index.first].empty()))
      {
        --index.first;
      }
      // we only found an empty scan
      if ((*map_)[index.first].empty()) throw FeatureFinderDefs::NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getPrevRt", index);

      // perform binary search to find the neighbour in mz dimension

      typename MapType::SpectrumType::ConstIterator it = lower_bound((*map_)[index.first].begin(),
                                                                     (*map_)[index.first].end(),
                                                                     (*map_)[index_first_tmp][index.second],
                                                                     typename PeakType::PositionLess());

      // if the found peak is at the end of the spectrum, there is not much we can do.
      if (it == (*map_)[index.first].end())
      {
        index.second = (*map_)[index.first].size() - 1;
      }
      // if the found peak is at the beginning of the spectrum, there is not much we can do.
      else if (it == (*map_)[index.first].begin())
      {
        index.second = 0;
      }
      // see if the next smaller one fits better
      else
      {
        // peak to the right is closer (in m/z dimension)
        if (it->getMZ() - mz_pos < mz_pos - (it - 1)->getMZ())
        {
          index.second = it - (*map_)[index.first].begin();
        }
        else
        {
          index.second = --it - (*map_)[index.first].begin();
        }
      }
    }

    ///Calculates the convex hull of a index @p set and adds it to the @p feature
    void addConvexHull(const FeatureFinderDefs::IndexSet & set, Feature & feature) const
    {
      std::vector<DPosition<2> > points;
      points.reserve(set.size());
      DPosition<2> tmp;
      for (FeatureFinderDefs::IndexSet::const_iterator it = set.begin(); it != set.end(); ++it)
      {
        tmp[Peak2D::MZ] = (*map_)[it->first][it->second].getMZ();
        tmp[Peak2D::RT] = (*map_)[it->first].getRT();
        points.push_back(tmp);
      }
      feature.getConvexHulls().resize(feature.getConvexHulls().size() + 1);
      // computes convex hull
      feature.getConvexHulls().back().addPoints(points);
    }

protected:
    ///Input data pointer
    const MapType * map_;
    ///Output data pointer
    FeatureMap * features_;
    ///Pointer to the calling FeatureFinder that is used to access the feature flags and report progress
    FeatureFinder * ff_;

private:
    /// Not implemented
    FeaFiModule();
    /// Not implemented
    FeaFiModule & operator=(const FeaFiModule &);
    /// Not implemented
    FeaFiModule(const FeaFiModule &);

  };   // class FeaFiModule

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H
