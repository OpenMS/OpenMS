// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost, Michał Startek, Mateusz Łącki $
// --------------------------------------------------------------------------

#pragma once

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>


// Override IsoSpec's use of mmap whenever it is available
#define ISOSPEC_GOT_SYSTEM_MMAN false
#define ISOSPEC_GOT_MMAN false

#include <OpenMS/../../thirdparty/IsoSpec/IsoSpec/isoSpec++.h>


namespace OpenMS
{

  /**
    * @brief Interface for the IsoSpec algorithm - a generator of infinitely-resolved theoretical spectra.
    *
    * Provides an interface to the IsoSpec algorithm. See IsoSpecWrapper for a
    * more convenient wrapper. Implements a generator pattern using the
    * nextConf function, which can be used to iterate through all
    * configurations:
    *
    * @code
    * IsoSpecGeneratorWrapperSubclass iso(...);
    *
    * while(iso.nextConf())
    * {
    *     Peak1D conf = iso.getConf(); // and/or getMass, getIntensity, etc.
    *     // do some computations on that conf;
    * }
    * @endcode
    *
    *
    * The computation is based on the IsoSpec algorithm
    *
    * @code
    * Łącki MK, Startek M, Valkenborg D, Gambin A.
    * IsoSpec: Hyperfast Fine Structure Calculator.
    * Anal Chem. 2017 Mar 21;89(6):3272-3277. doi: 10.1021/acs.analchem.6b01459.
    * @endcode
    *
    **/
  class OPENMS_DLLAPI IsoSpecGeneratorWrapper
  {

public:

    /**
     * @brief Move the generator to a next isotopologue
     *
     * Advance the internal generator to the next isotopologue. The value returned determines whether the
     * generator has been exhausted (that is, all eligible configurations have already been visited).
     * It is invalid to call any other generator methods before the first call to nextConf(), as well as
     * after this method returns false.
     *
     * @returns A boolean value stating whether the generator has been exhausted.
     */
    virtual bool nextConf() = 0;

    /**
     * @brief Obtain the current isotopologue
     *
     * @return The current isotopologue as a Peak1D
     *
     * @note It is invalid (undefined results) to call this method before the first call to nextConf(), or after it returns false
     */
    virtual Peak1D getConf() = 0;

    /**
     * @brief Obtain the mass of the current isotopologue
     *
     * @return The mass of the current isotopologue
     *
     * @note It is invalid (undefined results) to call this method before the first call to nextConf(), or after it returns false
     */
    virtual double getMass() = 0;

    /**
     * @brief Obtain the intensity (probability, relative peak height) of the current configuration
     *
     * @return The intensity (probability) of the current isotopologue
     *
     * @note It is invalid (undefined results) to call this method before the first call to nextConf(), or after it returns false
     */
    virtual double getIntensity() = 0;

    /**
     * @brief Obtain the natural logarithm of the intensity (probability, relative peak height) of the current configuration
     *
     * This will be more precise (and faster) than just calling std::log(getIntensity()) - it will produce correct results even
     * for configurations so unlikely that the double-precision floating point number returned from getIntensity() underflows to zero.
     *
     * @return The natural logarithm of intensity (probability) of the current isotopologue
     *
     * @note It is invalid (undefined results) to call this method before the first call to nextConf(), or after it returns false
     */
    virtual double getLogIntensity() = 0;

    /**
     * @brief Destructor
     */
    virtual inline ~IsoSpecGeneratorWrapper() {};
  };

  /** @brief A convenience class for the IsoSpec algorithm - easier to use than the IsoSpecGeneratorWrapper classes.
   *
   */
  class OPENMS_DLLAPI IsoSpecWrapper
  {
public:
    /**
      * @brief Run the algorithm
      *
      * This method will run the algorithm with parameters as set up by the constructor. It will return an
      * IsotopeDistribution containing the observed configurations. The configurations are explicitly stored
      * in memory, which may become a problem when considering some especially large distributions. If this,
      * or (a rather small) performance overhead is a concern, then the
      * generator methods (see IsoSpecGeneratorWrapper) should be used instead.
      *
      * This method is provided for convenience. As calling that method invalidates the object (the method should
      * not be called again, nor anything other than destroying the object should be done with it), the most common
      * usage pattern of IsoSpecGeneratorWrapper classes with the run method is:
      *
      * @code
      * IsotopeDistribution dist = IsoSpecGeneratorWrapperSubclass(...).run();
      * // do something with dist;
      * @endcode
      *
      * @note Calling this method invalidates the object! In future versions this limitation might be removed.
      *
      **/
    virtual IsotopeDistribution run() = 0;

    virtual inline ~IsoSpecWrapper() {};
  };

  //-------------------------------------------------------------------------- 
  // IsoSpecGeneratorWrapper classes
  //-------------------------------------------------------------------------- 

  /**
   * @brief Generate a p-set of configurations for a given p (that is, a set of configurations such that
   *        their probabilities sum up to p). The p in normal usage will usually be close to 1 (e.g. 0.99).
   *
   * An optimal p-set of isotopologues is the smallest set of isotopologues that, taken together, cover at
   * least p of the probability space (that is, their probabilities sum up to at least p). This means that
   * the computed spectrum is accurate to at least degree p, and that the L1 distance between the computed
   * spectrum and the true spectrum is less than 1-p. The optimality of the p-set means that it contains
   * the most probable configurations - any isotopologues outside of the returned p-set have lower intensity
   * than the configurations in the p-set.
   *
   * This is the method most users will want: the p parameter directly controls the accuracy of results.
   *
   * Advanced usage note: The algorithm works by computing an optimal p'-set for a p' slightly larger than
   * the requested p. By default these extra isotopologues are returned too (as they have to be computed
   * anyway). It is possible to request that the extra configurations be discarded, using the do_p_trim
   * parameter. This will *increase* the runtime and especially the memory usage of the algorithm, and
   * should not be done unless there is a good reason to.
   *
   * @note The eligible configurations are NOT guaranteed to be returned in any particular order.
   */
  class OPENMS_DLLAPI IsoSpecTotalProbGeneratorWrapper : public IsoSpecGeneratorWrapper
  {
public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param p Total coverage of probability space desired, usually close to 1 (e.g. 0.99)
      * @param do_p_trim Whether to discard extra configurations that have been computed
      *
      * @note This constructor is only useful if you need to define non-standard abundances
      *       of isotopes, for other uses the one accepting EmpiricalFormula is easier to use.
      *
      **/
    IsoSpecTotalProbGeneratorWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities,
             double p,
             bool do_p_trim = false);

    /**
      * @brief Setup the algorithm to run on an EmpiricalFormula
      *
      **/
    IsoSpecTotalProbGeneratorWrapper(const EmpiricalFormula& formula, double p, bool do_p_trim = false);

    virtual inline bool nextConf() override final { return ILG.advanceToNextConfiguration(); };
    virtual inline Peak1D getConf() override final { return Peak1D(ILG.mass(), ILG.prob()); };
    virtual inline double getMass() override final { return ILG.mass(); };
    virtual inline double getIntensity() override final { return ILG.prob(); };
    virtual inline double getLogIntensity() override final { return ILG.lprob(); };

protected:
    IsoSpec::IsoLayeredGenerator ILG;
  };

  /**
   * @brief Provides a threshold-based generator of isotopologues: generates all isotopologues
   *        more probable than a given probability threshold.
   *
   * This is the simplest generator - most users will however want to use IsoSpecTotalProbGeneratorWrapper.
   * The reason for it is that when thresholding by peak intensity one has no idea how far the obtained
   * spectrum is from a real spectrum. For example, consider human insulin: if the threshold is set at
   * 0.1 of the most intense peak, then the peaks above the threshold account for 99.7% of total probability
   * for water, 82% for substance P, only 60% for human insulin, and 23% for titin.
   * For a threshold of 0.01, the numbers will be: still 99.7% for water, 96% for substance P, 88% for human insulin
   * and 72% for titin (it also took 5 minutes on an average notebook computer to process the 17 billion configurations
   * involved).
   *
   * As you can see the threshold does not have a straightforward correlation to the accuracy of the final spectrum
   * obtained - and accuracy of final spectrum is often what the user is interested in. The IsoSpeTotalcProbGeneratorWrapper
   * provides a way to directly parametrise based on the desired accuracy of the final spectrum - and should be used
   * instead in most cases. The trade-off is that it's (slightly) slower than Threshold algorithm. This speed gap will
   * be dramatically improved with IsoSpec 2.0.
   *
   * @note The eligible isotopologues are NOT guaranteed to be generated in any particular order.
   */
  class OPENMS_DLLAPI IsoSpecThresholdGeneratorWrapper : public IsoSpecGeneratorWrapper
  {

public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param threshold Intensity threshold: will only compute peaks above this threshold
      * @param absolute Whether the threshold is absolute or relative (relative to the most intense peak)
      *
      * @note This constructor is only useful if you need to define non-standard abundances
      *       of isotopes, for other uses the one accepting EmpiricalFormula is easier to use.
      *
      **/
  IsoSpecThresholdGeneratorWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities,
             double threshold,
             bool absolute);

    /**
      * @brief Setup the algorithm to run on an EmpiricalFormula
      *
      **/
  IsoSpecThresholdGeneratorWrapper(const EmpiricalFormula& formula, double threshold, bool absolute);

  virtual inline bool nextConf() override final { return ITG.advanceToNextConfiguration(); };
  virtual inline Peak1D getConf() override final { return Peak1D(ITG.mass(), ITG.prob()); };
  virtual inline double getMass() override final { return ITG.mass(); };
  virtual inline double getIntensity() override final { return ITG.prob(); };
  virtual inline double getLogIntensity() override final { return ITG.lprob(); };


protected:
  IsoSpec::IsoThresholdGenerator ITG;
  };

  /**
   * @brief Generate the stream of configurations, ordered from most likely to least likely.
   *
   * This generator walks through the entire set of isotopologues of a given molecule, allowing
   * the user to terminate the search on the fly. The returned isotopologues are guaranteed to
   * be generated in order of descending probability (unlike the previous two generators which
   * make no such guarantees).
   *
   * This causes the algorithm to run in O(n*log(n)) and means that is it much slower than the
   * previous two.
   *
   * You should only use this generator if you don't know up-front when to stop the walk through
   * the configuration space, and need to make the decision on the fly. If you know the threshold
   * or the total probability needed, and only need the configurations sorted, it will be much
   * faster to generate them using one of the previous algorithms and sort them afterwards.
   */
  class OPENMS_DLLAPI IsoSpecOrderedGeneratorWrapper : public IsoSpecGeneratorWrapper
  {
public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      *
      **/
  IsoSpecOrderedGeneratorWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities);

    /**
      * @brief Setup the algorithm to run on an EmpiricalFormula
      *
      **/
  IsoSpecOrderedGeneratorWrapper(const EmpiricalFormula& formula);

  virtual inline bool nextConf() override final { return IOG.advanceToNextConfiguration(); };
  virtual inline Peak1D getConf() override final { return Peak1D(IOG.mass(), IOG.prob()); };
  virtual inline double getMass() override final { return IOG.mass(); };
  virtual inline double getIntensity() override final { return IOG.prob(); };
  virtual inline double getLogIntensity() override final { return IOG.lprob(); };

protected:
  IsoSpec::IsoOrderedGenerator IOG;
  };

  //-------------------------------------------------------------------------- 
  // IsoSpecWrapper classes
  //-------------------------------------------------------------------------- 

/**
  * @brief Create a p-set of configurations for a given p (that is, a set of configurations such that
  *        their probabilities sum up to p). The p in normal usage will usually be close to 1 (e.g. 0.99).
  *
  * An optimal p-set of isotopologues is the smallest set of isotopologues that, taken together, cover at
  * least p of the probability space (that is, their probabilities sum up to at least p). This means that
  * the computed spectrum is accurate to at least degree p, and that the L1 distance between the computed
  * spectrum and the true spectrum is less than 1-p. The optimality of the p-set means that it contains
  * the most probable configurations - any isotopologues outside of the returned p-set have lower intensity
  * than the configurations in the p-set.
  *
  * This is the method most users will want: the p parameter directly controls the accuracy of results.
  *
  * Advanced usage note: The algorithm works by computing an optimal p'-set for a p' slightly larger than
  * the requested p. By default these extra isotopologues are returned too (as they have to be computed
  * anyway). It is possible to request that the extra configurations be discarded, using the do_p_trim
  * parameter. This will *increase* the runtime and especially the memory usage of the algorithm, and
  * should not be done unless there is a good reason to.
  *
  * @note The eligible configurations are NOT guaranteed to be returned in any particular order.
  */
  class OPENMS_DLLAPI IsoSpecTotalProbWrapper : public IsoSpecWrapper
  {
public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param p Total coverage of probability space desired, usually close to 1 (e.g. 0.99)
      * @param do_p_trim Whether to discard extra configurations that have been computed
      *
      * @note This constructor is only useful if you need to define non-standard abundances
      *       of isotopes, for other uses the one accepting EmpiricalFormula is easier to use.
      *
      **/
    IsoSpecTotalProbWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities,
             double p,
             bool do_p_trim = false);

    /**
      * @brief Setup the algorithm to run on an EmpiricalFormula
      *
      **/
    IsoSpecTotalProbWrapper(const EmpiricalFormula& formula, double p, bool do_p_trim = false);

    virtual IsotopeDistribution run() override final;

protected:
    IsoSpec::IsoLayeredGenerator ILG;

  };

  /**
    * @brief A non-generator version of IsoSpecThresholdGeneratorWrapper
    *
    * This is the simplest algorithm - most users will however want to use IsoSpecTotalProbWrapper.
    * The reason for it is that when thresholding by peak intensity one has no idea how far the obtained
    * spectrum is from a real spectrum. For example, consider human insulin: if the threshold is set at
    * 0.1 of the most intense peak, then the peaks above the threshold account for 99.7% of total probability
    * for water, 82% for substance P, only 60% for human insulin, and 23% for titin.
    * For a threshold of 0.01, the numbers will be: still 99.7% for water, 96% for substance P, 88% for human insulin
    * and 72% for titin (it also took 5 minutes on an average notebook computer to process the 17 billion configurations
    * involved).
    *
    * As you can see the threshold does not have a straightforward correlation to the accuracy of the final spectrum
    * obtained - and accuracy of final spectrum is often what the user is interested in. The IsoSpeTotalcProbGeneratorWrapper
    * provides a way to directly parametrise based on the desired accuracy of the final spectrum - and should be used
    * instead in most cases. The trade-off is that it's (slightly) slower than Threshold algorithm. This speed gap will
    * be dramatically improved with IsoSpec 2.0.
    *
    * @note The eligible isotopologues are NOT guaranteed to be generated in any particular order.
    */
  class OPENMS_DLLAPI IsoSpecThresholdWrapper : public IsoSpecWrapper
  {

public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param threshold Intensity threshold: will only compute peaks above this threshold
      * @param absolute Whether the threshold is absolute or relative (relative to the most intense peak)
      *
      * @note This constructor is only useful if you need to define non-standard abundances
      *       of isotopes, for other uses the one accepting EmpiricalFormula is easier to use.
      *
      **/
    IsoSpecThresholdWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities,
             double threshold,
             bool absolute);

    /**
      * @brief Setup the algorithm to run on an EmpiricalFormula
      *
      **/
    IsoSpecThresholdWrapper(const EmpiricalFormula& formula, double threshold, bool absolute);

    virtual IsotopeDistribution run() override final;

protected:
    IsoSpec::IsoThresholdGenerator ITG;

  };

}

