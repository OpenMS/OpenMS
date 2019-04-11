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
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ISOTOPEDISTRIBUTION_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ISOTOPEDISTRIBUTION_H

#include <OpenMS/KERNEL/Peak1D.h>

#include <vector>

namespace OpenMS
{
  /**
        @ingroup Chemistry

        @brief Isotope distribution class

        A container that holds an isotope distribution. It consists from mass values and its 
        correspondent probabilities. The container mass values does not relate to the 
        atomic number of the of the molecule. It is a floating precision number and 
        it maps a specific atomic mass to an abundance of an isotope.

        To bypass this behavior the CoarseIsotopePatternGenerator can be used. It employs
        an Coarse Isotopic model so the calculations happen with atomic numbers instead of
        atomic masses.
        The CoarseIsotopePatternGenerator solver quantizes the atomic masses to integer 
        numbers that correspond to the atomic number. Then the calculation of the 
        IsotopeDistribution can produce nominal isotopes with accurate or rounded masses

    */
  class Element; 

  class OPENMS_DLLAPI IsotopeDistribution
  {
public:

    /// @name typedefs
    //@{
    /// container type, first holds the weight of the isotope, second the probability
    typedef Peak1D MassAbundance;
    typedef std::vector<MassAbundance> ContainerType;
    typedef ContainerType::iterator iterator;
    typedef ContainerType::iterator Iterator;
    typedef ContainerType::const_iterator const_iterator;
    typedef ContainerType::const_iterator ConstIterator;

    typedef ContainerType::reverse_iterator reverse_iterator;
    typedef ContainerType::reverse_iterator ReverseIterator;
    typedef ContainerType::const_reverse_iterator const_reverse_iterator;
    typedef ContainerType::const_reverse_iterator ConstReverseIterator;
    //@}

    enum Sorted {INTENSITY, MASS, UNDEFINED};

    /// @name Constructors and Destructors
    //@{
    /** Default constructor
    */
    IsotopeDistribution();

    /// Copy constructor
    IsotopeDistribution(const IsotopeDistribution&) = default;

    /// Move constructor
    IsotopeDistribution(IsotopeDistribution&&) noexcept = default;

    /// Destructor
    virtual ~IsotopeDistribution() = default;
    //@}

    /// @name Accessors
    //@{

    /// overwrites the container which holds the distribution using @p distribution
    void set(const ContainerType & distribution);

    void set(ContainerType && distribution);

    /// returns the container which holds the distribution
    const ContainerType & getContainer() const;

    /// returns the maximal weight isotope which is stored in the distribution
    Peak1D::CoordinateType getMax() const;

    /// returns the minimal weight isotope which is stored in the distribution
    Peak1D::CoordinateType getMin() const;

    /// returns the most abundant isotope which is stored in the distribution
    Peak1D getMostAbundant() const;

    /// returns the size of the distribution which is the number of isotopes in the distribution
    Size size() const;

    /// clears the distribution
    void clear();

    // resizes distribution container
    void resize(UInt size);

    /// remove intensities below the cutoff
    void trimIntensities(double cutoff);

    /// sort isotope distribution by intensity
    void sortByIntensity();

    /// sort isotope distribution by mass
    void sortByMass();

    /** @brief Re-normalizes the sum of the probabilities of all the isotopes to 1

        The re-normalisation may be needed as in some distributions with a lot
        of isotopes the calculations tend to be inexact.
    */
    void renormalize();

     /** @brief Merges distributions arbitrary data points with constant defined resolution.
        
        It creates a new IsotopeDistribution Container and assigns each isotope to the nearest bin.
        This function should be used to downsample the existing distribution.
        If the size of the new Container is larger this function throws an IllegalArgument Exception.
        
     */
    void merge(double resolution, double min_prob);

    /** @brief Trims the right side of the isotope distribution to isotopes with a significant contribution.

        If the isotope distribution is calculated for large masses, it might
        happen that many entries contain only small numbers. This function can
        be used to remove these entries.

        @note Consider normalising the distribution afterwards.
    */
    void trimRight(double cutoff);

    /** @brief Trims the left side of the isotope distribution to isotopes with a significant contribution.

        If the isotope distribution is calculated for large masses, it might
        happen that many entries contain only small numbers. This function can
        be used to remove these entries.

        @note Consider normalising the distribution afterwards.
    */
    void trimLeft(double cutoff);

    /// Compute average mass of isotope distribution
    double averageMass() const;
    //@}

    /// @name Operators
    //@{
    /// Assignment operator
    IsotopeDistribution & operator=(const IsotopeDistribution & isotope_distribution);

    /// equality operator, returns true if the @p isotope_distribution is identical to this, false else
    bool operator==(const IsotopeDistribution & isotope_distribution) const;

    /// inequality operator, returns true if the @p isotope_distribution differs from this, false else
    bool operator!=(const IsotopeDistribution & isotope_distribution) const;

    /// less operator
    bool operator<(const IsotopeDistribution & isotope_distribution) const;
    //@}

    /// @name Iterators
    //@{
    inline Iterator begin() { return distribution_.begin(); }

    inline Iterator end()   { return distribution_.end(); }

    inline ConstIterator begin() const { return distribution_.begin(); }

    inline ConstIterator end() const { return distribution_.end(); }

    inline ReverseIterator rbegin() { return distribution_.rbegin(); }

    inline ReverseIterator rend()   { return distribution_.rend(); }

    inline ConstReverseIterator rbegin() const { return distribution_.rbegin(); }

    inline ConstReverseIterator rend() const { return distribution_.rend(); }

    inline void insert(const Peak1D::CoordinateType& mass, const Peak1D::IntensityType& intensity)
    {
      distribution_.push_back(Peak1D(mass, intensity));
    }
    //@}

    /// @name Data Access Operators
    //@{
    /// operator which access a cell of the distribution and wraps it in SpectrumFragment struct
    Peak1D& operator[](const Size& index){ return distribution_[index];}

    //@}


protected:   

    /// sort wrapper of the distribution
    void sort_(std::function<bool(const MassAbundance& p1, const MassAbundance& p2)> sorter);
    /// takes a function as a parameter to transform the distribution
    void transform_(std::function<void(MassAbundance&)> lambda);

    /// stores the isotope distribution
    ContainerType distribution_;
  };


} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ISOTOPEDISTRIBUTION_H
