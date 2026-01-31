


    def __str__(self):
        """
        __str__(self: IsotopeDistribution) -> str
        
        Return a string representation of the IsotopeDistribution object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: IsotopeDistribution) -> str
        
        Return a string representation of the IsotopeDistribution object.

        Returns key properties in a readable format:
        IsotopeDistribution(num_isotopes=5, mass_range=[100.00, 105.00], most_abundant_mass=102.0500)
        """
        num_isotopes = self.size()

        parts = []
        parts.append(f"num_isotopes={num_isotopes}")

        if num_isotopes > 0:
            min_mass = self.getMin()
            max_mass = self.getMax()
            parts.append(f"mass_range=[{min_mass:.2f}, {max_mass:.2f}]")

            # Get most abundant isotope
            most_abundant = self.getMostAbundant()
            parts.append(f"most_abundant_mass={most_abundant.getMZ():.4f}")

        return f"IsotopeDistribution({', '.join(parts)})"
