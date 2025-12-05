


    def __str__(self):
        """
        Return the formula string (e.g., 'C6H12O6').
        This is the expected format for printing in equations/expressions.
        """
        return self.toString()

    def __repr__(self):
        """
        Return a string representation of the EmpiricalFormula object.

        Returns key properties in a readable format:
        EmpiricalFormula(formula='C6H12O6', mono_mass=180.063, charge=0)
        """
        # toString() already returns a Python string via autowrap
        formula_str = self.toString()
        cdef double mono_mass = self.getMonoWeight()
        cdef int charge = self.getCharge()

        parts = []
        if formula_str:
            parts.append(f"formula='{formula_str}'")
        parts.append(f"mono_mass={mono_mass:.4f}")
        if charge != 0:
            parts.append(f"charge={charge}")

        return f"EmpiricalFormula({', '.join(parts)})"
