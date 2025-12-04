


    def __str__(self):
        """
        Return a string representation of the EmpiricalFormula object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the EmpiricalFormula object.

        Returns key properties in a readable format:
        EmpiricalFormula(formula='C6H12O6', mono_mass=180.063, charge=0)
        """
        cdef String formula = self.toString()
        cdef double mono_mass = self.getMonoWeight()
        cdef int charge = self.getCharge()

        formula_str = formula.decode('utf-8') if formula.size() > 0 else ""

        parts = []
        if formula_str:
            parts.append(f"formula='{formula_str}'")
        parts.append(f"mono_mass={mono_mass:.4f}")
        if charge != 0:
            parts.append(f"charge={charge}")

        return f"EmpiricalFormula({', '.join(parts)})"
