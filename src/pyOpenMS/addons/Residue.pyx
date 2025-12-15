


    def __str__(self):
        """
        __str__(self: Residue) -> str
        
        Return a string representation of the Residue object.
        Returns the residue as a compact string (one letter code with optional modification).
        """
        one_letter = self.getOneLetterCode()
        if one_letter:
            if self.isModified():
                mod_name = self.getModificationName()
                if mod_name:
                    return f"{one_letter}({mod_name})"
            return one_letter
        return self.getName()

    def __repr__(self):
        """
        __repr__(self: Residue) -> str
        
        Return a detailed string representation of the Residue object.

        Returns key properties in a readable format:
        Residue(name='Methionine', one_letter='M', three_letter='Met', formula='C5H11NO2S', mono_mass=149.0510)
        """
        name = self.getName()
        one_letter = self.getOneLetterCode()
        three_letter = self.getThreeLetterCode()
        formula = self.getFormula()
        mono_mass = self.getMonoWeight()
        is_modified = self.isModified()

        parts = []
        if name:
            parts.append(f"name='{name}'")
        if one_letter:
            parts.append(f"one_letter='{one_letter}'")
        if three_letter:
            parts.append(f"three_letter='{three_letter}'")
        if formula:
            formula_str = formula.toString()
            if formula_str:
                parts.append(f"formula='{formula_str}'")
        parts.append(f"mono_mass={mono_mass:.4f}")
        if is_modified:
            mod_name = self.getModificationName()
            if mod_name:
                parts.append(f"modification='{mod_name}'")

        return f"Residue({', '.join(parts)})"
