


    def __str__(self):
        """
        __str__(self: ResidueModification) -> str

        Return a string representation of the ResidueModification object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: ResidueModification) -> str

        Return a string representation of the ResidueModification object.

        Returns key properties in a readable format:
        ResidueModification(id='Oxidation', name='Oxidation', mass_diff=15.9949, origin='M')
        """
        # autowrap methods return Python types, not C++ types
        mod_id_str = self.getId()
        name_str = self.getName()
        diff_mono = self.getDiffMonoMass()
        origin = self.getOrigin()

        parts = []
        if mod_id_str:
            parts.append(f"id='{mod_id_str}'")
        if name_str and name_str != mod_id_str:
            parts.append(f"name='{name_str}'")
        if diff_mono != 0.0:
            parts.append(f"mass_diff={diff_mono:.4f}")
        # origin is a string from autowrap, check if not empty and not 'X'
        if origin and origin != 'X':
            parts.append(f"origin='{origin}'")

        return f"ResidueModification({', '.join(parts)})"
