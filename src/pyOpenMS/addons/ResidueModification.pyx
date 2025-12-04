


    def __str__(self):
        """
        Return a string representation of the ResidueModification object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the ResidueModification object.

        Returns key properties in a readable format:
        ResidueModification(id='Oxidation', name='Oxidation', mass_diff=15.9949, origin='M')
        """
        cdef String mod_id = self.getId()
        cdef String name = self.getName()
        cdef double diff_mono = self.getDiffMonoMass()
        cdef char origin = self.getOrigin()

        mod_id_str = mod_id.decode('utf-8') if mod_id.size() > 0 else ""
        name_str = name.decode('utf-8') if name.size() > 0 else ""

        parts = []
        if mod_id_str:
            parts.append(f"id='{mod_id_str}'")
        if name_str and name_str != mod_id_str:
            parts.append(f"name='{name_str}'")
        if diff_mono != 0.0:
            parts.append(f"mass_diff={diff_mono:.4f}")
        if origin != 0 and origin != ord('X'):
            parts.append(f"origin='{chr(origin)}'")

        return f"ResidueModification({', '.join(parts)})"
