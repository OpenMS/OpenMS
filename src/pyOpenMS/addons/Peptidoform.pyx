


    @staticmethod
    def fromString(str proforma_string):
        """
        fromString(proforma_string: str) -> Peptidoform

        Parse a ProForma string into a Peptidoform object.

        This is a convenience method that wraps ProFormaParser.parse().

        Args:
            proforma_string: A ProForma notation string (e.g., "EM[UNIMOD:35]K")

        Returns:
            Peptidoform: The parsed peptidoform AST

        Raises:
            RuntimeError: If the ProForma string is invalid

        Example:
            >>> pf = Peptidoform.fromString("EM[UNIMOD:35]K")
            >>> str(pf)
            'EM[UNIMOD:35]K'
        """
        return ProFormaParser.parse(proforma_string)

    def toString(self, mode=None):
        """
        toString(self, mode: ProFormaWriteMode = None) -> str

        Convert this Peptidoform back to ProForma string notation.

        Args:
            mode: Write mode (LOSSLESS or CANONICAL). Defaults to LOSSLESS.

        Returns:
            str: The ProForma string representation

        Example:
            >>> pf = Peptidoform.fromString("EM[+15.99]K")
            >>> pf.toString()  # LOSSLESS - preserves original
            'EM[+15.99]K'
            >>> pf.toString(ProFormaWriteMode.CANONICAL)  # normalized
            'EM[+15.9900]K'
        """
        if mode is None:
            mode = ProFormaWriteMode.LOSSLESS
        return ProFormaParser.toString(self, mode)

    def getMonoWeight(self):
        """
        getMonoWeight(self) -> float

        Calculate the monoisotopic mass of this peptidoform in Daltons.

        Returns:
            float: The monoisotopic mass

        Raises:
            RuntimeError: If mass cannot be calculated (e.g., unresolved modifications)

        Example:
            >>> pf = Peptidoform.fromString("PEPTIDE")
            >>> pf.getMonoWeight()
            799.36...
        """
        return ProFormaParser.getMonoWeight(self)

    def canCalculateMass(self):
        """
        canCalculateMass(self) -> bool

        Check if the monoisotopic mass can be calculated for this peptidoform.

        Returns:
            bool: True if mass calculation is possible

        Example:
            >>> pf = Peptidoform.fromString("PEPTIDE")
            >>> pf.canCalculateMass()
            True
        """
        return ProFormaParser.canCalculateMass(self)

    def getMassCalculationIssues(self):
        """
        getMassCalculationIssues(self) -> list[ConversionIssue]

        Get a list of issues preventing mass calculation.

        Returns:
            list: List of ConversionIssue objects describing problems
        """
        return ProFormaParser.getMassCalculationIssues(self)

    def getMZ(self, int charge):
        """
        getMZ(self, charge: int) -> float

        Calculate m/z for this peptidoform at a given charge state.

        Args:
            charge: The charge state (e.g., 2 for [M+2H]2+)

        Returns:
            float: The m/z value
        """
        return ProFormaParser.getMZCharge(self, charge)

    def toAASequence(self, policy=None):
        """
        toAASequence(self, policy: AASequenceConversionPolicy = None) -> AASequence

        Convert this Peptidoform to an OpenMS AASequence.

        Args:
            policy: Conversion policy (STRICT, DROP_UNLOCALISED, or BEST_EFFORT).
                    Defaults to STRICT.

        Returns:
            AASequence: The converted amino acid sequence

        Raises:
            RuntimeError: If conversion fails (depends on policy)

        Example:
            >>> pf = Peptidoform.fromString("EM[UNIMOD:35]K")
            >>> seq = pf.toAASequence()
            >>> str(seq)
            'EM(Oxidation)K'
        """
        if policy is None:
            policy = AASequenceConversionPolicy.STRICT
        return ProFormaParser.toAASequence(self, policy)

    def isRepresentableAsAASequence(self):
        """
        isRepresentableAsAASequence(self) -> bool

        Check if this Peptidoform can be fully represented as an AASequence.

        Returns:
            bool: True if conversion to AASequence is lossless
        """
        return ProFormaParser.isRepresentableAsAASequence(self)

    def getAASequenceConversionIssues(self):
        """
        getAASequenceConversionIssues(self) -> list[ConversionIssue]

        Get a list of issues that would arise during AASequence conversion.

        Returns:
            list: List of ConversionIssue objects
        """
        return ProFormaParser.getAASequenceConversionIssues(self)

    def resolveModifications(self):
        """
        resolveModifications(self) -> None

        Resolve all modifications in this Peptidoform using ModificationsDB.

        Modifies this object in place.
        """
        ProFormaParser.resolveModifications(self)

    def toJSON(self):
        """
        toJSON(self) -> str

        Convert this Peptidoform to a JSON string representation.

        Returns:
            str: JSON string of the peptidoform AST
        """
        return ProFormaParser.peptidoformToJSON(self)

    @staticmethod
    def fromJSON(str json_str):
        """
        fromJSON(json_str: str) -> Peptidoform

        Construct a Peptidoform from a JSON string.

        Args:
            json_str: JSON representation of a Peptidoform

        Returns:
            Peptidoform: The deserialized peptidoform
        """
        return ProFormaParser.peptidoformFromJSON(json_str)

    @staticmethod
    def fromAASequence(seq):
        """
        fromAASequence(seq: AASequence) -> Peptidoform

        Create a Peptidoform from an OpenMS AASequence.

        Args:
            seq: An AASequence object

        Returns:
            Peptidoform: The converted peptidoform

        Example:
            >>> seq = AASequence.fromString("EM(Oxidation)K")
            >>> pf = Peptidoform.fromAASequence(seq)
            >>> str(pf)
            'EM[UNIMOD:35]K'
        """
        return ProFormaParser.fromAASequence(seq)

    def __str__(self):
        """
        Return the ProForma string representation (lossless mode).
        """
        return self.toString(ProFormaWriteMode.LOSSLESS)

    def __repr__(self):
        """
        Return a detailed string representation for debugging.
        """
        try:
            seq_str = self.toString(ProFormaWriteMode.LOSSLESS)
            if len(seq_str) > 50:
                seq_str = seq_str[:47] + "..."

            parts = [f"sequence='{seq_str}'"]

            if self.canCalculateMass():
                mass = self.getMonoWeight()
                parts.append(f"mono_mass={mass:.4f}")

            if len(self.unlocalised_mods) > 0:
                parts.append(f"unlocalised_mods={len(self.unlocalised_mods)}")
            if len(self.labile_mods) > 0:
                parts.append(f"labile_mods={len(self.labile_mods)}")
            if len(self.n_term_mods) > 0:
                parts.append("n_term_mod=True")
            if len(self.c_term_mods) > 0:
                parts.append("c_term_mod=True")

            return f"Peptidoform({', '.join(parts)})"
        except:
            return "Peptidoform()"
