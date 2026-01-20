


    @staticmethod
    def fromString(str proforma_string):
        """
        fromString(proforma_string: str) -> PeptidoformIon

        Parse a ProForma string into a PeptidoformIon object.

        This is a convenience method that wraps ProFormaParser.parseIon().

        Args:
            proforma_string: A ProForma notation string with optional charge
                            (e.g., "PEPTIDE/2", "PEPTIDE//SEQUENCE/3")

        Returns:
            PeptidoformIon: The parsed peptidoform ion AST

        Raises:
            RuntimeError: If the ProForma string is invalid

        Example:
            >>> pfi = PeptidoformIon.fromString("PEPTIDE/2")
            >>> str(pfi)
            'PEPTIDE/2'
        """
        return ProFormaParser.parseIon(proforma_string)

    def toString(self, mode=None):
        """
        toString(self, mode: ProFormaWriteMode = None) -> str

        Convert this PeptidoformIon back to ProForma string notation.

        Args:
            mode: Write mode (LOSSLESS or CANONICAL). Defaults to LOSSLESS.

        Returns:
            str: The ProForma string representation
        """
        if mode is None:
            mode = ProFormaWriteMode.LOSSLESS
        # PeptidoformIon needs to convert first chain for toString
        # Note: This assumes single chain; multi-chain needs different handling
        if len(self.chains) == 1:
            return ProFormaParser.toString(self.chains[0], mode)
        else:
            # For multi-chain, join with // separator
            return "//".join(ProFormaParser.toString(c, mode) for c in self.chains)

    def getMonoWeight(self):
        """
        getMonoWeight(self) -> float

        Calculate the monoisotopic mass of this peptidoform ion in Daltons.

        Returns:
            float: The monoisotopic mass
        """
        return ProFormaParser.getMonoWeightIon(self)

    def canCalculateMass(self):
        """
        canCalculateMass(self) -> bool

        Check if the monoisotopic mass can be calculated for this peptidoform ion.

        Returns:
            bool: True if mass calculation is possible
        """
        return ProFormaParser.canCalculateMassIon(self)

    def getMassCalculationIssues(self):
        """
        getMassCalculationIssues(self) -> list[ConversionIssue]

        Get a list of issues preventing mass calculation.

        Returns:
            list: List of ConversionIssue objects describing problems
        """
        return ProFormaParser.getMassCalculationIssuesIon(self)

    def getMZ(self):
        """
        getMZ(self) -> float

        Calculate m/z for this peptidoform ion using its charge state.

        Returns:
            float: The m/z value

        Raises:
            RuntimeError: If no charge state is specified
        """
        return ProFormaParser.getMZ(self)

    def toJSON(self):
        """
        toJSON(self) -> str

        Convert this PeptidoformIon to a JSON string representation.

        Returns:
            str: JSON string of the peptidoform ion AST
        """
        return ProFormaParser.peptidoformIonToJSON(self)

    @staticmethod
    def fromJSON(str json_str):
        """
        fromJSON(json_str: str) -> PeptidoformIon

        Construct a PeptidoformIon from a JSON string.

        Args:
            json_str: JSON representation of a PeptidoformIon

        Returns:
            PeptidoformIon: The deserialized peptidoform ion
        """
        return ProFormaParser.peptidoformIonFromJSON(json_str)

    def __len__(self):
        """
        Return the number of chains in this peptidoform ion.
        """
        return len(self.chains)

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
            parts.append(f"chains={len(self.chains)}")

            if self.canCalculateMass():
                mass = self.getMonoWeight()
                parts.append(f"mono_mass={mass:.4f}")
                try:
                    mz = self.getMZ()
                    parts.append(f"mz={mz:.4f}")
                except:
                    pass

            return f"PeptidoformIon({', '.join(parts)})"
        except:
            return "PeptidoformIon()"
