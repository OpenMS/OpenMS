


    @staticmethod
    def fromString(str proforma_string):
        """
        fromString(proforma_string: str) -> PeptidoformIon

        Parse a ProForma string into a PeptidoformIon object.

        This is a convenience method that wraps ProForma.parseIon().

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
        return ProForma.parseIon(proforma_string)

    def toString(self, mode=None):
        """
        toString(self, mode: WriteMode = None) -> str

        Convert this PeptidoformIon back to ProForma string notation.

        Args:
            mode: Write mode (LOSSLESS or CANONICAL). Defaults to LOSSLESS.

        Returns:
            str: The ProForma string representation
        """
        if mode is None:
            mode = ProForma.WriteMode.LOSSLESS
        # Use the C++ ProForma::toString(PeptidoformIon, mode) via toStringIon binding
        return ProForma.toStringIon(self, mode)

    def getMonoWeight(self):
        """
        getMonoWeight(self) -> float

        Calculate the monoisotopic mass of this peptidoform ion in Daltons.

        Returns:
            float: The monoisotopic mass
        """
        return ProForma.getMonoWeightIon(self)

    def canCalculateMass(self):
        """
        canCalculateMass(self) -> bool

        Check if the monoisotopic mass can be calculated for this peptidoform ion.

        Returns:
            bool: True if mass calculation is possible
        """
        return ProForma.canCalculateMassIon(self)

    def getMassCalculationIssues(self):
        """
        getMassCalculationIssues(self) -> list[ConversionIssue]

        Get a list of issues preventing mass calculation.

        Returns:
            list: List of ConversionIssue objects describing problems
        """
        return ProForma.getMassCalculationIssuesIon(self)

    def getMZ(self):
        """
        getMZ(self) -> float

        Calculate m/z for this peptidoform ion using its charge state.

        Returns:
            float: The m/z value

        Raises:
            RuntimeError: If no charge state is specified
        """
        return ProForma.getMZ(self)

    def toJSON(self):
        """
        toJSON(self) -> str

        Convert this PeptidoformIon to a JSON string representation.

        Returns:
            str: JSON string of the peptidoform ion AST
        """
        return ProForma.peptidoformIonToJSON(self)

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
        return ProForma.peptidoformIonFromJSON(json_str)

    def canGenerateSpectrum(self):
        """
        canGenerateSpectrum(self) -> bool

        Check if a theoretical spectrum can be generated for this peptidoform ion.

        Returns:
            bool: True if spectrum generation is possible

        Example:
            >>> pfi = PeptidoformIon.fromString("PEPTIDE/2")
            >>> pfi.canGenerateSpectrum()
            True
        """
        return ProForma.canGenerateSpectrumIon(self)

    def getSpectrumGenerationIssues(self):
        """
        getSpectrumGenerationIssues(self) -> list[ConversionIssue]

        Get a list of issues preventing spectrum generation.

        Returns:
            list: List of ConversionIssue objects (empty if spectrum can be generated)
        """
        return ProForma.getSpectrumGenerationIssuesIon(self)

    def generateSpectrum(self, int min_charge=1, int max_charge=1, str ion_types="by", bool add_losses=False, bool add_metainfo=True):
        """
        generateSpectrum(self, min_charge: int = 1, max_charge: int = 1, ion_types: str = "by", add_losses: bool = False, add_metainfo: bool = True) -> MSSpectrum

        Generate a theoretical MS/MS spectrum for this peptidoform ion.

        For single-chain peptidoforms, uses TheoreticalSpectrumGenerator.
        For cross-linked peptidoforms (// separator), uses TheoreticalSpectrumGeneratorXLMS.

        Args:
            min_charge: Minimum fragment ion charge state (default: 1)
            max_charge: Maximum fragment ion charge state (default: 1)
            ion_types: Ion types to generate, e.g. "by" (default) or "abyM". Use a/b/c/x/y/z for ion series, M for precursor, I for immonium.
            add_losses: If True, include neutral loss peaks (H2O, NH3) (default: False)
            add_metainfo: If True, include ion annotations in spectrum (default: True)

        Returns:
            MSSpectrum: Theoretical spectrum with fragment peaks

        Raises:
            RuntimeError: If spectrum generation fails

        Example:

            >>> # Single peptide
            >>> pfi = PeptidoformIon.fromString("PEPTIDE/2")
            >>> if pfi.canGenerateSpectrum():
            ...     spec = pfi.generateSpectrum(1, 2, "by")

            >>> # Cross-linked peptide
            >>> pfi = PeptidoformIon.fromString("PEPTK[-2.0#XL1]IDE//ANOTHERK[#XL1]PEPTIDE/3")
            >>> if pfi.canGenerateSpectrum():
            ...     spec = pfi.generateSpectrum(1, 2, "aby")
        """
        return ProForma.generateSpectrumIon(self, min_charge, max_charge, ion_types.encode('utf-8'), add_losses, add_metainfo)

    def __len__(self):
        """
        Return the number of chains in this peptidoform ion.
        """
        return len(self.chains)

    def __str__(self):
        """
        Return the ProForma string representation (lossless mode).
        """
        return self.toString(ProForma.WriteMode.LOSSLESS)

    def __repr__(self):
        """
        Return a detailed string representation for debugging.
        """
        try:
            seq_str = self.toString(ProForma.WriteMode.LOSSLESS)
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
                except Exception:
                    pass

            return f"PeptidoformIon({', '.join(parts)})"
        except Exception:
            return "PeptidoformIon()"
