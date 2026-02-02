


    # JSON serialization helpers with consistent naming
    @staticmethod
    def peptidoformToJSON(pf):
        """
        peptidoformToJSON(pf: Peptidoform) -> str

        Convert a Peptidoform to JSON string representation.
        Alias for toJSON(Peptidoform).

        Args:
            pf: The Peptidoform to serialize

        Returns:
            str: JSON string representation
        """
        return ProForma.toJSON(pf)

    @staticmethod
    def peptidoformIonToJSON(pfi):
        """
        peptidoformIonToJSON(pfi: PeptidoformIon) -> str

        Convert a PeptidoformIon to JSON string representation.
        Alias for toJSONIon(PeptidoformIon).

        Args:
            pfi: The PeptidoformIon to serialize

        Returns:
            str: JSON string representation
        """
        return ProForma.toJSONIon(pfi)

    def __repr__(self):
        return "ProForma()"
