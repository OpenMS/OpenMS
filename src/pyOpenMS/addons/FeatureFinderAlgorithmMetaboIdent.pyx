# two empty lines are important !


    @staticmethod
    def compounds_from_df(df):
        """
        compounds_from_df(df: pandas.DataFrame) -> List[FeatureFinderMetaboIdentCompound]

        Convert a pandas DataFrame to a list of FeatureFinderMetaboIdentCompound objects.

        This static method provides a convenient way to create compound libraries from
        tabular data. Column names are case-insensitive and support snake_case variants.

        Args:
            df: DataFrame with the following columns:
                - CompoundName (or compound_name): str, unique identifier (required)
                - SumFormula (or sum_formula): str, chemical formula (required, can be empty)
                - Mass (or mass): float, neutral mass; 0 = calculate from formula (required)
                - Charge (or charge): int, list[int], or comma-separated str (required)
                - RetentionTime (or retention_time): float, list[float], or comma-separated str (required)
                - RetentionTimeRange (or retention_time_range): float, list[float], or str; 0 = default window (required)
                - IsoDistribution (or iso_distribution): float, list[float], or str; 0 = calculate from formula (required)

        Returns:
            List of FeatureFinderMetaboIdentCompound objects

        Raises:
            ValueError: If required columns are missing, compound names are empty or duplicated,
                       or if charges/retention times are empty for any compound.

        Example:
            >>> import pandas as pd
            >>> df = pd.DataFrame({
            ...     'CompoundName': ['glucose', 'fructose'],
            ...     'SumFormula': ['C6H12O6', 'C6H12O6'],
            ...     'Mass': [0.0, 0.0],
            ...     'Charge': [-1, '-1,1'],
            ...     'RetentionTime': [123.4, '200.0,250.0'],
            ...     'RetentionTimeRange': [0.0, 0.0],
            ...     'IsoDistribution': [0.0, 0.0]
            ... })
            >>> compounds = FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        """
        try:
            import pandas as pd
            import numpy as np
        except ImportError as e:
            raise ImportError(
                "pandas and numpy are required for compounds_from_df. "
                "Install with: pip install pandas numpy"
            ) from e

        # Column name mapping: normalized lowercase -> canonical name
        column_map = {
            'compoundname': 'CompoundName',
            'compound_name': 'CompoundName',
            'sumformula': 'SumFormula',
            'sum_formula': 'SumFormula',
            'formula': 'SumFormula',
            'mass': 'Mass',
            'charge': 'Charge',
            'charges': 'Charge',
            'retentiontime': 'RetentionTime',
            'retention_time': 'RetentionTime',
            'rt': 'RetentionTime',
            'retentiontimerange': 'RetentionTimeRange',
            'retention_time_range': 'RetentionTimeRange',
            'rt_range': 'RetentionTimeRange',
            'rtrange': 'RetentionTimeRange',
            'isodistribution': 'IsoDistribution',
            'iso_distribution': 'IsoDistribution',
            'isotope_distribution': 'IsoDistribution',
        }

        # Build mapping from DataFrame columns to canonical names
        df_column_to_canonical = {}
        for col in df.columns:
            normalized = col.lower().strip()
            if normalized in column_map:
                df_column_to_canonical[col] = column_map[normalized]
            elif normalized in [v.lower() for v in column_map.values()]:
                # Direct match to canonical name (case-insensitive)
                for canonical in set(column_map.values()):
                    if normalized == canonical.lower():
                        df_column_to_canonical[col] = canonical
                        break

        # Rename columns to canonical names
        df_renamed = df.rename(columns=df_column_to_canonical)

        # Check for required columns
        required = ['CompoundName', 'SumFormula', 'Mass', 'Charge',
                    'RetentionTime', 'RetentionTimeRange', 'IsoDistribution']
        missing = [col for col in required if col not in df_renamed.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        def _to_list(value, dtype):
            """Convert scalar/list/string to list of dtype."""
            if value is None:
                return []
            # Handle numpy NaN
            if isinstance(value, float) and np.isnan(value):
                return []
            if isinstance(value, str):
                stripped = value.strip()
                if not stripped:
                    return []
                # Handle "0" or "0.0" as empty for optional fields
                return [dtype(x.strip()) for x in stripped.split(',') if x.strip()]
            if isinstance(value, (list, tuple, np.ndarray)):
                return [dtype(x) for x in value]
            # Single scalar value
            return [dtype(value)]

        compounds = []
        seen_names = set()

        for idx, row in df_renamed.iterrows():
            # Get compound name
            name = row['CompoundName']
            if pd.isna(name) or (isinstance(name, str) and not name.strip()):
                raise ValueError(f"Row {idx}: CompoundName cannot be empty")
            name = str(name).strip()

            # Check for duplicate names
            if name in seen_names:
                raise ValueError(f"Row {idx}: Duplicate CompoundName '{name}'")
            seen_names.add(name)

            # Get formula (can be empty string)
            formula = row['SumFormula']
            if pd.isna(formula):
                formula = ''
            else:
                formula = str(formula).strip()

            # Get mass
            mass = row['Mass']
            if pd.isna(mass):
                mass = 0.0
            else:
                mass = float(mass)

            # Get charges (must have at least one)
            charges = _to_list(row['Charge'], int)
            if not charges:
                raise ValueError(f"Row {idx} ('{name}'): Charge cannot be empty")

            # Get retention times (must have at least one)
            rts = _to_list(row['RetentionTime'], float)
            if not rts:
                raise ValueError(f"Row {idx} ('{name}'): RetentionTime cannot be empty")

            # Get retention time ranges (can be empty or zero for default window)
            rt_ranges = _to_list(row['RetentionTimeRange'], float)
            if not rt_ranges:
                rt_ranges = [0.0]

            # Get isotope distribution (can be empty or zero for automatic calculation)
            iso_dist = _to_list(row['IsoDistribution'], float)
            if not iso_dist:
                iso_dist = [0.0]

            compound = FeatureFinderMetaboIdentCompound(
                name, formula, mass, charges, rts, rt_ranges, iso_dist
            )
            compounds.append(compound)

        return compounds

    def run_from_df(self, df, features, spectra_path):
        """
        run_from_df(self: FeatureFinderAlgorithmMetaboIdent, df: pandas.DataFrame, features: FeatureMap, spectra_path: str) -> None

        Run feature extraction using compounds from a pandas DataFrame.

        This convenience method combines compounds_from_df() and run() into a single call.
        The spectra must be set beforehand using setMSData().

        Args:
            df: DataFrame with compound information (see compounds_from_df for format)
            features: FeatureMap to store detected features
            spectra_path: Path to the spectra file (annotated in the feature map)

        Example:
            >>> import pandas as pd
            >>> exp = MSExperiment()
            >>> MzMLFile().load("spectra.mzML", exp)
            >>> ff = FeatureFinderAlgorithmMetaboIdent()
            >>> ff.setMSData(exp)
            >>> df = pd.DataFrame({
            ...     'CompoundName': ['glucose'],
            ...     'SumFormula': ['C6H12O6'],
            ...     'Mass': [0.0],
            ...     'Charge': [-1],
            ...     'RetentionTime': [123.4],
            ...     'RetentionTimeRange': [0.0],
            ...     'IsoDistribution': [0.0]
            ... })
            >>> features = FeatureMap()
            >>> ff.run_from_df(df, features, "spectra.mzML")
        """
        compounds = FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.run(compounds, features, spectra_path)
