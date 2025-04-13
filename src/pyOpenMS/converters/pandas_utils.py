def to_dataframe(data_dict):
    """Convert OpenMS data dictionaries to pandas DataFrames (optional). 
    
    Args:
        data_dict (dict): Dictionary from OpenMS objects (MSSpectrum.get_data_dict())
    
    Returns:
        pandas.DataFrame: DataFrame representation
        
    Raises:
        ImportError: If pandas is not installed
    """
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "The 'pandas' package is required for DataFrame conversion.\n"
            "Install with: pip install pandas"
        ) from e
        
    return pd.DataFrame(data_dict)