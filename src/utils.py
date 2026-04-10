# src/utils.py
"""
Created on Fri Apr 10 14:22 2026

@author: Pioli Valeria
"""

def load_parameters(filepath):
    """
    Reads a text file with 'key = value' format and returns a dictionary.

    Parameters:
    -----------
    filepath : str
        Path to the .txt file.

    Returns:
    --------
    params : dict
        Dictionary containing parameter names as keys and floats as values.
    """
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            # Remove comments (everything after #) and whitespace
            clean_line = line.split('#')[0].strip()
            if not clean_line or '=' not in clean_line:
                continue

            # Split key and value
            key, value = clean_line.split('=')
            params[key.strip()] = float(value.strip())

    return params
