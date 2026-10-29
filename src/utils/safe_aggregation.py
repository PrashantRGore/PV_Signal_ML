"""
Safe Aggregation Utilities - Handle NaN, missing values, and edge cases
"""
import pandas as pd
import numpy as np


def safe_mean(series):
    """Calculate mean, handling NaN and empty series"""
    if len(series) == 0 or series.isna().all():
        return 0.0
    return float(series.mean())


def safe_median(series):
    """Calculate median, handling NaN and empty series"""
    if len(series) == 0 or series.isna().all():
        return 0.0
    return float(series.median())


def safe_std(series):
    """Calculate std, handling NaN and empty series"""
    if len(series) <= 1 or series.isna().all():
        return 0.0
    return float(series.std())


def safe_mode(series, default='Unknown'):
    """Get mode, handling NaN and empty series"""
    if len(series) == 0:
        return default
    
    clean = series.dropna()
    if len(clean) == 0:
        return default
    
    mode_vals = clean.mode()
    if len(mode_vals) == 0:
        return default
    
    return mode_vals[0]


def safe_percentage(series, condition):
    """Calculate percentage matching condition, handling NaN"""
    if len(series) == 0:
        return 0.0
    
    clean = series.dropna()
    if len(clean) == 0:
        return 0.0
    
    return float((clean == condition).mean())


def safe_threshold_percentage(series, threshold, operator='>'):
    """Calculate percentage above/below threshold"""
    if len(series) == 0:
        return 0.0
    
    clean = series.dropna()
    if len(clean) == 0:
        return 0.0
    
    if operator == '>':
        result = (clean > threshold).mean()
    elif operator == '>=':
        result = (clean >= threshold).mean()
    elif operator == '<':
        result = (clean < threshold).mean()
    elif operator == '<=':
        result = (clean <= threshold).mean()
    else:
        result = 0.0
    
    return float(result)


def safe_map_mean(series, mapping, default=0.0):
    """Map values and calculate mean"""
    if len(series) == 0:
        return default
    
    clean = series.dropna()
    if len(clean) == 0:
        return default
    
    mapped = clean.map(mapping)
    if mapped.isna().all():
        return default
    
    return float(mapped.mean())
