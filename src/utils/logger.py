"""
PV Signal ML - Logger Utility
Simple logger for production systems
"""
import logging
from pathlib import Path

def setup_logger(name: str = __name__):
    """Setup standardized logger for PV Signal ML."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    if not logger.handlers:
        # Console handler
        ch = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    
    return logger

logger = setup_logger(__name__)
