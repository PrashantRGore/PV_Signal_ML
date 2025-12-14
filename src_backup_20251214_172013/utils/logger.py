"""
Logging utility for PV_Signal_ML
Provides audit-compliant logging with timestamps
"""
import logging
import sys
from pathlib import Path
from datetime import datetime
import config


def setup_logger(name: str, level=logging.INFO) -> logging.Logger:
    """
    Setup logger with file and console handlers
    Args:
        name: Logger name (typically __name__)
        level: Logging level
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Avoid duplicate handlers
    if logger.handlers:
        return logger
    
    # Create formatters
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (if audit logging enabled)
    if config.ENABLE_AUDIT_LOGGING:
        file_handler = logging.FileHandler(
            config.AUDIT_LOG_FILE,
            encoding='utf-8'
        )
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger
