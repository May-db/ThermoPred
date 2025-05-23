# src/ThermoPred/__init__.py

"""
ThermoPred - A tool to predict thermodynamic feasibility of chemical reactions
"""

# You can optionally expose key modules or functions for easier access
from . import reaction_utils

# Define what gets imported with "from ThermoPred import *"
__all__ = ['reaction_utils']

# Package metadata
__version__ = '0.1.0'
__author__ = 'Maria, Solene, May'