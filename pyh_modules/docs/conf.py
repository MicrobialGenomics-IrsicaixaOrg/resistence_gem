# Configuration file for the Sphinx documentation builder.

import os
import sys
# Updated sys.path to point to the renamed package directory "pyh_modules"
sys.path.insert(0, os.path.abspath('../pyh_modules'))

project = 'resistance'
copyright = '2025, Oriol Careta'
author = 'Oriol Careta'

release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'alabaster'
html_static_path = ['_static']

