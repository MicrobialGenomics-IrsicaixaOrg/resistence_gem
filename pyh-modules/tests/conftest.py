import sys
import os

# Add the parent directory (pyh-modules) to the PYTHONPATH so that 'sample' can be imported
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

