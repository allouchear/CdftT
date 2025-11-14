import os
import sys

sys.path.insert(0, os.path.abspath("../source"))  # adjust as needed

project = "CdftT"
author = "Abdulrahman Allouche"
release = "0.1.0"

# File extensions that are regarded as sources.
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'restructuredtext',
    '.md': 'markdown',
}

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "sphinx_rtd_theme"

# Ensure autodoc works on Read the Docs
autodoc_mock_imports = [
    # add modules you cannot import at build time, e.g.
    # "numpy", "pandas"
]

# Options for html theme
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
}
