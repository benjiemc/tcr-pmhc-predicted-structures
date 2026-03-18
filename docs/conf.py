import os
import sys

sys.path.insert(0, os.path.abspath('..'))
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'tcr_pmhc_predicted_structures'
copyright = (
    '2026, Benjamin McMaster, Ali El Moselhy, Ilija Ilievski, Christopher J. Thorpe, Jamie Rossjohn, Charlotte M. Deane'
    ', and Hashem Koohy'
)
author = 'Benjamin McMaster'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'myst_parser', 'nbsphinx', 'sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
