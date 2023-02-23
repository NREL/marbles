import sys

extensions = [ 'sphinx.ext.mathjax']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = u'MARBLES'
copyright = u'2023, the Alliance for Sustainable Energy, LLC., through National Renewable Energy Laboratory. All rights reserved'
author = u'E. Young, M. T. Henry de Frahan, H. Sitaraman, R. Larsen, J. Rood'
version = u'0.1'
release = u'0.1'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
numfig = True
numfig_format = {'figure': '%s', 'table': '%s', 'code-block': '%s'}
html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']
htmlhelp_basename = 'marblesdoc'
latex_elements = { }
latex_documents = [
    (master_doc, 'marbles.tex', u'MARBLES Documentation',
     author, 'manual'),
]
texinfo_documents = [
    (master_doc, 'marbles', u'MARBLES Documentation',
     author, 'MARBLES', 'Lattice Boltzmann solver with adaptive mesh refinement.',
     'Miscellaneous'),
]
