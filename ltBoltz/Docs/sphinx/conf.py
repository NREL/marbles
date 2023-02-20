import sys

extensions = [ 'sphinx.ext.mathjax']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = u'ltBoltz'
copyright = u'2023, the Alliance for Sustainable Energy, LLC., through National Renewable Energy Laboratory. All rights reserved'
author = u'E. Young, M. T. Henry de Frahan, H. Sitaraman, R. Larsen'
version = u'0.1'
release = u'0.1'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
numfig = True
numfig_format = {'figure': '%s', 'table': '%s', 'code-block': '%s'}
html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']
htmlhelp_basename = 'ltBoltzdoc'
latex_elements = { }
latex_documents = [
    (master_doc, 'ltBoltz.tex', u'ltBoltz Documentation',
     author, 'manual'),
]
texinfo_documents = [
    (master_doc, 'ltBoltz', u'ltBoltz Documentation',
     author, 'ltBoltz', 'Lattice Boltzmann solver with adaptive mesh refinement.',
     'Miscellaneous'),
]
