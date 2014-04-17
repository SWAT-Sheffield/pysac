import sys, math, os

class Mock(object):
    __all__ = []
    def __init__(self, *args, **kwargs):
        for key, value in kwargs.iteritems():
            setattr(self, key, value)

    def __call__(self, *args, **kwargs):
        return Mock()

    __add__  = __mul__  = __getitem__ = __setitem__ = \
__delitem__ = __sub__ =  __floordiv__ = __mod__ = __divmod__ = \
__pow__ = __lshift__ = __rshift__ = __and__ = __xor__ = __or__ = \
__rmul__  = __rsub__  = __rfloordiv__ = __rmod__ = __rdivmod__ = \
__rpow__ = __rlshift__ = __rrshift__ = __rand__ = __rxor__ = __ror__ = \
__imul__  = __isub__  = __ifloordiv__ = __imod__ = __idivmod__ = \
__ipow__ = __ilshift__ = __irshift__ = __iand__ = __ixor__ = __ior__ = \
__neg__ = __pos__ = __abs__ = __invert__ = __call__

    def __getattr__(self, name):
        if name in ('__file__', '__path__'):
            return '/dev/null'
        if name == 'sqrt':
            return math.sqrt
        elif name[0] != '_' and name[0] == name[0].upper():
            return type(name, (), {})
        else:
            return Mock(**vars(self))

    def __lt__(self, *args, **kwargs):
        return True

    __nonzero__ = __le__ = __eq__ = __ne__ = __gt__ = __ge__ = __contains__ = \
__lt__


    def __repr__(self):
        # Use _mock_repr to fake the __repr__ call
        res = getattr(self, "_mock_repr")
        return res if isinstance(res, str) else "Mock"

    def __hash__(self):
        return 1

    __len__ = __int__ = __long__ = __index__ = __hash__

    def __oct__(self):
        return '01'

    def __hex__(self):
        return '0x1'

    def __float__(self):
        return 0.1

    def __complex__(self):
        return 1j

MOCK_MODULES = ['tvtk.api', 'h5py']

for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = Mock()


# Load all of the global Astropy configuration
from astropy.sphinx.conf import *

extensions += ['sphinx.ext.autosummary', 'sphinx.ext.mathjax']
numpydoc_show_class_members = False
# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'pySAC'
copyright = u'2013, Stuart Mumford'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '0.0.1a'
# The full version, including alpha/beta/rc tags.
release = '0.0.1a'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'default'
html_favicon = ''
