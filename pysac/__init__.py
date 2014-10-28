
__all__ = ['__version__', '__githash__']

#Import version and git hash
try:
    from .version import version as __version__
except ImportError:
    __version__ = ''
try:
    from .version import githash as __githash__
except ImportError:
    __githash__ = ''
