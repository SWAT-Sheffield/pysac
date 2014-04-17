# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 16:16:03 2013

@author: Stuart Mumford

pySAC: VAC / SAC routines.
"""

import setuptools

DOCLINES = __doc__.split("\n")

CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: GNU3 License',
    'Programming Language :: Python',
    'Topic :: Software Development',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Physics',
    'Operating System :: POSIX',
    'Operating System :: Unix',
]

setuptools.setup(
        author="Stuart Mumford",
        author_email="stuart@mumford.me.uk",
        classifiers=CLASSIFIERS,
        description=DOCLINES[0],
        install_requires=[
            'numpy',
            'scipy',
            'astropy>=0.3.0',
            'matplotlib>=1.2',
        ],
        license="GNU3",
        long_description="\n".join(DOCLINES[2:]),
        maintainer="Stuart Mumford",
        maintainer_email="stuart@mumford.me.uk",
        platforms=["Linux", "Solaris", "Unix"],
        provides=['pysac'],
        name="pysac",
        packages=setuptools.find_packages(),
        version="0.0.1",
    )
