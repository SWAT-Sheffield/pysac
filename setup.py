# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 16:16:03 2013

@author: Stuart Mumford

pySAC: VAC / SAC routines.
"""
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


def install(setup): #pylint: disable=W0621
    from setuptools import find_packages
    setup(
        author="Stuart Mumford",
        author_email="stuart@mumford.me.uk",
        classifiers=CLASSIFIERS,
        description=DOCLINES[0],
        install_requires=[
            'numpy',
            'h5py',
            'scipy',
            'matplotlib>=1.2',
        ],
        license="GNU3",
        long_description="\n".join(DOCLINES[2:]),
        maintainer="Stuart Mumford",
        maintainer_email="stuart@mumford.me.uk",
        platforms=["Linux", "Solaris", "Unix"],
        provides=['pysac'],
        name="pysac",
        packages=find_packages(),
        version="0.0.1",
    )

if __name__ == '__main__':
    from distribute_setup import use_setuptools
    use_setuptools()
    from setuptools import setup
    install(setup)