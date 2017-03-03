from setuptools import setup, find_packages, Extension
import numpy
import glob

exec(open('_version.py').read())

'''
filterfiles = glob.glob('data/filters/*.csv')
filterfiles.extend(glob.glob('data/filters/*.txt'))
templatefiles = glob.glob('data/templates/*.fits')
atmospherefiles = glob.glob('data/atmospheres/*.fits')
footprintfiles = glob.glob('data/footprints/*.txt')


scripts = ['scripts/combineAtmosphereFiles.py']

setup(
    name='fgcm_y3a1_tools',
    version=__version__,
    description='Tools for FGCM Y3A1',
    author='Eli Rykoff, Dave Burke',
    author_email='erykoff@slac.stanford.edu',
    packages=find_packages(),
    data_files=[('fgcm_y3a1_tools/data/filters', filterfiles),
                ('fgcm_y3a1_tools/data/templates',templatefiles),
                ('fgcm_y3a1_tools/data/atmospheres',atmospherefiles),
                ('fgcm_y3a1_tools/data/footprints',footprintfiles)],
    scripts=scripts
)
'''

