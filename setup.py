from setuptools import setup, find_packages

setup(
    name = 'QMDAnalysis',

    version = '0.1.0',

    description = '''A python package for analysising mixed quantum-classical
                     trajectory simulations, scattering calcualtions,
                     coordinate transforms, and PES generation.''',


    packages=find_packages("QMDAnalysis"),


    author = 'Kyle Acheson',
    author_email = 'kyle.acheson@warwick.ac.uk',

    long_description = open('README.md').read(),
    long_description_content_type = "text/markdown",

    url='https://github.com/KyleAcheson/QMDAnalysis',

    include_package_data=True,

    classifiers  = [
        'Development Status :: 1',
        'Programming Language :: Python :: 3',
        "License :: OSI Approved :: GNU General Public License",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],

    python_requires=">=3.10",


    install_requires = [

        'numpy ~= 1.26',
        'scipy ~= 1.12',
        'pyscf ~= 2.4.0'

    ],


    keywords = [
        'Quantum Dynamics',
        'Photochemistry',
        'Non-adiabatic Dynamics',
        'X-ray Scattering',
        'Electron Scattering',
        'Ultrafast Electron Diffraction'],

)
