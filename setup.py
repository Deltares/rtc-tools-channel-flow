"""Simple and hydraulic Modelica routing models for RTC-Tools 2.

RTC-Tools is the Deltares toolbox for control and optimization of water systems.
"""
from setuptools import setup, find_packages

import versioneer

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Information Technology
License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)
Programming Language :: Other
Topic :: Scientific/Engineering :: GIS
Topic :: Scientific/Engineering :: Mathematics
Topic :: Scientific/Engineering :: Physics
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    name="rtc-tools-channel-flow",
    version=versioneer.get_version(),
    author='Matthijs den Toom, Jorn Baayen, et al.',
    maintainer='Jack Vreeken',
    description=DOCLINES[0],
    url='http://www.deltares.nl/en/software/rtc-tools/',
    download_url='http://gitlab.com/deltares/rtc-tools-channel-flow/',
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    platforms=['Windows', 'Linux', 'Mac OS-X', 'Unix'],
    license="LGPL",
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=["rtc-tools"],
    include_package_data=True,
    cmdclass=versioneer.get_cmdclass(),
    entry_points={
        'rtctools.libraries.modelica': [
            'library_folder = rtctools_channel_flow:modelica',
        ]
    }
)
