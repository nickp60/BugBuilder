"""
"""

from setuptools import setup, find_packages
import re
from codecs import open
from os import path
import sys
from pip.req import parse_requirements

here = path.abspath(path.dirname(__file__))

VERSIONFILE = "BugBuilder/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: bugBuilder requires Python 3.5 " +
                     "or above...exiting.\n")
    sys.exit(1)

## parse requirements file
install_reqs = parse_requirements("requirements.txt",
                                  session=False)
requirements = [str(ir.req) for ir in install_reqs]
setup(
    name='BugBuilder',
    version=verstr,
    description='From reads to genome -- made easy!',
    long_description="BugBuilder is a pipeline for assembling bacterial " +
    "genomes with as little stress as possible.",
    url='https://github.com/jamesabbott/BugBuilder',
    author='James Abbott and Nick Waters',
    author_email='j.abbott@imperial.ac.uk and nickp60@gmail.com',
    license='Artistic License 2.0',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Artistic License 2.0',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='bioinformatics assembly genomics development',
    packages=['BugBuilder'],
    install_requires=requirements,
    include_package_data=True,
    package_data={
       '': [path.join(__name__, "BugBuilder", "config_data/*")],
    },
    entry_points={
       'console_scripts': [
           'BugBuilder=BugBuilder.BugBuilder:main',
       ],
    },
    # scripts=[
    #     "scripts/sis/multifasta.py",
    #     "scripts/sis/sis.py"
    # ],
)
