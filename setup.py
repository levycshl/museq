import setuptools
from distutils.core import setup, Extension
import distutils.sysconfig

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="museq",
    version="1.0.0",
    author="Dan Levy",
    author_email="levy@cshl.edu",
    description="Code for the analysis of MuSeq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/levycshl/museq/",
    packages=[],
    ext_modules=[Extension(
        "AlignmentExt",        
        sources=["alignment.cpp"],
        libraries=['npymath'],
        library_dirs = ['/usr/local/lib/python2.7/dist-packages/numpy/core/lib'],
        extra_compile_args=["-O3", "-std=c++11"]
    )],
    scripts=['museq_from_config_spurless.py',
             'museq_functions.py',
             'support_functions.py'],
    data_files=['museq.conf'],
    requires=['sys', 'os', 'time', 'numpy', 'itertools', 'operator',
              'gzip', 'shutil', 'collections', 'matplotlib'],
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=2.7, <3',
)
