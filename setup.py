#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
doclink = """
Documentation
-------------

The full documentation is at http://Absolute_Integrator.rtfd.org."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='Absolute_Integrator',
    version='0.1.0',
    description='STEM image analysis and quantification"',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Lewys Jones',
    author_email='lewysbjones@gmail.com',
    url='https://github.com/LewysJones/Absolute_Integrator',
    packages=[
        'Absolute_Integrator',
    ],
    package_dir={'Absolute_Integrator': 'Absolute_Integrator'},
    include_package_data=True,
    install_requires=["numpy",
    ],
    license='MIT',
    zip_safe=False,
    keywords='Absolute_Integrator',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)
