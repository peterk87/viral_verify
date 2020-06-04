#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=7.0',
                'biopython>=1.76',
                'pandas',
                'attrs']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Peter Kruczkiewicz",
    author_email='peter.kruczkiewicz@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="viralVerify rewrite/refactor for PyPI packaging and distribution",
    entry_points={
        'console_scripts': [
            'viral_verify=viral_verify.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='viral_verify',
    name='viral_verify',
    packages=find_packages(include=['viral_verify', 'viral_verify.*']),
    package_data={'viral_verify': ['data/classifier_table.txt', ]},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/peterk87/viral_verify',
    version='0.1.1',
    zip_safe=False,
)
