============
viral_verify
============


.. image:: https://img.shields.io/pypi/v/viral_verify.svg
        :target: https://pypi.python.org/pypi/viral_verify

.. image:: https://img.shields.io/travis/peterk87/viral_verify.svg
        :target: https://travis-ci.com/peterk87/viral_verify

.. image:: https://readthedocs.org/projects/viral-verify/badge/?version=latest
        :target: https://viral-verify.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


viralVerify_ rewrite/refactor for PyPI packaging and distribution, maintainability and clarity.

**NOTE:** BLAST+ search option has been removed. Results output table will be different than the original viralVerify_. Naive Bayes classifier training script has not been ported yet.


* Free software: MIT license
* Documentation: https://viral-verify.readthedocs.io.


Features
--------

* Gene prediction with Prodigal_ in metagenomic mode
* HMMer3_ ``hmmsearch`` for protein domains in predicted genes
* Naive Bayes classification of contigs as viral/not viral based on HMMer3_ results
* Output of detailed contig classification results table in CSV format
* Output of contigs based on classification into separate FASTA files

Requirements
------------

An HMMer3_ HMM database is required. For example, the latest version of Pfam-A HMM:

* http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz

**NOTE:** Please extract any compressed HMM DB (``$ gunzip Pfam-A.hmm.gz``)

Software dependencies:

* Prodigal_ v2.6.3
* HMMer3_ v3.3
* Python v3.6+

Python dependencies:

* attrs_
* Click_
* Pandas_
* Biopython_

Installation
------------

Conda
~~~~~

It's recommended that you use Conda_ to install the required software (Prodigal_ and HMMer3_) and Python dependencies.

.. code-block::

    $ conda env create -f environment.yml

Pip
~~~

If you have Prodigal_ and HMMer3_ installed in your ``$PATH``, and Python 3.6 or greater, you can use ``pip`` to install ``viral_verify``:

.. code-block::

    $ pip install viral_verify


Usage
-----

.. code-block::

    $ viral_verify --help
    Usage: viral_verify [OPTIONS]

      HMM and Naive Bayes classification of contig sequences as either viral,
      plasmid or chromosomal.

      Requires Prodigal for gene prediction and hmmsearch from HMMer3 for
      searching for Pfam HMM profiles.

    Options:
      -i, --input-fasta PATH          Input fasta file  [required]
      -o, --outdir PATH               Output directory  [required]
      -H, --hmm-db PATH               Path to Pfam-A HMM database  [required]
      -t, --threads INTEGER           Number of threads (default=16)
      -p, --output-plasmids-separately
                                      Output predicted plasmids separately?
      --prefix TEXT                   Output file prefix (default: None)
      --uncertainty-threshold FLOAT   Uncertainty threshold (Natural log
                                      probability) (default=3.0)

      --naive-bayes-classifier-table PATH
                                      Table of protein domain frequencies to use
                                      for Naive Bayes classification (default="/ho
                                      me/pkruczkiewicz/repos/viral_verify/viral_ve
                                      rify/data/classifier_table.txt")

      -v, --verbose                   Logging verbosity
      --version                       Show the version and exit.
      --help                          Show this message and exit.



Credits
-------

The original source code, design and conception can be found at viralVerify_. This is merely a rewrite for easier packaging via PyPI, adding some CI with Travis-CI and organizing the code for maintainability and clarity.

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _viralVerify: https://github.com/ablab/viralVerify/
.. _attrs: https://www.attrs.org/en/stable/
.. _Click: https://click.palletsprojects.com/en/7.x/
.. _Pandas: https://pandas.pydata.org/
.. _Biopython: https://github.com/biopython/biopython
.. _Conda: https://docs.conda.io/en/latest/
.. _HMMer3: http://hmmer.org/
.. _Prodigal: https://github.com/hyattpd/Prodigal
