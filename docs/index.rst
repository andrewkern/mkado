MKado Documentation
===================

**MKado** is a modern Python implementation of the McDonald-Kreitman test toolkit for detecting selection in molecular evolution.

.. note::

   MKado requires Python 3.12 or later.

Features
--------

- **Standard MK test**: Classic 2x2 contingency table with Fisher's exact test
- **Polarized MK test**: Uses a third outgroup to assign mutations to lineages
- **Asymptotic MK test**: Frequency-bin alpha estimates with exponential extrapolation (Messer & Petrov 2013)
- **Tarone-Greenland α_TG**: Weighted multi-gene estimator (Stoletzki & Eyre-Walker 2011)
- **Batch processing**: Process multiple genes with parallel execution
- **Volcano plots**: Visualize batch results with publication-ready volcano plots
- **Multiple output formats**: Pretty-print, TSV, and JSON

Quick Example
-------------

.. code-block:: bash

   # Standard MK test
   mkado test alignment.fa -i "dmel" -o "dsim"

   # Asymptotic MK test
   mkado test alignment.fa -i "dmel" -o "dsim" -a

   # Batch process a directory
   mkado batch alignments/ -i "dmel" -o "dsim"

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   batch-workflow
   alpha-tg
   file-formats
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
