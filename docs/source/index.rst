.. SafeBridge documentation master file, created by
   sphinx-quickstart on Fri Aug  1 14:32:04 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SafeBridge's documentation!
======================================

.. image:: ./image/safebridge_logo.png
   :width: 140px
   :align: left
   :alt: SafeBridge Logo

**SafeBridge** is an open-source Python package for estimating damage indicators for bridges, utilising Multi-Temporal InSAR (MT-InSAR) 
time-series data and geospatial operations to leverage bridge topologies within a Geographical Information System (GIS) environment. 
**SafeBridge** provides an efficient way to derive the damage indicators for bridges promptly in the era of big data.


This documentation serve as a comprehensive guide to using SafeBridge. 

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Background & Motivation
=======================

The growing use of satellite-based monitoring, particularly time-series data from MT-InSAR, has led to a significant increase in data volume and complexity, 
posing challenges for researchers and engineers in terms of storage, processing, and analysis. Handling such large datasets often requires access to high-performance 
computing (HPC) systems; however, access to these resources can be limited by institutional constraints, cost, or technical barriers. As an alternative, many professionals 
rely on in-house computing resources such as laptops, desktops, or workstations, which demand careful optimisation and the use of efficient, open-source tools to manage processing 
loads. Yet, selecting and effectively using these tools requires expertise across multiple domains, including geospatial analysis, remote sensing, and data science. This is 
where `SafeBridge` offers a practical solution—by providing a streamlined, tailored approach for analysing Multi-Temporal InSAR time-series data, it simplifies the complex 
workflow of bridge damage detection, enabling more users to extract meaningful insights from satellite observations without needing advanced computation infrastructures. 


Acknowledgements
================

This work is supported by Vidi project InStruct, project number 18912, financed by the Dutch Research Council (NWO).

