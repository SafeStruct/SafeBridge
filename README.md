# SafeBridge

## Overview

**SafeBridge** is a light-weight open-source tool for the Python programming language that calculates bridge damage indicators with low latency by combining remotely sensed and processed data (Persistent Scatterer InSAR time series) with algorithms to process topologies in the geospatial domain. 

## Installation

1. Create a virtual environment
```python
# via python
python -m venv safebridge
# or conda
conda create --name safebridge
```
2. Clone the reposity to your current directory
```bash
git clone https://github.com/shnmrt/SafeBridge.git
```
3. Install the library from direcoty using pip
```bash
pip install ./SafeBridge
```
4. Testing the installation using tutorial
```bash
cd SafeBridge/examples
python tutorial.py
```

## Acknowledgements & Funding

This work is supported by Vidi project InStruct, project number 18912, financed by the Dutch Research Council (NWO).
