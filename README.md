# SafeBridge

## Overview

**SafeBridge** is a light-weight open-source tool for the Python programming language that estimates bridge damage indicators with low latency by combining remotely sensed and processed data (Multi-Temporal InSAR time series) with algorithms to process topologies in the geospatial domain.

For furher information and underlying research, you can check the paper from [here](https://www.ndt.net/search/docs.php3?id=29731)

## Installation

1. Create a virtual environment and activate it

```bash
# via python note: safebridge requires at least version 3.11
python -m venv safebridge_env
# activate the environment
source safebridge_env/bin/activate # for unix
safebridge_env\Scripts\activate.bat # for windows

# or conda
conda create --name safebridge_env python=3.11
# activate environment
conda activate safebridge_env
```

2. Clone the reposity to your current directory

```bash
git clone https://github.com/SafeStruct/SafeBridge.git
```

3. Install the library from directory using pip

```bash
pip install -e ./SafeBridge
```

4. Testing the installation using tutorial

```bash
# before testing one extension needs to be installed
python -c "import duckdb; duckdb.install_extension('spatial')"
# change the directory to examples and run the tutorial
cd SafeBridge/examples
python tutorial.py
```

If everything works correctly, you should see a folder called **`SafeBridgeDB`** inside the example folder. Inside that folder, you will find the generated **PDF** report for the entire dataset provided.

## Disclaimer

The toy MT-InSAR data presented in the repository and used for the tutorial is a synthetic one and has no relation to real data. It was generated with the intent to present how the `SafeBridge` is working.

## Acknowledgements

This work is supported by Vidi project InStruct, project number 18912, financed by the Dutch Research Council (NWO).
