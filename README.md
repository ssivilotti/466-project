# Algae RNA Secondary Structure

Nussinov Algorithm Implementation for calculating the secondary structure of 23S Algae rRNA. 

## Installation

First clone the repository.

To run the algorithm, first install the required dependencies in your python environment by running  

```bash
pip install numpy matplotlib forgi
```

## To Run

The algorithm to generate secondary structures from a fasta file is implemented in nussinov.py. To generate structures for all the Algae rRNA sequences run the nussinov.py file in the terminal.

```bash
python nussinov.py
```

To generate structures for specific sequences, for example L42854, the below snippet of code can be used.

```python
from nussinov import compute_for_file

compute_for_file('microgreen_id_rna', 'data', sequences_to_read=['L42854'])
```



# RNA Secondary Structure Visualization Tool

## Overview
This tool provides an efficient method for visualizing RNA secondary structures. It is implemented in Python and utilizes the `forgi` library for generating graphical representations of RNA structures. The visualizations are displayed in an interactive web application built with Dash.

## Features
- Parsing RNA sequences and secondary structures from FASTA files.
- Generating visual plots of RNA secondary structures.
- Dynamic web interface for viewing and interacting with the generated plots.

## Prerequisites
- Python 3.x
- Jupyter Notebook or JupyterLab environment

## Installation
1. **Install Required Libraries**

   Create and activate a virtual environment (optional but recommended):

   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
   pip install -r requirements.txt
   ```

## Usage

1. **Open the Jupyter Notebook**

   Run Jupyter Notebook or JupyterLab:

   ```bash
   jupyter notebook
   ```

   or

   ```bash
   jupyter lab
   ```

2. **Navigate to the Notebook**

   In the Jupyter interface, open the `.ipynb` file containing the project.

3. **Run the Notebook**

   Execute the cells in the notebook to generate RNA structure visualizations and to start the Dash app.

4. **Viewing the Visualizations**

   Upon running the Dash app cell, a link will be provided in the output. Click on this link to view the RNA structure visualizations in your web browser.

