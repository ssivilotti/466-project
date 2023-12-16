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

## Analyzing Structures

The analysis of the secondary structures can be found in analyze.ipynb and creates 2D visualizations of the secondary structures produced. The fasta file used for the visualizations must only contain a single sequence.