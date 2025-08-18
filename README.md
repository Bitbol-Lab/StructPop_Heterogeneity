# StructPop_Heterogeneity
This code was used to generate figures in the paper _Environment heterogeneity creates fast amplifiers of natural selection in graph-structured populations_, by Cecilia Fruet, Arthur Alexandre, Alia Abbara, Claude Loverdo, and Anne-Florence Bitbol.

## Requirements

This project requires Python 3.11+ and the following packages:

```
ipython==8.15.0
matplotlib==3.7.2
numba==0.57.1
numpy==1.24.3
```

You can install all dependencies with:
```
pip install -r requirements.txt
```

## Usage
Clone the repository:
```
git clone https://github.com/yourusername/structured-population.git
cd structured-population
```

In your preferred Python:
```
from struct_pop_delta import *
```

And, after setting all parameters, the simulation is run with
```
et, ci95et, ft, ci95ft, fp = fixation_probability(
      in_numbers,
      migration_matrix_line,
      fitnesses,
      nb_sim,
      nb_cycles,
      t,
      K,
      size_follow_numbers=1000,
      print_frequency=1
)
```
