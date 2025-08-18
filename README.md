# Structured populations with heterogeneity
This code was used to generate figures in the paper _Environment heterogeneity creates fast amplifiers of natural selection in graph-structured populations_, by Cecilia Fruet, Arthur Alexandre, Alia Abbara, Claude Loverdo, and Anne-Florence Bitbol.

## Requirements

This project requires the following packages:

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
```bash
git clone https://github.com/yourusername/structured-population.git
cd structured-population
```

In Python, import the main module:
```python
from struct_pop_delta import *
```

And, after setting all parameters: 
- The inital number of individuals in each deme `in_numbers`
- The migration matrix `MigMat`
- Fitnesses `fs`
- Number of simulations `nb_sims`
- Number of cycles `nb_cycles`
- Growth time `t`
- Bottleneck size `k`

the simulation is run with:
```python
et, ci95et, ft, ci95ft, fp = fixation_probability(
      in_numbers,
      MigMat,
      fs,
      nb_sim,
      nb_cycles,
      t,
      k,
      size_follow_numbers=1000,
      print_frequency=1
)
```

`et` and `ft` are the extinction and fixation times for mutants, averaged over `nb_sim` simulations. `ci95et` and `ci95ft` are the 95% confidence intervals on those times.
`fp` is the mutant fixation probability. Demo notebook [here](Example_simulation.ipynb) for an example simulation.

## License

This project is licensed under the MIT License.

