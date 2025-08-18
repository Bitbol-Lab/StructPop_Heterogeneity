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
- The migration matrix `mig_mat`
- Fitnesses `fs`
- Number of simulations `nb_sims`
- Number of cycles `nb_cycles`
- Growth time `t`
- Bottleneck size `k`

the simulation is run with:
```python
et, ci95et, ft, ci95ft, fp = fixation_probability(
      in_numbers,
      mig_mat,
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
`fp` is the mutant fixation probability estimated over `nb_sim` simulations. Demo notebook [here](Example_simulation.ipynb) for an example simulation.

This code was inspired by and generalizes the implementation in the [StructuredPop Repo](https://github.com/Bitbol-Lab/Structured_pop), used in the paper [Frequent asymmetric migrations suppress natural selection in spatially structured populations](https://academic.oup.com/pnasnexus/article/2/11/pgad392/7420192) by Alia Abbara and Anne-Florence Bitbol.

## Citation
If you use this code in your research, please cite our preprint:
```bibtex
@article{fruet2025heterogeneity,
  title={Environment heterogeneity creates fast amplifiers of natural selection in graph-structured populations},
  author={Cecilia Fruet and Arthur Alexandre and Alia Abbara and Claude Loverdo and Anne-Florence Bitbol},
  journal={bioRxiv},
  pages={025.07.31.667961},
  year={2025},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2025.07.31.667961},
  url={https://www.biorxiv.org/content/10.1101/2025.07.31.667961v1}
}
```

## License

This project is licensed under the MIT License.

