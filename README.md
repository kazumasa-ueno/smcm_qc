# Quantum Algorithms for Stochastic Multicloud Model

Fortran and Python code for simulating Stochastic Multicloud Model (Khouider et al. 2010). 

## Files
`deterministic.f90` is the code for deterministic simulation.

`montecarlo.f90` is the code for classical Monte Carlo calculation.

`plot_for_article` is the jupyter notebook for plotting results.

`quantum.py` is the main quantum simulation algorithm working on classical computers using Qiskit Aer simulator.

## Install
```bash
git clone https://github.com/kazumasa-ueno/smcm_qc.git
cd smcm_qc
poetry install
```

## How to run quantum code
```bash
poetry run python quantum.py
```

## Required python libraries and tested versions
| Library | Version |
|--------------|------------|
| numpy        | 1.23.0     |
| scipy        | 1.13.1     |
| matplotlib   | 3.9.0      |
| qiskit       | 1.1.0      |
| qiskit-aer   | 0.14.2     |

It is important to note that versions of numpy other than 1.23.0 sometimes cause issues. For all dependencies, please refer to the `pyproject.toml` and `poetry.lock` files.

## Abstract
Quantum computers have attracted much attention in recent years. This may be because the development of the actual quantum machine is accelerating. Quantum computers use the properties of superposition, interference, and entanglement in quantum mechanics. Research on how to use quantum computers is active in the fields such as quantum chemistry and machine learning, where vast amounts of computation are required. However, in the field of weather and climate simulations, which also require vast amounts of computation, less research has been done on the use of quantum computers. In this study, a quantum computing algorithm was applied to a problem of the atmospheric science to demonstrate that it can achieve the same simulations as a conventional algorithm designed for the classical computers. The stochastic nature of a multi-cloud model fits well with quantum computers. One of the reasons is that probabilistic outputs can directly be obtained from computed quantum states. We have confirmed that stochastic fluctuations of cloud type ratios can be represented naturally on a quantum computer simulator. Our results demonstrate that quantum computers can suitably solve some types of the problems of the atmospheric and oceanic phenomena, in which stochasticity is widely inherent.