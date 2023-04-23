# Gravitools

This  package is a collection of tools to analyze star cluster simulations

# Prerequisites
- Python 3.x
- NumPy
- Pandas
- Numba
- math


# Installation
Clone the repository: 

> git clone https://github.com/FrancescopR/ML_bid_strategy_optimizer.git

Install the required packages:


> pip install -e .


# Usage
Usage example on how to compute core radius/density centre of a simulated star cluster:

```Python
import gravitools.starcluster_tools as sct

Rho, Xc, Vc        = density_center(pos, vel, mass)
rc, Mc, Rho_c, n_c = sct.core_radius_mass_density(R, M, Rho)
```

# https://www.youtube.com/watch?v=DhUpxWjOhME


# Data



# Authors

Francesco P. Rizzuto - Initial work
