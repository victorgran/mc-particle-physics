# Monte Carlo event generator for the simulation of electron-positron annihilation to quark-antiquark production 

## Set up

Clone the repository, create a virtual environment, activate it, and install the python package (requires Python >= 3.8, I didn't test with Python > 3.10, though).
```shell
git clone https://github.com/VictorG20/mc-particle-physics.git
cd mc-particle-physics  # Access the repository directory.
python -m venv .my_venv  # Create a virtual environment.
. .my_venv/bin/activate  # Activate it.
pip install ./simulator # Install the python package.
```

## Checking the exercises

Having installed the `simulator` package, the python scripts from [exercises/](exercises) can be executed from the terminal as follows
```bash
(my_venv) ... $ python exercises/exercise1d.py
```
By default, no images are saved and relevant data is loaded from [data/](data) when available to avoid computing it. This behavior can be changed within each script by setting the relevant variables to False.

**Note:** Every random process is seeded by default so, technically (hopefully), everything should be reproducible and the values should match those in the report.

## Relevant files

1. The (relevant) code required for parts a), c) and f) of exercise 1, is distributed along the files:
   * [distributions.py](simulator/src/simulator/integrator/distributions.py): contains the implementation of the Dirac, Uniform, and Breit-Wigner distributions.
   * [integrator.py](simulator/src/simulator/integrator/integrator.py): Monte Carlo integrator.
2. The use of `vegas` is contained in the corresponding exercise script [exercise1e.py](exercises/exercise1e.py) except for the plotting of the grids, which is in [grid.py](simulator/src/simulator/plotting/grid.py).
3. The code for exercise 2, part b), is split into two parts:
   1. A function which generates events from the Monte Carlo integrator sampling points, [exercise2b.py](exercises/exercise2b.py).
   2. The actual use of the parton shower class, [exercise2c.py](exercises/exercise2c.py).
4. The implementation of the Durham algorithm within the `Analysis` class is in the file [analysis.py](simulator/src/simulator/utils/analysis.py).

## Instructions

### How to start an instance of the Monte Carlo integrator

Since time is of essence, the Monte Carlo integrator has the squared matrix element hard-coded, so that's the only thing it can really integrate. Therefore, the only thing that remains is to specify distributions for the kinematic variables $s, \cos{\theta}$ and $\phi$, as well as the beam distribution $f(s)$ and the method to sum over light quark flavours. 

By default, an instance of the Integrator with no arguments yields the following:
* A fixed beam distribution at the mass of the Z boson, $f(s) = \delta(s - M_{Z}^{2})$
* Uniform distributions for $-1 \leq \cos{\theta} < 1$ and $0 \leq \phi < 2 \pi $.
* An explicit sum over light quark flavours.

```python3.10
from simulator import MonteCarloIntegrator

sample_size = 100_000  # Number of Monte Carlo points to generate.
integrator_b = MonteCarloIntegrator()  # Default integrator.
integrator_b.sampleDeltaSigma(sample_size)
sigma, mc_error = integrator_b.integrateCrossSection()
```

After integrating the cross-section, one can sample more the same integrator, in which case the total sample size is the sum of sample sizes passed, e.g.
```python3.10
from simulator import MonteCarloIntegrator

sample_size = 100_000  # Number of Monte Carlo points to generate.
integrator_b = MonteCarloIntegrator()  # Default integrator.
integrator_b.sampleDeltaSigma(sample_size)  # 100,000 samples are generated.
sigma, mc_error = integrator_b.integrateCrossSection()
integrator_b.sampleDeltaSigma(sample_size)  # 100,000 samples MORE are generated.
sigma2, mc_error2 = integrator_b.integrateCrossSection()  # Estimates for a sample size of 200,000.
```

Other distributions, such as a uniform or Breit-Wigner distributions can be defined and used to instantiate the Monte Carlo integrator as follows
```python3.10
from simulator import MonteCarloIntegrator, Uniform, BreitWigner, ZBoson

Z = ZBoson()
s_min = (Z.mass - 3 * Z.width) ** 2
s_max = (Z.mass + 3 * Z.width) ** 2
s_distro = BreitWigner(s_min=s_min, s_max=s_max, mass=Z.mass, decay_width=Z.width)
beam_distro = Uniform(lower=s_min, upper=s_max)
integrator = MonteCarloIntegrator(s_distro=s_distro, beam_distro=beam_distro, sum_quark_method="random")
```
