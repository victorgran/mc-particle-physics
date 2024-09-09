# Monte Carlo event generator for the simulation of electron-positron annihilation to quark-antiquark production 

## Instructions


### Set up

In order to have a not-so-painful experience:

```bash
$ git clone https://github.com/VictorG20/mc-particle-physics.git
$ cd mc-particle-physics  # Access the repository directory.
$ python -m venv .my_venv  # Create a virtual environment.
$ . .my_venv/bin/activate  # Activate it.
(my_venv) ... $ pip install -r requirements.txt # Install the python package with all its dependencies.
```

### Checking the exercises

Having installed the `simulator` package, the python scripts from `exercises/` can be executed from the terminal as follows
```bash
(my_venv) ... $ python exercise1/part_d.py
```

### How to start an instance of the Monte Carlo integrator

Since time is of essence, the Monte Carlo integrator has the squared matrix element hard-coded, so that's the only thing it can really integrate. Therefore, the only thing that remains is to specify distributions for the kinematic variables $s, \cos{\theta}$ and $\phi$, as well as the beam distribution $f(s)$ and the method to sum over light quark flavours. 

By default, an instance of the Integrator with no arguments yields the following:
* A fixed beam distribution at the mass of the Z boson, $f(s) = \delta(s - M_{Z}^{2})$
* Uniform distributions for $-1 \leq \cos{\theta} < 1$ and $0 \leq \phi < 2 \pi $.
* An explicit sum over light quark flavours.

```python
from simulator import MonteCarloIntegrator

sample_size = 100_000  # Number of Monte Carlo points to generate.
integrator_b = MonteCarloIntegrator()  # Default integrator.
integrator_b.sampleDeltaSigma(sample_size)
sigma, mc_error = integrator_b.integrateCrossSection()
```

After integrating the cross-section, one can sample more the same integrator, in which case the total sample size is the sum of sample sizes passed, e.g.
```python
from simulator import MonteCarloIntegrator

sample_size = 100_000  # Number of Monte Carlo points to generate.
integrator_b = MonteCarloIntegrator()  # Default integrator.
integrator_b.sampleDeltaSigma(sample_size)  # 100,000 samples are generated.
sigma, mc_error = integrator_b.integrateCrossSection()
integrator_b.sampleDeltaSigma(sample_size)  # 100,000 samples MORE are generated.
sigma2, mc_error2 = integrator_b.integrateCrossSection()  # Estimates for a sample size of 200,000.
```

Other distributions, such as a uniform or Breit-Wigner distributions can be defined and used to instantiate the Monte Carlo integrator as follows
```python
from simulator import MonteCarloIntegrator, Uniform, BreitWigner, ZBoson

Z = ZBoson()
s_min = (Z.mass - 3 * Z.width) ** 2
s_max = (Z.mass + 3 * Z.width) ** 2
s_distro = BreitWigner(s_min=s_min, s_max=s_max, mass=Z.mass, decay_width=Z.width)
beam_distro = Uniform(lower=s_min, upper=s_max)
integrator = MonteCarloIntegrator(s_distro=s_distro, beam_distro=beam_distro, sum_quark_method="random")
```

## Where are the relevant parts?

If you are the tutor reviewing the projects, you might be fed up with having to go through more or less the same things over and over again. A quick way to go through this repository, depending on what you want to have a look at, is the following:
* Where is the code for exercise 1, parts a), c) and f), Lebowski?
  * All the code for the integrator (or at least the one I could think you might want to see) is within these two files: [`integrator.py`](src/simulator/integrator/integrator.py) and [`distributions.py`](src/simulator/integrator/distributions.py) (the second one is relevant only because it has the definition of the Breit-Wigner distribution)
* How did you get your figures and results for exercise 1?
  * See (or just execute) the scripts in [`exercise1`](exercises). If executed, they will display all the relevant information one is supposed to get from them. By default, nothing is saved and, when available, the relevant data in [`data/`](data) is loaded to avoid computing stuff (although most of it should run in seconds at most).


**Note:** Every random process is seeded by default so, technically (hopefully), everything should be reproducible.