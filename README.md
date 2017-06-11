Python code to simulate the Ising model of a Ferromagnet. For a discussion of the theory, visit [my blog post](https://bdhammel.github.io/2017/06/10/ising-model.html).


The initial conditions of the ising lattice can be specifited by the `tempature`, `initial state`, and `size parameters` of the model.

Running the simulation will output a video of system as it changes through out the run steps.

__Example exicution:__

~~~python
from ising import IsingLattice
lattice = IsingLattice(tempature=.50, initial_state="r", size=(100,100))
lattice.run(video=True)
~~~
