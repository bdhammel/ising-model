import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from tqdm import tqdm
import click


class IsingLattice:

    def __init__(self, temperature, initial_state, size):
        self.size = size
        self.T = temperature
        self.system = self._build_system(initial_state)

    @property
    def sqr_size(self):
        return (self.size, self.size)

    def _build_system(self, initial_state):
        """Build the system

        Build either a randomly distributed system or a homogeneous system (for
        watching the deterioration of magnetization

        Parameters
        ----------
        initial_state : str: "r" or other
            Initial state of the lattice.  currently only random ("r") initial
            state, or uniformly magnetized, is supported
        """

        if initial_state == 'r':
            system = np.random.choice([-1, 1], self.sqr_size)
        elif initial_state == 'u':
            system = np.ones(self.sqr_size)
        else:
            raise ValueError(
                "Initial State must be 'r', random, or 'u', uniform"
            )

        return system

    def _bc(self, i):
        """Apply periodic boundary condition

        Check if a lattice site coordinate falls out of bounds. If it does,
        apply periodic boundary condition

        Assumes lattice is square

        Parameters
        ----------
        i : int
            lattice site coordinate

        Return
        ------
        int
            corrected lattice site coordinate
        """
        if i >= self.size:
            return 0
        if i < 0:
            return self.size - 1
        else:
            return i

    def energy(self, N, M):
        """Calculate the energy of spin interaction at a given lattice site
        i.e. the interaction of a Spin at lattice site n,m with its 4 neighbors

        - S_n,m*(S_n+1,m + Sn-1,m + S_n,m-1, + S_n,m+1)

        Parameters
        ----------
        N : int
            lattice site coordinate
        M : int
            lattice site coordinate

        Return
        ------
        float
            energy of the site
        """
        return -2*self.system[N, M]*(
            self.system[self._bc(N - 1), M] + self.system[self._bc(N + 1), M]
            + self.system[N, self._bc(M - 1)] + self.system[N, self._bc(M + 1)]
        )

    @property
    def internal_energy(self):
        e = 0
        E = 0
        E_2 = 0

        for i in range(self.size):
            for j in range(self.size):
                e = self.energy(i, j)
                E += e
                E_2 += e**2

        U = (1./self.size**2)*E
        U_2 = (1./self.size**2)*E_2

        return U, U_2

    @property
    def heat_capacity(self):
        U, U_2 = self.internal_energy
        return U_2 - U**2

    @property
    def magnetization(self):
        """Find the overall magnetization of the system
        """
        return np.abs(np.sum(self.system)/self.size**2)


def run(lattice, epochs, video=True):
    """Run the simulation
    """

    FFMpegWriter = manimation.writers['ffmpeg']
    writer = FFMpegWriter(fps=10)

    fig = plt.figure()

    with writer.saving(fig, "ising.mp4", 100):
        for epoch in tqdm(range(epochs)):
            # Randomly select a site on the lattice
            N, M = np.random.randint(0, lattice.size, 2)

            # Calculate energy of a flipped spin
            E = -1*lattice.energy(N, M)

            # "Roll the dice" to see if the spin is flipped
            if E <= 0.:
                lattice.system[N, M] *= -1
            elif np.exp(-E/lattice.T) > np.random.rand():
                lattice.system[N, M] *= -1

            if video and epoch % (epochs//75) == 0:
                img = plt.imshow(
                    lattice.system, interpolation='nearest', cmap='jet'
                )
                writer.grab_frame()
                img.remove()

    plt.close('all')


@click.command()
@click.option(
    '--temperature', '-t',
    default=0.5,
    show_default=True,
    help='temperature of the system'
)
@click.option(
    '--initial-state', '-i',
    default='r',
    type=click.Choice(['r', 'u'], case_sensitive=False),
    show_default=True,
    help='(R)andom or (U)niform initial state of the system'
)
@click.option(
    '--size', '-s',
    default=100,
    show_default=True,
    help='Number of sites, M, in the MxM lattice'
)
@click.option(
    '--epochs', '-e',
    default=1_000_000,
    type=int,
    show_default=True,
    help='Number of iterations to run the simulation for'
)
@click.option(
    '--video',
    is_flag=True,
    help='Record a video of the simulation progression'
)
def main(temperature, initial_state, size, epochs, video):
    lattice = IsingLattice(
        temperature=temperature, initial_state=initial_state, size=size
    )
    run(lattice, epochs, video)

    print(f"{'Net Magnetization [%]:':.<25}{lattice.magnetization:.2f}")
    print(f"{'Heat Capacity [AU]:':.<25}{lattice.heat_capacity:.2f}")


if __name__ == "__main__":
    plt.ion()
    main()
