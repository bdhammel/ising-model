import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from scipy.misc import factorial
import io


class IsingLattice:

    _EPOCHS = int(1e6)

    def __init__(self, tempature, initial_state='r', size=(100,100), show=False):
        """
        """
        self.sqr_size = size
        self.size = size[0]
        self.T = tempature
        self._build_system(initial_state)

    def _build_system(self, initial_state):
        """Build the system
        Build either a randomly distributed system or a homogeneous system (for
        watching the deterioration of magnetization

        Args
        ----
        initial_state (str: "r" or other) : Initial state of the lattice.
            currently only random ("r") initial state, or uniformly magnetized, 
            is supported 
        """

        if initial_state == 'r':
            system = np.random.randint(0, 1+1, self.sqr_size)
            system[system==0] = -1
        else:
            system = np.ones(self.sqr_size)

        self.system = system

    def _bc(self, i):
        """Apply periodic boundary condition

        Check if a lattice site coordinate falls out of bounds. If it does, 
        apply periodic boundary condition 
        
        Assumes lattice is square

        Args
        ----
        i (int) : lattice site coordinate

        Return
        ------
        (int) : corrected lattice site coordinate
        """
        if i+1 > self.size-1:
            return 0
        if i-1 < 0:
            return self.size-1
        else:
            return i

    def _energy(self, N, M):
        """Calculate the energy of spin interaction at a given lattice site
        i.e. the interaction of a Spin at lattice site n,m with its 4 neighbors

        - S_n,m*(S_n+1,m + Sn-1,m + S_n,m-1, + S_n,m+1)

        Args
        ----
        N (int) : lattice site coordinate
        M (int) : lattice site coordinate

        Return
        """
        return -self.system[N,M]*(
                    self.system[self._bc(N-1), M] 
                    + self.system[self._bc(N+1), M] 
                    + self.system[N, self._bc(M-1)] 
                    + self.system[N, self._bc(M+1)]
                    )

    @property
    def internal_energy(self):
        e=0; E=0; E_2=0

        for i in range(self.size):
            for j in range(self.size):   
                e = self._energy(i, j)
                E += e
                E_2 += e**2

        U = (1./self.size**2)*E
        U_2 = (1./self.size**2)*E_2 

        return U, U_2

    '''
    @staticmethod
    def _factorial(n):
        """Calculates N!

        If N is too large, approximate it using Sterling's approximation
        N! = N*log N - N
        """
        if n>100:
            return n*np.log(n)-n
        else:
            return np.log(factorial(n))
        
    @property
    def entropy(self):
        """Find the entropy of the system

        Calculates N!/(N_up! * Nb_down!)
        """

        nup = self.system[self.system == 1].size

        return self._factorial(self.size**2) / \
                (self._factorial(nup)*self._factorial(self.size**2-nup))
    '''

    @property
    def magnetization(self):
        """Find the overall magnetization of the system
        """
        return np.abs(np.sum(self.system)/self.size**2)

    def run(self, video=True):
        """Run the simulation
        """

        FFMpegWriter = manimation.writers['ffmpeg']
        writer = FFMpegWriter(fps=10)

        plt.ion()
        fig = plt.figure()

        with writer.saving(fig, "ising.mp4", 100):
            for epoch in range(self._EPOCHS):
                # Randomly select a site on the lattice
                N, M = np.random.randint(0, self.size, 2)

                # Calculate energy of a flipped spin
                E = -1*self._energy(N, M)
                
                # "Roll the dice" to see if the spin is flipped 
                if E <= 0.:
                    self.system[N,M]*=-1
                elif np.exp(-E/self.T) > np.random.rand():
                    self.system[N,M]*=-1

                if epoch % (self._EPOCHS//75) == 0: 
                    print("Epoch {:.2E} / {:.2E}".format(epoch, self._EPOCHS), end='\t')
                    print("Percent completed {:.2f}%".format(100*epoch/self._EPOCHS), end='\t')
                    print("Net Magnetization: {:.2f}".format(self.magnetization))
                    if video:
                        img = plt.imshow(self.system, interpolation='nearest')
                        writer.grab_frame()
                        img.remove()
            print("...done")

        plt.close('all')

if __name__ == "__main__":

    lattice = IsingLattice(tempature=.50, initial_state="r", size=(100,100))
    lattice.run()

