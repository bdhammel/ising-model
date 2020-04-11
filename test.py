from unittest import mock

import numpy as np
import pytest

from ising import IsingLattice


@pytest.mark.parametrize('i, t', [
    (100, 0),  # a value at 100 wraps around to 0
    (99, 99),  # a value at 99 is in the lattice
    (0, 0),    # a value at 0 is in the lattice
    (-1, 99)   # a value at -1 wraps around to 99
])
def test_boundary_conditions(i, t):
    mock_self = mock.Mock()
    mock_self.size = 100
    r = IsingLattice._bc(mock_self, i)
    assert r == t  # is the returned value the true value


@pytest.mark.parametrize('i, system, e', [
    ((1, 1), np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]), -8),           # most probable state
    ((1, 1), np.array([[-1, -1, -1], [-1, -1, -1], [-1, -1, -1]]), -8),  # most probable state
    ((1, 1), np.array([[-1, -1, -1], [-1, 1, -1], [-1, -1, -1]]), 8),    # least probable state
    ((1, 1), np.array([[1, -1, 1], [-1, 1, 1], [1, 1, 1]]), 0),          # 50/50 state
])
def test_energy_calculation_is_correct(i, system, e):
    lattice = IsingLattice(temperature=1, initial_state='u', size=3)
    lattice.system = system
    _e = lattice.energy(*i)
    assert _e == e  # is the returned value the true value


@pytest.mark.skip("TODO")
def check_internel_energy_is_correct():
    pass


@pytest.mark.parametrize('system, m', [
    (np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]), 1),           # uniform magnetization
    (np.array([[-1, -1, -1], [-1, -1, -1], [-1, -1, -1]]), 1),  # uniform magnetization
    (np.array([[-1, -1, -1], [-1, -1, 1], [1, 1, 1]]), 1/9),    # net magnetization (i.e. sum of system) 1/9
])
def test_magnitization_correct(system, m):
    lattice = IsingLattice(temperature=1, initial_state='u', size=3)
    lattice.system = system
    _m = lattice.magnetization
    assert _m == m
