from numpy import pi

BOHR = 0.52917721092  # Angstroms
BOHR_SI = BOHR * 1e-10
ATOMIC_MASS = 1e-3 / 6.022140857e23
HARTREE2J = 4.359744650e-18
HARTREE2EV = 27.21138602
LIGHT_SPEED_SI = 299792458

AU2Hz = ((HARTREE2J / (ATOMIC_MASS * BOHR_SI ** 2)) ** 0.5 / (2 * pi))

periodic_table = {
    'H': {'Z': 1.008, 'nelec': 1},
    'C': {'Z': 12.011, 'nelec': 6},
    'N': {'Z': 14.007, 'nelec': 7},
    'O': {'Z': 15.999, 'nelec': 8},
    'S': {'Z': 32.06, 'nelec': 16},
    'Br': {'Z': 79.904, 'nelec': 35},
    'I': {'Z': 126.90, 'nelec': 53}
}
