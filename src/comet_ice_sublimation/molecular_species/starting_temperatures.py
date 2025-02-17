from .molecular_species import *


def get_starting_temperature(species: MolecularSpecies) -> float:

    start_t = {
        MolecularSpecies.h2o: 190,
        MolecularSpecies.h2o_ch4: 190,
        MolecularSpecies.co2: 100,
        MolecularSpecies.co: 60,
    }

    return start_t[species]
