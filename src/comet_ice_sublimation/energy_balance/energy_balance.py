from dataclasses import dataclass
import math

import numpy as np

from comet_ice_sublimation.heat_of_sublimation import *
from comet_ice_sublimation.model_input import *
from comet_ice_sublimation.physical_constants import *


@dataclass
class SublimationRateIterationResult:
    # sublimation rate for the model, molecules per cm^2 per second
    z: float
    # if model converges, the resulting surface temperature
    # if not converged, suggested next temperature to try
    t_K: float
    # if energy balances succeeds, the model is converged and this is set to True
    converged: bool


def converge_energy_balance(
    smi: SublimationModelInput,
    average_projection_factor: float,
    t_init_K: float,
    num_iterations_max: int = 100000,
) -> SublimationRateIterationResult:

    niter = 0
    cur_temp_K = t_init_K

    if average_projection_factor > 0:
        while niter < num_iterations_max:
            niter += 1
            srir = energy_balance(smi, average_projection_factor, cur_temp_K)
            if srir.converged:
                break
            cur_temp_K = srir.t_K
        else:
            raise RuntimeError("Energy balance iteration did not converge.")
    else:
        return SublimationRateIterationResult(z=0.0, t_K=np.nan, converged=True)

    return srir


def energy_balance(
    smi: SublimationModelInput,
    frac: float,
    t_K: float,
) -> SublimationRateIterationResult:
    """Calculate temperature and sublimation rate.


    Parameters
    ----------
    species: MolecularSpecies
        Inputted species

    Av : float
        Visual albedo (Av > 0)

    Air : float
        Infrared albedo

    rh_au: float
        Heliocentric distance (au)

    frac: float
        Insolation scale factor of this obliquity and latitude

    t_K: float
        Estimated temperature at latitude (Kelvins).


    Returns
    -------
    z : float
        Sublimation rate.

    temperature : float
        Updated temperature estimate (Kelvins).

    converged : bool
        `True` if the energy equation is balanced.

    """

    hosd = heat_of_sublimation(species=smi.species, t_K=t_K)

    root = 1 / math.sqrt(hosd.mass_g * 2 * math.pi * boltz)
    root_t = math.sqrt(t_K)
    sun = f0 * frac * (1.0 - smi.visual_albedo) / smi.rh_au**2
    radiat = (1 - smi.infrared_albedo) * sigma * t_K**4
    evap = root / root_t * hosd.pressure * hosd.latent_heat_of_vaporization
    phi = radiat + evap - sun
    z = max(evap / hosd.latent_heat_of_vaporization, 1e-30)

    drad = 4 * radiat / t_K
    x1 = hosd.pressure_prime * hosd.latent_heat_of_vaporization
    x2 = hosd.pressure * hosd.latent_heat_of_vaporization_prime

    devap = root / root_t * (x1 + x2)
    phipri = drad + devap

    dt = math.copysign(min(10, abs(phi / phipri / 2)), phi / phipri)
    t_K -= dt

    convergence_threshold = 1e-6
    converged = (
        abs(phi / sun) < convergence_threshold or abs(phi) < convergence_threshold
    )

    return SublimationRateIterationResult(z=z, t_K=t_K, converged=converged)
