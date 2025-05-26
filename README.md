# Co-SIMP CRDM Module for DarkSUSY

This repository provides a custom implementation of Co-Scattering Strongly Interacting Massive Particles (Co-SIMPs) probed via Cosmic Ray Dark Matter (CRDM) within the [DarkSUSY](http://www.darksusy.org) framework. The modifications support relativistic $2\to3$ inelastic scattering, attenuation modeling, and two-loop elastic processes to derive realistic constraints on DM–SM interactions.

The code supports all numerical results from Edvard Rørnes' master's thesis (2025), and is built on top of DarkSUSY 6.4.1. No other versions have been tested.

## Repository Structure

All new programs are prefixed with `DDCR_`. Key computational routines are located in the `my_replaceable/` directory:

### Core Files
| File                              | Description                                                        |
|-----------------------------------|--------------------------------------------------------------------|
| `DDCR_CRflux.f`                  | Computes tabulated interstellar CR flux for $\psi \in \{p, \text{He}, \text{C}, \text{O}\}$. |
| `DDCR_flux.f`                    | Computes the $2\to3$ Co-SIMP CRDM flux.                                |
| `DDCR_target_recoil.f`           | Computes the nuclear recoil spectrum.                              |
| `DDCR_sigtot.f`                  | Computes integrated $2\to3$ cross section.                             |
| `DDCR_countrate.f`               | Computes $\Gamma_\text{CRDM}^\text{Co-SIMP} / \Gamma_\text{DM}$.   |
| `DDCR_limits.f`                  | Numerically inverts $c_\psi$ to find $\sigma_\text{NR}$.           |
| `DDCR_limits2.f`                 | Optimized version for $2\to3$ limits (no attenuation).                 |
| `DDCR_EnergyLoss.f`              | Computes average energy loss $\langle \omega_\chi \rangle$ for $\chi\psi \rightarrow \chi\chi\psi$. |
| `DDCR_EnergyLoss2t2.f`           | Same, but for $2\to2$ elastic scattering.                              |

### CRDM Module Extensions
| File                                    | Description |
|-----------------------------------------|-------------|
| `dsddDMCRflux.f`, `dsddDMCRflux_2loop.f` | $2\to3$ and 2-loop CRDM flux calculation. |
| `dsddDMCRsigCR.f`, `dsddDMCRsigCRff.f`  | Differential $2\to3$ cross sections (with and without form factor). |
| `dsddDMCRinvertcpsi.f`                  | Inverts $\sigma_\text{NR} \to C$ numerically via bisection. |
| `dsddDMCREnergyLoss.f`, `dsddDMCREnergyLoss2t2.f` | Computes $\langle \omega_\chi \rangle\sigma$ or $\sigma$ directly. |
| `dsddDMCRsigtot.f`, `dsddDMCRsigtot_2loop.f` | Total cross section integrators. |
| `dsddDMCRsigCRs35int.f`, `dsddDMCRsigCRt14int.f` | Kinematic integrals. |
| `dsddDMCRsigCR_2loop.f`                | 2-loop differential cross section. |

## Configuration Options

Edit `dsddcrdm_init.f` to change simulation behavior:

- `CRDM_form_factor` – Enable form factors $G(Q^2)$.
- `CRDM_cs` – Use $|\mathcal{M}|^2 = C/s$ if `.true.`, else $|\mathcal{M}|^2 = C$.
- `CRDM_high_acc` – Use high-accuracy (slower) integration.
- `CRDM_tab` – Use tabulated $2\to3$ CRDM flux.
- `CRDM_2loop` – Include two-loop elastic process in detector interaction.
- `CRDM_attenuation` – Enable attenuation (default on).
- `CRDM_both` – Include 2-loop elastic in CRDM flux.
- `CRDM_EnergyLoss` – Used internally for average energy loss. **Do not change.**

To support these options, modify `dsddcom.h` accordingly. A default effective depth `Deff = 5 kpc` is assumed.

## Notes

- The mapping $\sigma_\text{NR} \mapsto C$ is not analytically invertible and is solved via a bisection method.
- All physical results were validated using a specific subset of configuration flags.
- Changing multiple toggles simultaneously may require additional validation.

## Related References

- Master's Thesis: *Cosmic-Ray Upscattered Co-SIMPs* (Edvard Rørnes, 2025). See [duo.uio](https://www.duo.uio.no/handle/10852/14/discover?query=Cosmic-ray+boosted+Co-SIMPs&sort_by=dc.date.issued_dt&order=DESC&rpp=100&filtertype=type&filter_relational_operator=equals&filter=Masteroppgave) for more details. (Note it will be made available ~ 120 days after this update, exact link will be updated once available)
- DarkSUSY: [https://github.com/DarkSUSY/DarkSUSY](https://github.com/DarkSUSY/DarkSUSY)
- This module: [https://github.com/EdvardRornes/CoSIMP_CRDM](https://github.com/EdvardRornes/CoSIMP_CRDM)

## License

This project is distributed for academic purposes. Feel free to use without citation.
