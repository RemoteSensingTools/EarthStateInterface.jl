# Atmospheric-column design survey

This document records the initial design evidence. It is intentionally candid: the new
contract borrows useful patterns but does not preserve legacy APIs.

## AtmosTransport

Useful foundation:

- `AtmosGrid` separates horizontal mesh, vertical coordinate, architecture, and planetary
  parameters.
- `CellState` stores dry-air and tracer mass explicitly, making conservation testable.
- `MetState` keeps upstream meteorology separate from prognostic transport state.

Interface consequence: an AtmosTransport column adapter should be a view assembled from
grid, met, and tracer state. The shared package must not force the global transport state
into a radiative-transfer-shaped allocation.

## vSmartMOM

Current `AtmosphericProfile` combines source state and derived quantities (`p_full`,
`p_half`, `q`, water VMR, dry/wet VCD, gas VMR, and layer height). It uses concrete CPU
arrays, string-keyed dictionaries, and positional construction. Its two reduction paths
also hide materially different physics behind `reduce_profile(n; binavg=false)`.

Problems to remove:

- no positional constructor for a many-field physical state;
- no independent mutable copies of quantities that can be diagnosed from primary state;
- no unqualified `reduce_profile(n)`;
- no arithmetic averaging of dry VMR or aerosol state;
- no assumption that an interpolated coarse p/T state is optically equivalent to native
  layers.

## RRTMGP.jl

Useful hints from [CliMA/RRTMGP.jl](https://github.com/CliMA/RRTMGP.jl):

- atmospheric, cloud, and aerosol payloads are separate;
- arrays are device-generic and adaptable;
- host models exchange state through a getter contract;
- layer and interface quantities are recognized as distinct.

Patterns not adopted:

- packed `layerdata` whose first index implicitly means dry column, pressure,
  temperature, or RH;
- pressure documented as either Pa or mb;
- gas identity determined by lookup-table integer position;
- aerosol physical inputs and derived AOD stored in one `AerosolState`;
- a solver-owned atmospheric state as the ecosystem interchange object.

RRTMGP remains a solver adapter target, not the definition of the shared contract.

RRTMGP's dry-column kernel does provide one useful parity target: when latitude is
available it uses the Helmert approximation `g = 9.80665 - 0.02586 cos(2φ)` and otherwise
uses one constant gravity. `HelmertLatitudeGravity` reproduces that choice. The shared API
also exposes altitude-only and latitude-plus-altitude models instead of baking either
retrieval or reanalysis practice into `AtmosphericColumn`.

### Future radiation-intercomparison criterion

An RRTMGP-versus-vSmartMOM experiment is a future use case, not a current package
milestone. The interface should nevertheless make it straightforward for one physical
column to drive both solvers. In particular:

- neither solver's state container or spectral representation becomes the shared state;
- layer and interface quantities retain unambiguous orientation and units;
- gas, aerosol, and future cloud inputs describe physical state rather than derived
  two-stream or multi-stream optical properties;
- gravity, layer reduction, and concentration-basis choices remain explicit, allowing a
  comparison to hold them fixed;
- adapters may independently select correlated-k, lookup-table, or line-by-line optics
  without changing the source atmospheric column.

A future comparison should first hold spectral layer optics fixed while varying only the
radiative-transfer solver and stream count. More detailed spectroscopy, phase functions,
and cloud microphysics can then be introduced separately. Flux profiles, heating rates,
solver configuration, and spectral quadrature are radiative-transfer outputs or controls;
they do not belong in `EarthStateInterface` merely to support that experiment.

## Hydrostatic air-column amounts

For total pressure thickness `Δp` and specific humidity `q`, the diagnosed moist-air
mass per area is `Δp/g`. Water and dry-air mass are then `q Δp/g` and
`(1-q) Δp/g`, respectively. Dividing each by its own molar mass gives water-vapor and
dry-air moles; their sum is the total wet-air amount. This direct mass split avoids the
common ambiguity between wet- and dry-basis water VMR.

Gravity is a diagnostic strategy selected by multiple dispatch:

- `ConstantGravity` matches the usual reanalysis convention;
- `HelmertLatitudeGravity` isolates latitude and matches RRTMGP;
- `SphericalAltitudeGravity` isolates inverse-square altitude dependence;
- `SomiglianaAltitudeGravity` combines WGS 84 surface gravity with an explicitly
  documented inverse-square altitude continuation.

Altitude-aware calculations require geometric layer-center height. The interface does not
silently treat geopotential height as geometric altitude. The pressure state itself remains
independent of the selected gravity diagnostic, so the same column can be compared under
all four assumptions.

## GCHP / GEOS-Chem

A useful GCHP column can carry:

- `Met_DELP` / `Met_DELPDRY`, surface pressure, temperature, and specific humidity;
- `Met_AD` dry-air layer mass and `Met_AIRVOL` layer volume;
- dry pressure thickness, which gives dry-air column moles without requiring cell area;
- `SpeciesConcVV_<gas>` trace gases in mol mol⁻¹ dry;
- TOMAS bin number (`NK01...`) and per-bin component fields (`SF`, `SS`, `DUST`,
  `OCOB`, `OCIL`, `ECOB`, `ECIL`, `AW`, etc.).

This motivates a column with independent thermodynamic, gas, and sectional-aerosol
capabilities. TOMAS scheme metadata, refractive-index databases, mixing rules, and Mie
calculations remain outside the base state contract.

## Layer reduction and spectroscopy

There are two valid but different workflows.

### State then optics

1. Define an explicit partition of native layers.
2. Sum extensive quantities.
3. Average intensive quantities with their physically appropriate carrier amount.
4. Compute cross sections or aerosol optics at representative coarse-layer state.

This is efficient but approximate because cross sections are nonlinear in pressure and
temperature, while aerosol optics are nonlinear in humidity, size, and composition.

### Optics then merge

1. Compute cross sections and absorber amounts per native layer.
2. Form native optical depths.
3. Sum optical depths within each target layer.
4. If a coarse-layer cross section is required, diagnose
   `sigma_eff = sum(sigma_i * N_i) / sum(N_i)`.

This exactly preserves monochromatic absorption through the merge. Aerosol extinction
and scattering optical depths likewise add. Single-scattering albedo is reconstructed
from the summed scattering and extinction depths; phase moments must be weighted by
scattering optical depth.

The interface names these routes `StateThenOptics` and `OpticsThenMerge`; neither is an
undocumented Boolean option.

## Next adapters

1. Optional NCDatasets-based GCHP reader producing `AtmosphericColumn` plus
   `SectionalAerosolState`.
2. Zero-copy AtmosTransport column view.
3. vSmartMOM adapter replacing `AtmosphericProfile` at the model-construction boundary.
4. RRTMGP adapter filling its solver state through its getter/update contract.
