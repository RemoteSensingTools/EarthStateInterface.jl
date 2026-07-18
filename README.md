# EarthStateInterface.jl

`EarthStateInterface.jl` defines small, physically explicit contracts for exchanging
Earth-system state between Julia packages. The first implemented contract is a single
atmospheric column suitable for transport diagnostics, spectroscopy, radiative transfer,
and GCHP/GEOS-Chem ingestion.

It is intentionally not a model and not an I/O framework.

## Design principles

- Physical meaning is part of the type: dry VMR, wet VMR, layer mass, and number density
  cannot be silently interchanged.
- Pressure and all other core quantities use SI units.
- Vertical orientation is explicit.
- A physical state does not contain derived optical properties.
- Layer merging is an explicit `LayerPartition`, not an ambiguous target layer count.
- State-space and optics-space reduction are distinct APIs.
- External packages can implement the accessor contract for their own storage; they do
  not have to inherit from a package-owned state hierarchy.
- Hydrostatic air-column diagnostics dispatch on an explicit gravity model; using a
  reanalysis constant or retrieval-style latitude/altitude correction is never implicit.

## Atmospheric column

```julia
using EarthStateInterface

column = AtmosphericColumn(
    pressure = PressureCoordinate([0.0, 20_000, 60_000, 100_000]),
    temperature = [220.0, 260.0, 290.0],
    specific_humidity = [0.0, 0.003, 0.01],
    dry_air_mass = [1.0, 3.0, 5.0],
    dry_air_column_moles = [700.0, 1400.0, 2100.0],
    trace_gases = (
        CO2 = ConstituentField([410e-6, 412e-6, 415e-6], DryMoleFraction()),
    ),
)
```

## Gravity-aware dry and wet columns

Specific humidity and pressure thickness are sufficient to diagnose three distinct
amounts per layer:

```julia
amounts = air_column_amounts(
    column,
    SomiglianaAltitudeGravity();
    latitude = 34.2,                         # geodetic degrees
    altitude = [30_000.0, 10_000.0, 1_000.0], # geometric layer centers [m]
)

dry_air_moles(amounts)       # mol m⁻²
water_vapor_moles(amounts)   # mol m⁻²
wet_air_moles(amounts)       # their sum: the complete moist mixture
```

Available models isolate each assumption: `ConstantGravity`,
`HelmertLatitudeGravity` (the RRTMGP latitude formula),
`SphericalAltitudeGravity`, and `SomiglianaAltitudeGravity`. See
[`examples/compare_gravity_models.jl`](examples/compare_gravity_models.jl) for a direct
impact comparison.

## Two physically different layer-reduction routes

```julia
plan = LayerPartition([1:2, 3:3], 3)

# Approximate for nonlinear spectroscopy: merge p/T/composition first.
coarse_column = merge_column(column, plan)

# Exact for optical depth at the wavelengths already evaluated.
coarse_tau = merge_optical_depth(native_layer_tau, plan)

# Equivalent effective cross section for a merged layer.
sigma_eff = effective_cross_section(native_sigma, absorber_amount, plan)
```

See [the design survey](docs/design_survey.md) for how this maps to AtmosTransport,
vSmartMOM, RRTMGP, and GCHP. Concrete migration boundaries are recorded in the
[adapter notes](docs/adapter_notes.md).

## Optional GCHP reader

Loading `NCDatasets` activates a reader for one or more GCHP History collections:

```julia
using EarthStateInterface, NCDatasets

column = open_gchp([state_met_path, species_path]) do source
    read_column(source, GCHPLocation(x=12, y=8, face=3))
end
```

Sectional aerosol ingestion additionally requires an explicit `SectionalSizeGrid`.
The reader does not infer physical diameters from TOMAS mass-bin indices using an
arbitrary density.
