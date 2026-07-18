# Adapter notes

Adapters should live with the package that owns the source representation, or in a Julia
package extension. `EarthStateInterface` owns semantic functions and canonical exchange
containers; it should not depend on transport or radiative-transfer solvers.

## vSmartMOM

The transition can happen at `model_from_parameters`:

1. Parsers and GCHP readers produce an `AtmosphericColumn`.
2. Absorption consumes `pressure_interfaces`, `temperature`, and basis-tagged gas fields.
3. Dry-air column moles replace independently stored `vcd_dry`.
4. Water VMR, wet VCD, representative pressures, and layer height become diagnosed
   quantities or explicit cached diagnostics.
5. Profile reduction is removed from model construction. Callers choose an explicit
   `LayerPartition` and either `StateThenOptics` or `OpticsThenMerge`.

`examples/vsmartmom_profile_adapter.jl` is a temporary migration bridge for the current
`AtmosphericProfile`.

## AtmosTransport

An adapter must combine three owners of information:

- `AtmosGrid.vertical` supplies pressure interfaces from hybrid coefficients and local
  surface pressure;
- meteorological forcing supplies temperature and specific humidity;
- `CellState` supplies dry-air and tracer mass.

The correct implementation is a column view, not a copied global `AtmosphericColumn`
array. Trace gases represented as tracer mass must remain `LayerMass` until an explicit
conversion uses dry-air mass and molecular weights. This also permits signed transport
perturbations without pretending they are physical positive VMR fields.

## RRTMGP

RRTMGP already exposes a host getter contract and uses dry-air column molecules. Its
adapter should:

- convert molecules cm⁻² to mol m⁻²;
- attach gas names from the lookup table rather than preserve integer gas indices;
- declare whether its water VMR is dry- or wet-basis before converting to specific
  humidity;
- keep cloud and derived aerosol optics outside the physical aerosol state.

The unresolved water-basis declaration is intentionally not guessed in an executable
adapter.

## GCHP

Loading `NCDatasets` activates `open_gchp` and `read_column`. Multiple History files may
be opened together; variables are located across the collection. The reader:

- indexes cubed-sphere `(Xdim, Ydim, nf, lev, time)` variables by dimension name;
- converts pressure to Pa and humidity to kg kg⁻¹;
- normalizes surface-to-TOA GCHP layers to top-to-surface;
- reads gases as `DryMoleFraction`;
- reads `Met_AD`, `Met_AIRVOL`, and, when present, `Met_DELPDRY`;
- optionally converts TOMAS `NK` to number density and preserves component fields as
  dry-air mole fractions.

The TOMAS size grid is an explicit input. Converting a mass-bin boundary to physical
diameter requires density/composition assumptions and therefore does not belong in a
generic file reader.
