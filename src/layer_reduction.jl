"""
    LayerPartition(groups, source_layers)

Explicit, contiguous mapping from source layers to target layers. Every source layer must
occur exactly once, in storage order. Making this plan explicit avoids ambiguous APIs such
as `reduce_profile(n)` whose pressure grid and conservation behavior are implicit.
"""
struct LayerPartition{G}
    groups::G
    source_layers::Int
end

function LayerPartition(groups, source_layers::Integer)
    normalized = [Int(first(group)):Int(last(group)) for group in groups]
    isempty(normalized) && throw(ArgumentError("layer partition cannot be empty"))
    flattened = collect(Iterators.flatten(normalized))
    flattened == collect(1:Int(source_layers)) || throw(ArgumentError(
        "layer groups must cover 1:$source_layers exactly once in storage order"))
    return LayerPartition{typeof(normalized)}(normalized, Int(source_layers))
end

Base.length(partition::LayerPartition) = length(partition.groups)

function partition_by_count(source_layers::Integer, target_layers::Integer)
    1 <= target_layers <= source_layers || throw(ArgumentError(
        "target layer count must satisfy 1 ≤ target ≤ source"))
    groups = [begin
        first_layer = fld((target - 1) * source_layers, target_layers) + 1
        last_layer = fld(target * source_layers, target_layers)
        first_layer:last_layer
    end for target in 1:target_layers]
    return LayerPartition(groups, source_layers)
end

"""
    partition_at_interfaces(coordinate, target_interfaces; atol=0, rtol=sqrt(eps()))

Build a partition whose target interfaces are a subset of existing source interfaces.
This initial conservative contract does not split source layers.
"""
function partition_at_interfaces(coordinate::PressureCoordinate,
                                 target_interfaces;
                                 atol=zero(eltype(pressure_interfaces(coordinate))),
                                 rtol=sqrt(eps(eltype(pressure_interfaces(coordinate)))))
    source = pressure_interfaces(coordinate)
    indices = map(target_interfaces) do target
        index = findfirst(source) do value
            isapprox(value, target; atol=atol, rtol=rtol)
        end
        index === nothing && throw(ArgumentError(
            "target interface $target does not coincide with a source interface"))
        index
    end
    first(indices) == 1 && last(indices) == length(source) || throw(ArgumentError(
        "target interfaces must retain the full source column"))
    all(>(0), diff(indices)) || throw(ArgumentError(
        "target interfaces must follow source storage order"))
    groups = [indices[i]:(indices[i + 1] - 1) for i in 1:(length(indices) - 1)]
    return LayerPartition(groups, nlayers(coordinate))
end

"""Strategy for merging physical atmospheric state before computing diagnostics."""
abstract type StateReduction end

"""
    DryAirMassWeightedState()

Conservative for dry-air mass, dry-basis constituent amount, layer volume, aerosol
particle number, and water mass reconstructed from `(dry_air_mass, q)`. Temperature is a
dry-air-mass-weighted representative value; it is not claimed to conserve moist enthalpy.
"""
struct DryAirMassWeightedState <: StateReduction end

_sum_groups(values, partition) = [sum(@view values[group]) for group in partition.groups]

function _weighted_groups(values, weights, partition)
    return [sum(@view(values[group]) .* @view(weights[group])) / sum(@view weights[group])
            for group in partition.groups]
end

function _merge_constituent(field::ConstituentField, weights, volumes, partition)
    values = constituent_values(field)
    basis = concentration_basis(field)
    merged = if basis isa Union{LayerMass, LayerMoles, LayerParticleNumber}
        ndims(values) == 1 ? _sum_groups(values, partition) :
            hcat([sum(values[:, group]; dims=2) for group in partition.groups]...)
    elseif basis isa DryMoleFraction
        ndims(values) == 1 ? _weighted_groups(values, weights, partition) :
            hcat([sum(values[:, group] .* reshape(weights[group], 1, :); dims=2) ./
                  sum(weights[group]) for group in partition.groups]...)
    elseif basis isa NumberDensity || basis isa MassDensity
        isavailable(volumes) || throw(ArgumentError(
            "merging $(typeof(basis)) requires layer volume"))
        ndims(values) == 1 ?
            [sum(values[group] .* volumes[group]) / sum(volumes[group])
             for group in partition.groups] :
            hcat([sum(values[:, group] .* reshape(volumes[group], 1, :); dims=2) ./
                  sum(volumes[group]) for group in partition.groups]...)
    else
        throw(ArgumentError(
            "no physically unambiguous merge rule is defined for $(typeof(basis))"))
    end
    return ConstituentField(merged, basis; metadata=field.metadata)
end

function _merge_aerosols(::NoAerosols, weights, volumes, partition)
    return NoAerosols()
end

function _merge_aerosols(state::SectionalAerosolState, weights, volumes, partition)
    number = _merge_constituent(state.number, weights, volumes, partition)
    names = keys(state.components)
    values = map(names) do name
        _merge_constituent(getproperty(state.components, name), weights, volumes, partition)
    end
    components = NamedTuple{names}(Tuple(values))
    return SectionalAerosolState(state.size_grid, number, components;
                                 metadata=state.metadata)
end

function _dry_air_weights(column)
    isavailable(dry_air_mass(column)) && return dry_air_mass(column)
    isavailable(dry_air_column_moles(column)) && return dry_air_column_moles(column)
    throw(ArgumentError(
        "state-space layer merging requires dry-air mass [kg] or " *
        "dry-air column moles [mol m⁻²]"))
end

"""
    merge_column(column, partition, DryAirMassWeightedState())

Merge a column in physical-state space. Use [`merge_optical_depth`](@ref) instead when
optics have already been evaluated on the native layers.
"""
function merge_column(column,
                      partition::LayerPartition,
                      ::DryAirMassWeightedState=DryAirMassWeightedState())
    validate_column(column)
    partition.source_layers == nlayers(column) || throw(DimensionMismatch(
        "partition source layer count does not match the column"))
    weights = _dry_air_weights(column)
    volumes = layer_volume(column)

    merged_dry_mass = isavailable(dry_air_mass(column)) ?
        _sum_groups(dry_air_mass(column), partition) : Unavailable()
    merged_dry_column_moles = isavailable(dry_air_column_moles(column)) ?
        _sum_groups(dry_air_column_moles(column), partition) : Unavailable()
    merged_volume = isavailable(volumes) ? _sum_groups(volumes, partition) : Unavailable()
    merged_temperature = _weighted_groups(temperature(column), weights, partition)

    # q = water_mass / moist_air_mass and dry_mass = (1-q) * moist_air_mass.
    moist_mass = weights ./ (1 .- specific_humidity(column))
    water_mass = specific_humidity(column) .* moist_mass
    merged_moist_mass = _sum_groups(moist_mass, partition)
    merged_water_mass = _sum_groups(water_mass, partition)
    merged_humidity = merged_water_mass ./ merged_moist_mass

    gas_names = keys(trace_gases(column))
    gas_values = map(gas_names) do name
        _merge_constituent(getproperty(trace_gases(column), name),
                           weights, volumes, partition)
    end
    merged_gases = NamedTuple{gas_names}(Tuple(gas_values))
    merged_aerosols = _merge_aerosols(aerosols(column), weights, volumes, partition)

    source_interfaces = pressure_interfaces(column)
    interface_indices = [first(group) for group in partition.groups]
    push!(interface_indices, last(last(partition.groups)) + 1)
    merged_pressure = PressureCoordinate(source_interfaces[interface_indices];
        orientation=vertical_orientation(column))

    reduction_metadata = (method=:dry_air_mass_weighted,
                          source_layers=partition.source_layers,
                          groups=partition.groups)
    source_metadata = metadata(column)
    merged_metadata = source_metadata isa NamedTuple ?
        merge(source_metadata, (; reduction=reduction_metadata)) :
        (; source_metadata, reduction=reduction_metadata)
    return AtmosphericColumn(
        pressure=merged_pressure,
        temperature=merged_temperature,
        specific_humidity=merged_humidity,
        dry_air_mass=merged_dry_mass,
        dry_air_column_moles=merged_dry_column_moles,
        layer_volume=merged_volume,
        trace_gases=merged_gases,
        aerosols=merged_aerosols,
        metadata=merged_metadata,
    )
end
