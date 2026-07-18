"""Route chosen when reducing vertical resolution for an optical calculation."""
abstract type OpticalReductionRoute end

"""
    StateThenOptics()

Merge thermodynamic and composition state, then evaluate cross sections or particle
optics on the merged layers. Faster, but generally approximate because spectroscopy and
Mie optics are nonlinear in temperature, pressure, humidity, size, and composition.
"""
struct StateThenOptics <: OpticalReductionRoute end

"""
    OpticsThenMerge()

Evaluate native-layer optics first, then add optical depths into target layers. This
preserves monochromatic absorption/extinction through the merge. Effective cross
sections or phase properties may be diagnosed afterward with the appropriate weights.
"""
struct OpticsThenMerge <: OpticalReductionRoute end

function _check_optical_layers(values, partition)
    size(values, 1) == partition.source_layers || throw(DimensionMismatch(
        "the first dimension of optical data must be the source-layer dimension"))
end

"""
    merge_optical_depth(optical_depth, partition)

Add native-layer optical depth into each target layer. The first dimension is layers;
any trailing dimensions (for example wavelength or g-point) are preserved.
"""
function merge_optical_depth(optical_depth::AbstractArray,
                             partition::LayerPartition)
    _check_optical_layers(optical_depth, partition)
    output_size = (length(partition), Base.tail(size(optical_depth))...)
    merged = similar(optical_depth, output_size)
    trailing_indices = ntuple(_ -> Colon(), ndims(optical_depth) - 1)
    for (target, group) in enumerate(partition.groups)
        source = view(optical_depth, group, trailing_indices...)
        selectdim(merged, 1, target) .= dropdims(sum(source; dims=1); dims=1)
    end
    return merged
end

function merge_optical_depth(optical_depth::AbstractVector,
                             partition::LayerPartition)
    length(optical_depth) == partition.source_layers || throw(DimensionMismatch(
        "optical-depth vector must contain one value per source layer"))
    return _sum_groups(optical_depth, partition)
end

"""
    effective_cross_section(cross_section, absorber_amount, partition)

Compute a target-layer cross section that exactly reproduces the summed native-layer
optical depth:

`σ_eff = sum(σᵢ Nᵢ) / sum(Nᵢ)`.

The first dimension of `cross_section` is layers. `absorber_amount` is one extensive
amount per native layer and may be molecules, moles, or any consistent column amount.
"""
function effective_cross_section(cross_section::AbstractVector,
                                 absorber_amount::AbstractVector,
                                 partition::LayerPartition)
    length(cross_section) == partition.source_layers || throw(DimensionMismatch(
        "cross section must contain one value per source layer"))
    length(absorber_amount) == partition.source_layers || throw(DimensionMismatch(
        "absorber amount must contain one value per source layer"))
    return [begin
        total_amount = sum(absorber_amount[group])
        total_amount > zero(total_amount) || throw(ArgumentError(
            "effective cross section requires positive absorber amount in every target layer"))
        sum(cross_section[group] .* absorber_amount[group]) / total_amount
    end for group in partition.groups]
end

function effective_cross_section(cross_section::AbstractMatrix,
                                 absorber_amount::AbstractVector,
                                 partition::LayerPartition)
    _check_optical_layers(cross_section, partition)
    length(absorber_amount) == partition.source_layers || throw(DimensionMismatch(
        "absorber amount must contain one value per source layer"))
    merged = similar(cross_section, length(partition), size(cross_section, 2))
    for (target, group) in enumerate(partition.groups)
        weights = absorber_amount[group]
        total_amount = sum(weights)
        total_amount > zero(total_amount) || throw(ArgumentError(
            "effective cross section requires positive absorber amount in every target layer"))
        merged[target, :] .= vec(sum(cross_section[group, :] .*
                                     reshape(weights, :, 1); dims=1)) ./
                             total_amount
    end
    return merged
end
