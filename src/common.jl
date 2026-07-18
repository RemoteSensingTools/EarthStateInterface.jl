"""
    Unavailable()

Typed sentinel for an optional physical field that is not present in a source dataset.
Unlike `nothing`, it states that the quantity belongs to the contract but is unavailable.
"""
struct Unavailable end

isavailable(::Unavailable) = false
isavailable(::Any) = true

_require_available(value, name) = isavailable(value) ? value :
    throw(ArgumentError("$name is required for this operation but is unavailable"))

_allfinite(values) = all(isfinite, values)
