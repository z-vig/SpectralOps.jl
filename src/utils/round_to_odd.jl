# utils/round_to_odd.jl

"
    round_to_odd(x::Real) :: Int

Round a real number `x` to the nearest odd integer. If the rounded integer
is even, adjust it by 1 in the direction of `x` to get the nearest odd number.
"
function round_to_odd(x::Real) :: Int
    r = round(Int, x)  # Round to the nearest integer
    return isodd(r) ? r : r + sign(x - r)  # Adjust if even
end