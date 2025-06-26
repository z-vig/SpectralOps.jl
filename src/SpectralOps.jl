# SpectralOps.jl

module SpectralOps

include("utils/round_to_odd.jl")
include("utils/separate_tuple_matrix.jl")
include("utils/linearfit.jl")
include("utils/nanmean.jl")
include("utils/findwvl.jl")
include("utils/linear_interpolation.jl")
include("utils/make3d.jl")
include("utils/polyfit.jl")
export round_to_odd, separate_tuple_matrix, linearfit, nanmean,
       findλ, linear_interpolation, make3d, polyfit, vandermonde

include("smoothing/moving_average.jl")
include("smoothing/outlier_removal.jl")
export remove_outliers, moving_average

include("continuum_removal/lunar_double_line.jl")
include("continuum_removal/single_line.jl")
export lunar_double_line, single_line

include("spectrum.jl")
export Spectrum

include("spectral_cube.jl")
export SpectralCube, save_spectral_cube, read_spectral_cube

include("band_parameters/calculate_area.jl")
include("band_parameters/calculate_asymmetry.jl")
include("band_parameters/calculate_center.jl")
include("band_parameters/calculate_depth.jl")
export calculate_area, calculate_asymmetry, calculate_center, calculate_depth

include("band_parameters/absorption_band_stats.jl")
export AbsorptionBandStats

end  # SpectralOps