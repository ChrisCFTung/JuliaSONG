module JuliaSONG
using QuadGK, Interpolations
using Distributions
using Random
using DataFrames

include("Luminosity.jl")

#=
Evolution Models
=#
function SFR(z::Real)
    a = 0.015
    b = 2.7
    c = 2.9
    d = 5.6
    x = (1+z)
    density = a*(x^b) / (1.0 + (x/c)^d) 
    return density
end

function NoEvo(z::Real)
    return 1
end

#=
Parameters for the simulation
=#
# lambdaCDM parameter
struct Cosmology
    c::Real
    h::Real
    om0::Real
    ode::Real
    ok0::Real
    D_h::Real
    evol::Function
end

function Cosmology(;om0::Real=0.315, ode::Real=0.685, h::Real=0.674, 
                    evolution::Function=SFR)
    c = 299792458
    ok0 = 1-om0-ode
    D_h = c/1000/100/h
    return Cosmology(c, h, om0, ode, ok0, D_h, evolution)
end

# a default cosmology
cosmo = Cosmology()

#=
functions for calculating distance
=#

function E(z::Real; cosmo_par::Cosmology=cosmo)
    return sqrt(cosmo_par.om0*(1+z)^3+cosmo_par.ok0*(1+z)^2+cosmo_par.ode)
end

function comoving_distance(z::Real; z0::Real=0.0, cosmo_par::Cosmology=cosmo)
    return cosmo_par.D_h*QuadGK.quadgk(x->1/E(x), z0, z)[1]
end

function tranverse_comoving_distance(z::Real; cosmo_par::Cosmology=cosmo)
    if cosmo_par.ok0 == 0
        return comoving_distance(z, cosmo_par=cosmo_par)
    else
        ok0 = abs(cosmo_par.ok0)
        return cosmo_par.D_h/sqrt(ok0)*sinh(sqrt(ok0)/cosmo_par.D_h*
        comoving_distance(z, cosmo_par=cosmo_par))
    end
end

function angular_diameter_distance(z::Real; cosmo_par::Cosmology=cosmo)
    return tranverse_comoving_distance(z, cosmo_par=cosmo_par)/(1.0+z)
end

function luminosity_distance(z::Real; cosmo_par::Cosmology=cosmo)
    return (1.0+z)*tranverse_comoving_distance(z, cosmo_par=cosmo_par)
end

function dL_interpolation(; cosmo_par::Cosmology=cosmo)
    z_bin = 0.0:0.02:10.0
    dL = luminosity_distance.(z_bin, cosmo_par=cosmo_par)
    return Interpolations.LinearInterpolation(z_bin, dL)
end

function comoving_volume(z::Real; cosmo_par::Cosmology=cosmo)
    if cosmo_par.ok0 == 0
        return 4.0*pi/3.0 * comoving_distance(z, cosmo_par=cosmo_par)^3
    else
        D_m = tranverse_comoving_distance(z, cosmo_par=cosmo_par)
        D_ratio = D_m/cosmo_par.D_h
        ok0 = abs(cosmo_par.ok0)
        prefactor = 4.0*pi*cosmo_par.D_h^3/2.0/ok0
        return prefactor*(D_ratio*sqrt(1+ok0*D_ratio^2)
               -asinh(sqrt(ok0)*D_ratio)/sqrt(ok0))
    end
end

function diff_comoving_volume(z::Real; cosmo_par::Cosmology=cosmo)
    return cosmo_par.D_h*(1.0+z)^2*angular_diameter_distance(z, cosmo_par=cosmo_par)^2/
    E(z, cosmo_par=cosmo_par)
end

#=
Source Distribution
=#

# redshift up to which should be cartesian
__zlocal = 0.01 

function RedshiftDistribution(z::Real; cosmo_par::Cosmology=cosmo)
    return 4.0*pi*cosmo_par.evol(z)*diff_comoving_volume(z, cosmo_par=cosmo_par)
end

function RedshiftIntegral(z::Real; cosmo_par::Cosmology=cosmo)
    return QuadGK.quadgk(x -> RedshiftDistribution(x, cosmo_par=cosmo_par), 0, z)[1]
end

function NSources(local_density::Real, z_max::Real; cosmo_par::Cosmology=cosmo)
    vlocal = comoving_volume(__zlocal, cosmo_par=cosmo_par)
    Ntotal = local_density*vlocal/RedshiftIntegral(__zlocal, cosmo_par=cosmo_par)*
            RedshiftIntegral(z_max, cosmo_par=cosmo_par)
    return Ntotal
end

function StandardCandleFlux(flux_norm::Real, density::Real, index::Real;
                               cosmo_par::Cosmology=cosmo, z0=1.0, z_inf=10.0)
    all_sky_flux = flux_norm*4*pi
    Ntotal = NSources(density, z_inf, cosmo_par=cosmo_par)
    norm = RedshiftIntegral(z_inf, cosmo_par=cosmo_par)
    flux_in = QuadGK.quadgk(z->((1+z)/(1+z0))^(2-index)/luminosity_distance(z, cosmo_par=cosmo_par)^2*
                 RedshiftDistribution(z, cosmo_par=cosmo_par)/norm,
                 0, z_inf)[1]
    return all_sky_flux/Ntotal/flux_in/luminosity_distance(z0, cosmo_par=cosmo_par)^2
end

function Flux_to_NLuminosity(flux::Real, index::Real, z::Real; cosmo_par::Cosmology=cosmo)
    dN = flux .* 4pi .* luminosity_distance.(z, cosmo_par=cosmo_par) .^ 2 ./ (1 .+ z) .^ (2.0 - index)
    return dN
end

#=
Sampling
=#
function invCDF(pdf::Function, low::Real, high::Real;
                Nbin::Real=10000.0)
    bins = low:(high-low)/Nbin:high
    pdf_v = pdf.(bins)
    cdf_v = cumsum(pdf_v)
    cdf_v /= last(cdf_v)
    invcdf = LinearInterpolation(cdf_v, bins)
    return invcdf
end

function z_Sampling(nsources::Int; z_max::Real=10.0,
                    cosmo_par::Cosmology=cosmo)
    pdf(z) = RedshiftDistribution(z, cosmo_par=cosmo_par)
    icdf = invCDF(pdf, 0.0, z_max)
    return icdf(rand(nsources))
end

function generate_source(density::Real, flux_norm::Real, index::Real; 
             z_max::Real=10.0, cosmo_par::Cosmology=cosmo, lumi=StandardCandle)
    Ntotal = NSources(density, z_max)
    nsrc = rand(Poisson(Ntotal), 1)[1]
    z = z_Sampling(nsrc, z_max=z_max, cosmo_par=cosmo_par)
    candle_flux = StandardCandleFlux(flux_norm, density, index, 
                                        cosmo_par=cosmo_par)
    L0 = Flux_to_NLuminosity(candle_flux, index, 1.0)
    L = rand(lumi(L0), nsrc)
    dL = dL_interpolation(cosmo_par=cosmo_par)
    flux = L ./ (dL.(z) .^ 2) .* (1 .+ z) .^ (2.0-index)
    return DataFrame(Redshift=z, Flux=flux)
end

end
