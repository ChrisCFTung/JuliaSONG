using Distributions, Random
import Base.rand

struct StandardCandle{T<:Real} <: Distribution{Univariate, Continuous}
    value::T
    function StandardCandle{T}(value::T) where {T<:Real}
        return new{T}(value)
    end
end

function StandardCandle(value::Float64)
    return StandardCandle{Float64}(value)
end

function rand(d::StandardCandle, s::Int64)
    return d.value .* ones(s)
end

#= 
LogNormal
=#

function LogNormLumin(L0::Real, sigma::Real)
    mu = log(L0) - sigma^2/2
    return LogNormal(mu, sigma)
end

#= 
PowerLaw Distribution in the type of Distributions
=#

struct PowerLaw{T<:Real} <: Distribution{Univariate, Continuous}
    index::T
    x_mean::T
    width::T
    function PowerLaw{T}(index::T, x_mean::T, width::T; check_args=true) where {T<:Real}
        check_args && Distributions.@check_args(PowerLaw, index>1)
        return new{T}(index, x_mean, width)
    end
end

function PowerLaw(index::Float64, x_mean::Float64, width::Float64; check_args=true)
    return PowerLaw{Float64}(index, x_mean, width; check_args=check_args)
end
PowerLaw(index::Real, x_mean::Real, width::Real) = PowerLaw(promote(index, x_mean, width)...)

param(d::PowerLaw) = (d.index, d.x_mean, d.width)
#1 rand
function rand(d::PowerLaw, s::Int64)
    index, x_mean, width = param(d)
    beta = 1-index
    if index == 2
        x_min = x_mean / width / log(10)
    else
        x_min = (beta+1)/beta * x_mean * (10^(width*beta)-1) / (10^(width*(beta+1))-1)
    end
    x_max = x_min*10^width 
    u = Base.rand(s)
    return (x_min^beta .+ (x_max^beta-x_min^beta) .* u).^(1/beta)
end
