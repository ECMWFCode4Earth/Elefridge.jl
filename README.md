[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/milankl/Elefridge.jl.svg?branch=master)](https://travis-ci.com/milankl/Elefridge.jl.svg)
[![AppVeyor](https://img.shields.io/appveyor/ci/milankl/Sherlogs-jl?label=Windows&logo=appveyor&logoColor=white)](https://ci.appveyor.com/project/milankl/Sherlogs-jl)


# Elefridge.jl - Compressing atmospheric data into its real information

**Milan Klöwer**\
Atmospheric, Oceanic and Planetary Physics, University of Oxford\
*milan.kloewer@physics.ox.ac.uk*

This repository summarises the results on [ECMWF](https://www.ecmwf.int)'s [summer of weather code](https://esowc.ecmwf.int) [challege #14: Size, precision, speed - pick two](https://github.com/esowc/challenges_2020/issues/4). For comments and changes, please raise an [issue](https://github.com/esowc/Elefridge.jl/issues) or create a pull request. The original prosal is [here](https://github.com/esowc/Elefridge.jl/blob/master/proposal.md)

## Abstract


## 1. Linear and logarithmic quantisation

[]

## 2. Information entropy of quantisation

## 3. Rounding modes

Quantisation of floating-point numbers into a subset of floating-point numbers is achieved via rounding. Several rounding modes for floats have been proposed in the past. The IEEE-754 standard defines the round-to-nearest standard, in which a float `f` is round to the adjacent nearest quantised floats `f0` and `f1`, whichever is nearer in linear space. Special so-called tie-rules apply when `f` is exactly half-way between `f0` and `f1`, in which case the tie-to-even defines a rounding mode in which `f` gets round to the "even" (i.e. ending in a zero bit) float of `f0` and `f1`. 

```julia
julia> pi = Float32(π)
3.1415927f0

julia> bitstring(pi,:split)
"0 10000000 10010010000111111011011"

julia> bitstring(round(pi,7),:split)
"0 10000000 10010010000000000000000"
```

Alternatives to round-to-nearest have been proposed for data compression. Bit-shaving always sets the rounded bits to `0`, which effectivly rounds every `f` between `f0` and `f1` towards 0. Bit-shaving is results in a bias for data distributions that are not symmetrical around 0. Assuming a uniform distribution of floats between `f0` and `f1` yields that the expected absolute error of bit-shaving is ULP/2 where ULP (unit in the last place) is the distance between between `f0` and `f1`. In contrast, round-to-nearest introduces a rounding error that is ULP/4 in expectation, as individual absolute errors are ULP/2 at most. To reduce the bias introduced by bit-shaving, bit-grooming was proposed, which alternatingly sets the discarded bits to `0` and to `1`. The number `pi = Float32(π)` is round to 7 significant bits by setting the 16 least significant bits to 0 ("bit-shaving") or to 1 ("bit-setting") or 

```julia
julia> a = [pi,round(pi,7),shave(pi,7),set_one(pi,7)]

julia> bitstring.(a,:split)
4-element Array{String,1}:
 "0 10000000 10010010000111111011011"    # pi at full Float32 precision
 "0 10000000 10010010000000000000000"    # round-to-nearest for 7 significant bits
 "0 10000000 10010010000000000000000"    # bit-shaving
 "0 10000000 10010011111111111111111"    # bit-setting
```
In this example bit-shaving and round-to-nearest yield the same result. However, bit-setting introduces an error of ~ULP.

![](https://github.com/esowc/Elefridge.jl/blob/master/plots/groom_vs_round.png)
**Figure 3.** Mean, absolute and decimal error for different floating-point rounding modes: round-to-nearest, bit-grooming and bit-shaving. For each statistical distribution, rounding modes were applied to only retain the first 7 signficant bits. From each statistical distribution 10000 samples were drawn 10000 times, which result in the distribution of the error norms as shown. [[Creating script]](https://github.com/esowc/Elefridge.jl/blob/master/wip/groom_vs_round.jl)

The rounding mode round-to-nearest tie-to-even, as initially defined by the IEEE-754 standard, was found to perform best with respect to the error norms regarded here. We therefore do not recommend alternative rounding modes for data compression.

## 4. Bitwise information content of n-dimensional arrays

## 5. Rounding combined with lossless compression

## 6. 2-4D array floating-point compression
