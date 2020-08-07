[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/milankl/Elefridge.jl.svg?branch=master)](https://travis-ci.com/milankl/Elefridge.jl.svg)
[![AppVeyor](https://img.shields.io/appveyor/ci/milankl/Sherlogs-jl?label=Windows&logo=appveyor&logoColor=white)](https://ci.appveyor.com/project/milankl/Sherlogs-jl)


# Elefridge.jl - Compressing atmospheric data into its real information

**Milan Klöwer**\
Atmospheric, Oceanic and Planetary Physics, University of Oxford\
*milan.kloewer@physics.ox.ac.uk*

This repository summarises the results on [ECMWF](https://www.ecmwf.int)'s [summer of weather code](https://esowc.ecmwf.int) [challege #14: Size, precision, speed - pick two](https://github.com/esowc/challenges_2020/issues/4). For comments and changes, please raise an [issue](https://github.com/esowc/Elefridge.jl/issues) or create a pull request. The original prosal is [here](https://github.com/esowc/Elefridge.jl/blob/master/proposal.md).

## Abstract

Enormous amounts of data are produced at weather forecast centres worldwide. 
Compressing large data sets is inevitable to reduce storage requirements and can provide uncertainty estimates. 
Current compression techniques in weather forecast centres do not exploit the spatio-temporal correlation of many geo-physical and geo-chemical variables nor do they only compress the real information contained in floating-point numbers. 
Here, we find alternatives to the 24-bit linear quantisation compression in the Copernicus Atmospheric Monitoring (CAMS) data set and provide a perspective for climate data compression at high compression factors.
Logarithmic quantisation was found to be better suited for the data distributions in CAMS, allowing successful compression at 16-bit per value.
The bitwise information content is calculated, suggesting only 3-10 significant bits contain real information for most variables.
Floating-point quantisation with round-to-nearest can consequently set random bits that do not contain information to 0, making the data set suitable for available lossless compression algorithms.
The entire CAMS data set can be compressed into its real information with this technique by a factor of 13, relative to 32-bit floating-point numbers.
Most lossless compression algorithms work on one-dimensional arrays, but making use of the spatial correlation across three dimensions with zfp, an overall compression factor of 26 (1.2 bit per value) for the entire dataset is achieved.
This study provides evidence that climate and weather forecast data archives can be reduced by one to two orders of magnitude in size without losing information.

## 1. Linear and logarithmic quantisation

## 2. Information entropy of quantisation

## 3. Rounding modes

Quantisation of floating-point numbers into a subset of floating-point numbers is achieved via rounding. Several rounding modes for floats have been proposed in the past. The IEEE-754 standard defines the round-to-nearest standard, in which a float `f` is round to the adjacent nearest quantised floats `f0` and `f1`, whichever is nearer in linear space. Special so-called tie-rules apply when `f` is exactly half-way between `f0` and `f1`, in which case the tie-to-even defines a rounding mode in which `f` gets round to the "even" (i.e. ending in a zero bit) float of `f0` and `f1`. 

Alternatives to round-to-nearest have been proposed for data compression. Bit-shaving always sets the rounded bits to `0`, which effectivly rounds every `f` between `f0` and `f1` towards 0. Bit-shaving is results in a bias for data distributions that are not symmetrical around 0. Assuming a uniform distribution of floats between `f0` and `f1` yields that the expected absolute error of bit-shaving is ULP/2 where ULP (unit in the last place) is the distance between between `f0` and `f1`. In contrast, round-to-nearest introduces a rounding error that is ULP/4 in expectation, as individual absolute errors are ULP/2 at most. To reduce the bias introduced by bit-shaving, bit-grooming was proposed, which alternatingly sets the discarded bits to `0` and to `1`. The number `π` is round to 7 significant bits with the different rounding modes as

```julia
julia> pi = Float32(π)
3.1415927f0

julia> a = [pi,round(pi,7),shave(pi,7),set_one(pi,7)]

julia> bitstring.(a,:split)
4-element Array{String,1}:
 "0 10000000 10010010000111111011011"    # pi at full Float32 precision
 "0 10000000 10010010000000000000000"    # round-to-nearest for 7 significant bits
 "0 10000000 10010010000000000000000"    # bit-shaving
 "0 10000000 10010011111111111111111"    # bit-setting
```
In this example bit-shaving and round-to-nearest yield the same result. However, bit-setting introduces an error of ~ULP. To compare the different rounding modes quantitativly, the mean, absolute and decimal error is analysed for different uniform, normal and log-normal distributions in Fig. 3.

![](https://github.com/esowc/Elefridge.jl/blob/master/plots/groom_vs_round.png)
**Figure 3.** Mean, absolute and decimal error for different floating-point rounding modes: round-to-nearest, bit-grooming and bit-shaving. For each statistical distribution, rounding modes were applied to only retain the first 7 signficant bits. From each statistical distribution 10000 samples were drawn 10000 times, which result in the distribution of the error norms as shown. [[Creating script]](https://github.com/esowc/Elefridge.jl/blob/master/wip/groom_vs_round_plot.jl)

The rounding mode round-to-nearest tie-to-even, as initially defined by the IEEE-754 standard, was found to perform best with respect to the error norms regarded here. We therefore do not recommend alternative rounding modes for data compression.

## 4. Bitwise information content of n-dimensional arrays

The bitwise information content of a dataset has to be analysed to determine the number of bits that can be discarded in a rounding mode. For geo-physical and geo-chemical data, we expect the sign and the exponent bits to have a high real information content unless they are not used, e.g. the sign-bit does not contain information in non-negative data. The most significant bits presumably contain information as long as bits are not randomly occurring, which is assumed for the least-significant bits. We calculate the real bitwise information content for a dataset `A` based on unconditional and conditional entropies for a given bit in all values of `A`. All those bits form a bitstream `bi`, for which the information content `Ic` is calculated as
```
Ic(bi) = H - q0*H0 - q1*H1
```
with `H` being the unconditional entropy, `q0,q1` the probability of a bit being `0,1` and `H0,H1` are the conditional entropies. 
`H0` is the entropy calculated from the conditional probabilities that a bit is `0` or `1` given that the previous bit is `0`. 
Similarly for `H1`. 
Although the entropy `H` is 1 for random uniformly distributed bits (i.e. p(bi=`0`) = 0.5) the conditional probabilities p(0|0), p(1|0), p(0|1), p(1|1) are 0.5 too, such that the conditional entropy is high, reducing the information content to 0. 
In other words, knowing the state of a bit does not provide any further information to the state of the succeeding bit.
For correlated data, in contrast, the conditional entropy reduces (as the conditional probabilities are less uniform) increasing the information content. 
Bits with low information content are therefore either largely unused or independently distributed.

The information content calculation is repeated for every bit in a floating-point number across all elements in a 1-dimensional array `A`. 
For n-dimensional arrays, the conditional probabilities can be calculated in n directions by permuting the dimensions of `A` before unravelling into an 1-dimensional array. 
Summing the n information contents for n-dimensional arrays is the generalisation in which a bit's information can have predictive skill in any of the n dimensions. For a 3D-array `A` with dimensions (x,y,z) the information content is
```
Ic_xyz(A) = Ic_x + Ic_y + Ic_z = 3H - q0 * (H0x + H0y + H0z) - q1 * (H1x + H0y + H0z) 
```
where the subscript `x,y,z` denotes that the array `A` was first unravelled along that dimension. 
We normalise the n-dimensional information content by `1/n` to have a the maximum information content of 1 bit, meaning that this bit contains full information in all 3 dimensions.
To avoid a simulatenous bitflip of all exponent bits around 1 due to the biased-exponent formulation of floating-point numbers, we reinterpret the exponent bits in the sign-and-magnitude formulation. The first exponent bit is consequently the sign of the exponent, the only exponent bit flipping around 1. For the CAMS dataset this makes little difference as most variables are within the range [0,1).

![](https://github.com/esowc/Elefridge.jl/blob/master/plots/bitinformation_all.png)

**Figure 4.** Bitwise information content for all variables in the CAMS data set encoded as Float32. 
Bits that do not contain real information are grey-shaded. 
The total information is the sum of the real information bits.

Most variables in the CAMS dataset do not use the sign bit, nor the sign bit of the exponent as their information is 0. 
Exceptions are the variables derived from the wind velocities, divergence d, etadot, vorticity vo and vertical velocity w. 
Other exponent bits usually have a high information content as they are slowly varying throughout space. 
The information drops quickly to zero beyond the first significant bits and in most cases only the first 3-10 significant bits contain real information.
For some variables information re-emerges for the least significant bits, which is presumably caused by some unphysical quantisation artifacts in the underlying equations.
The total information per value, which is the sum of the information in the real information bits, rarely exceeds 7 bit.
Some variables like CO, CO2, CH4 (including its variants ch4_c, kch4) and temprature have a high share of information stored in the significant bits.

The number of significant bits that contain real information can be used to inform the compression algorithm about the required precision.

## 5. Rounding combined with lossless compression



![](https://github.com/esowc/Elefridge.jl/blob/master/plots/linlogroundzfp_all.png)

**Figure 5.** Compression factor versus absolute and decimal error for linear and logarithmic quantisation, round+lossless and zfp compression. 
Every circle represents the 90th percentile of the respective error norms for one variable. 
The geometric mean of compression factors over all variables is given as horizontal lines.


## 6. 2-4D array floating-point compression

![](https://github.com/esowc/Elefridge.jl/blob/master/maps/o3/round_o3_85.png)
![](https://github.com/esowc/Elefridge.jl/blob/master/maps/o3/zfp_precision3d_o3_85.png)


