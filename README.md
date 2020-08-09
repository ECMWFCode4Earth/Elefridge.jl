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

# Functionality of Elefridge.jl

### Linear quantisation

Linear quantisation of n-dimensional arrays (any number format that can be converted to `Float64` is supported, including `Float32, Float16, BFloat16`) into 8, 16 or 24 bit is achieved via

```julia
julia> A = rand(Float32,1000);

julia> L = LinQuant8Array(A)
1000-element LinQuantArray{UInt8,1}:
 0xc2
 0x19
 0x3e
 0x5b
    ⋮
```
and similarly with `LinQuant16Array, LinQuant24Array`. Decompression via
```julia
julia> Array(L)
1000-element Array{Float32,1}:
 0.76074356
 0.09858093
 0.24355145
 0.357177
    ⋮
```
`Array{T}()` optionally takes a type parameter `T` such that decompression to other number formats than the default `Float32` is possible (e.g. `Float64, BFloat16`).

### Logarithmic quantisation

In a similar way, `LogQuant8Array, LogQuant16Array, LogQuant24Array` compresses an n-dimensional array (non-negative elements only) via logarithmic quantisation.
```julia
julia> A = rand(Float32,100,100);

julia> A[1,1] = 0;

julia> L = LogQuant16Array(A)
100×100 LogQuantArray{UInt16,2}:
 0x0000  0xf22d  0xfdf6  0xf3e8  0xf775  …  
 0xe3dc  0xfdc0  0xedb5  0xed47  0xee5b     
 0xde3d  0xbe58  0xb541  0xf573  0x9885     
 0xf38b  0xfefe  0xea2f  0xfbb6  0xf0d2     
 0xd0d2  0xfe1f  0xff60  0xf6cd  0xec26        
 0xffa6  0xe621  0xf14d  0xfb2c  0xf50c  …  
 0xfcb7  0xe6fb  0xf237  0xecd5  0xfb0a     
 0xe4ed  0xf86f  0xf83d  0xff86  0xb686     
      ⋮                                  ⋱ 
```
Exception occurs for 0, which is mapped to `0x0`.
`Ox1` to `0xff, 0xffff, 0xffffff` are then the available bitpatterns to encode the range from `minimum(A)` to `maximum(A)` logarithmically.
Decompression as with linear quantisation via the `Array{T}()` function.

### Rounding modes

Elefridge.jl implements 4 rounding modes for floating-point numbers, `round` (which is round-to-nearest tie-to-even), `shave`, `halfshave`, `set_one` and `groom`. All rounding functions take an integer as the number of keepbits as a second argument and operate on scalars as well as arrays (`Float32` and `Float64` are supported).
```julia
julia> bitstring.(A,:split)             # bitwise representation (split in sign, exp, sig bits) of some random numbers
5-element Array{String,1}:
 "0 01111101 01001000111110101001000"
 "0 01111110 01010000000101001110110"
 "0 01111110 01011101110110001000110"
 "0 01111101 00010101010111011100000"
 "0 01111001 11110000000000000000101"

julia> bitstring.(round(A,3),:split)
5-element Array{String,1}:
 "0 01111101 01000000000000000000000"
 "0 01111110 01100000000000000000000"  # correct round-to-nearest by flipping the third significant bit
 "0 01111110 01100000000000000000000"  # same here
 "0 01111101 00100000000000000000000"  # and here
 "0 01111010 00000000000000000000000"  # note how the carry bits correctly carries into the exponent

julia> bitstring.(shave(A,3),:split)
5-element Array{String,1}:
 "0 01111101 01000000000000000000000"  # identical to round here
 "0 01111110 01000000000000000000000"
 "0 01111110 01000000000000000000000"
 "0 01111101 00000000000000000000000"
 "0 01111001 11100000000000000000000"

julia> bitstring.(set_one(A,3),:split)
5-element Array{String,1}:
 "0 01111101 01011111111111111111111"
 "0 01111110 01011111111111111111111"
 "0 01111110 01011111111111111111111"
 "0 01111101 00011111111111111111111"
 "0 01111001 11111111111111111111111"

julia> bitstring.(groom(A,3),:split)
5-element Array{String,1}:
 "0 01111101 01000000000000000000000"   # shave
 "0 01111110 01011111111111111111111"   # set to one
 "0 01111110 01000000000000000000000"   # shave
 "0 01111101 00011111111111111111111"   # etc.
 "0 01111001 11100000000000000000000"

julia> bitstring.(halfshave(A,3),:split)
5-element Array{String,1}:
 "0 01111101 01010000000000000000000"   # set all discarded bits to 1000...
 "0 01111110 01010000000000000000000"
 "0 01111110 01010000000000000000000"
 "0 01111101 00010000000000000000000"
 "0 01111001 11110000000000000000000"
```
### Bitpattern entropy

The bitpattern entropy (i.e. a measure for the effective use of a given bit-encoding/quantisation) can be calculated via the `bitentropy` function for n-dimensional arrays with types that fit into `8,16,24,32,40,48,56` or `64` bit. Unsigned integers of the respective size are automatically implemented via the [BitIntegers.jl](https://github.com/rfourquet/BitIntegers.jl) package. Consequently, all common number formats like `Float16/32/64` or even `BFloat16` or `Posit8/16/32` etc. are supported. The `bitentropy(A::Array)` function's bottleneck is the sorting function, which sorts the elements in `A` in the beginning by its corresponding bitpatterns. The entropy is by default calculated in units of bits, which can be changed with a second argument `::Real` if desired.
```julia
julia> A = rand(Float32,100000000);

julia> bitentropy(A)
22.938590744784577
```
Here, the entropy is about 23 bit, meaning that `9` bits are effectively unused.

### Information content

To calculate the information content of an n-dimensional array (any typ `T` is supported that can be reinterpreted as `8,16,24,32,40,48,56` or `64-bit` unsigned integer) the following functions are exported:

**bitcount**. The function `bitcount(A::Array)` counts all occurences of the 1-bit in every bit-position in every element of `A`. E.g.

```julia
julia> bitstring.(A)
5-element Array{String,1}:
 "10001111"
 "00010111"
 "11101000"
 "10100100"
 "11101011"

julia> bitcount(A)
8-element Array{Int64,1}:
 4
 2
 3
 1
 3
 3
 3
 3
 ```
The first bit of elements (here: `UInt8`) in `A` are 4 times `1` and 1 times `0`, etc. In contrast, elements drawn from a uniform distribution U(0,1)
 ```julia
julia> A = rand(Float32,100000);

julia> bitcount(A)
32-element Array{Int64,1}:
      0
      0
 100000
 100000
      ⋮
  37411
  25182
      0
```
have never a sign bit that is `0`, but the 2nd and third exponent bit is always `1`.

**bitcountentropy**. The `bitcountentropy(A::Array)` calculates the entropy of bit occurences in `A`. For random bits occuring at probabilities p(0) = 0.5, p(1) = 0.5 the entropy for every bit is maximised to 1 bit:
```julia
julia> A = rand(UInt8,100000);

julia> Elefridge.bitcountentropy(A)
8-element Array{Float64,1}:
 0.9999998727542938
 0.9999952725717266
 0.9999949724904816
 0.9999973408228667
 0.9999937649515901
 0.999992796900212
 0.9999970566115759
 0.9999998958374157
 ```
The converges to 1 for larger arrays.
 
**bitpaircount**. 

### Zfp Compression

Julia bindings to th [zfp compression library](https://computing.llnl.gov/projects/floating-point-compression/zfp-compression-ratio-and-quality) have been developed. This functionality is exported to a separate package: [ZfpCompression.jl](https://github.com/milankl/ZfpCompression.jl) and documentation can be found therein.

## Installation

Not yet registered in the Julia registry, hence do
```julia
(@v1.5) pkg> add https://github.com/esowc/Elefridge.jl
```
in the package manager (access via `]`).
