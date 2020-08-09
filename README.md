[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/milankl/Elefridge.jl.svg?branch=master)](https://travis-ci.com/milankl/Elefridge.jl.svg)
[![AppVeyor](https://img.shields.io/appveyor/ci/milankl/Sherlogs-jl?label=Windows&logo=appveyor&logoColor=white)](https://ci.appveyor.com/project/milankl/Sherlogs-jl)


# Elefridge.jl - Compressing atmospheric data into its real information

**Milan Klöwer**\
Atmospheric, Oceanic and Planetary Physics, University of Oxford\
*milan.kloewer@physics.ox.ac.uk*

This repository summarises the results on [ECMWF](https://www.ecmwf.int)'s [summer of weather code](https://esowc.ecmwf.int) [challege #14: Size, precision, speed - pick two](https://github.com/esowc/challenges_2020/issues/4). For comments and changes, please raise an [issue](https://github.com/esowc/Elefridge.jl/issues) or create a pull request. The original prosal is [here](https://github.com/esowc/Elefridge.jl/blob/master/proposal.md).

## Abstract

Enormous amounts of data are produced at weather and climate forecast centres
worldwide. Compressing large data sets is inevitable to reduce storage and data
sharing requirements. Current compression techniques in forecast centres do not
exploit the spatio-temporal correlation of many geophysical and geochemical
variables nor do they only compress the real information contained in 32 or
64-bit floating-point numbers. Here, we find alternatives to the default 24-bit
linear quantisation compression in the Copernicus Atmospheric Monitoring (CAMS)
data set and provide a perspective for climate data compression at large
compression factors. Logarithmic quantisation was found to better match the data
distributions in CAMS, allowing for successful compression at 16-bit per value.
The bitwise information content is calculated in the original 32-bit floating-point
number encoding, suggesting only 3-9 significant bits contain real information
for most variables. Floating-point quantisation can consequently set bits that
do not contain information to 0, facilitating available lossless compression
algorithms. The entire CAMS data set can be compressed into its real information
with this technique by a factor of 13, relative to 32-bit floating-point numbers.
Most lossless compression algorithms work on one-dimensional arrays, but making
use of the correlation in three spatial dimensions with zfp, an overall compression
factor of 26 (1.2 bit per value) for the entire dataset is achieved. This study
provides evidence that climate and weather forecast data archives can be reduced
by one to two orders of magnitude in size without losing valuable information.

## 1. Linear and logarithmic quantisation

24-bit linear quantisation is the current default compression method for CAMS data.
To compress an array `A`, the minimum and maximum is obtained
```julia
Amin = minimum(A)
Amax = maximum(A)
```
which allows the calculation of `Δ`, the inverse of the spacing between two
quantums
```julia
Δ = 2^(n-1)/(Amax-Amin)
```
where `n` is the number of bits used for quantisation, 24 in this case. For every
element `a` in `A` the corresponding quantum `q` which is closest in linear space
is calculated via
```julia
q = T(round((a-Amin)*Δ))
```
where `round` is the round-to-nearest function for integers and `T` the conversion
function to `UInt24` (or `UInt8, UInt16` for other choices of `n`). Consequently,
an array of all `q` and `Amin,Amax` have to be stored to allow for decompression,
which is obtained by reversing the conversion from `a` to `q`. Note that the rounding
error is introduced as the `round` function can only be approximately inverted.

Logarithmic quantisation distributes the quantums logarithmically, such that
more bitpatterns are reserved for values close to the minimum and fewer close to
the maximum in A. Logarithmic quantisation can be generalised to negative values
by introducing a sign-bit, however, we limit our application here to non-negative
values. We obtain the minimum and maximum value in `A` as follows
```julia
Alogmin = log(minpos(A))
Alogmax = log(maximum(A))
```
where zeros are ignored in the `minpos` function, which instead returns the smallest
positive value. The inverse spacing `Δ` is then
```julia
Δ = 2^(n-2)/(logmax-logmin)
```
Note, that only `2^(n-1)` (and not 2^n as for linear quantisation) bitpatterns
are used to resolve the range between minimum and maximum, as we want to reserve
the bitpattern 0x000000 for the value 0. The corresponding quantum `q` for `a`
`A` is then, unless `a=0` in which case `q=0x000000`,
```julia
q = T(round(c + Δ*log(a)))+0x1
```
`c` is a constant which can be set as `-Alogmin*Δ` such that we obtain essentially
the same compression function as for linear quantisation, except that every element
`a` in `A` is converted to their logarithm first. However, rounding to nearest
in logarithmic space will therefore be achieved, which is a biased rounding mode,
that has a bias away from zero. We can correct this round-to-nearest in logarithmic
space rounding mode with
```julia
c = 1/2 - Δ*log(minimum(A)*(exp(1/Δ)+1)/2)
```
which yields round-to-nearest in linear space.

The addition of `0x1` in the quantisation is to map the minpos-maximum range to
bitpatterns `0x000001` to `0xffffff` but to keep the `0x000000` free for encoding
0.

## 1a. Information entropy of quantisation

Choosing either linear or logarithmic quantisation imposes a linear or logarithmic
distribution of quantums on the data set. The distance between the data values
and the nearest quantums depend on the data distribution, which in turn affects
the overall rounding error. For a continuous data distribution, the quantisation
should be chosen such that more quantums are in the vicinity of data values
to minimise the distance between them.

Calculating the entropy of quantisation will quantify the amount of bitpatterns
that are effectively used by the quantisation, measured in bits, such that the
difference between the entropy and the available bits is ideally as small as possible.
The information entropy is defined as
```julia
H = - Σᵢ pᵢ * log2(pᵢ)
```
where `pᵢ` the chance of bitpattern `i` occurring in the quantisation. Theoretically,
a uniform distribution achieves maximum entropy with linear quantisation
as every bitpattern is used at the same frequency. Similarly, a log-uniform distribution
achieves maximum entropy with logarithmic quantisation. Measuring the entropy
of either linearly or logarithmically quantised arrays is therefore a way to
quantify whether the underlying data is rather linearly or logarithmically distributed.

![](https://github.com/esowc/Elefridge.jl/blob/master/plots/bitpattern_histogram.png)
**Figure 1.** Bitpattern histogram for (a) 24-bit linear quantisation and
(b) 16-bit logarithmic quantisation for NO2. The information entropy
of quantisation is given in the top-right of each panel.

The bitpattern histogram in Fig. 1a reveals that 24-bit linear quantisation of NO2 does
not use most of the available bitpatterns. Contrarily, high occurrences for bitpatterns
in the vicinity of `0x000000` are observed. The entropy is only about 16-bit,
such that 8-bit are basically unused, which translates to only 2^-8=0.4% of bitpatterns
are effectively used.

In contrast, the logarithmic quantisation reveals a histogram that
makes use of a wide range of available bitpatterns. The entropy is 15bit, such
that effectively only 1-bit is redundant. This shows that a logarithmic distribution
of quantums is much better suited for the variable NO2.

![](https://github.com/esowc/Elefridge.jl/blob/master/plots/entropy_linlog16.png)
**Figure 2.** Effectively used bit patterns when compressing the variables in
CAMS via linear or logarithmic quantisation as measured by the quantisation
entropy at 16 bit.

The quantisation entropy for most variables in the CAMS data set reveals that
almost all variables are rather logarithmically than linearly distributed (Fig. 2).
Variables like temperature, ozone or CO2 are an exception, but still show high
entropies for logarithmic quantisation. Consequently, no single variable in CAMS
would significantly benefit from a linear quantisation over a logarithmic.

## 1b. Error quantification

Although they are related, maximising the entropy does, in general, not guarantee
a minimisation of the rounding error. We therefore quantify the compression errors
of linear and logarithmic quantisation with the following error norms

### Normalised mean error

The mean error of an array A to its quantised array Q is
```julia
mean error = ∑ᵢ (Aᵢ - Qᵢ)
```
which quantifies a rounding bias between `A` and `Q`. The mean error can be normalised
to allow easier comparison between different data sets.
```julia
normalised mean error = ∑ᵢ (Aᵢ - Qᵢ) / ∑ᵢ Aᵢ
```
The normalisation does not change the qualitative results of comparing different
quantisation methods `Q1,Q2,...` as the division by the mean of A is for all identical.

### Normalised absolute error

The normalised absolute error is
```julia
normalised absolute error = abs(Aᵢ - Qᵢ) / ∑ᵢ Aᵢ
```
which measures the average distance of values in `A` to their respective quantums `Qᵢ`.
The normalised absolute error is an array of the same size as `A` and `Q`, such
that its mean, for example is the L1-norm of the linear error, which is invariant under
the addition of a constant to `A` and `Q`.

### Decimal error

The decimal error is a relative error, which is defined by
```julia
decimal error = abs(log10(Aᵢ/Qᵢ))
```
As with relative errors, the decimal error is invariant under multiplication
with a constant. It is therefore not necessary to normalise the decimal error,
as errors of different variables will be comparable by definition.

### Quantisation errors in CAMS

The mean, absolute and decimal error for 16 and 24-bit linear and logarithmic
quantisation are compared for variables in the CAMS data set. As shown by the tails
of the error distributions, logarithmic quantisation puts a much stronger bound
on the decimal error. Linear quantisation errors, however, can reach decimal errors
of 1 and more, meaning that the quantisation introduced an error on the order
of the magnitude of the value. Linear quantisation puts a slightly stronger bound
on the absolute error. The mean error is small for all quantisation methods,
only 16-bit linear quantisation may pose an intolerable error on the mean.

Comparing the quantisation for all variables, we conclude
that logarithmic quantisation into 16-bit can replace the current 24-bit
linear quantisation method safely, in even reduce the error for many variables.
Due to the 16-bit word length, the entire CAMS dataset can therefore be archived
at 67% of the current archive size.

![](https://github.com/esowc/Elefridge.jl/blob/master/plots/linvslog_all.png)
**Figure 3.** Error comparison for linear and logarithmic quantisation for variables
in the CAMS dataset. Variables are sorted by the absolute error of LinQuant16.

## 3. Rounding modes

Quantisation of floating-point numbers into a subset of floating-point numbers
is achieved via rounding. Several rounding modes for floats have been proposed
in the past. The IEEE-754 standard defines the round-to-nearest standard, in which
a float `f` is round to the adjacent nearest quantised floats `f0` and `f1`,
whichever is nearer in linear space. Special so-called tie-rules apply when
`f` is exactly half-way between `f0` and `f1`, in which case the tie-to-even
defines a rounding mode in which `f` gets round to the "even" (i.e. ending in
a zero bit) float of `f0` and `f1`.

Alternatives to round-to-nearest have been proposed for data compression.
Bit-shaving always sets the rounded bits to `0`, which effectively rounds every
`f` between `f0` and `f1` towards 0. Bit-shaving is results in a bias for data
distributions that are not symmetrical around 0. Assuming a uniform distribution
of floats between `f0` and `f1` yields that the expected absolute error of bit-shaving
is ULP/2 where ULP (unit in the last place) is the distance between between `f0`
and `f1`. In contrast, round-to-nearest introduces a rounding error that is ULP/4
in expectation, as individual absolute errors are ULP/2 at most. To reduce the
bias introduced by bit-shaving, bit-grooming was proposed, which alternatingly
sets the discarded bits to `0` and to `1`. The number `π` is round to 7 significant
bits with the different rounding modes as

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
In this example bit-shaving and round-to-nearest yield the same result. However,
bit-setting introduces an error of ~ULP. To compare the different rounding modes
quantitatively, the mean, absolute and decimal error is analysed for different
uniform, normal and log-normal distributions in Fig. 3.

![](https://github.com/esowc/Elefridge.jl/blob/master/plots/groom_vs_round.png)
**Figure 3.** Mean, absolute and decimal error for different floating-point rounding modes: round-to-nearest, bit-grooming and bit-shaving. For each statistical distribution, rounding modes were applied to only retain the first 7 significant bits. From each statistical distribution 10000 samples were drawn 10000 times, which result in the distribution of the error norms as shown. [[Creating script]](https://github.com/esowc/Elefridge.jl/blob/master/wip/groom_vs_round_plot.jl)

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

**Figure 5.** Compression factor versus absolute and decimal error for linear and
logarithmic quantisation, round+lossless and zfp compression. Every circle represents
the 90th percentile of the respective error norms for one variable. Lossless compression
is the LZ4HC algorithm (level 5). The geometric mean of compression factors over
all variables is given as horizontal lines.


## 6. 2-4D array floating-point compression

![](https://github.com/esowc/Elefridge.jl/blob/master/maps/o3/round_o3_85.png)
![](https://github.com/esowc/Elefridge.jl/blob/master/maps/o3/zfp_precision3d_o3_85.png)
**Figure 6.** Compression of ozone (O3) at different levels of precision with
(a) round+lossless and (b) zfp compression. The precision decreases from left to
right with increasing compression factors. Only one vertical level in the high
troposphere (model level 85) is shown, but compression factors include all
vertical levels.

# Conclusion: A roadmap for atmospheric data compression

## 1. Short-term: Logarithmic quantisation

## 2. Medium-term: Round+lossless

## 3. Long-term: Zfp multi-dimensional compression


# Functionality of Elefridge.jl

### Linear quantisation

Linear quantisation of n-dimensional arrays (any number format that can be
converted to `Float64` is supported, including `Float32, Float16, BFloat16`)
into 8, 16 or 24 bit is achieved via

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
`Array{T}()` optionally takes a type parameter `T` such that decompression to
other number formats than the default `Float32` is possible (e.g. `Float64, BFloat16`).
The compression involves a rounding error which follows a round-to-nearest in
linear space (see logarithmic quantisation for alternatives).

### Logarithmic quantisation

In a similar way, `LogQuant8Array, LogQuant16Array, LogQuant24Array` compresses
an n-dimensional array (non-negative elements only) via logarithmic quantisation.
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
`Ox1` to `0xff, 0xffff, 0xffffff` are then the available bitpatterns to encode
the range from `minimum(A)` to `maximum(A)` logarithmically.
By default the rounding mode for logarithmic quantisation is round-to-nearest
in linear space. Alternatively, a second argument can be either `:linspace` or
`:logspace`, which allows for round-to-nearest in logarithmic space.
Decompression as with linear quantisation via the `Array{T}()` function.

### Rounding modes

Elefridge.jl implements 4 rounding modes for floating-point numbers, `round`
(which is round-to-nearest tie-to-even), `shave`, `halfshave`, `set_one` and
`groom`. All rounding functions take an integer as the number of bits kept as a
second argument and operate on scalars as well as arrays (`Float32` and `Float64`
are supported).
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

The bitpattern entropy (i.e. a measure for the effective use of a given bit-encoding/quantisation)
can be calculated via the `bitentropy` function for n-dimensional arrays with types
that fit into `8,16,24,32,40,48,56` or `64` bit. Unsigned integers of the respective
size are automatically implemented via the [BitIntegers.jl](https://github.com/rfourquet/BitIntegers.jl) package.
Consequently, all common number formats like `Float16/32/64` or even `BFloat16`
or `Posit8/16/32` etc. are supported. The `bitentropy(A::Array)` function's
bottleneck is the sorting function, which sorts the elements in `A` in the
beginning by its corresponding bitpatterns. The entropy is by default calculated
in units of bits, which can be changed with a second argument `::Real` if desired.
```julia
julia> A = rand(Float32,100000000);

julia> bitentropy(A)
22.938590744784577
```
Here, the entropy is about 23 bit, meaning that `9` bits are effectively unused.

### Information content

To calculate the information content of an n-dimensional array (any typ `T` is
   supported that can be reinterpreted as `8,16,24,32,40,48,56` or `64-bit`
   unsigned integer) the following functions are exported:

**bitcount**. The function `bitcount(A::Array)` counts all occurences of the
1-bit in every bit-position in every element of `A`. E.g.

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
The first bit of elements (here: `UInt8`) in `A` are 4 times `1` and 1 times
`0`, etc. In contrast, elements drawn from a uniform distribution U(0,1)
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

**bitcountentropy**. The `bitcountentropy(A::Array)` calculates the entropy of
bit occurences in `A`. For random bits occuring at probabilities p(0) = 0.5, p(1) = 0.5
the entropy for every bit is maximised to 1 bit:
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

**bitpaircount**. The `bitpaircount(A::Array)` function returns a `4xn` (with `n`
being the number of bits in every element of `A`) array, the counts the occurrences
of `00`,`01`,`10`,`11` for all bit-positions in `a in A` across all elements `a` in `A`.
For a length `N` of array `A` (one or multi-dimensional) the maximum occurrence
is `N-1`. E.g.
```julia
julia> A = rand(UInt8,5);
julia> bitstring.(A)
5-element Array{String,1}:
 "01000010"
 "11110110"
 "01010110"
 "01111111"
 "00010100"

julia> bitpaircount(A)
4×8 Array{Int64,2}:
 2  0  0  0  2  0  0  2    # occurences of `00` in the n-bits of UInt8
 1  0  2  1  1  1  0  1    # occurences of `01`
 1  1  2  0  1  0  1  1    # occurences of `10`
 0  3  0  3  0  3  3  0    # occurences of `11`
 ```
The first bit of elements in `A` is as a sequence `01000`. Consequently,
`00` occurs 2x, `01` and `10` once each, and `11` does not occur.
Multi-dimensional arrays are unravelled into a vector, following Julia's
memory layout (column-major).

**bitcondprobability**. Based on `bitpaircount` we can calculate the conditional
entropy of the state of one bit given the state of the previous bit (previous
meaning in the same bit position but in the previous element in the array `A`).
In the previous example we obtain
```julia
julia> Elefridge.bitcondprobability(A)
4×8 Array{Float64,2}:
 0.666667  NaN     0.0  0.0  0.666667  0.0  NaN     0.666667
 0.333333  NaN     1.0  1.0  0.333333  1.0  NaN     0.333333
 1.0         0.25  1.0  0.0  1.0       0.0    0.25  1.0
 0.0         0.75  0.0  1.0  0.0       1.0    0.75  0.0
```
Given the previous bit being `0` there is a 2/3 chance that th next bit is a `0`
too, and a 1/3 change that the next bit is a `1`, i.e. p(next=0|previous=0) = 2/3,
p(1|0), such that p(0|0)+p(1|0)=1 always, if not NaN, and similarly for p(0|1)
and p(1|1).

**bitinformation**. Base on the previous functions, the bitwise information
content defined as
```
Ic = H - q0*H0 - q1*H1
```
for any sequence of bits. `H` is the uncoditional entropy (calculated similarly
to `bitcountentropy` on `A[1:end-1]`), `q0/q1` is the probability of the `0/1`
bit in the sequence and `H0/H1` are the conditional entropies.
`H0 = entropy(p(0|0),p(1|0))` and `H1 = entropy(p(0|1),p(1|1))`.
The bitwise inforamtion content can be calculated with `bitinformation(A::Array)`, e.g.
```julia
julia> A = rand(UInt8,1000000);

julia> bitinformation(A)
8-element Array{Float64,1}:
 2.2513712005789444e-9
 4.3346579969849586e-7
 1.2269593584468552e-6
 1.71376803870249e-6
 1.035394191328809e-6
 2.0511669801548393e-6
 3.1941966260884413e-7
 1.1631417273783029e-7
```
The information of random uniform bits is 0 as the knowledge of a given bit does
not provide any information for the succeeding bits. However, correlated arrays
(which we achieve here by sorting)
```julia
julia> A = rand(Float32,1000000);

julia> sort!(A);

julia> bitinformation(A)
32-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
 0.0005324203865138733
 0.0654345833904318
 0.5219048602771518
 0.9708027628747311
 0.9180796560925731
 0.9993719214848099
 0.9989466303997855
 0.9980955861347737
 0.996610924254921
 0.9941803365302689
 0.989869550776399
 0.9828792401091413
 0.9710813685504847
 0.9515753574258816
 0.9209285439473764
 0.8726093489385611
 0.8010696825274235
 0.6983323983355056
 0.5597023618984831
 0.3923518757327869
 0.21652879399986952
 0.08053181875714122
 0.019144066375869684
 0.00712349142391111
 0.010814913217907396
 0.02885101299472209
 0.0728553098975398
 0.0
 ```
have only zero information in the sign (unused for random uniform distribution
U(0,1)), and in the first exponent bits (also unused due to limited range) and
in the last significant bit (flips randomly). The information is maximised to
1bit for the last exponent and the first significant bits, as knowing the state
of such a bit one can expect the next (or previous) bit to be the same due to
the correlation.

### Round+lossless

The compression technique which involves rounding in combination with a lossless
compression algorithm can be achieved as follows (using the highly correlated
array `A` from previous example).

```julia
julia> using Blosc                    # Blosc library incl blosc & LZ4HC compression

julia> Blosc.set_compressor("lz4hc")  # use LZ4HC

julia> Ar = round(A,5);               # round A to 5 significant bits

julia> Ac = compress(Ar)              # compress into UInt8 array
18562-element Array{UInt8,1}:
 0x02
 0x01
 0x21
 0x04
    ⋮

julia> sizeof(A)/sizeof(Ac)           # compression factor
215.49402004094387                    # very high due to unrealisitically high correlation in A
```
Decompression via
```julia
julia> Ad = decompress(Float32,Ac)

julia> Ad == Ar
true                # indeed lossless
```
It is recommended to round to `k` number of keepbits as informed by `bitinformation`.

### Zfp Compression

Julia bindings to th [zfp compression library](https://computing.llnl.gov/projects/floating-point-compression/zfp-compression-ratio-and-quality) have been developed. This functionality is exported to a separate package: [ZfpCompression.jl](https://github.com/milankl/ZfpCompression.jl) and documentation can be found therein.

## Installation

Not yet registered in the Julia registry, hence do
```julia
(@v1.5) pkg> add https://github.com/esowc/Elefridge.jl
```
in the package manager (access via `]`).
