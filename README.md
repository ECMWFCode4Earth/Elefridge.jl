[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/milankl/Elefridge.jl.svg?branch=master)](https://travis-ci.com/milankl/Elefridge.jl.svg)
[![AppVeyor](https://img.shields.io/appveyor/ci/milankl/Sherlogs-jl?label=Windows&logo=appveyor&logoColor=white)](https://ci.appveyor.com/project/milankl/Sherlogs-jl)


# Elefridge.jl
***How to put an elephant into a fridge?***  
***Size**: You may buy a very large fridge, but that’s a very large fridge.*  
***Precision**: You may store a picture, but then you don’t know how it looks from behind.*  
***Speed**: You may store its DNA, but then it takes years to recreate an elephant.*  

This is a work-in-progress repository on fast and precise but lossy data compression for climate data.
This repository is inteded for the work on [ECMWF](https://www.ecmwf.int)'s [summer of weather code](https://esowc.ecmwf.int) [challege #14: Size, precision, speed - pick two](https://github.com/esowc/challenges_2020/issues/4).

# Proposal

### Challenges of climate data compression

Weather and climate models produce enormous amounts of data that has to be archived in a fast way, without losing valuable information and, ideally, consuming the smallest storage space possible. This usually poses a dilemma in which not all of the three requirements size, precision, and speed can be satisfied (1). Depending on the data this trade-off has to be optimised for a useful data compression. Climate model data is predominantly calculated with 64-bit floating-point numbers to minimize rounding errors in the calculations, although the majority of bits does not contain real information (2). Simple data compression approaches include the conversion to floats with 32 bit or 16 bit to reduce file sizes, while simultaneously removing many bits from the significand that are presumably meaningless. Bit-grooming (or shaving, resetting) are other approaches (3), which reset meaningless bits such that conventional data compression algorithms will identify them as statistically redundant bit patterns. Removing or grooming the least significant random bits is important as they will not compress with a lossless compression algorithm.

Quantization is a standard technique in signal processing that can also be used for climate data (4,5). Bit patterns will be redefined such that a mapping function is required for conversion. The GRIB data format uses fixed-point numbers to span the range between a minimum and a maximum value per data set to encode all values with 16 bit. Most quantization methods have an equi-distant spacing on a linear-scale, such that fixed-point numbers (including integers) are the natural binary format to be used. However, geophysical data is rarely uniformly distributed on a linear-scale which means that these quantization methods are sub-optimal to minimise the rounding error in compression. Floating-point encoding is approximately equi-distant on a log-scale, such that only log-uniform distributed data will maximise the information entropy.

It is therefore important to understand the data distribution to maximise the precision while reducing the file size. Simplifying the mapping function as well as its inverse between floating-point and the compressed data type is essential for speed. This ECMWF summer of weather code proposal will address these issues on climate data compression as described in the following.

### Maximising the information entropy for climate data

Posits are a recently proposed alternative to floating-point numbers (6), which were found to be better suited for weather and climate applications (7). Posits have a pyramid-shaped precision in log-space, such that log-normal distributed data can be encoded with smaller precision losses compared to floats or fixed-point numbers (Fig. 1). Posits have the advantage that the conversion to floats is relatively simple, as both formats encode a real number in terms of a power of two multiplied with a significand that is encoded as a fixed-point number. Using posit numbers for data compression is presumably a step towards full entropy encoding while maintaining a simple conversion function.

Further flexibility regarding the distribution of the data can be achieved with entropy encoding. A recently developed number format are self-organizing numbers, so called Sonums (8), which can be set-up to maximize the entropy or, often closely related, minimize the decimal rounding error (minimizing the L1 or L2-error is another possibility). Sonums require an alphabet to be stored with the data, which is, however, for large data sets a negligible overhead. The sonum conversion is asymmetric in speed and complexity. While the creation of an alphabet can be a time consuming operation, the back conversion is simple array access.

Data compression can be applied not just on a number-by-number basis, as discussed, but also on a set of numbers. This comes with the advantage that interdependencies of number occurrences can be exploited, especially important due to the typical power-law spectral variance distribution in geophysical data. While this is an attractive approach to climate data compression4,5, it interlinks the compressed data in a way which may slow down the decompression of subsets, as individual numbers cannot be decompressed independently.

### Assessment of precision: Error norms

Errors of lossy compressions can be measured in many ways. A natural way to think about errors on a linear-scale is the L2-error, which quantifies the distance in an Euclidean space of a compressed array to its reference. However, using floating-point numbers this distance changes not just with the number of bits used but also with values itself. Quantifying the L1-log error, which is identical to the decimal (or binary) error, returning the number of correct decimal places (or significant bits), is another way to measure a scale-invariant error. Reducing the number of bits in a floating-point encoded number creates an error that is approximately proportional to the L1-log error, such that an investigation of various error norms to assess the precision retained is proposed.

### Timeline of key developments

This proposal is divided into smaller projects as follows

#### 1. Information entropy of various numbers formats for CAMS data

The first two months will be spent to understand the data distributions of the various chemical tracers in the CAMS data. This is an important step to understand the inefficiency in the current GRIB compression.

#### 2. A posit-based quantization method

Presumably, one bottleneck of the current GRIB compression is the fixed-point number quantization, which encodes non-uniformly distributed with low entropy. A possible improvement would be to replace the fixed-point numbers with posits.

#### 3. Sonums for CAMS data

A next step towards a higher entropy encoding of the CAMS data is the use of sonums, which would learn the optimal quantization from the CAMS data itself. We will investigate whether this learning process has to be repeated for different variables of the CAMS data, but hopefully a single sonum alphabet is sufficient for large parts of the data set.

#### 4. Error norm validation

To validate the precision retained in the compressed data sets we will use various error norms and investigate which norms are important to minimize to maximize the bitwise information. Additional validations can be performed through typical analyses that are based on the CAMS data.

#### 5. Benchmarks for read and write performance

To find the optimal compression method in the size-precision-speed space for CAMS data we will benchmark the different methods for compression and decompression speeds. Based on the requirements and limitations within the ECMWF computing environment the different compression will be identified as more or less suitable. Ideally, this project will find various compression methods, such that the user can specify a compression level with a trade-off between speed and size. Frequently used data would be stored at larger size but with higher decompression speeds, and vice versa for infrequently used data.

#### 6. Documentation and tests

The last month of the project will focus on the documentation to make the resulting compression methods easily accessible. This includes the publication in open-source packages that come with tests and benchmark for reproducibility.

The following timeline is suggested for the project covering the months May to August.

|               | May           | June  | July  | August|
| ------------- |:-------------:|:-----:|:-----:|:-----:|
| Information entropy for CAMS data with floats, fixed-points, logarithmic fixed-points, and posits      | O | O |   |   |
| Development of a posit-based quantization             | O | O | O |   |
| Development of sonums for CAMS data                   |   | O | O |   |
| Validation with different error norms                 |   | O | O |   |
| Benchmarks for read and write performance             |   |   | O | O |
| Documentation and tests for repository                |   |   |   | O |


### References

1.    [Silver, J. D. & Zender, C. S. The compression–error trade-off for large gridded data sets. Geosci. Model Dev. 10, 413–423 (2017).](http://dust.ess.uci.edu/ppr/ppr_SiZ17.pdf)
2.    [Jeffress, S., Düben, P. & Palmer, T. Bitwise efficiency in chaotic models. Proc. R. Soc. Math. Phys. Eng. Sci. 473, 20170144 (2017).](https://royalsocietypublishing.org/doi/10.1098/rspa.2017.0144)
3.    [Zender, C. S. Bit Grooming: statistically accurate precision-preserving quantization with compression, evaluated in the netCDF Operators (NCO, v4.4.8+). Geosci. Model Dev. 9, 3199–3211 (2016).](https://www.geosci-model-dev.net/9/3199/2016/)
4.    [Düben, P. D., Leutbecher, M. & Bauer, P. New Methods for Data Storage of Model Output from Ensemble Simulations. Mon. Weather Rev. 147, 677–689 (2018).](https://journals.ametsoc.org/doi/full/10.1175/MWR-D-18-0170.1)
5.    [Düben, P. D. A New Number Format for Ensemble Simulations. J. Adv. Model. Earth Syst. 10, 2983–2991 (2018).](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001420)
6.    [Gustafson, J. L. & Yonemoto, I. Beating Floating Point at its Own Game: Posit Arithmetic. Supercomput. Front. Innov. 4, 16 (2017).](http://www.johngustafson.net/pdfs/BeatingFloatingPoint.pdf)
7.    [Klöwer, M., Düben, P. D. & Palmer, T. N. Posits as an alternative to floats for weather and climate models. in Proceedings of the Conference for Next Generation Arithmetic 2019 on   - CoNGA’19 1–8 (ACM Press, 2019). doi:10.1145/3316279.3316281.](https://dl.acm.org/doi/abs/10.1145/3316279.3316281)
8.    [Klöwer, M. Sonums - a maximum entropy number format. (Zenodo, 2019). doi:10.5281/zenodo.3531887.](https://github.com/milankl/Sonums.jl)
