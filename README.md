[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

# Elefridge.jl - Compressing atmospheric data into its real information

**Milan Klöwer**\
Atmospheric, Oceanic and Planetary Physics, University of Oxford\
*milan.kloewer@physics.ox.ac.uk*

This repository contains analysis and plotting scripts for an upcoming publication, which is submitted as 

M Klöwer, M Razinger, JJ Dominguez, PD Düben, TN Palmer, 2021. *Compressing atmospheric data into its real information content.*

Analysis notebooks can be found in [/nb](https://github.com/esowc/Elefridge.jl/blob/master/nb/bitinformation.ipynb). This repository
also summarises the results on [ECMWF](https://www.ecmwf.int)'s [summer of weather code](https://esowc.ecmwf.int)
[challege #14: Size, precision, speed - pick two](https://github.com/esowc/challenges_2020/issues/4)
in [summary.md](https://github.com/esowc/Elefridge.jl/blob/master/summary.md). The original prosal is in 
[proposal.md](https://github.com/esowc/Elefridge.jl/blob/master/proposal.md).

As part of this project, the following Julia packages have been developed
- [BitInformation.jl](https://github.com/milankl/BitInformation.jl)
- [LinLogQuantization.jl](https://github.com/milankl/LinLogQuantization.jl)
- [ZfpCompression.jl](https://github.com/milankl/ZfpCompression.jl)

## Abstract

Hundreds of petabytes of data are produced annually at weather and climate forecast centres worldwide.
Compression is inevitable to reduce storage and to facilitate data sharing. Current techniques do not
distinguish the real from the false information in data. We define the bitwise real information content
from information theory for data from the Copernicus Atmospheric Monitoring Service (CAMS). Most variables
contain less than 7 bits of real information per value, which are also highly compressible due to
spatio-temporal correlation. Rounding bits without real information to zero facilitates lossless
compression algorithms and encodes the uncertainty within the data itself. The entire CAMS data
is compressed by a factor of 17x, relative to 64-bit floats, while preserving 99% of real information.
Combined with 4-dimensional compression to exploit the spatio-temporal correlation, factors beyond 60x
are achieved without an increase in forecast errors. A data compression Turing test is proposed to
optimize compressibility while minimizing information loss for the end use of weather and climate
forecast data. 

