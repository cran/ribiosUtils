*ribiosUtils*: Swiss-army knife for computational biology tasks in drug discovery
===

![R-CMD-check badge](https://github.com/bedapub/ribiosUtils/workflows/R-CMD-check/badge.svg)

## What is *ribiosUtils*?

*ribiosUtils* is a R package that performs various routine tasks for bioinformatics and computational biology research in drug discovery. It is distributed under the GPL-3 license.

## Installation

Run following commands in the R console.

```{R}
library(devtools)
devtools::install_github("bedapub/ribiosUtils")
```

## History

*ribiosUtils* is part of the *ribios* software suite (*ribios* stands for *R* *i*nterface to *BIOS*). The package was previously a sub-directory of the *ribios* project. It started as a R implementation of the BIOS library, publicly known as Bioinfo-C library, developed by Clemens Broger *et al.* at F. Hoffmann-La Roche for many years. Since its inception in 2012, many functionalities have been added, which are implemented in either C (or C++) or R. 

[Jitao David Zhang](mailto:jitao_david.zhang@roche.com) maintains and develops *ribiosUtils* and other ribios packages in memory of Clemens Broger, a pioneer of bioinformatics and cheminformatics in drug discovery, a man true to himself.
