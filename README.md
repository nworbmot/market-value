# Market Value

This repository provides the code for the simulations in [Decreasing
market value of variable renewables is a result of policy, not
variability](https://arxiv.org/abs/2002.05209) by Tom Brown and Lina
Reichenberg (2020).


# Data Requirements

The script `solve_network.py` requires the input files from the EMMA
model provided in the supplementary material of the paper [The Market
Value of Variable Renewables: The Effect of Solar and Wind Power
Variability on their Relative
Price](https://doi.org/10.1016/j.eneco.2013.02.004) in a sub-folder
`emma/`. Some version of `pandas` required empty columns to be deleted
in the Excel files.


# Running the code

The code uses the
[snakemake](https://snakemake.readthedocs.io/en/stable/) workflow
manager. It calculates scenarios defined in `config.yaml`.

# License

Copyright 2019-2020 Tom Brown (KIT)

This program is free software: you can redistribute it and/or modify
it under the terms of the [GNU General Public
License](http://www.gnu.org/licenses/gpl-3.0.en.html) as published by
the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
