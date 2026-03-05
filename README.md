# Model A 2D

Forked from https://github.com/vskokov/modelA.

Usage

`julia -t N modelA2D.jl [seed] [L] [iter]`

where 'iter' represents what checkpoint file to load, if the checkpoint file doesn't already exist, it will generate it. If you have no checkpoint files already saved, it will generate them all up to the value of `iter`. For example if you have no checkpoint files previously, `julia -t 4 modelA2D.jl 398124 16 5` will generate checkpoint files from 1 through 5.

## To do:
Suppose you have file 4 but no other files before it and want to create a 5th file, currently the program will unnecessarily thermalize and create files 1-3 before creating file 5. This will be fixed. May add an optional argument to enable/disable the bulk creation of checkpoint files.