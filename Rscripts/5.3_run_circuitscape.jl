#-------------------------------------------------------------------------------
## 5.3 Run Circuitscape
## 2020
## Inputs: INI files, resistance rasters
## Outputs: connectivity rasters
#-------------------------------------------------------------------------------

using Pkg ; Pkg.activate(".") ; Pkg.instantiate()

using Distributed
using LinearAlgebra.BLAS

# Set BLAS threads to 1 to avoid oversubscription
BLAS.set_num_threads(1)

# Read the ini files
searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
dir = "config/all/"
ext = "ini"
ini_list = dir .* searchdir(dir, ext)
# Get cores
#cores = 2
# Add cores with prokect flag
#addprocs(cores, exeflags="--project")
# Still need to declare Circuitscape everywhere
using Circuitscape
# META-PARALLELIZATION => Call to pmap, batch_size size in question
compute(ini_list[1])
#pmap(compute, ini_list)
