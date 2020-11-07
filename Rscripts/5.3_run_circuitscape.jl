#-------------------------------------------------------------------------------
## 5.3 Run Circuitscape
## 2020
## Inputs: INI files, resistance rasters
## Outputs: connectivity rasters
#-------------------------------------------------------------------------------

using Pkg ; Pkg.activate(".") ; Pkg.instantiate()

using Distributed
using LinearAlgebra.BLAS
using Circuitscape

# Set BLAS threads to 1 to avoid oversubscription
BLAS.set_num_threads(1)

# Read the ini files
searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
dir = "config/all/"
ext = "ini"
ini_list = dir .* searchdir(dir, ext)

for file in ini_list
    compute(file)
end
#pmap(compute, ini_list)
