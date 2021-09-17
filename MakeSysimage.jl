# run this to make a precompiled system image with the packages used
# in Econometrics pre-compiled, and with some commonly used functions
# also pre-compiled.
#
# NOTE: this system image does not work with VS Code.
using Pkg, PackageCompiler
Pkg.activate(".")
PackageCompiler.create_sysimage(:Econometrics; sysimage_path="EconometricsSysimage.so",
                                                precompile_execution_file="warmup.jl")

