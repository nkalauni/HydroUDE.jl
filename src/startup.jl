
# Change current directory to parent folder.
cd(@__DIR__)
cd("..")
pwd()

using Pkg
#Pkg.instantiate()

# Load all external packages
using Revise

using DataFrames, Dates, Statistics
using CSV

using OrdinaryDiffEq
# using Lux, DiffEqFlux
using Lux
using ComponentArrays
using SciMLSensitivity

using Optimization, OptimizationBBO
using OptimizationOptimisers
using OptimizationPolyalgorithms
using OptimizationNLopt
using Zygote

using Interpolations

using Plots

# using Pkg
# Pkg.add("JSON")
# using JSON

# using Pkg
# Pkg.add("JLD2")
using JLD2  

# using Pkg
# Pkg.add("OptimizationOptimJL")
using OptimizationOptimJL

using Random
Random.seed!(300)   #Seed for reproducibility

