# Change current directory to parent folder.
cd(@__DIR__)
cd("..")
pwd()

# Load all external packages
using Revise

using DataFrames, Dates, Statistics
using CSV

using OrdinaryDiffEq
using Lux
using ComponentArrays
using SciMLSensitivity

using Optimization, OptimizationBBO
using OptimizationOptimisers
using OptimizationPolyalgorithms
using OptimizationNLopt
using OptimizationOptimJL
using Zygote

using Interpolations
using Plots
using JLD2
using BenchmarkTools

using Random
Random.seed!(300)   #Seed for reproducibility