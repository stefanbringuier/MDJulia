#!/usr/bin/julia

# ----------------------------------------------------------------------
#   Part of code to calculate Bond-Angle Distribution from MD simulation.
#
#   Copyright (2015) Stefan Bringuier
#   This software is distributed under the GNU General Public License.
#   See LICENSE file for more details.
#
#   See the README file for program description.
#------------------------------------------------------------------------- 

#TODO - add command-line parser
#include("neighbor.jl")
#include("angle.jl")
#Load in current working directory
push!(LOAD_PATH,pwd())
using Base.Test
using neighbor
using angle
using readlammps
using writelammps
using lattice
using orderparam


function testneighbor(func::Function)
    # run test on neighbor list cubic diamond Silicon
    # func::Function - should be a function that builds a neighborlist and returns
    # the neighborlist and distance array.

    println("Running Test 1:")
    println("Neighbor list build test for cubic diamond Silicon (PBC)")

    boxd =  [5.43 0. 0.;
             0. 5.43 0.; 
             0. 0. 5.43]
    natoms = 8
    frac = [0.0 0.0 0.0; 
            0.5 0.5 0.0; 
            0.5 0.0 0.5; 
            0.0 0.5 0.5; 
            0.25 0.25 0.25; 
            0.75 0.75 0.25; 
            0.75 0.25 0.75; 
            0.25 0.75 0.75]
    pos = frac*boxd
    types = Int[1,1,1,1,1,1,1,1] 
    rcut = 2.50 # Si-Si Angstroms

    neighlist,rijarry = func(pos,types,natoms,boxd,rcut)

    for i=1:natoms
        @test length(neighlist[i]) == 4 
    end
    println("Test 1 Passed!")
end
   

function testangle(func1::Function,func2::Function)
    # Run bond angle test on  cubic diamond Silicon
    # func1::Function - should be a function that builds a neighborlist and returns
    # the neighborlist and distance array.
    # func2::Function - a function which calculates the bond angles and returns an array.

    println("Running Test 2:")
    println("Bond Angle calculation Routine")

    boxd =  [5.43 0. 0.;
             0. 5.43 0.;
             0. 0. 5.43]
    natoms = 8
    frac = [0.0 0.0 0.0; 
            0.5 0.5 0.0; 
            0.5 0.0 0.5; 
            0.0 0.5 0.5; 
            0.25 0.25 0.25; 
            0.75 0.75 0.25; 
            0.75 0.25 0.75; 
            0.25 0.75 0.75]
    pos = frac*boxd
    types = Int[1,1,1,1,1,1,1,1] 
    rcut = 2.50 # Si-Si Angstroms

    neighlist,rijarry = func1(pos,types,natoms,boxd,rcut)
    
    angleall = func2(natoms,rijarry,types,neighlist)
    
    #TODO - use test macro 
    #for i=1:natoms
    #    @test_approx_eq_eps angle[i,end] 109.430 1.0e-2
    #end

    if any((angleall[:,end] .< 106.00 ) | (angleall[:,end] .> 112.00))
        throw(error("Failed angle unit test"))
    end

    println("Test 2 Passed!")
end
   
function testread(readfunction::Function,file)
    println("Running Test 3:")
    println("Read LAMMPS dump style format file")

    try
        readfunction(file,1)
    catch err
        println("Could not read file: $err")
    end
end

function testwrite(writefunction::Function,f)
    println("Running Test 4:")
    println("Write LAMMPS dump style format file")

    #Creat lattice
    typ,box,pos = lattice.fccubic((4,4,4),3.80)
    lx,ly,lz = maximum(pos[:,1]),maximum(pos[:,2]),maximum(pos[:,3])
    try
        writefunction(1,box,typ,pos,file=f)
    catch err
        println("Could not write file: $err")
    end
   
end

function testorderparam(f)
        println("Running Test 5:")
        println("NiTi order parameter")
        n,b,d = readlammps.readdump(f,1,1,5)
        neigh,rija = neighbor.neighborlist(d[:,3:end],d[:,2],n,b,5)
        angles = angle.anglecalc(n,rija,d[:,2],neigh)
        niti = orderparam.orderniti(n,d[:,2],neigh,rija)
        for i=1:n
            @test_approx_eq_eps niti[i] -1.0000 1.0e-2
        end
        println("Test 5 Passed!")
end

#Main
flag=ARGS[1] 

if flag == "1"
    @profile begin
        testneighbor(neighbor.neighborlist)
    end
    @profile Profile.print(format=:flat)
elseif flag == "2"
    @time testangle(neighbor.neighborlist,angle.anglecalc)
elseif flag == "3"
    @time testread(readlammps.readdump,"../testfiles/dump.NiTi")
    @time testwrite(writelammps.writedump,"../testfiles/dump.fccubic.testcase")
elseif flag == "4"
    @time testorderparam("../testfiles/dump.NiTi")
else
    println("Testing all unit test functions!")
    @time testneighbor(neighbor.neighborlist)
    @time testangle(neighbor.neighborlist,angle.anglecalc)
    @time testread(readlammps.readdump,"../testfiles/dump.NiTi")
    @time testwrite(writelammps.writedump,"../testfiles/dump.fccubic.testcase")
    @time testorderparam("../testfiles/dump.NiTi")
end

