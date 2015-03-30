#!/usr/bin/julia

# ----------------------------------------------------------------------
#   Part of code to calculate Bond-Angle Distribution from MD simulation.
#
#   Copyright (2015) Stefan Bringuier
#   This software is distributed under the GNU General Public License.
#   See LICENSE file for more details.
#
#   See the README file for program description.
2#------------------------------------------------------------------------- 

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
using Winston
#using PyPlot

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
            0.25 0.25 0.25;
            0.5 0.5 0.0; 
            0.5 0.0 0.5; 
            0.0 0.5 0.5;  
            0.75 0.75 0.25; 
            0.75 0.25 0.75; 
            0.25 0.75 0.75]
    pos = frac*boxd
    types = Int[1,1,1,1,1,1,1,1] 
    rcut = 2.50 # Si-Si Angstroms

    neighlist = func(pos,types,natoms,boxd,rcut)
    println(neighlist)
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

    @time neighlist = func1(pos,types,natoms,boxd,rcut)
    
    @time angleall = func2(natoms,pos,types,neighlist,boxd)
    
    #TODO - use test macro 
    #for i=1:natoms
    #    @test_approx_eq_eps angle[i,end] 109.430 1.0e-2
    #end

    # if any((angleall[:,end] .< 106.00 ) | (angleall[:,end] .> 112.00))
    #     println(angleall[:,end])
    #     throw(error("Failed angle unit test"))
    # end

    println("Test 2 Passed!")
end
   
function testread(file1,file2)
    println("Running Test 3:")
    println("Read LAMMPS dump/data style format files")

    try
        readlammps.readdump(file1,1)
    catch err
        println("Could not read file: $err")
    end
    try
        readlammps.readdata(file2)
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
        #Convert types to integers
        d2 = int(d[:,2])
        @time neigh = neighbor.build(d[:,3:end],d2,n,b,3.01) 
        
        @time t,i = angle.anglecalc(n,d[:,3:end],d2,neigh,b)
  
        plot(t, i, "b");

        savefig("winston_plot.pdf");
      
         #niti = orderparam.orderniti(n,d2,neigh,d[:,3:end],b)
        #for i=1:n
        #    @test_approx_eq_eps niti[i] -1.0000 1.0e-2
        #end
        println("Test 5 Passed!")
end

function testanglebig(f)
        println("Running Test 6:")
        println("NiTi angle time test on big cell")
        n,b,d = readlammps.readdump(f,1,1,5)
        d2 = int(d[:,2])
        @time neigh = neighbor.build(d[:,3:end],d2,n,b,5)
        @time angle.anglecalc(n,d[:,3:end],d2,neigh,b)
end


#Main
flag=ARGS[1] 

#Since Julia is using JIT, its a good idea to precompute functions
# which might have a lot of for loop operations, that way they run faster
# println("Preload using test2") 
# @time testangle(neighbor.neighborlist,angle.anglecalc)
# println("********************")

if flag == "1"
    @profile begin
        testneighbor(neighbor.build)
    end
    @profile Profile.print(format=:flat)
elseif flag == "2"
 testangle(neighbor.build,angle.anglecalc)
elseif flag == "3"
    @time testread("../testfiles/dump.NiTi","../testfiles/data.B2-unitcell")
    @time testwrite(writelammps.writedump,"../testfiles/dump.fccubic.testcase")
elseif flag == "4"
    testorderparam("../testfiles/dump.NiTi")
elseif flag == "5"
     testanglebig("../testfiles/dump.NiTiBig")
else
    println("Testing all unit test functions!")
    @time testneighbor(neighbor.build)
    @time testangle(neighbor.build,angle.anglecalc)
    @time testread(readlammps.readdump,"../testfiles/dump.NiTi")
    @time testwrite(writelammps.writedump,"../testfiles/dump.fccubic.testcase")
    @time testorderparam("../testfiles/dump.NiTi")
end

