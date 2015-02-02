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
include("neighbor.jl")
include("angle.jl")

#Tester Functions
function tester(flag=1)
    if flag == 1
        tester1()
    elseif flag == 2
        tester2()
    elseif flag == 3
        tester3()
    end
end

function tester1()
    # run test on neighbor list cubic diamond Silicon
    #include("neighbor.jl")
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

    neighlist,rijarry = neighborlist(pos,types,natoms,boxd,rcut)

    for i=1:natoms
        if length(neighlist[i]) != 4
            println("Atom: $i fail")
        else
            println("Atom: $i pass")
        end
    end
end
   

function tester2()
    println("Running Test 2:")
    println("Angle calculation routine test")

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

    neighlist,rijarry = neighborlist(pos,types,natoms,boxd,rcut)
    
    angleall = angle(natoms,rijarry,types,neighlist)

    if any((angleall[:,end] .< 106.00 ) | (angleall[:,end] .> 112.00))
        println("Bond angle failed!")
    else
        println("Bond angle passed!")
    end
    
end
   
function tester3()
    println("Not implemented!")
end