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
module neighbor
using types

#Docstring function should be supported in Julia v0.4
#@doc """ Build neighbor list for atomic configuration
# Inputs:
# pos - atomic positions in cartesian coordinates
# types - array of atomic types
# natoms - number of atoms
# boxd - simulation cell dimensions (array)
# rcut - cutoff distance to build neighbor list in units of position.
# 
# Outputs:
# neighlist - neighbor list for atomic configuration
# rij_arry - multidimension array of seperation vectors between 
#             atoms i and j
# 
#TODO - Memory gc is fairly high try to optimize.
#       Profiling shows high usage of abstract_eval()
#       Reduce number of function inputs
#""" ->
function neighborlist(pos::Array,types::Array,natoms::Integer,boxd::Array,rcut::Real)
    ## Function to Build neighbor list of atoms

    boxdinv = inv(boxd) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator
    
    #Initialize,list comprehensions
    index = Int[ i for i=1:natoms ] # or [1:natoms]
    neighlist = [ Int32[] for i=1:natoms ] #Array of type Int32 Arrays to store atom IDs
    
    for i=1:natoms

        #sizehint(neighlist[i],convert(Int,round(natoms/2)))
        
        #PBC - Min. Image Conv., origin 0,0,0
        si = spos[i,:]
        sij = spos .- si
        sij -= round(sij)
        rij = sij * boxd 

        #Fill distance array based on cutoff
        #Note: . before a operator provide bit/element-wise behavior
        dij = √sum(rij.^2,2) #.^ element wise
        slice = (dij .< rcut) & (dij .> 0.0e-12)

        #find function returns indices where true
        #neighlist is somewhat like a list of empty lists
        keep = find(slice) 
        neighlist[i] = index[keep]
        #rij_arry[i,keep,:] = rij[keep,:]

    end
    return neighlist 
end

function neighborlist2(atoms::Atoms,boxd::Cell,rcut::Real)
    ## Function to Build neighbor list of atoms
    
    boxdinv  = inv(boxd.cell) #Invert matrix to scale coordinates
    spos = atoms.pos * boxdinv #Here the * operator is an array-array operator
    
    #Initialize,list comprehensions
    neighlist = [ Int32[] for i=1:atoms.num ] #Array of type Int32 Arrays to store atom IDs

    # NOTE - Causes issues with memory.
    # NOTE - this is probably a poor programming choice and I'm not sure it actually saves time
    # and is detrimental towards memory use 
    # rij_arry = zeros(Float32,natoms,natoms,3) # Distance array, Float32 type to save on memory
    
    for i=1:atoms.num

        #sizehint(neighlist[i],convert(Int,round(natoms/2)))
        
        #PBC - Min. Image Conv., origin 0,0,0
        si = spos[i,:]
        sij = spos .- si
        sij -= round(sij)
        rij = sij * boxd.cell 

        #Fill distance array based on cutoff
        #Note: . before a operator provide bit/element-wise behavior
        dij = √sum(rij.^2,2) #.^ element wise
        slice = (dij .< rcut) & (dij .> 0.0e-12)

        #find function returns indices where true
        #neighlist is somewhat like a list of empty lists
        keep = find(slice) 
        neighlist[i] = atoms.id[keep]
        
    end
    return neighlist 
end


function parallelneighborlist()
    println("Not yet implemented!")
end

end