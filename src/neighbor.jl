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
#TODO - Memory usage is significant, try to optimize
function neighborlist(pos,types,natoms,boxd,rcut)
    ## Function to Build neighbor list of atoms

    boxdinv = inv(boxd) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator
    
    #Initialize,list comprehensions
    index = [ i for i=1:natoms ] # or [1:natoms]
    neighlist = [ Int32[] for i=1:natoms ] #Array of type Int32 Arrays to store atom IDs
    rij_arry = zeros(Float32,natoms,natoms,3) # Distance array, Float32 type to save on memory
    
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
        rij_arry[i,keep,:] = rij[keep,:]

    end
    return neighlist, rij_arry
end

#TODO - Untested
function neighborlist_bf(pos,types,natoms,boxd,rcut)
    ## Function to Build neighbor list of atoms: Bruteforce
       
    boxdinv = inv(boxd) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator
    
    #Initialize,list comprehensions
    index = [ i for i=1:natoms ] # or [1:natoms]
    neighlist = [ Int32[] for i=1:natoms ] #Array of type Any to store atom IDs
    
    for i=1:natoms-1
        for j=i+1:natoms
        #PBC - Min. Image Conv., origin 0,0,0
        si = spos[i,:]
        sj = spos[j,:]
        sij = sj .- si
        sij -= round(sij)
        rij = sij * boxd 

        #Fill distance array based on cutoff
        #Note: . before a operator provide bit/element-wise behavior
        dij = √sum(rij.^2) #.^ element wise
        keep = (dij < rcut) & (dij > 0.0e-12)
            if slice == true
                push!(neighlist[i], j)
                push!(neighlist[j], i)
            end
        end
    end
    return neighlist
end

end