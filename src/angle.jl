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

module angle

#Docstring function should be supported in Julia v0.4
#@doc """ Calculate Bond Angle Distribution
# Inputs:
# natoms - number of atoms
# rdist - multidimension array of seperation vectors between 
#         atoms i and j. Same as rij_arry in neighbor.jl
# type -  array of atomic types
# neighlist - neighbor list for atomic configuration (see neighbor.jl)
#
# Outputs:
# anglearry - 2D array of bond angles format: i j k theta.
#             Atom i is the central atom.
#
#TODO - Optimize if possible.
#       Reduce number of function inputs.
#""" ->
function anglecalc(natoms,rdist,types,neighlist)
    anglelist = Vector{Float64}[] #Empty Vector containg Array of float64 type
    for i=1:natoms
        itype = types[i]
        for jpair in enumerate(neighlist[i][1:end-1]) #jpair is a tuple (index,entry)
            jtype = types[jpair[2]]
            j = jpair[2]
            for k in neighlist[i][jpair[1]+1:end]
                ktype = types[k]

                rij = reshape(rdist[i,j,:],3) #reshape to column vector
                rik = reshape(rdist[i,k,:],3)
               
                dotprod = dot(rij,rik) 
                
                magv1 = √(sum(rij.^2))
                magv2 = √(sum(rik.^2))
                
                costheta = dotprod / (magv1*magv2)
                
                #TODO - see if there is a better way to do this
                #Add Array to Array anglelist, the data type is
                #an Array of size n,1 where the entry corresponds to
                #an array. What I really want is an Array of Float64 that is
                #of the size n by 4
                push!(anglelist,[itype,jtype,ktype,costheta])                
             end
        end
    end

    
    #NOTE - anglelist is not an n x 4 Array{Float64,2}
    # however I want an array of n x 4
    n = length(anglelist) 
    m = length(anglelist[1])
    anglearry = zeros(n,m)
    for i=1:n
        entry = anglelist[i]
        theta = rad2deg(acos(entry[end]))
        anglearry[i,1:end-1] = entry[end-1]
        anglearry[i,end] = theta
    end

    return anglearry
end

end