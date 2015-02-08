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
# pos - position array of atoms
# (legacy)rdist - multidimension array of seperation vectors between 
#         atoms i and j. Same as rij_arry in neighbor.jl
# type -  array of atomic types
# neighlist - neighbor list for atomic configuration (see neighbor.jl)
#
# Outputs:
# anglearry - 2D array of bond angles format: i j k theta.
#             Atom i is the central atom.
#
#TODO - Optimize if possible.
#       REDUCE ARGUEMENTS REQUIRED BY FUNCTION
#""" ->
function anglecalc(natoms,pos,types,neighlist,boxd)
    
    boxdinv = inv(boxd) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator

    anglelist = Vector{Float64}[] #Empty Vector containg Array of float64 type
    for i=1:natoms
        itype = types[i]
        for jpair in enumerate(neighlist[i][1:end-1]) #jpair is a tuple (index,entry)
            jtype = types[jpair[2]]
            j = jpair[2]
            for k in neighlist[i][jpair[1]+1:end]
                ktype = types[k]

                #TODO - How can I improve this?
                #NOTE: spos[i,:] is 1x3 2D-Array
                #PBC - Min. Image Conv., origin 0,0,0
                sij = spos[j,:] - spos[i,:]
                sij -= round(sij)
                rij = vec(sij * boxd) 
                sik = spos[k,:] - spos[i,:]
                sik -= round(sik)
                rik = vec(sik * boxd) 
                
                #Legacy approach
                #rij = reshape(rdist[i,j,:],3) #reshape to column vector
                #rik = reshape(rdist[i,k,:],3)
               
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
        #TODO - Dumb! better handle DomainError()
        if entry[end] > 1.0000000 
            theta=0.000
        elseif  entry[end] < -1.0000000
            theta=180.000
        else
            theta = rad2deg(acos(entry[end]))
        end
        anglearry[i,1:end-1] = entry[end-1]
        anglearry[i,end] = theta
    end

    return anglearry
end


#Different Approach
function angletesting(natoms,pos,types,neighlist,boxd)
    
    boxdinv = inv(boxd) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator

    anglelist = Vector{Float64}[] #Empty Vector containg Array of float64 type
    for i=1:natoms
        itype = types[i]
        slice = neighlist[i]
        s = spos[slice,:] - spos[i,:]
        s -= round(s)
        r = s * boxd 
        for j in neighlist[i] #jpair is a tuple (index,entry)                      
            jtype = types[j]
            for k in neighlist[i][jpair[1]+1:end]
                ktype = types[k]
               
                dotprod = dot(r[i,j],r[i,k]) 
                
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
        #TODO - Dumb! better handle DomainError()
        if entry[end] > 1.0000000 
            theta=0.000
        elseif  entry[end] < -1.0000000
            theta=180.000
        else
            theta = rad2deg(acos(entry[end]))
        end
        anglearry[i,1:end-1] = entry[end-1]
        anglearry[i,end] = theta
    end

    return anglearry
end

end