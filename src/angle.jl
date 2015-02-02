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


function angle(natoms,dist,types,neighlist)
    anglelist = Vector{Float64}[] #Empty Vector containg Array of float64 type
    for i=1:natoms
        itype = types[i]
        for jpair in enumerate(neighlist[i])
            jtype = types[pair[2]]
            for kpair in enumerate(neighlist[j])
                ktype = types[kpair[2]]
                rij = rdist[i,j,:]
                rik = rdist[i,k,:]
                dotprod = rij * rik
                
                magv1 = √(∑(rij.^2))
                magv2 = √(∑(rik.^2))
                
                costheta = dotprod / (magv1*magv2)
                
                #TODO - see if there is a better way to do this
                #Add Array to Array anglelist, the data type is
                #an Array of size n,1 where the entry corresponds to
                #an array. What I really want is an Array of Float64 that is
                #of the size n by 4
                push!(anglelist,[i,pairj[1],pairk[1],costheta])

             end
        end
    end

    #NOTE - anglelist is not an n x 4 Array{Float64,2}
    # however I want an array of n x 4
    n = length(anglelist) 
    m = length(anglelist[1])
    anglearry = zeros(n,m)
    for i=1:n
        entry = b[n]
        theta = rad2deg(acos(entry[end]))
        anglearry[i,:end-1] = entry[end-1]
        anglearry[i,end] = theta
    end

    return anglearry
