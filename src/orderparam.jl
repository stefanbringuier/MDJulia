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
module orderparam

function orderniti(natoms,types,neighlist,pos,boxd)
    println("Calculating order parameter based on NiTi Mutter et al.")

    #Heuristic constants from Mutter et al. paper
    #Characteristic distances between NiTi
    const d0B19 = 2.460232 #Angstrom
    const d1B19 = 2.646524 #ditto
    const d0B2 = 2.611067  
    const d1B2 = 2.611067
    orderp = zeros(natoms)
    
    boxdinv = inv(boxd) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator

    # Bruteforce Approach to get distance 
    # between  i-j of different TYPES 
    dist = [ Float64[] for i=1:natoms ]
    for i=1:natoms
        itype = types[i]
        #Empty array to store distances
        tmpdij = Float64[]
        for j in neighlist[i]
            jtype = types[j]
            if itype == jtype #Only Ni-Ti bonds
                continue
            end
            
            #TODO - How can I improve this?
            #NOTE: spos[i,:] is 1x3 2D-Array
            #PBC - Min. Image Conv., origin 0,0,0
            sij = spos[j,:] - spos[i,:]
            sij -= round(sij)
            rij = vec(sij * boxd) 

            # rijvec size n x n x3            
            #rij = reshape(rijvec[i,j,:],(1,3))
            dij = √(sum(rij.^2))
            push!(tmpdij,dij)
        end
        dist[i] = tmpdij
    end

    #Now calculate order parameter, no direction dependenece used
    for i=1:natoms
        sort!(dist[i])
        #Order Parameter from D. Mutter et al. 
        d0 = 0.00
        d1 = 0.00
        
        #Ensure 8 neighbors
        if size(dist[i],1) < 8 
            throw(DomainError())
        end
        #Only use first 8 neighbors
        for k=1:8
            if k < 7
                d0 += dist[i][k]
            else
                d1 += dist[i][k]
            end
        end   
        d0 /= 6.00
        d1 /= 2.00
        chi = (d0*(d1B2+d1B19) - d1*(d0B2+d0B19)) / (d0B2*(d0B19-d1B19))
        orderp[i] = chi
    end

    return orderp
end

end