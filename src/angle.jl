#!/usr/bin/julia
function angle(natoms,dist,types,neighlist)
    anglelist = {} #Empty any type array
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
                
                push!(anglelist,[i,pairj[1],pairk[1],costheta])
             end
        end
    end
    anglelist[:            

def angle(natoms,types,neighlist,rdist):
    ''' TODO - Make more general
    still requires edits for a specific system'''
    
    print 'Starting Angle Calculation'
    anglelist = []
    for i in xrange(natoms):
        itype = types[i]
        for jj,j in enumerate(neighlist[i][:-1]):
            jtype = types[j]

            if itype != jtype:
                rij = rdist[i,j,:]
                
            for k in neighlist[i][jj+1:]:              
                ktype = types[k]

                rij = rdist[i,j,:]
                rik = rdist[i,k,:]
                dotprod = np.dot(rij, rik) 
                
                # sqrt(vx*vx + vy*vy + vz*vz)
                magv1 = np.sqrt(np.sum(np.power(rij,2)))
                magv2 = np.sqrt(np.sum(np.power(rik,2)))
                magnitude = magv1 * magv2

                costheta = dotprod / magnitude

                anglelist.append([itype,jtype,ktype,costheta])

                
    #Numpy arry of angle list
    anglearry = np.array(anglelist)
    anglearry[:,3] = np.rad2deg(np.arccos(anglearry[:,3]))
    print "*************************"
    return anglearry

function neighborlist(pos,types,natoms,boxd,rcut)
    ## Function to Build neighbor list of atoms

    boxdinv = inv(boxd) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator
    
    #Initialize,list comprehensions
    index = [ i for i=1:natoms ] # or [1:natoms]
    neighlist = [ Int32[] for i=1:natoms ] #Array of type Int32 to store atom IDs
    rij_arry = zeros(Float32,natoms,natoms,3) # Distance array, Float32 type to save on memory
    
    for i=1:natoms
        #
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
