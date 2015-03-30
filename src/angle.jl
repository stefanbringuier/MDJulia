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
using types
const radtodeg = 180.000/pi

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
# thetas - angles in degrees
# intensity - frequency of an angle occurance
#
#TODO - Handle function call so that specific types are used
#       REDUCE ARGUEMENTS REQUIRED BY FUNCTION
#""" ->
function anglecalc(natoms::Integer,pos::Array,types::Array,neighlist::Array,boxd::Array,nbins=90)
    
    binsize = pi/nbins 
    degbinsize = rad2deg(binsize) 
    thetas = Float64[ (i*binsize) for i=1:nbins ]
    thetas = rad2deg(thetas)
    intensity = zeros(Int,nbins)

    boxdinv = inv(boxd) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator

    for i=1:natoms
        itype = types[i] :: Int  
        for jpair in enumerate(neighlist[i][1:end-1]) #jpair is a tuple (index,entry)
            jtype = types[jpair[2]] :: Int 
            j = jpair[2]
            sij = spos[j,:] - spos[i,:]
            sij -= round(sij)
            rij = vec(sij * boxd)
            magv1 = √(sum(rij.^2)) ::Float64
            for k in neighlist[i][jpair[1]+1:end]
                ktype = types[k] :: Int

                #NOTE: spos[i,:] is 1x3 2D-Array
                #PBC - Min. Image Conv., origin 0,0,0
                sik = spos[k,:] - spos[i,:]
                sik -= round(sik)
                rik = vec(sik * boxd) 
                magv2 = √(sum(rik.^2)) ::Float64
                dotprod = dot(rij,rik) ::Float64 
                
                costheta = dotprod / (magv1*magv2)::Float64
                
                #Might want to rework this so that binning is done
                #with costheta, since acos() can give domainerror()
                #when round-off occurs so we have to check this,
                #for example when costheta=-1.00001
                #Used defined constant radtodeg instead of rad2deg()
                if costheta <= -1.000000 
                    theta = 180.000000
                else
                    theta = acos(costheta) * radtodeg
                end

                ibin = int(theta/degbinsize) 
                intensity[ibin] = intensity[ibin] + 1
               
             end
        end
    end
    return thetas,intensity
end


#Multiple Dispatch - use types
function anglecalc(atoms::Atoms,boxd::Cell,neighlist::Array,nbins=90)
    
    binsize = pi/nbins 
    degbinsize = rad2deg(binsize) 
    thetas = Float64[ (i*binsize) for i=1:nbins ]
    thetas = rad2deg(thetas)
    intensity = zeros(Int,nbins)

    #Type Atoms does this already
    #boxdinv = inv(boxd) #Invert matrix to scale coordinates
    #spos = pos * boxdinv #Here the * operator is an array-array operator

    for i=1:atoms.num
        itype = atoms.typ[i]   
        for jpair in enumerate(neighlist[i][1:end-1]) #jpair is a tuple (index,entry)
            jtype = atoms.typ[jpair[2]]  
            j = jpair[2]
            sij = atoms.spos[j,:] - atoms.spos[i,:]
            sij -= round(sij)
            rij = vec(sij * boxd.cell)
            magv1 = √(sum(rij.^2)) ::Float64
            for k in neighlist[i][jpair[1]+1:end]
                ktype = atoms.typ[k] 

                #NOTE: spos[i,:] is 1x3 2D-Array
                #PBC - Min. Image Conv., origin 0,0,0
                sik = atoms.spos[k,:] - atoms.spos[i,:]
                sik -= round(sik)
                rik = vec(sik * boxd.cell) 
                magv2 = √(sum(rik.^2)) ::Float64
                dotprod = dot(rij,rik) ::Float64 
                
                costheta = dotprod / (magv1*magv2)::Float64
                
                #Might want to rework this so that binning is done
                #with costheta, since acos() can give domainerror()
                #when round-off occurs so we have to check this,
                #for example when costheta=-1.00001
                #Used defined constant radtodeg instead of rad2deg()
                if costheta <= -1.000000 
                    theta = 180.000000
                else
                    theta = acos(costheta) * radtodeg
                end

                ibin = int(theta/degbinsize) 
                intensity[ibin] = intensity[ibin] + 1
               
             end
        end
    end
    return thetas,intensity
end


end