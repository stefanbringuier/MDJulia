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
#TODO - Rewrite this entire module. Julia can "shine" here by using multiple dispatch approach
# i.e. one function for any structure
module lattice

#Define C struc like object
type Atom
    id::Int
    typ::Int
    x::Float64
    y::Float64
    z::Float64
end

function scubic(multiple,latparam)
    #Generate simple cubic structure 
    # multiple should be a tuple like object
    
    basis = [0.0 0.0 0.0] #Row vector 1x3
    N = multiple[1]*multiple[2]*multiple[3]
    
    #Try out Atoms type
    #lattice = zeros(N,3)
    #id = [i for i=1:N*length(basis)]
    types = [ 1 for i=1:N*length(basis)]
    atoms = [Atom(i,1,0.,0.,0.) for i=1:N*length(basis)]
    
    n = 1 #Atomic index
    #Origin: 0,0,0
    for i=0:multiple[1]-1
        for j=0:multiple[2]-1
            for k=0:multiple[3]-1
                #lattice[n,1] = (i+basis[1])*latparam
                #lattice[n,2] = (j+basis[2])*latparam
                #lattice[n,3] = (k+basis[3])*latparam
                atoms[n].x = (i+basis[1])*latparam
                atoms[n].y = (j+basis[2])*latparam
                atoms[n].z = (k+basis[3])*latparam
                n += 1
            end
        end
    end
    return atoms
    #return types,lattice
end

function fccubic(multiple,latparam)
    # Generate face-centered cubic structure 
    # multiple should be a tuple like object
    
    basis = [0.0 0.0 0.0;
             0.5 0.5 0.0;
             0.5 0.0 0.5;
             0.0 0.5 0.5] #4x3 Array
    N = size(basis,1)*multiple[1]*multiple[2]*multiple[3]
    lx = (multiple[1]+1)*latparam # plus one is for PBC
    ly = (multiple[2]+1)*latparam 
    lz = (multiple[3]+1)*latparam 
    cell = [lx 0.0 0.0;
            0.0 ly 0.0;
            0.0 0.0 lz]
    lattice = zeros(N,3)
    types = [ 1 for i=1:N]
    
    #TODO - reduce loops (i.e. vector operations)
    n = 1 #Atomic index
    #Origin: 0,0,0
    for i=0:multiple[1]-1
        for j=0:multiple[2]-1
            for k=0:multiple[3]-1
                for a=1:size(basis,1)
                lattice[n,1] = (i+basis[a,1])*latparam
                lattice[n,2] = (j+basis[a,2])*latparam
                lattice[n,3] = (k+basis[a,3])*latparam
                n += 1
                end
            end
        end
    end
    return types,cell,lattice
end


function bccubic(multiple,latparam)
    # Generate face-centered cubic structure 
    # multiple should be a tuple like object
    
    basis = [0.0 0.0 0.0;
             0.5 0.5 0.5]
                        
    N = size(basis,1)*multiple[1]*multiple[2]*multiple[3]
    
    lattice = zeros(N,3)
    types = [ 1 for i=1:N]
    
    #TODO - reduce loops (i.e. vector operations)
    n = 1 #Atomic index
    #Origin: 0,0,0
    for i=0:multiple[1]-1
        for j=0:multiple[2]-1
            for k=0:multiple[3]-1
                for a=1:size(basis,1)
                lattice[n,1] = (i+basis[a,1])*latparam
                lattice[n,2] = (j+basis[a,2])*latparam
                lattice[n,3] = (k+basis[a,3])*latparam
                n += 1
                end
            end
        end
    end
    return types,lattice
end

function diamond(multiple,latparam)
    # Generate diamond cubic structure 
    # multiple should be a tuple like object
    
    basis = [0.0 0.0 0.0;
             0.5 0.0 0.5;
             0.0 0.5 0.5;
             0.5 0.5 0.0;
             0.25 0.25 0.25;
             0.75 0.75 0.25;
             0.25 0.75 0.75;
             0.75 0.25 0.75]
             
                        
    N = size(basis,1)*multiple[1]*multiple[2]*multiple[3]
    
    lattice = zeros(N,3)
    types = [ 1 for i=1:N]
    
    #TODO - reduce loops (i.e. vector operations)
    n = 1 #Atomic index
    #Origin: 0,0,0
    for i=0:multiple[1]-1
        for j=0:multiple[2]-1
            for k=0:multiple[3]-1
                for a=1:size(basis,1)
                lattice[n,1] = (i+basis[a,1])*latparam
                lattice[n,2] = (j+basis[a,2])*latparam
                lattice[n,3] = (k+basis[a,3])*latparam
                n += 1
                end
            end
        end
    end
    return types,lattice
end

end 