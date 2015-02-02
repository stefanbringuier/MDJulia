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


function simplecubic(multiple,latparam)
    #Generate simple cubic structure 
    # multiple should be a tuple like object
    
    basis = [0.0 0.0 0.0] #Row vector 1x3
    N = multiple[1]*multiple[2]*multiple[3]
    
    lattice = zeros(N,3)
    types = [ 1 for i=1:N*length(basis)]
    
    n = 1 #Atomic index
    #Origin: 0,0,0
    for i=0:multiple[1]-1
        for j=0:multiple[2]-1
            for k=0:multiple[3]-1
                lattice[n,1] = (i+basis[1])*latparam
                lattice[n,2] = (j+basis[2])*latparam
                lattice[n,3] = (k+basis[3])*latparam
                n += 1
            end
        end
    end
    return types,lattice
end