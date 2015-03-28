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

module RDF
using types

function rdf(atoms::Atoms,box::Cell,nbins)

    binsize = box.cut/nbins
    distance = zeros(nbins,1)
    distribution = zeros(nbins,1)
    for i=1:atoms.num-1
        for j=i+1:atoms.num
            sij = atoms.spos[j,:] - atoms.spos[i,:]
            sij -= round(sij)
            rij = unscale(box,sij)
            dij = √sum(rij.^2) #.^ element wise
            if  dij <= box.cut
                ibin = int(dij/binsize) + 1 #place in ith bin, use int to round TODO -check this
                distribution[ibin]  = distribution[ibin] + 1
            end
        end 
    end
    # Normalize and create proper distances
    factor = 2.0*box.vol/(4.0*pi*atoms.num^2*binsize)
    for i=1:nbins
        distance[i] = (i+0.5)*binsize;
        distribution[i] = factor*distribution[i]/(i*binsize)^2
    end

    return distance,distribution
end

function write_rdf(distance,distribution,tag)
    println("Not yet implemented!")
end

end
