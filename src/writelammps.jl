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
module writelammps

function writedump(savestp,box,types,xyzcoord;param=false,file="dump.lammps")
    ### box - LAMMPS format box
    ###    lx,0,0
    ###    xy,ly,0
    ###    xz,yz,lz

    # Number of atoms
    natoms = (size(xyzcoord,1))
    # Sequential IDs
    ids = [1:natoms]

    #Always origin 0,0,0
    xlo,ylo,zlo = 0.00,0.00,0.00
    xhi,yhi,zhi = box[1,1],box[2,2],box[3,3]
    xy,xz,yz = box[2,1],box[3,1],box[3,2]
    xlo_bound  = xlo + min(0.0,xy,xz,xy+xz)
    xhi_bound  = xhi + max(0.0,xy,xz,xy+xz)
    ylo_bound  = ylo + min(0.0,yz)
    yhi_bound  = yhi + max(0.0,yz)
    zlo_bound  = zlo
    zhi_bound  = zhi
    

    out = open(file,"w")
    write(out,"ITEM: TIMESTEP \n")
    write(out,"$savestp \n")
    write(out,"ITEM: NUMBER OF ATOMS \n")
    write(out,"$natoms"," \n")
    write(out,"ITEM: BOX BOUNDS xy xz yz pp pp pp \n")
    write(out,"$xlo_bound $xhi_bound $xy \n")
    write(out,"$ylo_bound $yhi_bound $xz \n")
    write(out,"$zlo_bound $zhi_bound $yz \n")
    write(out,"ITEM: ATOMS id type x y z \n")
    
    #NOTE - this can be done in one line with writedlm 
    #however this would require hcat which  might be expensive
    #when the number of atoms is large
    if natoms <= 1000
        writedlm(out,hcat(ids,types,xyzcoord))
    else
        for i=1:natoms
            x,y,z = xyzcoord[i,1],xyzcoord[i,2],xyzcoord[i,3]
            write(out,"$(ids[i]) $(types[i]) x y z \n")
        end
    end
    close(out)
end

end