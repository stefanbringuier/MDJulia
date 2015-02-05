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
module readlammps

#TODO - see if you can rewrite this in a better fashion
function readdump(FileName,SnapShots,SaveSnap,Nfield,scaled=false)
    ##read lammps dumpfile with header saved for specified file
    ##number of atoms cannot change
    ##format: id type x y z .....
    ##"""

    file = open(FileName)

    #Get number of atoms and box then rewind
    garb1 = readline(file)
    garb2 = readline(file)
    garb3 = readline(file)
    numatoms= int(readline(file))
    garb4 = readline(file)
    #Exception to see if tilt factor is included, only 1 line
    flag = false
    try
        testxlo,testxhi,testxy = float(split(readline(f)))
    catch
        flag = true
    end
    seek(file,0)
    
    #Start reading dumpfile
    data = zeros(numatoms,Nfield,SnapShots)
    t = 1
    while t <= SnapShots
        #read header
        h1 = readline(file)
        tmptime = readline(file)
        h2 = readline(file)
        tmpatoms = readline(file)
        h3 = readline(file)
        #Based on parsing of box
        if flag != true
            tmpxlo,tmpxhi,tmpxy = split(readline(file))
            tmpylo,tmpyhi,tmpxz = split(readline(file))
            tmpzlo,tmpzhi,tmpyz = split(readline(file))
        else
            tmpxlo,tmpxhi = split(readline(file))
            tmpylo,tmpyhi = split(readline(file))
            tmpzlo,tmpzhi = split(readline(file))
            tmpxy,tmpxz,tmpyz = 0.00,0.00,0.00
        end
        h4 = readline(file)

        #Save desired configuration box (tuple)
        if t == SaveSnap
            xlo_bound,xhi_bound,xy = tmpxlo,tmpxhi,tmpxy
            ylo_bound,yhi_bound,xz = tmpylo,tmpyhi,tmpxz
            zlo_bound,zhi_bound,yz = tmpzlo,tmpzhi,tmpyz
        end
        #Read atomic positions etc.
        for a = 1:numatoms
            #Read string -> strip '\n' char -> split into new list
            line = readline(file)
            data[a,:,t] = float(split(line))
        end
        t += 1
    end
    close(file)
    
    #only return data for SaveSnap
    datasave = data[:,:,SaveSnap]
    
    #Convert data types & LAMMPS Dump bounding box
    xlo_bound,xhi_bound = float(xlo_bound),float(xhi_bound)
    ylo_bound,yhi_bound = float(ylo_bound),float(yhi_bound)
    zlo_bound,zhi_bound = float(zlo_bound),float(zhi_bound)
    xy,xz,yz = float(xy),float(xz),float(yz)
 
    xlo  = xlo_bound - min(0.0,xy,xz,xy+xz)
    xhi  = xhi_bound - max(0.0,xy,xz,xy+xz)
    ylo  = ylo_bound - min(0.0,yz)
    yhi  = yhi_bound - max(0.0,yz)
    zlo  = zlo_bound
    zhi  = zhi_bound
    
    #    Shift box to origin 0,0,0
    if xlo < 0.00
       datasave[:,2] += abs(xlo)
    elseif ylo < 0.00
       datasave[:,3] += abs(ylo)
    elseif zlo < 0.00
       datasave[:,4] += abs(zlo)
    end
    
    #Define box lengths and tilts
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo

    #lbox = [lx,ly,lz]
    #ltilt = [xy,xz,yz]
    boxarry = [lx 0.0 0.0;
              0.0 ly 0.0;
              0.0 0.0 lz]
    #boxarry = boxarry.T
    println("Simulation Box Dimensions: ")
    println(boxarry)
    
    #If scaled coordinates unscale
    if scaled == true
        datasave[:,2] *= lx
        datasave[:,3] *= ly
        datasave[:,4] *= lz
    end

    return numatoms,boxarry,datasave
end

function readdata()
    println("Not implemented yet!")
end

end