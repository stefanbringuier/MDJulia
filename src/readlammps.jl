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

#FIXME: Need to check if I'm correctly shifting to origin.

#Docstring function should be supported in Julia v0.4
#@doc """ Read LAMMPS dumpfile format
# Inputs:
# FileName - file
# SnapShots - number of trajectories/frames (default=1)
# SaveSnap - which trajectory to save for processing (default=1)
# Scaled - units scaled (default=false)
#
# Outputs:
# numatoms - number of atoms
# boxarry - box dimensions as an array
# datasave - atom id type x y z for specified trajectory
#
#TODO - see if this can be rewritten so that SnapShots is not required
#""" ->
function readdump(FileName,SnapShots=1,SaveSnap=1,scaled=false)
    ##read lammps dumpfile with header saved for specified file
    ##number of atoms cannot change
    ##format: id type x y z .....
    ##"""

    file = open(FileName)

    #Parse first snapshot header get simulation info
    #i.e. number of atoms, box dimensions, number of columns
    #then rewind file 
    garb1 = readline(file)
    garb2 = readline(file)
    garb3 = readline(file)
    numatoms= int(readline(file))
    garb4 = readline(file)
    #Exception to see if tilt factor is included, only 1 line
    flag = false
    vec1 = readline(file)
    try
        testxlo,testxhi,testxy = float(split(vec1))
    catch
        flag = true
    end
    #Go to Column header line
    #skip(file,2)
    vec2 = readline(file)
    vec3 = readline(file)
    #ITEM ATOM: id type x y z .....
    nfield = size(split(readline(file)),1)-2
    #Go back to beginning
    seek(file,0)
    
    #Start reading dumpfile
    data = zeros(numatoms,nfield,SnapShots)
    t = 1
    #while t <= SnapShots
     while !eof(file)
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
              xy ly 0.0;
              xz yz lz]
    #boxarry = boxarry.T
    #println("Simulation Box Dimensions: ")
    #println(boxarry)
    
    #If scaled coordinates unscale
    if scaled == true
        datasave[:,2] *= lx
        datasave[:,3] *= ly
        datasave[:,4] *= lz
    end

    return numatoms,boxarry,datasave
end

#@doc """ Read lammps data structure file
#   Input:
#    FileName - file name
#   Output:
#   numatoms - number of atoms
#   boxarry - simulation box
#   data - atomic coordinates
#   NOTE: this is not a robust function to read
#   lammps data file and thus you must strictly 
#   adhere to the format of empty lines.
#   TODO: implement keyword parser            
# """ ->
function readdata(FileName)
    #format: id type x y z .....

    file = open(FileName)

    header = readline(file)
    empty1 = readline(file) #emtpy line
    numatoms = int(split(readline(file))[1])
    atomtypes = int(split(readline(file))[1])
    empty2 = readline(file) #emtpy line
    xlo,xhi = float(split(readline(file))[1:2])
    ylo,yhi = float(split(readline(file))[1:2])
    zlo,zhi = float(split(readline(file))[1:2])
    flag = false
    try
        xy,xz,yz = float(split(readline(file)))
    catch
        flag = true
    end
    if flag == true
        xy,xz,yz = 0.00 , 0.00 , 0.00
    else
        empty3 = readline(file)
    end

    #empty3 = readline(file)
    header1 = readline(file)
    empty4 =  readline(file)
    
    #read coordinates
    data = zeros(numatoms,5)
    for i=1:numatoms
        line = readline(file)
        data[i,:] = float(split(line))
    end 
    close(file)
 
    #FIXME: Forgot if data file needs this conversion
    #Convert to vectors a,b,c
    xlo  = xlo - min(0.0,xy,xz,xy+xz)
    xhi  = xhi - max(0.0,xy,xz,xy+xz)
    ylo  = ylo - min(0.0,yz)
    yhi  = yhi - max(0.0,yz)
    
    #FIXME: This  isn't correct for nonorthogonal boxes
    # Shift box to origin 0,0,0
    if xlo < 0.00
       data[:,2] += abs(xlo)
    elseif ylo < 0.00
       data[:,3] += abs(ylo)
    elseif zlo < 0.00
       data[:,4] += abs(zlo)
    end
    
    #Define box lengths and tilts
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo

    boxarry = [lx 0.0 0.0;
              xy ly 0.0;
              xz yz lz]
    
   
    return numatoms,boxarry,data
    
end

end