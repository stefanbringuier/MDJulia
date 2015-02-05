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
module read

#TODO - see if you can rewrite this in a better fashion
function readdump(FileName,SnapShots,SaveSnap,Nfield,gzipf=False,scaled=False):
    ##read lammps dumpfile with header saved for specified file
    ##number of atoms cannot change
    ##format: id type x y z .....
    ##"""

    #if gzipf == True
    #    #File = gzip.open(FileName,'r')
    #else
        #File = open(FileName,'r')
   # end
    file = open(FileName)
    #Get number of atoms then rewind
    garb1 = readline(file)
    garb2 = readline(file)
    garb3 = readline(file)
    numatoms= int(readline(file))
    garb4 = readline(file)
    try
        testxlo,testxhi,testxy = float(split(readline(f)))
        Flag = 0
    catch
        Flag = 1
    end

    seek(file,0)

    data = zeros(numatoms,Nfield,SnapShots)
    
    t = 0
    while t < SnapShots
        #read header
        h1 = readline(file)
        tmptime = readline(file)
        h2 = readline(file)
        tmpatoms = readline(file)
        h3 = readline(file)
        if Flag != 1:
            tmpxlo,tmpxhi,tmpxy = readline(file).strip('\n').split()
            tmpylo,tmpyhi,tmpxz = readline(file).strip('\n').split()
            tmpzlo,tmpzhi,tmpyz = readline(file).strip('\n').split()        
        else:
            tmpxlo,tmpxhi = File.readline().strip('\n').split()
            tmpylo,tmpyhi = File.readline().strip('\n').split()
            tmpzlo,tmpzhi = File.readline().strip('\n').split()
            tmpxy,tmpxz,tmpyz = 0.00,0.00,0.00

        h4 = readline(file)

        if t == SaveSnap:
            xlo_bound,xhi_bound,xy = tmpxlo,tmpxhi,tmpxy
            ylo_bound,yhi_bound,xz = tmpylo,tmpyhi,tmpxz
            zlo_bound,zhi_bound,yz = tmpzlo,tmpzhi,tmpyz

        for a in range(numatoms):
            #Read string -> strip '\n' char -> split into new list
            line = File.readline().strip('\n').split()
            data[a,:,t] = line

        t += 1
    end
    File.close()
    
    #only use data for SaveSnap
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
    if xlo < 0.00:
       datasave[:,2] += abs(xlo)
    if ylo < 0.00:
       datasave[:,3] += abs(ylo)
    if zlo < 0.00:
       datasave[:,4] += abs(zlo)

    
    #Define box lengths and tilts
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo

    #lbox = [lx,ly,lz]
    #ltilt = [xy,xz,yz]
    boxarry = np.array([[lx,0.,0.],
                      [xy,ly,0.],
                      [xz,yz,lz]])
    #boxarry = boxarry.T
    print 'Simulation Box: '
    print boxarry
    
    #If scaled coordinates unscale
    if scaled == True:
        datasave[:,2] *= lx
        datasave[:,3] *= ly
        datasave[:,4] *= lz


    return numatoms,boxarry,datasave
end

function readdata()
    println("Not implemented yet!")
end

end