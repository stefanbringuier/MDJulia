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

# This module is going to be the main driver of the program
# it will parse an input file and run analysis 

module main

function readinput()
    println("Not implemented yet!")
end

function process()
       #Flow:
       #flags is a dictionary
       #flags = readinput()
       # readfile(flags[Name],flags[Type])
       #neighborlist build
       # if flags["RDF"]
            #RDF()
             #Take atoms of types (x,y)
       #elseif flags["Angle"]
            #Angle()
        #elseif flags["Order"]
             #NiTi()
          #  end
         
     
end

#Run
#process(sys.ARGV[1])