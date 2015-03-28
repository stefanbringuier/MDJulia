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

module types

export Atoms, Cell, Sim, scale, unscale


#Pass vectors of fields
type Atoms 
    num::Int
    id::Vector{Int}
    typ::Vector{Int}
    pos::Array{Float64,2}
    spos::Array{Float64,2}
    vel::Array{Float64,2}
    forc::Array{Float64,2}
    Atoms(num,id,typ,pos,box) = new(num,id,typ,pos,
                                    scale(box,pos),zeros(num,3),zeros(num,3))
end

type Cell
    cell::Array{Float64,2}
    cut::Float64
    pbc::Bool
    vol::Float64
    Cell(cell,cut,pbc) = new(cell,cut,pbc,det(cell))
end

type Sim
    dt::Real
    dim::Int
end    

function scale(box::Cell,pos::Array{Float64,2})
    boxdinv  = inv(box.cell) #Invert matrix to scale coordinates
    spos = pos * boxdinv #Here the * operator is an array-array operator
    return spos
end

function unscale(box::Cell,spos::Array{Float64,2})
    pos = spos * box.cell #Here the * operator is an array-array operator
    return pos
end


end