"""
This code is based on the original work on the Scilab DoE package by the
following individuals
    Copyright (C) 2012 - 2013 - Michael Baudin
    Copyright (C) 2012 - Maria Christopoulou
    Copyright (C) 2010 - 2011 - INRIA - Michael Baudin
    Copyright (C) 2009 - Yann Collette
    Copyright (C) 2009 - CEA - Jean-Marc Martinez
and
    Copyright (C) 2014, Abraham D. Lee
who converted the package to Python.
"""

"""
$(TYPEDSIGNATURES)

Constructs a central composite design.

```jldoctest
julia> ccdesign(3)
22×3 Matrix{Float64}:
 -1.0      -1.0      -1.0
  1.0      -1.0      -1.0
 -1.0       1.0      -1.0
  1.0       1.0      -1.0
 -1.0      -1.0       1.0
  1.0      -1.0       1.0
 -1.0       1.0       1.0
  1.0       1.0       1.0
  0.0       0.0       0.0
  0.0       0.0       0.0
  ⋮
  1.82574   0.0       0.0
  0.0      -1.82574   0.0
  0.0       1.82574   0.0
  0.0       0.0      -1.82574
  0.0       0.0       1.82574
  0.0       0.0       0.0
  0.0       0.0       0.0
  0.0       0.0       0.0
  0.0       0.0       0.0

```
"""
function ccdesign(n::Int, center::Array{Int}=[4, 4], alpha::Symbol=:orthogonal, face::Symbol=:circumscribed)
    @assert n>1
    @assert alpha in [:orthogonal, :rotatable]
    @assert face in [:circumscribed, :inscribed, :faced]
    @assert length(center) == 2

    H2, a = star(n, alpha, center)
    if face == :inscribed
        H1 = explicit_fullfactorial(fullfactorial(Tuple(fill([-1, 1], n))))
        H1 = H1/a
        H2, a = star(n, :faced, [1,1])
    elseif face == :faced
        H2, a = star(n, :faced, [1,1])
        H1 = explicit_fullfactorial(fullfactorial(Tuple(fill([-1, 1], n))))
    elseif face == :circumscribed
        H1 = explicit_fullfactorial(fullfactorial(Tuple(fill([-1, 1], n))))
    else
        print("Invalid value of face")
    end

    C1 = zeros(center[1], n)
    C2 = zeros(center[2], n)

    H1 = vcat(H1, C1)
    H2 = vcat(H2, C2)
    H = vcat(H1, H2)
end

function star(n, alpha, center)
    if alpha  == :faced
        a = 1
    elseif alpha == :orthogonal
        nc = 2.0^n
        nco = center[1]
        na = 2.0*n
        nao = center[2]
        a = (n*(1 + nao/na)/(1+nco/nc))^0.5
    elseif alpha == :rotatable
        nc = 2.0^n
        a = nc^0.25
    else
        print("Invalid value of alpha")
    end

    H = zeros(2*n, n)
    for i in 1:n
         H[2*i-1:2*i, i] = [-1, 1]
    end
    H *= a
    return H, a
end
