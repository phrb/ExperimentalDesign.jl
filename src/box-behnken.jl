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

Constructs a Box-Behnken design with size `matrix_size`.

```jldoctest
julia> boxbehnken(4)
27×4 transpose(::Matrix{Float64}) with eltype Float64:
 -1.0  -1.0   0.0   0.0
  1.0  -1.0   0.0   0.0
 -1.0   1.0   0.0   0.0
  1.0   1.0   0.0   0.0
 -1.0   0.0  -1.0   0.0
  1.0   0.0  -1.0   0.0
 -1.0   0.0   1.0   0.0
  1.0   0.0   1.0   0.0
 -1.0   0.0   0.0  -1.0
  1.0   0.0   0.0  -1.0
 -1.0   0.0   0.0   1.0
  1.0   0.0   0.0   1.0
  0.0  -1.0  -1.0   0.0
  0.0   1.0  -1.0   0.0
  0.0  -1.0   1.0   0.0
  0.0   1.0   1.0   0.0
  0.0  -1.0   0.0  -1.0
  0.0   1.0   0.0  -1.0
  0.0  -1.0   0.0   1.0
  0.0   1.0   0.0   1.0
  0.0   0.0  -1.0  -1.0
  0.0   0.0   1.0  -1.0
  0.0   0.0  -1.0   1.0
  0.0   0.0   1.0   1.0
  0.0   0.0   0.0   0.0
  0.0   0.0   0.0   0.0
  0.0   0.0   0.0   0.0

```
"""
function boxbehnken(matrix_size::Int)
    boxbehnken(matrix_size, matrix_size)
end

"""
$(TYPEDSIGNATURES)

Constructs a Box-Behnken design with size `matrix_size` with specified number of center points.

```jldoctest
julia> boxbehnken(4,0)
  24×4 transpose(::Matrix{Float64}) with eltype Float64:
  -1.0  -1.0   0.0   0.0
   1.0  -1.0   0.0   0.0
  -1.0   1.0   0.0   0.0
   1.0   1.0   0.0   0.0
  -1.0   0.0  -1.0   0.0
   1.0   0.0  -1.0   0.0
  -1.0   0.0   1.0   0.0
   1.0   0.0   1.0   0.0
  -1.0   0.0   0.0  -1.0
   1.0   0.0   0.0  -1.0
  -1.0   0.0   0.0   1.0
   1.0   0.0   0.0   1.0
   0.0  -1.0  -1.0   0.0
   0.0   1.0  -1.0   0.0
   0.0  -1.0   1.0   0.0
   0.0   1.0   1.0   0.0
   0.0  -1.0   0.0  -1.0
   0.0   1.0   0.0  -1.0
   0.0  -1.0   0.0   1.0
   0.0   1.0   0.0   1.0
   0.0   0.0  -1.0  -1.0
   0.0   0.0   1.0  -1.0
   0.0   0.0  -1.0   1.0
   0.0   0.0   1.0   1.0

```
"""
function boxbehnken(matrix_size::Int, center::Int)
    @assert matrix_size>=3

    A_fact = explicit_fullfactorial(Tuple([-1,1] for i = 1:2))

    rows = floor(Int, (0.5*matrix_size*(matrix_size-1))*size(A_fact)[1])
    A = zeros(rows, matrix_size)

    l = 0
    for i in 1:matrix_size-1
        for j in i+1:matrix_size
            l = l +1
            A[max(0, (l - 1)*size(A_fact)[1])+1:l*size(A_fact)[1], i] = A_fact[:, 1]
            A[max(0, (l - 1)*size(A_fact)[1])+1:l*size(A_fact)[1], j] = A_fact[:, 2]
        end
    end

    if center == matrix_size
        if matrix_size <= 16
            points = [0, 0, 3, 3, 6, 6, 6, 8, 9, 10, 12, 12, 13, 14, 15, 16]
            center = points[matrix_size]
        end
    end

    A = transpose(hcat(transpose(A), transpose(zeros(center, matrix_size))))
end
