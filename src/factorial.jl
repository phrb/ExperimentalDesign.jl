"""
$(TYPEDSIGNATURES)

Receives a tuple of arrays representing categorical factor levels, and returns a
`Base.Iterators.ProductIterator`.   This allows  full  factorial  designs to  be
arbitrarily large and  only be computed as needed.  To  compute an explicit full
factorial design, use [`explicit_fullfactorial`](@ref).

```jldoctest
julia> fullfactorial(Tuple([-1, 1] for i = 1:10))
Base.Iterators.ProductIterator{NTuple{10,Array{Int64,1}}}(([-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1]))

```
"""
function fullfactorial(factors::Tuple)
    Base.Iterators.product(factors...)
end

"""
$(TYPEDSIGNATURES)

Receives  a  `Base.Iterators.ProductIterator`  and  computes  an  explicit  full
factorial design. The generated array is exponentially large.

```jldoctest
julia> explicit_fullfactorial(fullfactorial(([-1, 1], [:a, :b, :c])))
6×2 Array{Any,2}:
 -1  :a
  1  :a
 -1  :b
  1  :b
 -1  :c
  1  :c

```
"""
function explicit_fullfactorial(iterator::Base.Iterators.ProductIterator)
    permutedims(foldl((x, y) -> hcat(collect(x), collect(y)), iterator))
end

"""
$(TYPEDSIGNATURES)

Receives a tuple of arrays  representing categorical factor levels, and computes
an explicit full factorial design. The generated array is exponentially large.

```jldoctest
julia> explicit_fullfactorial(([-1, 1], [:a, :b, :c], [1, 2]))
12×3 Array{Any,2}:
 -1  :a  1
  1  :a  1
 -1  :b  1
  1  :b  1
 -1  :c  1
  1  :c  1
 -1  :a  2
  1  :a  2
 -1  :b  2
  1  :b  2
 -1  :c  2
  1  :c  2

```
"""
function explicit_fullfactorial(factors::Tuple)
    explicit_fullfactorial(fullfactorial(factors))
end
