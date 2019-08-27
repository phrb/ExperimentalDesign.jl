"""
$(TYPEDSIGNATURES)

Receives a tuple of arrays representing categorical factor levels, and returns a
`Base.Iterators.ProductIterator`.   This allows  full  factorial  designs to  be
arbitrarily large and  only be computed as needed.  To  compute an explicit full
factorial design, use [`explicit_fullfactorial`](@ref).

"""
function fullfactorial(factors::Tuple)
    Base.Iterators.product(factors...)
end

"""
$(TYPEDSIGNATURES)

Receives  a  `Base.Iterators.ProductIterator`  and  computes  an  explicit  full
factorial design. The generated array is exponentially large.

"""
function explicit_fullfactorial(iterator::Base.Iterators.ProductIterator)
    permutedims(foldl((x, y) -> hcat(collect(x), collect(y)), iterator))
end

"""
$(TYPEDSIGNATURES)

Receives a tuple of arrays  representing categorical factor levels, and computes
an explicit full factorial design. The generated array is exponentially large.

"""
function explicit_fullfactorial(factors::Tuple)
    explicit_fullfactorial(fullfactorial(factors))
end
