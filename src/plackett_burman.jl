"""
$(TYPEDSIGNATURES)

Gets the next prime `p`, starting from `n`, for which `(p + offset) % divisor ==
0` holds.

```jldoctest
julia> next_offset_divisible_prime(3, 1, 4, 1000)
3

julia> next_offset_divisible_prime(5, 1, 4, 1000)
7

julia> next_offset_divisible_prime(4, 1, 4, 1000)
7

```
"""
function next_offset_divisible_prime(n::Int, offset::Int, divisor::Int, tries::Int)
    for i = 1:tries
        prime = nextprime(n, i)

        if (prime + offset) % divisor == 0
            return prime
        end
    end

    error("There are no primes from $n, within $tries tries, " *
          "divisible by $divisor with offset $offset")
end

"""
$(TYPEDSIGNATURES)

To check if  a given design is  a Plackett-Burman design, we must check for the
following properties, obtained in the original Plackett-Burman paper:

1. Each component is replicated at each  of its values the same number of times,
   that is, the sum of elements in each column is zero

2. Each  pair of components  occur together at  every combination of  values the
   same number of times, that is, the sum of each pair of columns will produce a
   column with the  same number of occurrences  of ``2`` and ``-2``,  and twice that
   number of occurrences of ``0``

> Plackett,  R.L. and  Burman, J.P.,  1946. The  design of  optimum multifactorial
> experiments. Biometrika, 33(4), pp.305-325.

```jldoctest
julia> isplackettburman(plackettburman(2))
true

julia> isplackettburman(plackettburman(4))
true

julia> isplackettburman(plackettburman(16))
true

julia> isplackettburman(rand(4,4))
ERROR: MethodError: no method matching isplackettburman(::Array{Float64,2})
[...]

julia> isplackettburman(rand(Int, 4,4))
false

```
"""
function isplackettburman(d::Matrix{Int})
    sum(sum(d, dims = 1)) == 0 && sum(sum(d' * d, dims = 1) / size(d, 1)) == size(d, 2)
end

"""
$(TYPEDSIGNATURES)

The [Paley  construction](https://en.wikipedia.org/wiki/Paley_construction) is a
method for constructing Hadamard matrices using finite fields.

```jldoctest
julia> paley(Matrix{Int}(undef, 8, 8))
8×8 Array{Int64,2}:
 -1   1   1  -1   1   1   1  -1
  1   1  -1   1   1   1  -1  -1
  1  -1   1   1   1  -1  -1   1
 -1   1   1   1  -1  -1   1   1
  1   1   1  -1  -1   1   1  -1
  1   1  -1  -1   1   1  -1   1
  1  -1  -1   1   1  -1   1   1
 -1  -1   1   1  -1   1   1   1

```
"""
function paley(matrix::Matrix{Int})
    dimension::Tuple{Int, Int} = size(matrix)

    starting_row::Int = (dimension[1] == dimension[2] ? 1 : 2)
    ending_row::Int = dimension[1]

    prime::Int = dimension[2]

    residues_length::Int = 1 + floor((prime - 1) / 2)
    residues::Array{Int, 1} = Array{Int, 1}(undef, residues_length)

    residues[1] = 0
    residues[2] = 1

    @inbounds for i = 2:(residues_length - 1)
        residues[i + 1] = (i * i) % prime
    end

    sort!(residues)

    @inbounds for i = starting_row:ending_row
       @inbounds for j = 1:prime
           offset_value::Int = (i + j - 1) % prime
           if offset_value in residues
               matrix[i, j] = -1
           else
               matrix[i, j] = 1
           end
       end
    end

    return matrix
end

"""
$(TYPEDSIGNATURES)

Constructs a Plackett-Burman  design with size `matrix_size` if  possible, or to
the closest, largest, number for which it is possible.

```jldoctest
julia> plackettburman(4)
8×7 Array{Int64,2}:
  1   1   1   1   1   1   1
 -1   1  -1   1   1  -1  -1
  1  -1   1   1  -1  -1  -1
 -1   1   1  -1  -1  -1   1
  1   1  -1  -1  -1   1  -1
  1  -1  -1  -1   1  -1   1
 -1  -1  -1   1  -1   1   1
 -1  -1   1  -1   1   1  -1

```
"""
function plackettburman(matrix_size::Int)
    p = next_offset_divisible_prime(matrix_size, 1, 4, 1000)

    A = Matrix{Int}(undef, p + 1, p)
    A[1, :] = ones(p)
    paley(A)
end
