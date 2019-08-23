"""
$(TYPEDSIGNATURES)
"""
function primes_divisible_offset(factor::Int,
                                 offset::Int,
                                 max_prime::Int)
    all_primes = primes(max_prime)
    selected_primes = Int[]

    sizehint!(selected_primes, length(all_primes))

    @inbounds for i = 1:length(all_primes)
        (all_primes[i] + offset) % factor == 0 && push!(selected_primes, all_primes[i])
    end

    return selected_primes
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

"""
function isplackettburman(d::Matrix{Int})
    sum(sum(d, dims = 1)) == 0 && sum(sum(d'*d, dims = 1) / size(d, 1)) == size(d, 2)
end

"""
$(TYPEDSIGNATURES)
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
the closest, smallest, number for which it is possible.

"""
function plackettburman(matrix_size::Int)
    # Paley construction works for primes where `(p +  1) % 4 == 0` Therefore, we must
    # use    `factor    =    4`    and    `offset   =    1`    in    the    call    to
    # [`primes_divisible_offset`](@ref).
    p = primes_divisible_offset(4, 1, matrix_size)[end]

    A = Matrix{Int}(undef, p + 1, p)
    A[1, :] = ones(p)
    paley(A)
end
