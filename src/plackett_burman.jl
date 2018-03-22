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

function is_plackett_burman(design::Matrix{Int})
    design_size = size(design, 2)

    for i = 1:design_size
        if sum(design[:, i]) != 0
            return false
        end
    end

    for i = 1:design_size
        for j = 1:design_size
            if i != j
                values = counts(design[:, i] + design[:, j])
                if length(values) != 5 ||
                   values[1] != values[5] ||
                   values[1] != values[3] / 2
                    return false
                end
            end
        end
    end

    return true
end

function binary_search(array::Array{Int, 1}, target::Int)
    n::Int = length(array)

    left::Int = 1
    right::Int = n

    @inbounds while left <= right
        middle = (left + right) >>> 1

        if array[middle] < target
            left = middle + 1
        elseif array[middle] > target
            right = middle - 1
        else
            return true
        end
    end

    return false
end

function paley(matrix::Matrix{Int})
    dimension::Tuple{Int, Int} = size(matrix)

    starting_row::Int = (dimension[1] == dimension[2] ? 1 : 2)
    ending_row::Int = dimension[1]

    prime::Int = dimension[2]

    residues_length::Int = 1 + floor((prime - 1) / 2)
    residues::Array{Int, 1} = Array{Int, 1}(residues_length)

    residues[1] = 0
    residues[2] = 1

    @inbounds for i = 2:(residues_length - 1)
        residues[i + 1] = (i * i) % prime
    end

    sort!(residues)

    @inbounds for i = starting_row:ending_row
       @inbounds for j = 1:prime
           offset_value::Int = (i + j - 1) % prime
           if binary_search(residues, offset_value)
               matrix[i, j] = -1
           else
               matrix[i, j] = 1
           end
       end
    end

    return matrix
end

function plackett_burman(matrix_size::Int)
    # Paley construction works for primes where (p + 1) % 4 == 0
    # Therefore, 'factor::Int = 4' and 'offset::Int = 1' in the call
    # to 'primes_divisible_offset':
    p = primes_divisible_offset(4, 1, matrix_size)[end]

    # The matrix for the Plackett-Burman design must have
    # 'p + 1' rows and 'p' columns, corresponding to
    # 'p + 1' experiments for 'p' factors.
    A = Matrix{Int}(p + 1, p)
    A[1, :] = ones(p)
    paley(A)
end
