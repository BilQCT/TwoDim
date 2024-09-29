# Utility functions for TwoDim code:

using Polymake

function nerve_face_map(s,d)
    """
    nerve_face_map(s, d)

    Compute the faces of an n-tuple of dit-strings using Equation (12) of SimpCont.

    # Arguments
    - `s`: A dit-string representing the input data.
    - `d`: An integer representing the dimension.

    # Output
    - `ds`: An array of faces of the n-tuple of dit-strings.

    # Description
    The `nerve_face_map` function computes the faces of an n-tuple of dit-strings using Equation (12) of SimpCont. It iterates over each dimension of the n-tuple and calculates the corresponding face based on the input dit-string `s` and the dimension `d`.

    # Example
    ```julia
    s = [0, 1, 0, 1]
    d = 2
    faces = nerve_face_map(s, d)
    """
    ds = [];
    for i in 1:(length(s)+1)
        di_s = [];
        # d0 face:
        if i == 1
            di_s = s[2:end]
        # d1 face:
        elseif i==2
            a = (s[1]+s[2]) % d;
            di_s = vcat([a],[s[j] for j in (i+1):length(s)]);
        # di face: 2 < i < n+1
        elseif (i>2) && (i < (length(s)+1))
            a = [s[j] for j in 1:(i-2)];
            b = [(s[i-1]+s[i])%d];
            if length(vcat(a,b)) == length(s)-1
                di_s = vcat(a,b)
            else
                c = [s[j] for j in (i+1):length(s)];
                di_s = vcat(vcat(a,b),c);
            end
        # dn face:
        else
            di_s = s[1:(end-1)]
        end
        #println(di_s)
        push!(ds,di_s);
    end
    return ds
end



function faces_in_simplex(X1,X2)
    """
    Determine all faces (edges) occurring in each higher simplex (triangle).

    # Arguments
    - `X1`: Two-dimensional simplicial set representing Delta_{n-1}.
    - `X2`: Two-dimensional simplicial set representing Delta_n.

    # Output
    An array of [face, simplex] pairs for each f in Delta_{n-1}.

    """
    F = [];
    for x1 in X1
        # isolate edge identifier:
        tau = x1[1]; f = [];

        # check if tau is a face of each triangle sigma:
        for x2 in X2
            sigma = x2[1]; faces = x2[2]; tau_in_sigma = findall(y->y==tau,faces);

            # if tau is a face, register it:
            if isempty(tau_in_sigma) == false
                for k in tau_in_sigma; push!(f,[k,sigma]) end;
            end
        end
        push!(F,f);
    end
    return F
end




function generate_dit_strings(d,N,n,s)
    """
    Recursive loop to construct all dit strings.

    This function generates all dit strings of length `N` with alphabet size `d` recursively.

    # Arguments
    - `d`: Number of outcomes (alphabet size).
    - `N`: Length of the dit strings.
    - `n`: Current depth of recursion.
    - `s`: Current list of dit strings.

    # Output
    An array containing all dit strings of length `N` with alphabet size `d`.

    """
    if n < N
        Zd = [(i-1) for i in 1:d] ; dit_strings = [];
        for i in 1:length(s)
            for x in Zd
                #println(i)
                push!(dit_strings,push!(copy(s[i]),x));
            end
        end
        s = dit_strings; n = n+1;
        generate_dit_strings(d,N,n,s)
    else
        return s
    end
end



function all_dit_strings(d,N)
    """
    Generate all dit strings.

    This function generates all dit strings of length `2N` with an alphabet size `d`.

    # Arguments
    - `d`: Number of outcomes (alphabet size).
    - `N`: Number of generators (dimension of the simplex).

    # Output
    An array containing all dit strings of length `2N` with an alphabet size `d`.

    """

    s = [[i-1] for i in 1:d]; n = 1;
    return generate_dit_strings(d,N,n,s)
end



function compatible_face_maps(d,s_k,i)

    """
    Find all dit-strings of length k+1 with face maps coinciding on target dit-string s_k.

    This function finds all dit-strings of length `k+1` with face maps coinciding on target dit-string `s_k` for the given `d` and `i`.

    # Arguments
    - `d`: Number of outcomes (alphabet size).
    - `s_k`: Target dit-string.
    - `i`: Index of the face map.

    # Output
    An array containing all dit-strings of length `k+1` with face maps coinciding on target dit-string `s_k`.

    """

    k = length(s_k); higher_bitstrings = all_dit_strings(d,k+1);

    # simplices sn of (NZd)_n such that di(sn) = sn_1
    S = [];
    for s in higher_bitstrings
        ds = nerve_face_map(s,d);
        if ds[i] == s_k
            push!(S,s)
        end
    end
    return S
end




function indicator_function(bit_string,faces)
    if bit_string in faces; return 1; else; return 0; end;
end



function polytope_from_inequalities(A::Matrix{Int64})
    """
    inequalities_to_polytope(A::Matrix{Int64}) -> Polytope

    Create a polytope object from a matrix of inequalities.

    # Arguments
    - `A::Matrix{Int64}`: A matrix representing a system of linear inequalities, where each row
    corresponds to an inequality.

    # Output
    - `Polytope`: A polymake polytope object constructed from the input inequalities.

    # Example
    julia
    A = [1 0 0;
    0 1 0;
    0 0 1]
    polytope = inequalities_to_polytope(A)
    """
    return pm.polytope.Polytope(INEQUALITIES = A)
end

function polytope_from_inequalities_and_equations(A::Matrix{Int64},E::Matrix{Int64})
    """
    polytope_from_inequalities_and_equations(A::Matrix{Int64},E::Matrix{Int64}) -> Polytope

    Create a polytope object from a matrix of inequalities.

    # Arguments
    - `A::Matrix{Int64}`: A matrix representing a system of linear inequalities, where each row
    corresponds to an inequality.

    - `E::Matrix{Int64}`: A matrix representing a system of linear equalities.

    # Output
    - `Polytope`: A polymake polytope object constructed from the input inequalities and equations.

    """
    return pm.polytope.Polytope(INEQUALITIES = A, EQUATIONS = E)
end



