########################################################
using Polymake
const pm = Polymake

path = joinpath(dirname(@__DIR__),"libs/utilities.jl"); include(path)

# Methods for edge coordinates

function edge_outcome_assignments(ab,T = [])

    """
    edge_outcome_assignments(ab, [T])

    Compute the outcome assignments for the given two-bit string (ab) and optional twisting parameters.

    # Arguments
    - `ab`: An array representing the two-bit string (ab) with elements 0 or 1.
    - `T`: Optional twisting parameters as an array `[beta, di]`, where `beta` is either 0 or 1, and `di` is one of 1, 2, or 3. Default is an empty array.

    # Output
    - An array of length three giving the outcome assignments for `[b, a + b + beta, a]`.

    # Description
    The `edge_outcome_assignments` function computes the outcome assignments for the given two-bit string `ab` and optional twisting parameters `T`. It calculates the assignments based on the simplicial distributions on a 2-simplex. If twisting parameters are provided, it applies the twisting to the outcome assignments accordingly.

    # Example
    ```julia
    ab = [0, 1]
    T = [1, 2]
    outcomes = edge_outcome_assignments(ab, T)
    """

    # default input no twisting, otherwise T=[beta,di] determines twisting:
    if isempty(T) == false; beta = T[1]; di = T[2]; else;    beta = 0; di = 1; end;
    
    # edge assignments with no twisting: (simplicial distributions on 2-simplex)
    edge_assignments = [ab[2],((ab[1]+ab[2]) % 2),ab[1]];

    # edge assignments with twisting:
    edge_assignments[di] = ((edge_assignments[di]+beta) % 2);

    return edge_assignments
end


########################################################



function edge_inequalities(X::Vector{Vector},T = [[]])

    """
    edges_to_inequalities(X::Vector{Vector{Int}}, T::Vector{Vector{Int}} = []) -> A

    Generate a matrix of inequalities from simplicial set and twists.

    # Arguments
    - `X::Vector{Vector{Int}}`: A vector of vectors representing vertices, edges, and triangles.
        - `X[1]`: Vector of vertices.
        - `X[2]`: Vector of edges.
        - `X[3]`: Vector of triangles.
    - `T::Vector{Vector{Int}}`: A vector of vectors representing twists.
        - Each twist is represented as `[twist_index, beta_value, face_index]`.

    # Output
    - `A::Matrix{Int64}`: A matrix of inequalities where each row corresponds to an inequality
    representing a linear combination of projectors over edges and their coefficients.
    """

    X1 = X[2]; X2 = X[3];
    # number of edges:
    N = length(X1); d = 2;
    # edges
    E = [X1[i][1] for i in 1:N]; dict = Dict(zip(E,[i for i in 1:N]));

    # extract twisted edges:
    if isempty(T[1]) == false; twisted_idx = [T[i][1] for i in 1:length(T)]; else; twisted_idx = []; end;
    
    # bit strings
    Z2_2 = [[i,j] for i in [0,1] for j in [0,1]];

    # initialize A matrix for P(A,b)
    A = zeros(Int64, 4*length(X2), N+1);

    # For each 2-simplex
    for i in 1:length(X2)
        # extract beta value and face to be applied to:
        if X2[i][1] in twisted_idx
            # find unique identifier: twist = [beta,di]:
            idx = (findall(x->x==X2[i][1],twisted_idx))[1]; twist = [T[idx][2],T[idx][3]];
        else
            twist = [];
        end
        
        
        # assignment of outcomes to edges: (ab) -> [b,a+b,a]:
        edge_assignments = [edge_outcome_assignments(ab,twist) for ab in Z2_2];
        # assign coefficients of edges:
        edge_coefficients = [[(-1)^s[1],(-1)^s[2],(-1)^s[3]] for s in edge_assignments];
        
        
        for j in 1:4
            # initialize row vector of zeros
            a = zeros(Int8,N+1);
            
            
            # row vector defined by sum of projectors (-1)^a*P1+(-1)^b*P2+(-1)^a+b*P3
            for k in 1:3
                P = zeros(Int64,N+1);
                e = X2[i][2][k];
                # if degenerate modify constant term:
                if e < 0
                    P[1] = P[1]+edge_coefficients[j][k];
                elseif e > 0
                    x = dict[e];
                    P[x+1] = edge_coefficients[j][k];
                end
                a = a + P;
            end
            
            # add constant term:
            a[1] = a[1]+1;
            
            # determine ordering or inequalities:
            row = 4*(i-1)+j;
            # append to matrix
            A[row,:] = a;
            
        end
    end
    
    return A
end


########################################################

mutable struct Edge
    """
    Structure to store a pauli strings, bit strings, and dictionary mapping between the two 

    Edge(X::Vector{Vector{Int}}, T::Vector{Vector{Int}} = [[]]) -> Edge

    Create an `Edge` object from a given set of inputs.

    # Arguments
    - `X::Vector{Vector{Int}}`: A vector of vectors representing the input data for constructing the `Edge`.
    - `T::Vector{Vector{Int}} = [[]]`: An optional vector of vectors representing additional input data for constructing the `Edge`.

    # Output
    - `Edge`: An `Edge` object containing the constructed inequalities and polytope.

    # Description
    The `Edge` constructor takes a set of input data and constructs an `Edge` object. It first generates inequalities from the input data using the `edges_to_inequalities` function, and then constructs a polytope from these inequalities using the `inequalities_to_polytope` function.

    # Example
    ```julia
    X = [[], [[1, 2], [2, 3], [3, 1]], [[1, 2, 3], [2, 3, 4]]]
    T = []
    edge = Edge(X, T)
    """

    EdgeInequalities::Matrix{Int64}
    EdgePolytope::Polymake.LibPolymake.BigObjectAllocated

    function Edge(X::Vector{Vector}, T = [[]])
        EdgeInequalities = edge_inequalities(X, T)
        EdgePolytope = polytope_from_inequalities(EdgeInequalities)
        new(EdgeInequalities, EdgePolytope)
    end
end
