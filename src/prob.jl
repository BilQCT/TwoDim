########################################################
using Polymake
const pm = Polymake

path = joinpath(dirname(@__DIR__),"libs/utilities.jl"); include(path)

# Methods for prob coordinates

# call Linear Algebra package:
using LinearAlgebra


function nonnegative(Xn,d)

    """
    Determine a matrix of inequalities for nonnegative constraints.

    This function generates a matrix of inequalities for nonnegative constraints based on the given simplices and alphabet size.

    # Arguments
    - `Xn`: Array of nondegenerate simplices.
    - `d`: Alphabet size.

    # Output
    A matrix of inequalities for nonnegative constraints.

    """

    # extract dimension of simplices:
    dim = (length(Xn[1][2])-1);

    # Determine number of probability parameters to marginalize over:
    L = length(Xn); N_prob = d^(dim); N = N_prob*L;

    # create matrix of inequalities:
    Id = Matrix{Int}(I,N,N); Z = (zeros(Int8,N));
    #concatenate
    A = hcat(Z,Id);             
    return A
end




function normalization_matrix(Xn,d)

    """
    Generate a matrix of normalization constraints.

    This function generates a matrix of normalization constraints based on the given two-dimensional simplicial set.

    # Arguments
    - `Xn`: Two-dimensional simplicial set.
    - `d`: Alphabet size.

    # Output
    A matrix of normalization constraints.

    """

    L = length(Xn); simplex_dim = [(length(Xn[i][2])-1) for i in 1:L];
    N_simplex = [d^x for x in simplex_dim]; N = sum(N_simplex);

    # initialize matrix of constraints:
    init = zeros(Int8, L, N); col = -ones(Int8,L);
    NORM = hcat(col,init);

    # set counter:
    m = 0; for i in range(1,L) n=m+N_simplex[i]; NORM[i,(m+2):n+1] .= 1; m=n;  end
    return NORM
end



function collapse_conditions(d::Int64,X2,T = [[]])

    """
    Generate a matrix of collapse conditions.

    This function generates a matrix of collapse conditions based on the given two-dimensional simplicial set and optional twisting.

    # Arguments
    - `d`: Alphabet size.
    - `X2`: Array of nondegenerate simplices whose dimension is given by the length of the array.
    - `T`: Array of twisted triangles. Defaults to an empty array.

    # Output
    A matrix of collapse conditions.

    """

    # Determine number of probability parameters to marginalize over:
    L = length(X2); N_prob = d^(2); N = N_prob*L;
    
    # Triangle indices:
    S = [X2[i][1] for i in 1:L]
    
    # Twisted triangles:
    if T == [[]]
        Ti = []; twist_dict = Dict(zip([0],[0]))
    else
        Ti = [T[i][1] for i in 1:length(T)];
        twist_dict = Dict(zip(Ti,[i for i in 1:length(T)]))
    end
    
    # initialize constraint matrix
    C = Array{Int64}(undef, 0, N+1);
    
    # generate all dit-strings of length 2:
    Zd_2 = all_dit_strings(d,2);
    # dictionary: map z in Zd_2 to {1,..,d^2}:
    dict = Dict(zip(Zd_2,[i for i in 1:(d^2)]))
    
    # search for degenerate edges:
    for i in 1:L
        E = X2[i][2]; s = X2[i][1];
        D = findall(x->x < 0,E);
        
        # Twisted faces:
        if s in Ti
            beta = T[twist_dict[s]][2];
            dk = T[twist_dict[s]][3]
        else
            dk = []; beta = 0;
        end
        
        
        # elements of Deg are faces of 2-simplex
        for dd in D
            if dd == dk; z = beta; else; z = 0; end;
            f = compatible_face_maps(d,[z],dd);
            f_perp = [z for z in Zd_2 if z ∉ f];
            for k in f_perp
                n = dict[k]; column = (d^2)*(i-1)+1+n;
                
                # create row
                c = zeros(Rational{Int64},(1,N+1));
                c[column] = 1; C = vcat(C,c)
            end
        end
    end
    return C
end




using Combinatorics

function two_dimensional_compatibility(X1, X2, d, T = [[]])


    """
    Generate a matrix of compatibility relationships for a two-dimensional simplicial set.

    This function generates a matrix of compatibility relationships based on the given two-dimensional simplicial set and optional twisting.

    # Arguments
    - `X1`: Array of nondegenerate edges.
    - `X2`: Array of nondegenerate triangles.
    - `d`: Alphabet size.
    - `T`: Array of twisted triangles. Defaults to an empty array.

    # Output
    A matrix of compatibility relationships.

    """

    # extract dimension of simplices:
    dim1 = (length(X1[1][2])-1); dim2 = (length(X2[1][2])-1);

    # Determine number of probability parameters to marginalize over:
    L = length(X2); N_prob = d^(dim2); N = N_prob*L;
    
    # Array of edges occuring in each triangle: [f_i,sigma_k]
    F = faces_in_simplex(X1,X2);

    # map simplex identifiers to position index in array (X2):
    dict = Dict(zip([x[1] for x in X2],[i for i in 1:length(X2)]));

    # NS condition whenever f in X1 occurs twice:
    # for all edges that appear more than once, we impose compatibility:
    NS = Array{Int64}(undef, 0, N+1); # initialize NS matrix
    for i in 1:length(F)
        # n-simplices (f[2]) for which the (n-1)-simplex (edge) (f[1]) is a face: 
        if length(F[i]) > 1
            # generate all pairs of simplices glued along same edge:
            F_combinations = combinations(F[i],2);
            for f in F_combinations
                # initialize row vector for constraint:
                init1 = zeros(Int64, 1, N+1); init2 = zeros(Int64, 1, N+1);
                
                # specify all n (n-1) simplices of NZ_d:
                lower_outcomes = all_dit_strings(d,dim1); higher_outcomes = all_dit_strings(d,dim2);

                # extract simplex idx:
                simplex_1 = f[1][2]; simplex_2= f[2][2];
                # extract face maps:
                di_1 = f[1][1]; di_2 = f[2][1];

                # twisted edges:
                if T == [[]]
                    idx1 = []; idx2 = [];
                else
                    idx1 = findall(x->x[1]==simplex_1,T); idx2 = findall(x->x[1]==simplex_2,T);
                end

                if isempty(idx1) == false
                    beta1 = T[idx1[1]][2]; face1 = T[idx1[1]][3];
                else
                    beta1 = 0; face1 = copy(di_1);
                end

                if isempty(idx2) == false
                    beta2 = T[idx2[1]][2]; face2 = T[idx2[1]][3];
                else
                    beta2 = 0; face2 = copy(di_2);
                end


                for s in lower_outcomes
                    # Determine parameters that enter into marginalization:
                    if (isempty(idx1)==false) && (face1 == di_1)
                        s_prime = [(a + beta1) % 2 for a in s];
                        s1 = compatible_face_maps(d,s_prime,di_1);
                    else
                        s1 = compatible_face_maps(d,s,di_1);
                    end
                    
                    if (isempty(idx2)==false) && (face2 == di_2)
                        s_prime = [(a + beta2) % 2 for a in s];
                        s2 = compatible_face_maps(d,s_prime,di_2);
                    else
                        s2 = compatible_face_maps(d,s,di_2);
                    end

                    eta1 = [indicator_function(x,s1) for x in higher_outcomes];
                    eta2 = [-indicator_function(x,s2) for x in higher_outcomes];

                    # define map of d^n length array into N x d^n length array:
                    sigma1 = f[1][2]; i1 = dict[sigma1]; sigma2 = f[2][2]; i2 = dict[sigma2];
                    END1 = i1*(N_prob)+1; INT1 = END1-(N_prob)+1;
                    END2 = i2*(N_prob)+1; INT2 = END2-(N_prob)+1;#

                    # define interval where constraint is non-trivial:
                    p1 = [i for i in INT1:END1]; p2 = [i for i in INT2:END2];
                    init1[p1] = eta1; init2[p2] = eta2; init = init1+init2;
                    NS = vcat(NS,init);
                end
            end
        end
    end
    return NS
end

########################################################


function prob_conditions(X,d, T = [[]])
    X1 = X[2]; X2 = X[3]
    # define matrix of nonnegativity inequalities:
    A = nonnegative(X2,d);
    N = normalization_matrix(X2,d);
    C = two_dimensional_compatibility(X1,X2,d,T);
    Co = collapse_conditions(d,X2,T)
    return A, vcat(N,C,Co)
end


########################################################
mutable struct Prob
    """
    Mutable struct for representing a probabilistic polytope.

    # Fields
    - `ProbInequalities`: Matrix of inequalities.
    - `ProbEquations`: Matrix of equations.
    - `ProbPolytope`: Big object allocated for the polytope.

    # Constructor
    - `Prob(X::Vector{Vector}, d::Int64, T::Vector{Vector{Int}} = [])`: Constructs a `Prob` object based on the given data.
        - `X`: Array of vectors representing the simplicial set.
        - `d`: Alphabet size.
        - `T`: Optional array of twisted triangles. Defaults to an empty array.
    """

    ProbInequalities::Matrix{Int64}
    ProbEquations::Matrix{Int64}
    ProbPolytope::Polymake.LibPolymake.BigObjectAllocated

    function Prob(X::Vector{Vector}, d::Int64 = 2, T = [[]])
        ProbInequalities, ProbEquations = prob_conditions(X, d, T)
        ProbPolytope = polytope_from_inequalities_and_equations(ProbInequalities, ProbEquations)
        new(ProbInequalities, ProbEquations, ProbPolytope)
    end
end

