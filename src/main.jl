################################################################################
################################################################################
################################ Main function #################################
################################################################################
################################################################################

include("edge.jl"); include("prob.jl")


using Polymake
const pm = Polymake


mutable struct TwoDimDist

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

    EdgeInequalities::Matrix{Int64}
    EdgePolytope::Polymake.LibPolymake.BigObjectAllocated

    function TwoDimDist(X::Vector{Vector}, d::Int64, T = [[]])
        if d == 2
            ProbInequalities = Prob(X,d,T).ProbInequalities
            ProbEquations = Prob(X,d,T).ProbEquations
            ProbPolytope = Prob(X,d,T).ProbPolytope
    
            EdgeInequalities = Edge(X,T).EdgeInequalities
            EdgePolytope = Edge(X,T).EdgePolytope
    
            new(ProbInequalities, ProbEquations,ProbPolytope,EdgeInequalities,EdgePolytope)
        elseif d > 2
            ProbInequalities = Prob(X,d,T).ProbInequalities
            ProbEquations = Prob(X,d,T).ProbEquations
            ProbPolytope = Prob(X,d,T).ProbPolytope
    
            new(ProbInequalities,ProbEquations,ProbPolytope)
        end
    end
    
end