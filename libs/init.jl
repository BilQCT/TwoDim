# initialize project:
current_dir = @__DIR__
src_path = joinpath(dirname(current_dir),"src");
#proj_path = joinpath(dirname(dirname(current_dir)),"MyProject")

#using Pkg; Pkg.activate(proj_path);;

# include auxiliary Julia scripts:
for file in readdir(src_path)
    # Check if the file ends with ".jl" extension
    if endswith(file, ".jl")
        include(joinpath(src_path, file))
    end
end
