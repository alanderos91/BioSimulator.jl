using Weave

repo_directory = joinpath(@__DIR__,"..")

dir = joinpath(repo_directory,"notes")
files = [
  "sample_path.jmd",
  "plotting_configurations.jmd",
  "profiling.jmd"
]

for file in files
  tmp = joinpath(dir, file)
  weave(tmp, doctype = "github", out_path=dir)
end
