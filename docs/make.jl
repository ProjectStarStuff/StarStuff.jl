push!(LOAD_PATH,"../src")

using Documenter, StarStuff

makedocs(sitename="StarStuff.jl",
	  pages = ["Home" => "home.md"
		   ]
	  )

