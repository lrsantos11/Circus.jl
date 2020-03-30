using Documenter, Circus

makedocs(;
    modules=[Circus],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/lrsantos11/Circus.jl/blob/{commit}{path}#L{line}",
    sitename="Circus.jl",
    authors="["Roger Behling", "Yunier Bello-Cruz", "Luiz-Rafael Santos", "Hugo Lara Urdaneta", "Davi Sales Barreira"]",
    assets=String[],
)

deploydocs(;
    repo="github.com/lrsantos11/Circus.jl",
)
