include("JF2012.jl")
include("Cash-Karp.jl")
include("Propa_JF2012.jl")
i = readdlm("./input.txt", ' ' , String);
alpha = Propa_JF2012(float(i[1]),float(i[2]),float(i[3]),float(i[4]),float(i[5]),1,true,false,false,20000.);