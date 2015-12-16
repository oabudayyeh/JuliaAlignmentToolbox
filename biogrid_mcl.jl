
f = open("biogrid.txt")
c = 0
d=Dict()
id = 1
for ln in eachline(f)
    entries = split(ln,"\t")
    a = entries[2]
    b = entries[3]
    s = entries[19]
    if c > 0
        if get(d,a, 0) == 0
            d[a] = id
            id +=1
        end
        if get(d,b, 0) == 0
            d[b] = id
            id +=1
        end
    end
    c+=1
end
close(f)

println(length(keys(d)))

f = open("biogrid.txt")
n = length(keys(d))
m = Array(Float64, (n,n))
c = 0
for ln in eachline(f)
    entries = split(ln,"\t")
    a = entries[2]
    b = entries[3] 
    s = entries[19]
    if c > 0 && length(s) > 1
        m[d[a],d[b]] = float(s)*-1
        m[d[b],d[a]] = float(s)*-1
    end
    c+=1
end
close(f)

addprocs(1)
println(nworkers())

include("parallel_library.jl")

M = add_diag(m,0)
M = normalize(M)
t0 = time()
Y = mcl(M)
t1 = time()
println("elapsed time for regular MCL: ", t1-t0, " seconds")




#c = get_clusters(Y)
#println(length(c))

for testing purposes to make it fit the varied # of processors.
s = 4000
m2 = hcat(m,fill(0.0000000001,(n,s-n)))
m2 = vcat(m2,fill(0.000000001,(s-n,s)))

println(size(m2))

cores = [40]
for p in cores
    println(p)
    rmprocs(workers());
    addprocs(p)
    include("parallel_library.jl")

    M = add_diag(m2,0)
    M = normalize(M)
    t0 = time()
    Y = mclp(M,p,2,2,25,1)
    t1 = time()
    println("elapsed time for parallel MM: ", t1-t0, " seconds")
    c = get_clusters(Y)
    println(length(c))


    #sa = SharedArray(Float64,(n,n),init = s -> s[localindexes(s)] = rand(length(localindexes(s))))
    M = add_diag(m2,0)
    M = normalize(M)
    sa = convert(SharedArray,M)
    t0 = time()
    Y = mclpshared(sa,p,2,2,25,1)
    t1 = time()
    println("elapsed time for shared parallel MM: ", t1-t0, " seconds")
    c = get_clusters(Y)
    println(length(c))

end


println(c)