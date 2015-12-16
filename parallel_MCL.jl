
#rmprocs(workers());
#p=1
#addprocs(p)
#println(nworkers())

include("parallel_library.jl")

#@everywhere using DistributedArrays

#println(Base.CPU_CORES)



sizes = [1600,2400,3200]

for size in sizes
    n = size
    println(n)
    rfull = rand(size,size)

    M = add_diag(rfull,1)
    M = normalize(M)
    t0 = time()
    Y = mcl(M)
    t1 = time()
    println("elapsed time for regular MCL: ", t1-t0, " seconds")
end

cores = [2,4,8,16,32, 40]
for p in cores
    println(p)
    rmprocs(workers());
    addprocs(p)
    include("parallel_library.jl")

    for size in sizes
        n = size
        println(size)
        rfull = rand(size,size)

        M = add_diag(rfull,1)
        M = normalize(M)
        t0 = time()
        Mout = mclp(M,p,2,2,25,1)
        t1 = time()
        println("elapsed time for parallel MM: ", t1-t0, " seconds")

        #sa = SharedArray(Float64,(n,n),init = s -> s[localindexes(s)] = rand(length(localindexes(s))))
        sa = convert(SharedArray,rfull)
        t0 = time()
        Mout = mclpshared(sa,p,2,2,25,1)
        t1 = time()
        println("elapsed time for shared parallel MM: ", t1-t0, " seconds")
    end

end



