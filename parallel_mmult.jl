
#This is for testing parallel matrix multiplication 

# rmprocs(workers());
# p=32
 addprocs(1)
println(nworkers())

include("parallel_library.jl")

#@everywhere using DistributedArrays

#println(Base.CPU_CORES)


 sizes = [80]
for size in sizes
    n = size
    println(n)
    rfull = rand(size,size)


    t0 = time()
    b = rfull*rfull
    t1 = time()
    println("elapsed time for regular MM: ", t1-t0, " seconds")
end

M = add_diag(rfull2,1)
M = normalize(M)
M = distribute(M)


t0 = time()
Y = mcl(M)
t1 = time()
println("elapsed time: ", t1-t0, " seconds")

cores = [1,2,4,8,16,32,40]
cores = [40]
for p in cores
    println(p)
    rmprocs(workers());
    addprocs(p)
    include("parallel_library.jl")
    
    rfull = 

    for size in sizes
        n = size
        println(size)
        rfull = rand(size,size)

        t0 = time()
        Mout = p_m_mult(rfull,p)
        t1 = time()
        println("elapsed time for parallel MM: ", t1-t0, " seconds")

        sa = SharedArray(Float64,(n,n),init = s -> s[localindexes(s)] = rand(length(localindexes(s))))
        sb = SharedArray(Float64,(n,n),init = s -> s[localindexes(s)] = rand(length(localindexes(s))))


        sc = SharedArray(Float64,(n,n));
        t0 = time()
        Sout = sharedmult(n,p,sa,sb,sc)
        t1 = time()
        println("elapsed time for shared parallel MM: ", t1-t0, " seconds")
    end

end
# t0 = time()
# test = sa*sb
# t1 = time()
# println("elapsed time: ", t1-t0, " seconds")
