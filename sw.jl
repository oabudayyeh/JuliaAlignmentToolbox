# rmprocs(workers());
# addprocs(1)
# println(nworkers())

include("sw_library.jl")

#seq1 = "ATGCATGCATGC"
#seq2 = "ATGGGCATG"

bases = ["T","C","G","A"]

# t0 = time()
# matrix1,path = SW(seq1,seq2)
# t1 = time()
# println("elapsed time: ", t1-t0, " seconds")

# t0 = time()
# a1, a2 = traceback(matrix1,path,seq1,seq2)
# t1 = time()
# println("elapsed time: ", t1-t0, " seconds")


# print_matrix(matrix1,seq1,seq2)
# print_matrix(path,seq1,seq2)
# println(a1)
# println(a2)

n=2000
sizes = [40,200,1000,3000,5000,7500]
#sizes = [15000]
procs = [2,4,5,10,20,40]
#procs = [2,4]
for s in sizes
    println(s)
    seq1 = join(rand(bases,s))
    seq2 = join(rand(bases,s))

    t0 = time()
    matrix1,path = SW(seq1,seq2)
    t1 = time()
    println("elapsed time: ", t1-t0, " seconds")

    for p in procs
        rmprocs(workers());
        addprocs(p)
        println(nworkers())
        include("sw_library.jl")

        #println(p)

        t0 = time()
        matrix2,path = gridSW(seq1,seq2,p)
        t1 = time()
        println("elapsed time: ", t1-t0, " seconds")

        #print_matrix(matrix2,seq1,seq2)

        if matrix1 == matrix2
            println("YES")
        else
            println("NO")
            #print_matrix(matrix1,seq1,seq2)
            #print_matrix(matrix2,seq1,seq2)
        end
    end
end