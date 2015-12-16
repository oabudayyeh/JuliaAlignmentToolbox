
#MCL non-parallelized 
function normalize(A)
    return(A./sum(A,1))
end

function inflate(A,r)
    return(normalize(A.^r))
end

function expand(A,p)
    return(A^p)
end

function add_diag(A,mult)
    return(A+mult*ones(size(A)[1],size(A)[2]))
end

function stop(M,ind)
    if (ind % 5) == 4
        m = maximum(M^2 -M) - minimum(M^2-M)
        if m < 10.0^-200
            return(true)
        end
    end
    return(false)
end

function mcl(M,p=2,r=2,maxl=25, mult=1)
    count = 0
    for i = 1:maxl
        println(i)
        M = inflate(M,r)
        M = expand(M,p)
        if stop(M,i)
            break
        end
        count = count + 1
    end
    println(string("# iterations ", count))
    return(M)
end

function get_clusters(A)
    row = size(A)[1]
    clusters = Array[]
    for i =1:row
        ind = find(A[i,:] .> 0)
        if !isempty(ind)
            ind = unshift!(ind,i)
            push!(clusters,ind)
        end
    end

    meanlen = 0
    count = 0
    singlets = 0
    for i = 1:length(clusters)
        meanlen += length(clusters[i])
        if length(clusters[i]) > 5
            count+= 1
        end
        if length(clusters[i]) ==2
            singlets+= 1
        end
    end

    meanlen /= length(clusters)

    println("Average cluster length: ", meanlen)
    println("# Clusters with >5 members: ",count)
    println("# Clusters with singlets: ",singlets)
    return clusters
end   

#MCL implemented with parallel functions
@everywhere function mult(M,v,i)
    return(M*v,i)
end

@everywhere function hammult(M,pow,i)
    return(M.*pow,i)
end

function p_m_mult(M,procs)
    threads = []
    a = Array(Float64, size(M)[1], size(M)[2])
    step = round(Int,size(M)[2]/procs)
    for i=1:step:size(M)[2]
        push!(threads, @spawn mult(M,M[:,i:(i+step-1)],i))
    end

    for thread in threads
        y, ind = fetch(thread)
        a[:,ind:(ind+step-1)] = y
    end

    return(a)
end

function p_m_hammult(M,procs,pow)
    threads = []
    a = Array(Float64, size(M)[1], size(M)[2])
    step = round(Int,size(M)[2]/procs)
    for i=1:step:size(M)[2]
        push!(threads, @spawn hammult(M[:,i:(i+step-1)],pow,i))
    end

    for thread in threads
        y, ind = fetch(thread)
        a[:,ind:(ind+step-1)] = y
    end

    return(a)
end


function expand2(A,p,procs)
    for i =1:(p-1)
        A = p_m_mult(A,procs)
    end
    return(A)
end


function pstop(M,ind,procs)
    if (ind % 5) == 4
        Mp = expand2(M,2,procs)
        m = maximum(Mp -M) - minimum(Mp-M)
        if m < 10.0^-200
            return(true)
        end
    end
    return(false)
end

function mclp(M,procs=7,p=2,r=2,maxl=25, mult=1)
    count = 0
    for i = 1:maxl
        println(i)
        
        M = inflate(M,r)
        M = expand2(M,p,procs)
        if pstop(M,i,procs)
            break
        end
        count = count + 1
    end
    println(string("# iterations ", count))
    return(M)
end

function expand3(A,p,procs)
    n = size(A)[1]
    sc = SharedArray(Float64,(n,n))
    B = sharedmult(n,procs,A,A,sc)
    
    for i =1:(p-2)
        B = sharedmult(n,procs,A,B,sc)
    end
    return(B)
end


function pstop2(M,ind,procs)
    if (ind % 5) == 4
        Mp = expand3(M,2,procs)
        m = maximum(Mp -M) - minimum(Mp-M)
        if m < 10.0^-200
            return(true)
        end
    end
    return(false)
end

function mclpshared(M,procs=7,p=2,r=2,maxl=25, mult=1)
    count = 0
    for i = 1:maxl
        #println(i)
        
        M = inflate(M,r)
        M = convert(SharedArray,M)
        M = expand3(M,p,procs)
        #println(M[1:10,1:10])
        if pstop2(M,i,procs)
            break
        end
        count = count + 1
    end
    println(string("# iterations ", count))
    return(M)
end





#Matrix multiplication for shared arrays
@everywhere function mymatmul!(n,w,sa,sb,sc,p)
    range = 1+(w-1) * div(n,p) : (w) * div(n,p)
    sc[:,range] = sa[:,:] * sb[:,range]
end


function sharedmult(n,p,sa,sb,sc)
	@sync begin
	    for (i,w) in enumerate(workers())
	        @async remotecall_wait(w, mymatmul!, n, i, sa, sb, sc,p)
	    end
	end
	return sc
end


@everywhere function mymathamm!(n,w,sa,sc,p,pow)
    range = 1+(w-1) * div(n,p) : (w) * div(n,p)
    sc[:,range] = sa[:,range].^pow
end


function sharedhamm(n,p,sa,sc,pow)
	@sync begin
	    for (i,w) in enumerate(workers())
	        @async remotecall_wait(w, mymathamm!, n, i, sa, sc,p,pow)
	    end
	end
	return sc
end





