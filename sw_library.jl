function print_matrix(matrix,seq1,seq2)
    seq1 = "^$seq1"
    seq2 = "^$seq2"
    println('\t',join(collect(seq1),'\t'))
    for i =1:length(seq2)
        println(seq2[i],'\t',join(matrix[i,:],'\t'))    
    end
end
    
function initialize_path(seq1,seq2)
    col = length(seq1)
    row = length(seq2)
    matrix = fill(0,(row,col)) 
    path = fill("N",(row,col))
    return matrix,path
end


#For Serial implementation
function SW(seq1,seq2)
    indel = -1
    match = 2
    
    seq1 = "^$seq1"
    seq2 = "^$seq2"
    
    col = length(seq1)
    row = length(seq2)
    
    matrix,path = initialize_path(seq1,seq2)
    
    pathvals = ["-","|","M"]
    
    for i = 2:row
        for j = 2:col
            scores = []
            
            push!(scores,matrix[i,j-1]+indel)
            push!(scores,matrix[i-1,j]+indel)
            if seq1[j] == seq2[i]
                push!(scores,matrix[i-1,j-1]+match)
            else
                push!(scores,matrix[i-1,j-1]+indel)
            end
            
            val,ind = findmax(scores)
            if val < 0
                matrix[i,j] = 0
            else
                matrix[i,j] = val
            end
            
            path[i,j] = pathvals[ind]
            
        end
    end
    return matrix,path
end

function traceback(matrix,path,seq1,seq2)
    seq1 = "^$seq1"
    seq2 = "^$seq2"
    
    irow = size(matrix)[1]
    jcol = size(matrix)[2]
    
    align1 = []
    align2 = []
    
    maxval,ind = findmax(matrix[irow,:])
    
    jcol = ind
    check = false
    while !check
        if irow == 1 && jcol == 1
            break
        elseif path[irow,jcol] == "M"
            push!(align1,seq1[jcol])
            push!(align2,seq2[irow])
            irow -= 1
            jcol -= 1
        elseif path[irow,jcol] == "-"
            push!(align1,seq1[jcol])
            push!(align2,"-")
            jcol -= 1
        elseif path[irow,jcol] == "|"
            push!(align1,"-")
            push!(align2,seq2[irow])
            irow -= 1
        elseif path[irow,jcol] == "N"
            if irow > 1
                push!(align2,seq2[irow])
                irow -=1
            end
            if jcol > 1
                push!(align1,seq1[jcol])
                jcol -=1
            end
        end
    end
    
    return join(align1[end:-1:1]),join(align2[end:-1:1])
end


#For shared parallel implementation

@everywhere function shared_get_score!(matrix,path,equal,indel,match,i,j,scores,val,ind)
    @inbounds scores[1] = matrix[i,j-1]+indel
    @inbounds scores[2] = matrix[i-1,j]+indel
    if equal
        @inbounds scores[3] = matrix[i-1,j-1]+match
    else
        @inbounds scores[3] = matrix[i-1,j-1]+indel
    end
    
    val,ind = findmax(scores)
    if val < 0
        @inbounds matrix[i,j] = 0
    else
        @inbounds matrix[i,j] = val
    end
    
    @inbounds path[i,j] = ind
end

function shared_initialize_path(seq1,seq2)
    col = length(seq1)
    row = length(seq2)
    matrix = convert(SharedArray,fill(0,(row,col))) 
    path = convert(SharedArray,fill(0,(row,col)))
    return matrix,path
end


function spSW(seq1,seq2,p)
    indel = -1
    match = 2
    
    pathvalscode = ["-","|","M"]
    pathvals = [1,2,3]

    val = 0
    ind = 0

    seq1 = "^$seq1"
    seq2 = "^$seq2"
    
    scores = convert(SharedArray,fill(0,3))

    col = length(seq1)
    row = length(seq2)
    
    wl = workers()

    matrix,path = shared_initialize_path(seq1,seq2)

    for j = 2:col
        jcol = j
        irow = 2
        @sync begin
            count = 1
            w = workers()
            while jcol > 1 && irow < row + 1
                    #println(j," ",irow," ",jcol)
                    @inbounds if seq1[jcol] == seq2[irow]
                        equal = true
                    else
                        equal = false
                    end
                    
                    if count > length(w)
                        count = 1
                    end

                    @inbounds @async remotecall_wait(w[count],shared_get_score!,matrix,path,equal,indel,match,irow,jcol,scores,val,ind)
                    jcol -= 1
                    irow += 1
                    count += 1
            end
        end
    end
        
        
    for i = 3:row
        jcol = col
        irow = i
        @sync begin
            count = 1
            w=workers()
            while irow < row+1 && jcol > 1
                    #println(j," ",irow," ",jcol)
                    @inbounds if seq1[jcol] == seq2[irow]
                        equal = true
                    else
                        equal = false
                    end
                    if count > length(w)
                        count = 1
                    end

                    @inbounds @async remotecall_wait(w[count],shared_get_score!,matrix,path,equal,indel,match,irow,jcol,scores,val,ind)
                    jcol -= 1
                    irow += 1
                    count += 1
            end
        end
    end
    return matrix,path
end






#For the parallel grid implementation of SW
@everywhere function grid_get_score!(matrix,path,seq1,seq2,indel,match,scores,startx,starty,endx,endy)
    for i = startx:endx
        for j = starty:endy            
            scores[1] = matrix[i,j-1]+indel
            scores[2] = matrix[i-1,j]+indel
            if seq1[j] == seq2[i]
                scores[3] = matrix[i-1,j-1]+match
            else
                scores[3] = matrix[i-1,j-1]+indel
            end
            
            val,ind = findmax(scores)
            if val < 0
                matrix[i,j] = 0
            else
                matrix[i,j] = val
            end
            
            path[i,j] = ind
            
        end
    end
end


function gridSW(seq1,seq2,p)
    indel = -1
    match = 2
    
    pathvalscode = ["-","|","M"]
    pathvals = [1,2,3]

    n = length(seq1)

    val = 0
    ind = 0

    seq1 = "^$seq1"
    seq2 = "^$seq2"
    
    scores = convert(SharedArray,fill(0,3))

    col = length(seq1)
    row = length(seq2)
    
    wl = workers()

    matrix,path = shared_initialize_path(seq1,seq2)

    lwl = length(wl)
    for j = 1:lwl
        jcol = j
        irow = 1
        @sync begin
            count = 1
            w = workers()
            while jcol > 0 && irow < lwl+1                 
                    if count > length(w)
                        count = 1
                    end

                    startx = 1+(irow-1)*div(n,p)+1
                    starty = 1+(jcol-1)*div(n,p)+1

                    endx = irow*div(n,p)+1 
                    endy = jcol*div(n,p)+1

                    @async remotecall_wait(w[count],grid_get_score!,matrix,path,seq1,seq2,indel,match,scores,startx,starty,endx,endy)
                    jcol -= 1
                    irow += 1
                    count += 1
            end
        end
    end
        
        
    for i = 2:lwl
        jcol = lwl
        irow = i
        @sync begin
            count = 1
            w=workers()
            while irow < lwl+1 && jcol > 0
                    if count > length(w)
                        count = 1
                    end
                    startx = 1+(irow-1)*div(n,p)+1
                    starty = 1+(jcol-1)*div(n,p)+1

                    endx = irow*div(n,p)+1 
                    endy = jcol*div(n,p)+1

                    @async remotecall_wait(w[count],grid_get_score!,matrix,path,seq1,seq2,indel,match,scores,startx,starty,endx,endy)
                    jcol -= 1
                    irow += 1
                    count += 1
            end
        end
    end
    return matrix,path
end
