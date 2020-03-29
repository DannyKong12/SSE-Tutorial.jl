struct Lattice{D, T}
    dims::Dims{D}
    sites::Array{CartesianIndex{D}}
    neighbours::Array{CartesianIndex{D}, T}
    spins::Array{Int64, D}
    function Lattice(N::Int64)
        dims = (N,N)
        z = length(dims)
        sites = CartesianIndices(dims)
        neighbours = cat(circshift(sites, (-1,  0)),
                         circshift(sites, ( 1,  0)),
                         circshift(sites, ( 0, -1)),
                         circshift(sites, ( 0,  1)),
                         dims = z + 1)
        spins = rand(-1:2:1, dims)
        new{z, z+1}(dims, sites, neighbours, spins)
    end
end

mutable struct Operator
    op::Symbol
    bond::Int64
    # TODO: change struct wrapping to avoid mutable
end

struct SSE{D}
    N::Int64             # System size Lx=Ly=N
    Nb::Int64            # # of bonds, in 1D N=Nb but doesn't generalize to higher dimensions
    M::Base.RefValue{Int64}   # Taylor Expansion cutoff point (can we do better than ref/refvalue?)
    n::Base.RefValue{Int64}             # # of operators in list
    β::Float64           # thermodynamic beta
    ops::Array{Operator}
    lattice::Lattice
    vertexlist::Array{Int64}
    function SSE(N::Int64, β::Float64)
        lattice = Lattice(N)
        Nb = 2*N*N
        M = Int64(round(Nb*1.5*β))
        z = N*N
        ops = [Operator(:null,0) for i=1:M]
        vertexlist = fill(-1, M*4)
        new{z}(N, Nb, Ref(M), Ref(0), β, ops, lattice, vertexlist)
    end
end

function getBonds(index::Int64, sites::Array{CartesianIndex{2}, 2})
    (x,y) = size(sites)
    if index > x*y
        return (index, sites[index-x^2], sites[mod1(index-x^2+x,x*x)])
    else
        return (index, sites[index], sites[x*div(index-1,x)+mod1(index+1,x)])
    end
end

function next(v::Int64)
    return xor(v-1,1)+1
end

function diagonalUpdate!(state::SSE)
    for i in state.ops
        if i.op == :diag
            P = state.Nb*state.β/2/(state.M.x-state.n.x)
            if rand() < P
                i.op = :null
                i.bond = 0
                state.n.x -= 1
            end
        elseif i.op == :null
            index = rand(1:state.Nb)
            bond = getBonds(index, state.lattice.sites)
            if state.lattice.spins[bond[2]] != state.lattice.spins[bond[3]]
                P = 2*(state.M.x-state.n.x+1)/state.Nb/state.β
                if rand() < P
                    i.op = :diag
                    i.bond = index
                    state.n.x += 1
                end
            end
        elseif i.op == :offdiag
            # propogate spins
            bond = getBonds(i.bond, state.lattice.sites)
            state.lattice.spins[bond[2]] *= -1
            state.lattice.spins[bond[3]] *= -1
        end
    end
end

function buildVertexList(state::SSE)
    lastVertex = fill(-1,state.N,state.N)
    firstVertex = fill(-1,state.N,state.N)
    fill!(state.vertexlist, -1)
    for (i,op) in enumerate(state.ops)
        if op.op != :null
            v0 = 4*(i-1)+1
            b = getBonds(op.bond, state.lattice.sites)
            v1 = lastVertex[b[2]]
            v2 = lastVertex[b[3]]
            if v1 == -1
                firstVertex[b[2]] = v0
            else
                state.vertexlist[v0] = v1
                state.vertexlist[v1] = v0
            end
            if v2 == -1
                firstVertex[b[3]] = v0 + 1
            else
                state.vertexlist[v0+1] = v2
                state.vertexlist[v2] = v0+1
            end
            lastVertex[b[2]] = v0+2
            lastVertex[b[3]] = v0+3
        end
    end
    for site in 1:state.N*state.N
        f = firstVertex[site]
        l = lastVertex[site]
        if f != -1
            state.vertexlist[ f ] = l
            state.vertexlist[ l ] = f
        end
    end
end


function loopUpdate!(state::SSE)
    v0 = 1
    for (i,v) in enumerate(state.vertexlist)
        if v != -1
            v0 = i
            println(state.ops[div(i-1,4)+1])
            break
        end
    end
    flip = true # rand(Bool)
    v1 = next(v0)
    v1 = state.vertexlist[v1]
    while v1 != v0
        if flip
            if state.ops[div(v1-1,4)+1].op == :diag
                state.ops[div(v1-1,4)+1].op = :offdiag
            else
                state.ops[div(v1-1,4)+1].op = :diag
            end
        end
        v1 = next(v1)
        v1 = state.vertexlist[v1]
    end
    if flip
        if state.ops[div(v0-1,4)+1].op == :diag
            state.ops[div(v0-1,4)+1].op = :offdiag
        else
            state.ops[div(v0-1,4)+1].op = :diag
        end
    end
end


function MCSweep!(s::SSE)
    diagonalUpdate!(s)
    buildVertexList(s)
    loopUpdate!(s)
end

s = SSE(10, 10.0)

for i in 1:10
    MCSweep!(s)
end

for (i,op) in enumerate(s.ops)
    if op.op != :asdf
        print(i, "\t")
        for j in 1:4
            print(s.vertexlist[4*(i-1)+j], "\t")
        end
        println("\t", op.op)
    end
end
