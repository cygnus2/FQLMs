using Logging

@info "Working with $(Threads.nthreads()) threads."

a = fill(0,(10))
Threads.@threads for i = 1:10
    a[i] = Threads.threadid()
end
println(a)
