type UpdateManager
    overseer::Overseer

    dt::Float64
    next::Float64
    blocksize::Int
    block_end::Int
    maxindex::Int
    index::Int

    function UpdateManager(overseer, dt, blocksize, itr)
        new(overseer, dt, 0.0, blocksize, blocksize, blocksize*itr, 1)
    end
end

overseer(u::UpdateManager)  = u.overseer
dt(u::UpdateManager)        = u.dt
next(u::UpdateManager)      = u.next
blocksize(u::UpdateManager) = u.blocksize
block_end(u::UpdateManager) = u.block_end
maxindex(u::UpdateManager)  = u.maxindex
index(u::UpdateManager)     = u.index

init_updater(::Explicit,  overseer, dt, n, itr) = UpdateManager(overseer, dt, 0, itr)
init_updater(::Uniform,   overseer, dt, n, itr) = UpdateManager(overseer, dt, n, itr)
init_updater(::Histogram, overseer, dt, n, itr) = UpdateManager(overseer, dt, itr, 1)

function update!(::Explicit, u, t)
    notify!(overseer(u), index(u), t)
    u.index += 1
end

function update!(::Uniform, u, t)
    while t >= next(u)
        i = index(u)
        if i > block_end(u); return; end
        notify!(overseer(u), i, next(u))
        u.next  += dt(u)
        u.index += 1
    end
end

function update!(::Histogram, u, t)
    return;
end

function final_update!(::Explicit, u, t)
    notify!(overseer(u), index(u), t)
    u.index += 1
end

function final_update!(::Uniform, u, t)
    while index(u) <= block_end(u)
        i = index(u)

        #if i > block_end(u); return; end
        notify!(overseer(u), i, next(u))
        u.next  += dt(u)
        u.index += 1
    end

    u.next = 0.0
    if block_end(u) < maxindex(u)
        u.block_end = block_end(u) + blocksize(u)
    end
end

# Safeguard on indices?
function final_update!(::Histogram, u, t)
    i = index(u)
    notify!(overseer(u), i)
    u.index += 1
end
