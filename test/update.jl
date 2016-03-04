# Setup
m = kendall(); m.reactions[:Birth].istracked = true

model    = Simulation(m)
spcs     = model.state
rxns     = model.rxns
sname    = model.sname
stracked = model.stracked
rname    = model.rname
rtracked = model.rtracked

# Test Observers
o1 = make_observers(Explicit(), sname, stracked, rname, rtracked, spcs, rxns, n, itr)
o2 = make_observers(Uniform(), sname, stracked, rname, rtracked, spcs, rxns, n, itr)
o3 = make_observers(Histogram(), sname, stracked, rname, rtracked, spcs, rxns, n, itr)

# Explicit output
u = init_updater(Explicit(), Observer[], 0.0, 0, 0)
u = init_updater(Explicit(), O)
# Uniform output

# Histrogram output
