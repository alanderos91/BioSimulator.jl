# Setup
init = [Species("X", 1, false)]
curr = [Species("Y", 0, true)]
# Copy init into curr, then change the population of curr[1]
reset!(curr,init)
curr[1].pop = 0
@test init[1].pop != curr[1].pop
