# Functions for use with the sybil toolbox

#read in models
# must have three files: description, reactions, metabolites
readTSVmod(prefix = "")
# prints exchange reactions
findExchReact()
# change bounds of reactions
changeBounds(model, model["reaction_id"], lb = c(0), ub = c(1000))
# uptake reactions- didn't work with examples
uptReact(model)

# run FBA
opt <- optimizeProb(model, algorithm = "fba")
# get the value for the obj function
lp_obj(opt)
# to run minimize total flux (mtf), new obj function
optimizeProb(model, algorithm = "mtf", wtobj = mod_obj(opt))
# get the flux distribution of the mtf solution
getFluxDist(mtf)
# get flux distribution of exchange reactions
ex <- findExchReact(model)
fd <- getFluxDist(mtf, ex)
getNetFlux(fd)
# access the objective function of the model
mod_obj(mtf)

# gene knockouts
ko <- optimizeProb(model, gene = "gene_name", lb = 0, ub = 0)

# MOMA
# change algorithm = "lmoma" (linear) or "moma" (quadratic)

# ROOM
# regulatory on/off minimization, change algorithm = "room"

# simulating multiple knock-outs
oneGeneDel() # all single gene knockouts
doubleGeneDel() # all pairwise gene knockouts
geneDeletion() # simultaneous deletion of n genes

# FVA
fluxVar(model, percentage = 80)

# Summary of results from simulations
summaryOptsol(result, model)

# Class modelorg - store model framework
# Class optcol - store optimization solution
