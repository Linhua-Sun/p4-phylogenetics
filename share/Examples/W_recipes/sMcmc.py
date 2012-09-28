# Do an MCMC.
read("d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=1, val=0.2)
m = Mcmc(t, nChains=4, runNum=0, sampleInterval=1000, checkPointInterval=250000)
m.run(1000000)
func.summarizeMcmcPrams(skip=500)


