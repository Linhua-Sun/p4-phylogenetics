#os.system("rm -f mcmc*")
read('inTrees.phy')
stm = STMcmc(var.trees, runNum=0, nChains=1, sampleInterval=8, checkPointInterval=4000, defaultBeta=2.0)
stm.run(8000)
stm = STMcmc(var.trees, runNum=1, nChains=1, sampleInterval=8, checkPointInterval=4000, defaultBeta=2.0)
stm.run(8000)