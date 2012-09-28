# This is STMcmc, for super tree mcmc.
# Started 18 March 2011, first commit 22 March 2011.

import pf,func
from Var import var
import math,random,string,sys,time,copy,os,cPickle,types,numpy,glob
from Glitch import Glitch
from TreePartitions import TreePartitions
from Constraints import Constraints
import datetime

haveLoadedP4stm = False
try:
    import p4stm
    import pyublas # not explicitly used--but makes converters available
    haveLoadedP4stm = True
except ImportError:
    pass


def bitReduce(bk, txBits, lLen, sLen, allOnes):
    #print "bitReduce: bk %i, txBits %i, lLen %i, sLen %i, allOnes %i" % (bk, txBits, lLen, sLen, allOnes)
    newBk = 0L
    counter = 0
    pops = 0
    for pos in range(lLen):
        tester = 1L << pos
        #print "pos %2i, tester: %3i" % (pos, tester)
        if tester & txBits:
            #print "    tester & txBits -- True"
            if tester & bk:
                adder = 1L << counter
                #print "        adding:", adder
                newBk += adder
                pops += 1
            else:
                #print "        not adding"
                pass
            counter += 1
    if (1 & newBk):
        #print "flipping"
        newBk = allOnes ^ newBk
        pops = sLen - pops
    #print "returning newBk %i, pops %i" % (newBk, pops)
    return newBk, pops

if 0: # test bitReduce
    sk = 6   # always at least 2 bits, even
    txBits = 30
    lLen = 5
    sLen = 4
    allOnes = 15
    print "     sk: %3i  %s" % (sk, func.getSplitStringFromKey(sk, lLen))
    print "taxBits: %3i  %s" % (txBits, func.getSplitStringFromKey(txBits, lLen))
                          
    rsk, popcount = bitReduce(sk, txBits, lLen, sLen, allOnes)
    print "    rsk: %3i  %s" % (rsk, func.getSplitStringFromKey(rsk, sLen))

def maskedSymmetricDifference(skk, skkSet, taxBits, longLen, shortLen, allOnes):
    if 0:
        print "skk (from the current tree)"
        for sk in skk:
            print func.getSplitStringFromKey(sk, longLen)
        print "skkSet (from input tree)"
        for sk in skkSet:
            print func.getSplitStringFromKey(sk, shortLen)
        print "taxBits:", taxBits, func.getSplitStringFromKey(taxBits, longLen)
            
    newSkk = []
    for sk in skk:
        reducedSk, popcount = bitReduce(sk, taxBits, longLen, shortLen, allOnes)
        #print "taxBits: %s  " % func.getSplitStringFromKey(taxBits, longLen),
        #print "%4i %s  " % (sk, func.getSplitStringFromKey(sk, longLen)),
        #print "%4i %s , %i" % (reducedSk, func.getSplitStringFromKey(reducedSk, shortLen), popcount)
        if popcount in [0, 1, shortLen - 1]:
            pass
        else:
            newSkk.append(reducedSk)
    newSkkSet = set(newSkk)
    #print newSkkSet, skkSet
    return len(newSkkSet.symmetric_difference(skkSet))

def slowQuartetDistance(st, inputTree):
    dst = st.dupe()
    toRemove = []
    for n in dst.iterLeavesNoRoot():
        if n.name not in inputTree.taxNames:
            toRemove.append(n)
    for n in toRemove:
        dst.removeNode(n)
    qd = dst.topologyDistance(inputTree, metric='scqdist')
    return qd
    

class STChain(object):

    def __init__(self, aSTMcmc):
        self.stMcmc = aSTMcmc
        self.tempNum = -1 # 'temp'erature, not 'temp'orary

        self.curTree = aSTMcmc.tree.dupe()
        self.propTree = aSTMcmc.tree.dupe()

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

        self.getTreeLogLike(self.curTree)  
        self.getTreeLogLike(self.propTree)

        self.stm = None
        if self.stMcmc.usingP4stmModule:
            self.startStm()

        #print "STChain init()"
        #self.curTree.draw()
        #print "logLike is %f" % self.curTree.logLike

    def startStm(self):
        #print "usingP4stmModule"
        self.stm = p4stm.Stm(len(self.stMcmc.taxNames))
        self.bigTr = self.stm.setBigT(len(self.propTree.nodes), self.propTree.nTax, self.propTree.postOrder)

        for n in self.propTree.nodes:
            if n.parent:
                self.bigTr.setParent(n.nodeNum, n.parent.nodeNum)
            if n.leftChild:
                self.bigTr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
            else:
                self.bigTr.setNodeTaxNum(n.nodeNum, self.stMcmc.taxNames.index(n.name))
            if n.sibling:
                self.bigTr.setSibling(n.nodeNum, n.sibling.nodeNum)

        #stm.dump()
        if 1:
            for t in self.stMcmc.trees:
                tr = self.stm.appendInTree(len(t.nodes), t.nTax, t.postOrder, t.beta)
                for n in t.nodes:
                    if n.parent:
                        tr.setParent(n.nodeNum, n.parent.nodeNum)
                    if n.leftChild:
                        tr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
                    else:
                        tr.setNodeTaxNum(n.nodeNum, self.stMcmc.taxNames.index(n.name))
                    if n.sibling:
                        tr.setSibling(n.nodeNum, n.sibling.nodeNum)
        self.stm.setInTreeTaxBits()
        self.stm.setInTreeInternalBits()
        self.stm.maybeFlipInTreeBits()
        self.stm.setBigTInternalBits()
        

    def getTreeLogLike(self, aTree):
        aTree.makeSplitKeys()
        aTree.skk = [n.br.splitKey for n in aTree.iterInternalsNoRoot()]
        aTree.logLike = 0.0
        for t in self.stMcmc.trees:
            if self.stMcmc.dMetric == 'sd':
                thisDist = t.beta *  maskedSymmetricDifference(aTree.skk, t.skSet, t.taxBits, self.stMcmc.nTax, t.nTax, t.allOnes)
            elif self.stMcmc.dMetric == 'scqdist':
                thisDist = t.beta * slowQuartetDistance(self.propTree, t)
            else:
                raise Glitch, "STChain.getTreeLogLike() unknown metric '%s'" % self.stMcmc.dMetric
            aTree.logLike -= thisDist

    def p4stm_getTreeLogLike(self):
        self.stm.wipeBigTPointers()
        for n in self.propTree.nodes:
            if n.parent:
                self.bigTr.setParent(n.nodeNum, n.parent.nodeNum)
            if n.leftChild:
                self.bigTr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
            #else:
            #    bigTr.setNodeTaxNum(n.nodeNum, tNames.index(n.name))
            if n.sibling:
                self.bigTr.setSibling(n.nodeNum, n.sibling.nodeNum)
        self.stm.setBigTInternalBits()
        sd = self.stm.getSymmDiff()
        self.propTree.logLike = -sd
        
    def propose(self, theProposal):
        gm = ['STChain.propose()']
        #print "propose() About to propose %s" % theProposal.name


        if theProposal.name == 'nni':
            #self.proposeNni(theProposal)
            self.propTree.nni()
            if theProposal.doAbort:
                pass
            else:
                if not self.propTree.preAndPostOrderAreValid:
                    self.propTree.setPreAndPostOrder()
        elif theProposal.name == 'spr':
            #self.proposeNni(theProposal)
            self.propTree.randomSpr()
            if theProposal.doAbort:
                pass
            else:
                if not self.propTree.preAndPostOrderAreValid:
                    self.propTree.setPreAndPostOrder()

        # elif theProposal.name == 'polytomy':
        #     self.proposePolytomy(theProposal)
        #     if not self.propTree.preAndPostOrderAreValid:
        #         self.propTree.setPreAndPostOrder()
        #     self.propTree.setCStuff()
        #     pf.p4_setPrams(self.propTree.cTree, -1)
            
        else:
            gm.append('Unlisted proposal.name=%s  Fix me.' % theProposal.name)
            raise Glitch, gm

        if theProposal.doAbort:
            return 0.0
        else:
            #print "...about to calculate the likelihood of the propTree."
            if self.stMcmc.usingP4stmModule:
                self.p4stm_getTreeLogLike()
            else:
                
                self.getTreeLogLike(self.propTree)
            #print "propTree logLike is", self.propTree.logLike

            logLikeRatio = self.propTree.logLike - self.curTree.logLike

            theSum = logLikeRatio + self.logProposalRatio + self.logPriorRatio
            return theSum


    def gen(self, aProposal):
        gm = ['STChain.gen()']

        # doAborts means that it was not a valid generation,
        # neither accepted or rejected.  Give up, by returning True.

        acceptMove = False

        pRet = self.propose(aProposal)

        #print "pRet = %.6f" % pRet,
        if not aProposal.doAbort:
            if pRet < -100.0:  # math.exp(-100.) is 3.7200759760208361e-44
                r = 0.0
            elif pRet >= 0.0:
                r = 1.0
            else:
                r = math.exp(pRet)

            if r == 1.0:
                acceptMove = True
            elif random.random() < r:
                acceptMove = True

        #print "acceptMove = %s" % acceptMove

        aProposal.nProposals[self.tempNum] += 1
        if acceptMove:
            aProposal.accepted = True
            aProposal.nAcceptances[self.tempNum]  += 1
            

        #if not aProposal.doAbort:
        if acceptMove:
            a = self.propTree
            b = self.curTree
        else:
            a = self.curTree
            b = self.propTree


        if aProposal.name in ['nni', 'spr']:
            b.logLike = a.logLike
            a.copyToTree(b)

        else:
            gm.append('Unlisted proposal.name = %s  Fix me.' % aProposal.name)
            raise Glitch, gm




# for proposal probs
fudgeFactor = {}
fudgeFactor['local'] = 1.5
fudgeFactor['brLen'] = 1.0
fudgeFactor['polytomy'] = 1.0
fudgeFactor['root3'] = 0.02           
fudgeFactor['compLocation'] = 0.01    
fudgeFactor['rMatrixLocation'] = 0.01 
fudgeFactor['gdasrvLocation'] = 0.01  



class STMcmcTunings(object):
    def __init__(self):
        object.__setattr__(self, 'chainTemp', 0.15)  # was 0.2
        object.__setattr__(self, 'nni', None)  
        object.__setattr__(self, 'spr', None)  


    def __setattr__(self, item, val):
        #print "Got request to set %s to %s" % (item, val)
        if item in self.__dict__.keys():
            # Here is where I should do the sanity checking of the new vals.  Some day.
            #print "    Setting tuning '%s' to %s" % (item, val)
            object.__setattr__(self, item, val)
        else:
            print self.dump()
            gm = ["\nSTMcmcTunings.__setattr__()"]
            gm.append("Can't set tuning '%s'-- no such tuning." % item)
            raise Glitch, gm
    
    def reprString(self, advice=True):
        lst = ["\nSTMcmc.tunings:"]
        spacer = ' ' * 4
        lst.append("%s%15s: %s" % (spacer, 'chainTemp', self.chainTemp))
        lst.append("%s%15s: %s" % (spacer, 'nni', self.nni))
        lst.append("%s%15s: %s" % (spacer, 'spr', self.spr))
        return string.join(lst, '\n')

    def dump(self):
        print self.reprString()

    def __repr__(self):
        return  self.reprString()


class STMcmcProposalProbs(dict):
    """User-settable relative proposal probabilities.

    An instance of this class is made as STMcmc.prob, where you can
    do, for example,
        yourSTMcmc.prob.nni = 2.0

    These are relative proposal probs, that do not sum to 1.0, and
    affect the calculation of the final proposal probabilities (ie the
    kind that do sum to 1).  It is a relative setting, and the default
    is 1.0.  Setting it to 0 turns it off.  For small
    probabilities, setting it to 2.0 doubles it.  For bigger
    probabilities, setting it to 2.0 makes it somewhat bigger.

    Check the effect that it has by doing a
        yourSTMcmc.writeProposalIntendedProbs()
    which prints out the final calculated probabilities. 
    """
    
    def __init__(self):
        object.__setattr__(self, 'nni', 1.0)
        object.__setattr__(self, 'spr', 1.0)


    def __setattr__(self, item, val):
        # complaintHead = "\nSTMcmcProposalProbs.__setattr__()"
        gm = ["\nSTMcmcProposalProbs(). (set %s to %s)" % (item, val)]
        theKeys = self.__dict__.keys()
        if item in theKeys:
            try:
                val = float(val)
                if val < 1e-9:
                    val = 0
                object.__setattr__(self, item, val)
            except:
                gm.append("Should be a float.  Got '%s'" % val)
                raise Glitch, gm
                
        else:
            self.dump()
            gm.append("    Can't set '%s'-- no such proposal." % item)
            raise Glitch, gm

    def reprString(self):
        stuff = ["\nUser-settable relative proposal probabilities, from yourMcmc.prob"]
        stuff.append("  To change it, do eg ")
        stuff.append("    yourMcmc.prob.comp = 0.0 # turns comp proposals off")
        stuff.append("  Current settings:")
        theKeys = self.__dict__.keys()
        theKeys.sort()
        for k in theKeys:
            stuff.append("        %15s: %s" % (k, getattr(self, k)))
        return string.join(stuff, '\n')

    def dump(self):
        print self.reprString()

    def __repr__(self):
        return  self.reprString()


                           
class STProposal(object):
    def __init__(self, theSTMcmc=None):
        self.name = None
        self.variant = 'gtr'  # only for rMatrix.  2p or gtr
        self.stMcmc = theSTMcmc            # reference loop!
        self.nChains = theSTMcmc.nChains
        self.pNum = -1
        self.mtNum = -1
        self.weight = 1.0
        #self.tuning = None
        self.nProposals = [0] * self.nChains
        self.nAcceptances = [0] * self.nChains
        self.accepted = 0
        #self.topologyChanged = 0
        #self.nTopologyChangeAttempts = [0] * self.nChains
        #self.nTopologyChanges = [0] * self.nChains
        self.doAbort = False
        self.nAborts = [0] * self.nChains

    def dump(self):
        print "proposal name=%-10s pNum=%2i, mtNum=%2i, weight=%5.1f, tuning=%7.2f" % (
            '%s,' % self.name, self.pNum, self.mtNum, self.weight, self.tuning)
        print "    nProposals   by temperature:  %s" % self.nProposals
        print "    nAcceptances by temperature:  %s" % self.nAcceptances
        
    def _getTuning(self):
        if self.name in ['nni', 'spr']:
            #print "getting tuning for %s, returning %f" % (self.name, getattr(self.mcmc.tunings, self.name))
            #print self.stMcmc.tunings
            return getattr(self.stMcmc.tunings, self.name)
        elif self.name in ['comp', 'rMatrix', 'gdasrv', 'pInvar']:
            #print "getting tuning for %s, partNum %i, returning %f" % (
            #    self.name, self.pNum, getattr(self.mcmc.tunings.parts[self.pNum], self.name))
            # the variant attribute is new, and can mess up reading older pickles.
            if self.name == 'rMatrix' and hasattr(self, 'variant') and self.variant == '2p':
                return getattr(self.mcmc.tunings.parts[self.pNum], 'twoP')
            else:
                return getattr(self.mcmc.tunings.parts[self.pNum], self.name)
        elif self.name in ['compLocation', 'rMatrixLocation']:  # new, and causes older checkpoints to gag.
            #print "getting tuning for %s, partNum %i, returning %f" % (
            #    self.name, self.pNum, getattr(self.mcmc.tunings.parts[self.pNum], self.name))
            if hasattr(self.mcmc.tunings.parts[self.pNum], self.name):
                return getattr(self.mcmc.tunings.parts[self.pNum], self.name)
            else:
                return None
        else:
            return None
        
    def _setTuning(self, whatever):
        raise Glitch, "Can't set tuning this way."
    def _delTuning(self):
        raise Glitch, "Can't del tuning."
    
    tuning = property(_getTuning, _setTuning, _delTuning) 
            


class STMcmc(object):
    """An MCMC for making supertrees from a set of input trees.

    inTrees
        A list of p4 tree objects.

    nChains
        The number of chains in the MCMCMC, default 1

    runNum
        You may want to do more than one 'run' in the same directory,
        to facilitate convergence testing (another idea stolen from
        MrBayes, so thanks to the authors).  The first runNum would be
        0, and samples, likelihoods, and checkPoints are written to
        files with that number.

    sampleInterval
        Interval at which the (cold) chain is sampled, including
        writing a tree, and the logLike.

    checkPointInterval
        Intervals at which the MCMC is checkpointed, meaning that the
        whole thing is written to a pickle file.  You can re-start
        from a checkpoint, eg in the event of a crash, or if you just
        want to make the MCMC run longer.  You can turn off
        checkPointing by setting it to zero or None.

    defaultBeta
        By default, the defaultBeta is 1.0.  The beta is the weight as
        given in Steel and Rodrigo 2008.  Each input tree needs one
        (although it may be that each $X_i$ needs one ...), and they
        can be assigned before being passed given to this class.  If
        theInputTree.beta is already assigned, it is left alone, but
        if it has not been assigned, then a default value as given
        here is assigned to it.

    To prepare for a run, instantiate an Mcmc object::

        m = STMcmc(treeList, sampleInterval=10, checkpointInterval=2000)

    To start it running, do this::

        m.run(10000) # Tell it the number of generations to do

    You will want to make the last generation land on a
    checkPointInterval.  So in this case 10000 is ok, but 9000 is not.

    As it runs, it saves trees and likelihoods at sampleInterval
    intervals (actually whenever the current generation number is
    evenly divisible by the sampleInterval).

    Whenever the current generation number is evenly divisible by the
    checkPointInterval it will write a checkPoint file.  A checkPoint
    file is the whole MCMC, pickled.  Using a checkPoint, you can
    re-start an STMcmc from the point you left off.  Or, in the event
    of a crash, you can restart from the latest checkPoint.  But the
    most useful thing about them is that you can query checkPoints to
    get information about how the chain has been running, and about
    convergence diagnostics.

    In order to restart the MCMC from the end of a previous run:: 

        # read the last checkPoint file
        m = func.unPickleStMcmc(0)  # runNum 0
        m.run(20000)

    Its that easy if your previous run finished properly.  However, if
    your previous run has crashed and you want to restart it from a
    checkPoint, then you will need to repair the sample output files
    to remove samples that were taken after the last checkPoint, but
    before the crash.  Fix the trees, likelihoods, prams, and sims.
    (You probably do not need to beware of confusing gen (eg 9999) and
    gen+1 (eg 10000) issues.)  When you remove trees from the tree
    files be sure to leave the 'end;' at the end-- p4 needs it, and
    will deal with it.

    The checkPoints can help with convergence testing.  To help with
    that, you can use the STMcmcCheckPointReader class.  It will print
    out a table of average standard deviations of split supports
    between 2 runs, or between 2 checkPoints from the same run.  It
    will print out tables of proposal acceptances to show whether they
    change over the course of the MCMC.

    To make a consensus from the trees from more than one run, you can
    add trees to an existing TreePartitions object, like this::

        tp = TreePartitions('mcmc_trees_0.nex', skip=500)
        tp.read('mcmc_trees_1.nex', skip=500)
        t = tp.consensus()
        # and perhaps something like ...
        for n in t.iterInternalsNoRoot():
            n.name = '%.0f' % (100. * n.br.support)
        t.writeNexus('cons.nex')
    """

    
    #def __init__(self, inTrees, nChains=1, runNum=0, sampleInterval=100, checkPointInterval=10000, simulate=None, writePrams=True, constraints=None, verbose=True):
    def __init__(self, inTrees, nChains=1, runNum=0, sampleInterval=100, checkPointInterval=None, verbose=True, defaultBeta=1.0, allowPolytomyInTrees=False, dMetric='sd', useP4stm=False):
        gm = ['STMcmc.__init__()']

        self.verbose = verbose

        # Each inTree needs a beta, default 1.0
        for t in inTrees:
            if not hasattr(t, 'beta'):
                t.beta = defaultBeta
            else:
                try:
                    t.beta = float(t.beta)
                except ValueError:
                    if t.name:
                        gm.append("input tree %s" % t.name)
                    gm.append("I was unable to convert beta %s to a float" % t.beta)

        if allowPolytomyInTrees:
            pass
        else:
            for t in inTrees:
                if t.isFullyBifurcating():
                    pass
                else:
                    gm.append("At the moment STMcmc wants trees that are fully bifurcating.")
                    raise Glitch, gm
        
        try:
            nChains = int(nChains)
        except (ValueError,TypeError):
            gm.append("nChains should be an int, 1 or more.  Got %s" % nChains)
            raise Glitch, gm
        if nChains < 1:
            gm.append("nChains should be an int, 1 or more.  Got %s" % nChains)
            raise Glitch, gm
        self.nChains = nChains
        self.chains = []
        self.gen = -1
        self.startMinusOne = -1
        self.constraints = None
        self.simulate = None

        goodDMetrics = ['sd', 'scqdist']
        if dMetric not in goodDMetrics:
            gm.append("dMetric should be one of %s" % goodDMetrics)
            gm.append("got '%s'" % dMetric)
            raise Glitch, gm
        self.dMetric = dMetric

        try:
            runNum = int(runNum)
        except (ValueError, TypeError):
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise Glitch, gm
        if runNum < 0:
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise Glitch, gm
        self.runNum = runNum

        # Check that we are not going to over-write good stuff
        ff = os.listdir(os.getcwd())
        hasPickle = False
        for fName in ff:
            if fName.startswith("mcmc_checkPoint_%i." % self.runNum):
                hasPickle = True
                break
        if hasPickle:
            gm.append("runNum is set to %i" % self.runNum)
            gm.append("There is at least one mcmc_checkPoint_%i.xxx file in this directory." % self.runNum)
            gm.append("This is a new STMcmc, and I am refusing to over-write exisiting files.")
            gm.append("Maybe you want to re-start from the latest mcmc_checkPoint_%i file?" % self.runNum)
            gm.append("Otherwise, get rid of the existing mcmc_xxx_%i.xxx files and start again." % self.runNum)
            raise Glitch, gm

        if var.strictRunNumberChecking:
            # We want to start runs with number 0, so if runNum is more than that, check that there are other runs.
            if self.runNum > 0:
                for runNum2 in range(self.runNum):
                    hasTrees = False
                    for fName in ff:
                        if fName.startswith("mcmc_trees_%i" % runNum2):
                            hasTrees = True
                            break
                    if not hasTrees:
                        gm.append("runNum is set to %i" % self.runNum)
                        gm.append("runNums should go from zero up.")
                        gm.append("There are no mcmc_trees_%i.nex files to show that run %i has been done." % (runNum2, runNum2))
                        gm.append("Set the runNum to that, first.")
                        raise Glitch, gm
                              
        
        
        self.sampleInterval = sampleInterval
        self.checkPointInterval = checkPointInterval

        self.proposals = []
        self.proposalsHash = {}
        self.propWeights = []
        self.cumPropWeights = []
        self.totalPropWeights = 0.0

        self.treePartitions = None
        self.likesFileName = "mcmc_likes_%i" % runNum
        self.treeFileName = "mcmc_trees_%i.nex" % runNum
        #self.simFileName = "mcmc_sims_%i" % runNum
        #self.pramsFileName = "mcmc_prams_%i" % runNum
        #self.writePrams = writePrams

        self.lastTimeCheck = None

        # if simulate:
        #     try:
        #         simulate = int(simulate)
        #     except (ValueError, TypeError):
        #         gm.append("Arg 'simulate' should be an int, 1-31, inclusive.")
        #         raise Glitch, gm
        #     if simulate <= 0 or simulate > 31:
        #         gm.append("Arg 'simulate' should be an int, 1-31, inclusive.")
        #         raise Glitch, gm
        # self.simulate = simulate
        # if self.simulate:
        #     self.simTree = self.tree.dupe()
        #     self.simTree.data = self.tree.data.dupe()
        #     self.simTree.calcLogLike(verbose=False)
        #else:
        #    self.simTree = None

        if self.nChains > 1:
            self.swapMatrix = []
            for i in range(self.nChains):
                self.swapMatrix.append([0] * self.nChains)
        else:
            self.swapMatrix = None

        self.tunings = STMcmcTunings() 
        self.prob = STMcmcProposalProbs()                
                
        # Zap internal node names
        # for n in aTree.root.iterInternals():
        #     if n.name:
        #         n.name = None

        allNames = []
        for t in inTrees:
            t.unsorted_taxNames = [n.name for n in t.iterLeavesNoRoot()]
            allNames += t.unsorted_taxNames
        self.taxNames = list(set(allNames))
        self.nTax = len(self.taxNames)
        for t in inTrees:
            sorted_taxNames = []
            t.taxBits = 0L
            for tNum in range(self.nTax):
                tN = self.taxNames[tNum]
                if tN in t.unsorted_taxNames:
                    sorted_taxNames.append(tN)
                    adder = 1L << tNum
                    t.taxBits += adder
            t.taxNames = sorted_taxNames
            t.allOnes = 2L**(t.nTax) - 1
            t.makeSplitKeys()
            t.skSet = set([n.br.splitKey for n in t.iterInternalsNoRoot()])

        self.trees = inTrees
        self.tree = func.randomTree(taxNames=self.taxNames, name='stTree', randomBrLens=False)
        self.tree.makeSplitKeys()

        self.usingP4stmModule = False
        if useP4stm:
            if haveLoadedP4stm:
               self.usingP4stmModule = True
            else:
                gm.append("useP4stm is set, but I could not load the p4stm or the pyusblas modules.")
                raise Glitch, gm
                



    def _makeProposals(self):
        """Make proposals for the stmcmc."""

        gm = ['STMcmc._makeProposals()']


        # nni
        if self.prob.nni:
            p = STProposal(self)
            p.name = 'nni'
            p.weight = self.prob.nni # * (len(self.tree.nodes) - 1) * fudgeFactor['nni']
            #p.tuning = self.tunings.local
            self.proposals.append(p)
            #object.__setattr__(self.tuningsUsage, 'local', p)

        if self.prob.spr:
            p = STProposal(self)
            p.name = 'spr'
            p.weight = self.prob.spr # * (len(self.tree.nodes) - 1) * fudgeFactor['nni']
            #p.tuning = self.tunings.local
            self.proposals.append(p)
            #object.__setattr__(self.tuningsUsage, 'local', p)


        if not self.proposals:
            gm.append("No proposals?")
            raise Glitch, gm
        self.propWeights = []
        for p in self.proposals:
            self.propWeights.append(p.weight)
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise Glitch, gm
        for p in self.proposals:
            self.proposalsHash[p.name] = p

    def _refreshProposalProbsAndTunings(self):
        """Adjust proposals after a restart."""

        gm = ['STMcmc._refreshProposalProbsAndTunings()']

        for p in self.proposals:
            # nni
            if p.name == 'nni':
                #p.weight = self.prob.local * (len(self.tree.nodes) - 1) * fudgeFactor['local']
                p.weight = self.prob.nni


        self.propWeights = []
        for p in self.proposals:
            self.propWeights.append(p.weight)
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise Glitch, gm




    def writeProposalAcceptances(self):
        """Pretty-print the proposal acceptances."""

        if (self.gen - self.startMinusOne) <= 0:
            print "\nSTMcmc.writeProposalAcceptances()  There is no info in memory. "
            print " Maybe it was just emptied after writing to a checkpoint?  "
            print "If so, read the checkPoint and get the proposalAcceptances from there."
        else:

            spacer = ' ' * 8
            print "\nProposal acceptances, run %i, for %i gens, from gens %i to %i, inclusive." % (
                self.runNum, (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen)
            print "%s %15s %10s %13s%8s" % (spacer, 'proposal', 'nProposals', 'acceptance(%)', 'tuning')
            for p in self.proposals:
                print "%s" % spacer,
                print "%15s" % p.name,
                print "%10i" % p.nProposals[0],

                if p.nProposals[0]: # Don't divide by zero
                    print "       %5.1f " % (100.0 * float(p.nAcceptances[0]) / float(p.nProposals[0])),
                else:
                    print "           - ",

                if p.tuning == None:
                    print "      -",
                elif p.tuning < 2.0:
                    print "  %5.3f" % p.tuning,
                else:
                    print "%7.1f" % p.tuning,
                print

            # # Tabulate topology changes, if any were attempted.
            # doTopol = 0
            # p = None
            # try:
            #     p = self.proposalsHash['local']
            # except KeyError:
            #     pass
            # if p:
            #     for tNum in range(self.nChains):
            #         if p.nTopologyChangeAttempts[tNum]:
            #             doTopol = 1
            #             break
            #     if doTopol:
            #         p = self.proposalsHash['local']
            #         print "'Local' proposal-- attempted topology changes"
            #         print "%s tempNum   nProps nAccepts percent nTopolChangeAttempts nTopolChanges percent" % spacer
            #         for tNum in range(self.nChains):
            #             print "%s" % spacer,
            #             print "%4i " % tNum,
            #             print "%9i" % p.nProposals[tNum],
            #             print "%8i" % p.nAcceptances[tNum],
            #             print "  %5.1f" % (100.0 * float(p.nAcceptances[tNum]) / float(p.nProposals[tNum])),
            #             print "%20i" % p.nTopologyChangeAttempts[tNum],
            #             print "%13i" % p.nTopologyChanges[tNum],
            #             print "  %5.1f" % (100.0 * float(p.nTopologyChanges[tNum])/float(p.nTopologyChangeAttempts[tNum]))
            #     else:
            #         print "%sFor the 'local' proposals, there were no attempted" % spacer
            #         print "%stopology changes in any of the chains." % spacer


            # Check for aborts.
            # p = None
            # try:
            #     p = self.proposalsHash['local']
            # except KeyError:
            #     pass
            # if p:
            #     if hasattr(p, 'nAborts'):
            #         if p.nAborts[0]:
            #             print "The 'local' proposal had %i aborts." % p.nAborts[0]
            #             print "(Aborts might be due to brLen proposals too big or too small)"
            #             if self.constraints:
            #                 print "(Or, more likely, due to violated constraints.)"
            #         else:
            #             print "The 'local' proposal had no aborts (either due to brLen proposals"
            #             print "too big or too small, or due to violated constraints)."
            # for pN in ['polytomy', 'compLocation', 'rMatrixLocation', 'gdasrvLocation']:
            #     p = None
            #     try:
            #         p = self.proposalsHash[pN]
            #     except KeyError:
            #         pass
            #     if p:
            #         if hasattr(p, 'nAborts'):
            #             print "The %15s proposal had %5i aborts." % (p.name, p.nAborts[0])
            
            
            

    def writeSwapMatrix(self):
        print "\nChain swapping, for %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen)
        print "    Swaps are presented as a square matrix, nChains * nChains."
        print "    Upper triangle is the number of swaps proposed between two chains."
        print "    Lower triangle is the percent swaps accepted."
        print "    The current tunings.chainTemp is %5.3f\n" % self.tunings.chainTemp
        print " " * 10,
        for i in range(self.nChains):
            print "%7i" % i,
        print
        print " " * 10,
        for i in range(self.nChains):
            print "   ----",
        print
        for i in range(self.nChains):
            print " " * 7, "%2i" % i,
            for j in range(self.nChains):
                if i < j: # upper triangle
                    print "%7i" % self.swapMatrix[i][j],
                elif i == j:
                    print "      -",
                else:
                    if self.swapMatrix[j][i] == 0: # no proposals
                        print "      -",
                    else:
                        print "  %5.1f" % (100.0 * float(self.swapMatrix[i][j]) / float(self.swapMatrix[j][i])),
            print
        

    def _makeChainsAndProposals(self):
        """Make chains and proposals."""

        gm = ['STMcmc._makeChainsAndProposals()']

        #random.seed(0)

        # Make chains, if needed
        if not self.chains:
            self.chains = []
            for chNum in range(self.nChains):
                aChain = STChain(self)
                aChain.tempNum = chNum  # Temperature.  Set this way to start, but it changes.
                self.chains.append(aChain)
        if not self.proposals:
            self._makeProposals()

            # If we are going to be doing the resolution class prior
            # in the polytomy move, we want to pre-compute the logs of
            # T_{n,m}.  Its a vector with indices (ie m) from zero to
            # nTax-2 inclusive.
            # if self.proposalsHash.has_key('polytomy') and self.tunings.doPolytomyResolutionClassPrior:
            #     p = self.proposalsHash['polytomy']
            #     bigT = func.nUnrootedTreesWithMultifurcations(self.tree.nTax)
            #     p.logBigT = [0.0] * (self.tree.nTax - 1)
            #     for i in range(1, self.tree.nTax - 1):
            #         p.logBigT[i] = math.log(bigT[i])
            #     #print p.logBigT

    def _setOutputTreeFile(self):
        """Setup the (output) tree file for the stmcmc."""

        gm = ['STMcmc._setOutputTreeFile()']

        # Write the preamble for the trees outfile.
        self.treeFile = file(self.treeFileName, 'w')
        self.treeFile.write('#nexus\n\n')
        self.treeFile.write('begin taxa;\n')
        self.treeFile.write('  dimensions ntax=%s;\n' % self.tree.nTax)
        self.treeFile.write('  taxlabels')
        for tN in self.tree.taxNames:
            self.treeFile.write(' %s' % func.nexusFixNameIfQuotesAreNeeded(tN))
        self.treeFile.write(';\nend;\n\n')

        self.treeFile.write('begin trees;\n')
        self.translationHash = {}
        i = 1
        for tName in self.tree.taxNames:
            self.translationHash[tName] = i
            i += 1

        self.treeFile.write('  translate\n')
        for i in range(self.tree.nTax - 1):
            self.treeFile.write('    %3i %s,\n' % (
                i + 1, func.nexusFixNameIfQuotesAreNeeded(self.tree.taxNames[i])))
        self.treeFile.write('    %3i %s\n' % (
            self.tree.nTax, func.nexusFixNameIfQuotesAreNeeded(self.tree.taxNames[-1])))
        self.treeFile.write('  ;\n')
        self.treeFile.write('  [Tree numbers are gen+1]\n')
        self.treeFile.close()


    def run(self, nGensToDo, verbose=True):
        """Start the STMcmc running."""

        gm = ['STMcmc.run()']

        #Keep track of the first gen of this call to run(), maybe restart
        firstGen = self.gen + 1

        if self.checkPointInterval:
            # We want a couple of things:
            #  1.  The last gen should be on checkPointInterval.  For
            #      example, if the checkPointInterval is 200, then doing
            #      100 or 300 generations will not be allowed cuz the
            #      chain would continue past the checkPoint-- bad.  Or if
            #      you re-start after 500 gens and change to a
            #      checkPointInterval of 200, then you won't be allowed to
            #      do 500 gens.
            #if ((self.gen + 1) + nGensToDo) % self.checkPointInterval == 0:
            if nGensToDo % self.checkPointInterval == 0:
                pass
            else:
                gm.append("With the current settings, the last generation won't be on a checkPointInterval.")
                gm.append("self.gen+1=%i, nGensToDo=%i, checkPointInterval=%i" % ((self.gen + 1),
                                                                                  nGensToDo, self.checkPointInterval)) 
                raise Glitch, gm
            #  2.  We also want the checkPointInterval to be evenly
            #      divisible by the sampleInterval.
            if self.checkPointInterval % self.sampleInterval == 0:
                pass
            else:
                gm.append("The checkPointInterval (%i) should be evenly divisible" % self.checkPointInterval)
                gm.append("by the sampleInterval (%i)." % self.sampleInterval)
                raise Glitch, gm


        if self.proposals:
            # Its either a re-start, or it has been thru autoTune().
            # I can tell the difference by self.gen, which is -1 after
            # autoTune()
            if self.gen == -1:
                self._makeChainsAndProposals()
                self._setOutputTreeFile()
                #if self.simulate:
                #    self.writeSimFileHeader(self.tree)
            # The probs and tunings may have been changed by the user.
            self._refreshProposalProbsAndTunings()

            # This stuff below should be the same as is done after pickling, see below.
            self.startMinusOne = self.gen

            # Start the tree partitions over.
            self.treePartitions = None
            # Zero the proposal counts
            for p in self.proposals:
                p.nProposals = [0] * self.nChains
                p.nAcceptances = [0] * self.nChains
                p.nTopologyChangeAttempts = [0] * self.nChains
                p.nTopologyChanges = [0] * self.nChains
            # Zero the swap matrix
            if self.nChains > 1:
                self.swapMatrix = []
                for i in range(self.nChains):
                    self.swapMatrix.append([0] * self.nChains)
            
        else:
            self._makeChainsAndProposals()
            self._setOutputTreeFile()
            #if self.simulate:
            #    self.writeSimFileHeader(self.tree)
        if verbose:
            self.writeProposalIntendedProbs()
            sys.stdout.flush()

        coldChainNum = 0
        
        # If polytomy is turned on, then it is possible to get a star
        # tree, in which case local will not work.  So if we have both
        # polytomy and local proposals, we should also have brLen.
        # if self.proposalsHash.has_key("polytomy") and self.proposalsHash.has_key("local"):
        #     if not self.proposalsHash.has_key('brLen'):
        #         gm.append("If you have polytomy and local proposals, you should have a brLen proposal as well.")
        #         gm.append("It can have a low proposal probability, but it needs to be there.")
        #         gm.append("Turn it on by eg yourMcmc.prob.brLen = 0.001")
        #         raise Glitch, gm

        
        if self.gen > -1:
            # it is a re-start, so we need to back over the "end;" in the tree files.
            f2 = file(self.treeFileName, 'a+')
            pos = -1
            while 1:
                f2.seek(pos, 2)
                c = f2.read(1)
                if c == ';':
                    break
                pos -= 1
            #print "pos now %i" % pos
            pos -= 3 # end;
            f2.seek(pos, 2)
            c = f2.read(4)
            #print "got c = '%s'" % c
            if c != "end;":
                gm.append("Mcmc.run().  Failed to find and remove the 'end;' at the end of the tree file.")
                raise Glitch, gm
            else:
                f2.seek(pos, 2)
                f2.truncate()
            f2.close()
                
            if verbose:
                print
                print "Re-starting the MCMC run %i from gen=%i" % (self.runNum, self.gen)
                print "Set to do %i more generations." % nGensToDo
                #if self.writePrams:
                #    if self.chains[0].curTree.model.nFreePrams == 0:
                #        print "There are no free prams in the model, so I am turning writePrams off."
                #        self.writePrams = False
                sys.stdout.flush()

            self.startMinusOne = self.gen
        else:
            if verbose:
                print "Starting the MCMC %s run %i" % ((self.constraints and "(with constraints)" or ""), self.runNum)
                print "Set to do %i generations." % nGensToDo
                # if self.writePrams:
                #     if self.chains[0].curTree.model.nFreePrams == 0:
                #         print "There are no free prams in the model, so I am turning writePrams off."
                #         self.writePrams = False
                #     else:
                #         pramsFile = file(self.pramsFileName, 'a')
                #         self.chains[0].curTree.model.writePramsProfile(pramsFile)
                #         #pramsFile.write("     genPlus1")
                #         pramsFile.write("genPlus1")
                #         self.chains[0].curTree.model.writePramsHeaderLine(pramsFile)
                #         pramsFile.close()
                sys.stdout.flush()
            
        if verbose:
            print "Sampling every %i." % self.sampleInterval
            if self.checkPointInterval:
                print "CheckPoints written every %i." % self.checkPointInterval
            if nGensToDo <= 20000:
                print "One dot is 100 generations."
            else:
                print "One dot is 1000 generations."
            sys.stdout.flush()

        self.treePartitions = None
        realTimeStart = time.time()
        self.lastTimeCheck = time.time()

        abortableProposals = ['nni', 'spr', 'polytomy']

        for gNum in range(nGensToDo):
            self.gen += 1
            #Do an initial time estimate based on 100 gens 
            if nGensToDo > 100 and self.gen-firstGen == 100:
                diff_secs = time.time() - realTimeStart
                total_secs = (float(nGensToDo)/float(100))*float(diff_secs)
                deltaTime = datetime.timedelta(seconds = int(round(total_secs)))
                print "Estimated completion time: %s days, %s" % (
                    deltaTime.days, time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds)))

            # Above is a list of proposals where it is possible to abort.
            # When a gen(aProposal) is made, below, aProposal.doAbort
            # might be set, in which case we want to skip it for this
            # gen.  But we want to start each 'for chNum' loop with
            # doAborts all turned off.
            
            for chNum in range(self.nChains):
                failure = True
                nAttempts = 0
                while failure:
                    # Get the next proposal
                    gotIt = False
                    safety = 0
                    while not gotIt:
                        theRan = random.uniform(0.0, self.totalPropWeights)
                        for i in range(len(self.cumPropWeights)):
                            if theRan < self.cumPropWeights[i]:
                                break
                        aProposal = self.proposals[i]
                        gotIt = True
                        if aProposal.name == 'local':
                            if self.chains[chNum].curTree.nInternalNodes == 1:  # Can't do local on a star tree.
                                aProposal = self.proposalsHash['brLen']
                        elif aProposal.name == 'root3':
                            if self.chains[chNum].curTree.nInternalNodes == 1:  # Can't do root3 on a star tree.
                                gotIt = False
                        if aProposal.doAbort:
                            gotIt = False
                        safety += 1
                        if safety > 1000:
                            gm.append("Could not find a proposal after %i attempts." % safety)
                            gm.append("Possibly a programming error.")
                            gm.append("Or possibly it is just a pathologically frustrating Mcmc.")
                            raise Glitch, gm

                    #if gNum % 2:
                    #    aProposal = self.proposalsHash['brLen']
                    #else:
                    #    aProposal = self.proposalsHash['comp']

                    if 0:
                        print "==== gNum=%i, chNum=%i, aProposal=%s (part %i)" % (
                            gNum, chNum, aProposal.name, aProposal.pNum),
                        sys.stdout.flush()
                        #print gNum,

                    failure = self.chains[chNum].gen(aProposal)  # success returns None

                    if 0:
                        if failure:
                            print "    failure"
                        else:
                            print
                        
                    nAttempts += 1
                    if nAttempts > 1000:
                        gm.append("Was not able to do a successful generation after %i attempts." % nAttempts)
                        raise Glitch, gm
                #print "   Mcmc.run(). finished a gen on chain %i" % (chNum)
                for pr in abortableProposals:
                    if self.proposalsHash.has_key(pr):
                        self.proposalsHash[pr].doAbort = False
                        
            # Do swap, if there is more than 1 chain.
            if self.nChains == 1:
                coldChain = 0
            else:
                # Chain swapping stuff was lifted from MrBayes.  Thanks again.
                chain1,chain2 = random.sample(self.chains, 2)

                # Use the upper triangle of swapMatrix for nProposed's
                if chain1.tempNum < chain2.tempNum:
                    self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1
                else:
                    self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1

                lnR =  (1.0 / (1.0 + (self.tunings.chainTemp * chain1.tempNum))) * chain2.curTree.logLike
                lnR += (1.0 / (1.0 + (self.tunings.chainTemp * chain2.tempNum))) * chain1.curTree.logLike
                lnR -= (1.0 / (1.0 + (self.tunings.chainTemp * chain1.tempNum))) * chain1.curTree.logLike
                lnR -= (1.0 / (1.0 + (self.tunings.chainTemp * chain2.tempNum))) * chain2.curTree.logLike


                if lnR < -100.0:
                    r = 0.0
                elif lnR >= 0.0:
                    r = 1.0
                else:
                    r = math.exp(lnR)

                acceptSwap = 0
                if random.random() < r:
                    acceptSwap = 1

                if acceptSwap:
                    # Use the lower triangle of swapMatrix to keep track of nAccepted's
                    if chain1.tempNum < chain2.tempNum:
                        self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1
                    else:
                        self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1

                    # Do the swap
                    chain1.tempNum, chain2.tempNum = chain2.tempNum, chain1.tempNum

                # Find the cold chain, the one where tempNum is 0
                coldChainNum = -1
                for i in range(len(self.chains)):
                    if self.chains[i].tempNum == 0:
                        coldChainNum = i
                        break
                if coldChainNum == -1:
                    gm.append("Unable to find which chain is the cold chain.  Bad.")
                    raise Glitch, gm

            # If it is a writeInterval, write stuff
            if (self.gen + 1) % self.sampleInterval == 0:
                if 1:
                    likesFile = file(self.likesFileName, 'a')
                    likesFile.write('%11i %f\n' % (self.gen + 1, self.chains[coldChainNum].curTree.logLike))
                    likesFile.close()
                    treeFile = file(self.treeFileName, 'a')
                    treeFile.write("  tree t_%i = [&U] " % (self.gen + 1))
                    self.chains[coldChainNum].curTree.writeNewick(treeFile,
                                                               withTranslation=1,
                                                               translationHash=self.translationHash,
                                                               doMcmcCommandComments=False)
                    treeFile.close()

                # if self.writePrams:
                #     pramsFile = file(self.pramsFileName, 'a')
                #     #pramsFile.write("%12i " % (self.gen + 1))
                #     pramsFile.write("%12i" % (self.gen + 1))
                #     self.chains[coldChainNum].curTree.model.writePramsLine(pramsFile)
                #     pramsFile.close()

                # Do a simulation
                if self.simulate:
                    #print "about to simulate..."
                    self.doSimulate(self.chains[coldChainNum].curTree)
                    #print "...finished simulate."
                
                # Do other stuff.
                if hasattr(self, 'hook'):
                    self.hook(self.chains[coldChainNum].curTree)

                if 0 and self.constraints:
                    print "Mcmc x1c"
                    print self.chains[0].verifyIdentityOfTwoTreesInChain()
                    print "b checking curTree .."
                    self.chains[0].curTree.checkSplitKeys()
                    print "b checking propTree ..."
                    self.chains[0].propTree.checkSplitKeys()
                    print "Mcmc xxx"

                # Add curTree to treePartitions
                if self.treePartitions:
                    self.treePartitions._getSplitsFromTree(self.chains[coldChainNum].curTree)
                else:
                    self.treePartitions = TreePartitions(self.chains[coldChainNum].curTree)
                # After _getSplitsFromTree, need to follow, at some point,
                # with _finishSplits().  Do that when it is pickled, or at the end of the run.

                # Checking and debugging constraints
                if 0 and self.constraints:
                    print "Mcmc x1d"
                    print self.chains[coldChainNum].verifyIdentityOfTwoTreesInChain()
                    print "c checking curTree ..."
                    self.chains[coldChainNum].curTree.checkSplitKeys()
                    print "c checking propTree ..."
                    self.chains[coldChainNum].propTree.checkSplitKeys()
                    #print "c checking that all constraints are present"
                    #theSplits = [n.br.splitKey for n in self.chains[0].curTree.iterNodesNoRoot()]
                    #for sk in self.constraints.constraints:
                    #    if sk not in theSplits:
                    #        gm.append("split %i is not present in the curTree." % sk)
                    #        raise Glitch, gm
                    print "Mcmc zzz"

                # Check that the curTree has all the constraints
                if self.constraints:
                    splitsInCurTree = [n.br.splitKey for n in self.chains[coldChainNum].curTree.iterInternalsNoRoot()]
                    for sk in self.constraints.constraints:
                        if sk not in splitsInCurTree:
                            gm.append("Programming error.")
                            gm.append("The current tree (the last tree sampled) does not contain constraint")
                            gm.append("%s" % func.getSplitStringFromKey(sk, self.tree.nTax))
                            raise Glitch, gm

                # If it is a checkPointInterval, pickle
                if self.checkPointInterval and (self.gen + 1) % self.checkPointInterval == 0:
                    self.checkPoint()

                    # The stuff below needs to be done in a re-start as well.  See above "if self.proposals:"
                    self.startMinusOne = self.gen

                    # Start the tree partitions over.
                    self.treePartitions = None
                    # Zero the proposal counts
                    for p in self.proposals:
                        p.nProposals = [0] * self.nChains
                        p.nAcceptances = [0] * self.nChains
                        p.nTopologyChangeAttempts = [0] * self.nChains
                        p.nTopologyChanges = [0] * self.nChains
                        p.nAborts = [0] * self.nChains
                    # Zero the swap matrix
                    if self.nChains > 1:
                        self.swapMatrix = []
                        for i in range(self.nChains):
                            self.swapMatrix.append([0] * self.nChains)
                

            # Reassuring pips ...
            if firstGen != self.gen: #We want to skip the first gen of every call to run()
                if nGensToDo <= 20000:
                    if (self.gen-firstGen) % 1000 == 0:
                        if verbose:
                            deltaTime = self._doTimeCheck(nGensToDo, firstGen, 1000)
                            if deltaTime.days:
                                timeString = "%s days, %s" % (
                                    deltaTime.days, time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds)))
                            else:
                                timeString = time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds))
                            print "%10i - %s" % (self.gen, timeString)

                        else:
                            sys.stdout.write(".")
                            sys.stdout.flush()
                    elif (self.gen-firstGen) % 100 == 0:
                        sys.stdout.write(".")
                        sys.stdout.flush()
                else:
                    if (self.gen-firstGen) % 50000 == 0:
                        if verbose:
                            deltaTime = self._doTimeCheck(nGensToDo, firstGen, 50000)
                            if deltaTime.days:
                                timeString = "%s days, %s" % (
                                    deltaTime.days, time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds)))
                            else:
                                timeString = time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds))
                            print "%10i - %s" % (self.gen, timeString)
                        else:
                            sys.stdout.write(".")
                            sys.stdout.flush()
                    elif (self.gen-firstGen) % 1000 == 0:
                        sys.stdout.write(".")
                        sys.stdout.flush()

        # Gens finished.  Clean up.
        print
        if verbose:
            print "Finished %s generations." % nGensToDo

        treeFile = file(self.treeFileName, 'a')
        treeFile.write('end;\n\n')
        treeFile.close()
        
    def _doTimeCheck(self, nGensToDo, firstGen, genInterval):
        """Time check 
        
        firstGen is the first generation of this call to Mcmc.run() else
        timing fails on restart"""
        nowTime = time.time()
        diff_secs = nowTime - self.lastTimeCheck
        total_secs = (float(nGensToDo-(self.gen-firstGen))/float(genInterval))*float(diff_secs)
        deltaTime = datetime.timedelta(seconds = int(round(total_secs)))
        self.lastTimeCheck = nowTime
        return deltaTime

    # def writeSimFileHeader(self, curTree):
    #     simFile = file(self.simFileName, 'a')
    #     simFile.write(" genPlus1")
    #     if 1 & self.simulate:  # If self.simulate contains a 1, do unconstrained log like
    #         for pNum in range(self.simTree.data.nParts):
    #             simFile.write(' uncLike%i' % pNum)
    #     if 2 & self.simulate:  # If self.simulate contains a 2, do bigX^2
    #         for pNum in range(self.simTree.model.nParts):
    #             simFile.write(' bigXSq%i' % pNum)
    #     if 4 & self.simulate:  # If self.simulate contains a 4, do meanNCharPerSite
    #         for pNum in range(self.simTree.model.nParts):
    #             simFile.write(' meanCharsPerSite%i' % pNum)
    #     if 8 & self.simulate:  # If self.simulate contains an 8, do c_m, the compStatFromCharFreqs
    #         for pNum in range(self.simTree.model.nParts):
    #             simFile.write(' c_mSim%i  c_mOrig%i' % (pNum, pNum))
    #     if 16 & self.simulate:  # If self.simulate contains a 16, do constant sites count
    #         for pNum in range(self.simTree.model.nParts):
    #             simFile.write(' nConstSites%i' % pNum)
    #     simFile.write('\n')
    #     simFile.close()
        

    # def doSimulate(self, curTree):
    #     curTree.copyToTree(self.simTree)
    #     curTree.model.copyValsTo(self.simTree.model)
    #     self.simTree.simulate()
    #     simFile = file(self.simFileName, 'a')
    #     simFile.write(" %11i" % (self.gen + 1))
    #     if 1 & self.simulate:  # If self.simulate contains a 1, do unconstrained log like
    #         for p in self.simTree.data.parts:
    #             simFile.write(' %f' % pf.getUnconstrainedLogLike(p.cPart))
    #     if 2 & self.simulate:  # If self.simulate contains a 2, do bigX^2
    #         #ret2 = self.simTree.data.compoChiSquaredTest(verbose=0, skipColumnZeros=True)
    #         #for pNum in range(self.simTree.model.nParts):
    #         #    simFile.write(' %f' % ret[pNum][0])
    #         ret = self.simTree.data.simpleBigXSquared()
    #         for pNum in range(self.simTree.model.nParts):
    #             simFile.write(' %f' % ret[pNum])
    #         #for i in range(len(ret)):
    #         #    if math.fabs(ret[i][0] - ret2[i]) > 0.000001:
    #         #        print "The two methods of bigXSquared calculations differed.  %f and %f" % (ret[i], ret2[i])
    #     if 4 & self.simulate:  # If self.simulate contains a 4, do meanNCharPerSite
    #         ret = self.simTree.data.meanNCharsPerSite()
    #         # ret is a list, one number per part
    #         for pNum in range(self.simTree.model.nParts):
    #             simFile.write(' %f' % ret[pNum])
    #     if 8 & self.simulate:  # If self.simulate contains an 8, do c_m, the compStatFromCharFreqs
    #         ret = self.simTree.compStatFromCharFreqs()
    #         ret2 = curTree.compStatFromCharFreqs()
    #         # ret is a list, one number per part
    #         for pNum in range(self.simTree.model.nParts):
    #             simFile.write(' %f  %f' % (ret[pNum], ret2[pNum]))
    #             #print ' compStatFromCharFreqs: %f  %f' % (ret[pNum], ret2[pNum])
    #     if 16 & self.simulate:  # If self.simulate contains a 16, do constant sites count
    #         ret = self.simTree.data.simpleConstantSitesCount()
    #         # ret is a list, one number per part
    #         for pNum in range(self.simTree.model.nParts):
    #             simFile.write(' %i' % ret[pNum])
    #     simFile.write('\n')
    #     simFile.close()

    def checkPoint(self):

        # Maybe we should not save the inTrees? -- would make it more lightweight.
        if 0:
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                print "chain %i ==================" % chNum
                ch.curTree.summarizeModelThingsNNodes()
        # the Stm object does not pickle, and this commented-out bit does not fix it.  
        savedStms = []
        savedBigTrs = []
        if self.usingP4stmModule:
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                savedStms.append(ch.stm)
                ch.stm = None
                savedBigTrs.append(ch.bigTr)
                ch.bigTr = None
            
        theCopy = copy.deepcopy(self)

        theCopy.treePartitions._finishSplits()
        theCopy.likesFile = None
        theCopy.treeFile = None
        #theCopy.treePartitions = None    # this can be the biggest part of the pickle.

        # Pickle it.
        fName = "mcmc_checkPoint_%i.%i" % (self.runNum, self.gen + 1)
        f = file(fName, 'w')
        cPickle.dump(theCopy, f, 1)
        f.close()

        if self.usingP4stmModule:
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                ch.stm = savedStms[chNum]
                ch.bigTr = savedBigTrs[chNum]
       

    

    def writeProposalProbs(self):
        """(Another) Pretty-print the proposal probabilities.

        See also STMcmc.writeProposalAcceptances().
        """

        nProposals = len(self.proposals)
        if not nProposals:
            print "STMcmc.writeProposalProbs().  No proposals (yet?)."
            return
        #intended = self.propWeights[:]
        #for i in range(len(intended)):
        #    intended[i] /= self.totalPropWeights
        #if math.fabs(sum(intended) - 1.0 > 1e-15):
        #    raise Glitch, 'bad sum of intended proposal probs. %s' % sum(intended)

        nAttained = [0] * nProposals
        nAccepted = [0] * nProposals
        for i in range(nProposals):
            nAttained[i] = self.proposals[i].nProposals[0]
            nAccepted[i] = self.proposals[i].nAcceptances[0]
        sumAttained = float(sum(nAttained)) # should be zero or nGen
        if not sumAttained:
            print "STMcmc.writeProposalProbs().  No proposals have been made."
            print "Possibly, due to it being a checkPoint interval, nProposals have all been set to zero."
            return
        #assert int(sumAttained) == self.gen + 1, "sumAttained is %i, should be gen+1, %i." % (
        #    int(sumAttained), self.gen + 1)
        probAttained = []
        for i in range(len(nAttained)):
            probAttained.append(100.0 * float(nAttained[i]) / sumAttained)
        if math.fabs(sum(probAttained) - 100.0 > 1e-13):
            raise Glitch, 'bad sum of attained proposal probs. %s' % sum(probAttained)

        spacer = ' ' * 4
        print "\nProposal probabilities (%)"
        #print "There are %i proposals" % len(self.proposals)
        print "For %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen)
        print "%2s %11s %11s  %11s %10s %18s %5s %5s" % ('', 'nProposals', 'proposed(%)',
                                                    'accepted(%)', 'tuning', 'proposal', 'part', 'num')
        for i in range(len(self.proposals)):
            print "%2i" % i,
            p = self.proposals[i]
            print "   %7i " % self.proposals[i].nProposals[0],
            print "   %5.1f    " % probAttained[i],
            if nAttained[i]:
                print "   %5.1f   " % (100.0 * float(nAccepted[i]) / float(nAttained[i])),
            else:
                print "       -   ",

            if p.tuning == None:
                print "       -    ",
            elif p.tuning < 2.0:
                print "   %7.3f  " % p.tuning,
            else:
                print "   %7.1f  " % p.tuning,

            print " %15s" % p.name,
            if p.pNum != -1:
                print " %3i " % p.pNum,
            else:
                print "   - ",
            if p.mtNum != -1:
                print " %3i " % p.mtNum,
            else:
                print "   - ",
            print


    def writeProposalIntendedProbs(self):
        """Tabulate the intended proposal probabilities.
        """

        nProposals = len(self.proposals)
        if not nProposals:
            print "STMcmc.writeProposalIntendedProbs().  No proposals (yet?)."
            return
        intended = self.propWeights[:]
        for i in range(len(intended)):
            intended[i] /= self.totalPropWeights
        if math.fabs(sum(intended) - 1.0 > 1e-14):
            raise Glitch, 'bad sum of intended proposal probs. %s' % sum(intended)

        spacer = ' ' * 4
        print "\nIntended proposal probabilities (%)"
        #print "There are %i proposals" % len(self.proposals)
        print "%2s %11s %18s %5s %5s" % ('', 'intended(%)', 'proposal', 'part', 'num')
        for i in range(len(self.proposals)):
            print "%2i" % i,
            p = self.proposals[i]
            print "   %6.2f    " % (100. * intended[i]),

            print " %15s" % p.name,
            if p.pNum != -1:
                print " %3i " % p.pNum,
            else:
                print "   - ",
            if p.mtNum != -1:
                print " %3i " % p.mtNum,
            else:
                print "   - ",
            print

class STMcmcCheckPointReader(object):
    """Read in and display mcmc_checkPoint files.

    Three options--
    
    To read in a specific checkpoint file, specify the file name by
    fName=whatever

    To read in the most recent (by os.path.getmtime()) checkpoint
    file, say last=True
    
    If you specify neither of the above, it will read in all the
    checkPoint files that it finds.

    Where it looks is determined by theGlob, which by default is '*',
    ie everything in the current directory.  If you want to look
    somewhere else, you can specify eg

        theGlob='SomeWhereElse/*' 

    or, if it is unambiguous, just

        theGlob='S*/*' 

    So you might say

        cpr = STMcmcCheckPointReader(theGlob='*_0.*')

    to get all the checkpoints from the first run, run 0.  Then, you
    can tell the cpr object to do various things.  Eg

        cpr.writeProposalAcceptances()

    But perhaps the most powerful thing about it is that it allows
    easy access to the checkpointed Mcmc objects, in the list mm.  Eg
    to get the first one, ask for

        m = cpr.mm[0]

    and m is an STMcmc object, complete with all its records of
    proposals and acceptances and so on.  And the TreePartitions
    object.  

    (Sorry!  -- Lazy documentation.  See the source code for more that it can do.)
    """
    
    def __init__(self, fName=None, theGlob='*', last=False, verbose=True):
        self.mm = []
        if not fName:
            #fList = [fName for fName in os.listdir(os.getcwd()) if fName.startswith("mcmc_checkPoint")]
            #fList = glob.glob(theGlob)
            #print "Full glob = %s" % fList
            fList = [fName for fName in glob.glob(theGlob) if 
                     os.path.basename(fName).startswith("mcmc_checkPoint")]
            #print fList
            if not fList:
                raise Glitch, "No checkpoints found in this directory."
            if last:
                # Find the most recent
                mostRecent = os.path.getmtime(fList[0])
                mostRecentFileName = fList[0]
                if len(fList) > 1:
                    for fName in fList[1:]:
                        mtime = os.path.getmtime(fName)
                        if mtime > mostRecent:
                            mostRecent = mtime
                            mostRecentFileName = fName
                f = file(mostRecentFileName)
                m = cPickle.load(f)
                f.close()
                self.mm.append(m)
                
            else:
                # get all the files
                for fName in fList:
                    f = file(fName)
                    m = cPickle.load(f)
                    f.close()
                    self.mm.append(m)

                self.mm = func.sortListOfObjectsOn2Attributes(self.mm, "gen", 'runNum')
        else:
            # get the file by name
            f = file(fName)
            m = cPickle.load(f)
            f.close()
            self.mm.append(m)
        if verbose:
            self.dump()

    def dump(self):
        print "STMcmcCheckPoints (%i checkPoints read)" % len(self.mm)
        print "%12s %12s %12s %12s" % (" ", "index", "run", "gen+1")
        print "%12s %12s %12s %12s" % (" ", "-----", "---", "-----")
        for i in range(len(self.mm)):
            m = self.mm[i]
            #print "    %2i    run %2i,  gen+1 %11i" % (i, m.runNum, m.gen+1)
            print "%12s %12s %12s %12s" % (" ", i, m.runNum, m.gen+1)
            

    def compareSplits(self, mNum1, mNum2, verbose=True, minimumProportion=0.1):
        """Should we be only looking at splits within the 95% ci of the topologies?"""
        m1 = self.mm[mNum1]
        m2 = self.mm[mNum2]
        tp1 = m1.treePartitions
        tp2 = m2.treePartitions
        
        if verbose:
            print "\nSTMcmcCheckPointReader.compareSplits(%i,%i)" % (mNum1, mNum2)
            print "%12s %12s %12s %12s %12s" % ("mNum", "runNum", "start", "gen+1", "nTrees")
            for i in range(5):
                print "   ---------",
            print
            for mNum in [mNum1, mNum2]:
                print " %10i " % mNum,
                m = self.mm[mNum]            
                print " %10i " % m.runNum,
                print " %10i " % (m.startMinusOne + 1),
                print " %10i " % (m.gen + 1),
                #for i in m.splitCompares:
                #    print i
                print " %10i " % m.treePartitions.nTrees

        asdos = self.compareSplitsBetweenTwoTreePartitions(tp1, tp2, minimumProportion, verbose=verbose)
        if asdos == None and verbose:
                print "No splits > %s" % minimumProportion
        return asdos

        
    def compareSplitsBetweenTwoTreePartitions(tp1, tp2, minimumProportion, verbose=False):
        ret = tp1.compareSplits(tp2, minimumProportion=minimumProportion)
        if ret != []:
            sumOfStdDevs = 0.0
            if ret and len(ret):
                nSplits = len(ret)
                for i in ret:
                    #print "            %.3f  %.3f    " % (i[2][0], i[2][1]),
                    stdDev = math.sqrt(func.variance(i[2]))
                    #print "%.5f" % stdDev
                    sumOfStdDevs += stdDev
                if verbose:
                    #print "  %f " % sumOfStdDevs,
                    print "     nSplits=%i, average of std devs of splits %.4f " % (nSplits, sumOfStdDevs/nSplits)
            return sumOfStdDevs/nSplits
        else:
            return None

    compareSplitsBetweenTwoTreePartitions = staticmethod(compareSplitsBetweenTwoTreePartitions)

    def compareSplitsAll(self):
        nM = len(self.mm)
        nItems = ((nM * nM) - nM)/2
        results = numpy.zeros((nM, nM), numpy.float)
        vect = numpy.zeros(nItems, numpy.float)
        vCounter = 0
        for mNum1 in range(1, nM):
            for mNum2 in range(mNum1):
                ret = self.compareSplits(mNum1, mNum2, verbose=False)
                #print "+++ ret = %s" % ret
                if ret == None:
                    ret = 0.0
                results[mNum1][mNum2] = ret
                results[mNum2][mNum1] = ret
                vect[vCounter] = ret
                vCounter += 1
                if 0:
                    print " %10i " % mNum1,
                    print " %10i " % mNum2,
                    print "%.3f" % ret
        print results
        
        print "For the %i values in one triangle," % nItems
        print "max =  ", vect.max()
        print "min =  ", vect.min()
        print "mean = ", vect.mean()
        print "var =  ", vect.var()
            
        
    def writeProposalAcceptances(self):
        for m in self.mm:
            m.writeProposalAcceptances()

    def writeSwapMatrices(self):
        for m in self.mm:
            if m.nChains > 1:
                m.writeSwapMatrix()
            
    def writeProposalProbs(self):
        for m in self.mm:
            m.writeProposalProbs()
