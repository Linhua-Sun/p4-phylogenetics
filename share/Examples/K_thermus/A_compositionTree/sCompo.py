read('../noTRuberNoGapsNoAmbiguities.nex')
a = var.alignments[0]
d = Data()
print 'Using all sites ...'
d.compoChiSquaredTest(verbose=1)
a.setNexusSets()
a.excludeCharSet('constant')
d = Data()
print 'After constant sites removal ...'
d.compoChiSquaredTest(verbose=1)
