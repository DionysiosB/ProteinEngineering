#!/usr/bin/env python
from CHOMP import *
import sys,os,time
from numpy import log as ln

inputpdb = '1edmBH.pdb'

# load global data
atomtypes = AtomTypes()
amino_dir = AminoAcidDirectory(atomtypes)
if os.path.isfile('data_dir.binary'):
    data_dir = DataDirectory('data_dir.binary')
else:
    data_dir = DataDirectory(os.environ['CHOMPDATA'])
if os.path.isfile('dunbrack.binary'):
    dunbrack_lib = DunbrackLibrary(amino_dir,'dunbrack.binary')
else:
    dunbrack_lib = DunbrackLibrary(amino_dir,os.environ['DUNBRACK'])

EnergyFunction.SetPluginPath(os.environ['CHOMPDIST'])
energy_function = EnergyFunction(os.environ['CHOMP_EF'], atomtypes, data_dir, dunbrack_lib)
 
# load a system and score it
s = System(atomtypes, amino_dir)
s.Parse(open(inputpdb).read())

print 'generating rotamers...'
n = NeighborMap(s, data_dir.PairCutoffTable)
rl = RotamerLibrary(s)
rg = RotamerGenerator(n, dunbrack_lib)
rl.GenerateRotamers(rg)
ri = RotamerIterationMap(rl)

print 'filling energy graph'
eg = EnergyGraph()
energy_function.FillEnergyGraph(rl, n, ri, eg)


class BeliefPropagator:
        "Contains the IterateSumProduct method and a bunch of useful lookup lists, dictionaries"

	def normalize(self,vec):
	    tot = sum(vec)
	    if tot>0:
	        return list( numpy.array(vec)/tot )
	    else:
	        print 'Cannot normalize',vec
	        sys.exit()

	def P(self,i,j=-1):
	    if j>=0:
	        try:
	            return numpy.exp(-self.eg.GetEdgeEnergy(i,j))
	        except RuntimeError:
	            #print 'Did not find GetEdgeEnergy(%d,%d)' % (i,j)
	            return 1.0
	    else:
	        try:
	            return numpy.exp(-self.eg.GetVertexEnergy(i))
	        except RuntimeError:
	            print 'Failed to GetVertexEnergy(%d)' % i
	            sys.exit()

	def max_rot_rot_interaction(self,resA,resB):
	    """Find the interaction between rotamers of resA and rotamers of resB with the largest amplitude"""
	    A,B = self.res2rots[resA],self.res2rots[resB]
	    max = -1
	    for a in A:
	        #print 'Checking rotamer',a
	        partners = [x.first for x in self.eg.GetVertexEdges(a)]
	        for p in partners:
	            if p in B:
	                #print 'against',p
	                if abs(self.P(a,p)) > max:
	                    max = abs(self.P(a,p))
	    return max


        def __init__(self,rl,n,eg):
                self.rl = rl
                self.n = n
                self.eg = eg

                self.rotids = eg.GetVertexIDs()
                self.resids = [r.ID for r in rl.GetSystems()[0]]
		self.residues = [r for r in rl.GetSystems()[0]]
		self.res2partners = dict( (r.ID,[]) for r in self.residues )
		self.res2rots = dict((r.ID,[x for x in range(rl.GetRotamerIDRange(r).first,rl.GetRotamerIDRange(r).second+1)]) for r in self.residues)
		self.res2neighbors = dict((r.ID,[x for x in n.GetNeighborIDs(r)]) for r in self.residues)
		self.rot2res = {}

		for rid in [r.ID for r in self.residues]:
		    for rot in self.res2rots[rid]:
		        self.rot2res[rot]=rid

		self.resedges = []
		self.rotedges = []
		for r in self.res2neighbors.keys():
		    for neighbor in self.res2neighbors[r]:
		        if self.max_rot_rot_interaction(r,neighbor)>0:
		            if (r,neighbor) in self.resedges: continue
		            self.resedges.append((r,neighbor))
		            self.res2partners[r].extend([neighbor])
		 
		            if (neighbor,r) in self.resedges: continue
		            self.resedges.append((neighbor,r))
		            self.res2partners[neighbor].extend([r])
	 
		            for rot1 in self.res2rots[r]:
		               for rot2 in self.res2rots[neighbor]:
		                    self.rotedges.append((rot1,rot2))
		self.res2q = dict( (resid,len(self.res2partners[resid])) for resid in self.resids)

	def sharpen_beliefs(self,beliefs):
	    maxbelief = max(beliefs)
	    nummax = beliefs.count(maxbelief)
	    newbeliefs = []
	    for belief in beliefs:
	        if belief==maxbelief:
	            newbeliefs.append(1.0/nummax)
	        else:
	            newbeliefs.append(0)
	    return newbeliefs
	 
	def beliefs2rotamer(self,resid,beliefs):
	    assignment = self.sharpen_beliefs(beliefs)
	    somerotids = self.res2rots[resid]
	    ## Break ties in favor of the lower rotamer index. Dangerous!
	    for index,rotid in enumerate(somerotids):
	        if assignment[index]>0:
	            return rotid
	    print 'failed to find an assignment for',resid,beliefs,assignment

	def calc_beliefs(self,i,messages):
	    """Calculate the node beliefs given the potential for each choice and the current messages"""
	    beliefs = []
	    incoming_messages = [messages[(j,i)] for j in self.res2partners[i]]
	    incoming_product  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages)
	    for xii,xi in enumerate(self.res2rots[i]):
	        beliefs.append( self.P(xi)*incoming_product[xii] )
	    return self.normalize(beliefs)
	 
	def sumproduct_message(self,i,j,messages):
	    """Calculate an updated message from node i to node j, given the graph potentials and the current messages"""
	    new_message = []
	    incoming_messages = [messages[(k,i)] for k in self.res2partners[i] if k != j]
	    incoming_product  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages)
	    for xj in self.res2rots[j]:
	        sumproduct = sum([self.P(xi)*self.P(xi,xj)*incoming_product[xii]  for xii,xi in enumerate(self.res2rots[i])])
	        new_message.append(sumproduct)
	    return self.normalize(new_message)
	
	def sumproduct_round(self,messages):
	    """Inefficient, random flood of edge updates"""
	    newmessages = {}
	    for i,j in messages.keys():
	        newmessages[(i,j)] = self.sumproduct_message(i,j,messages)
	    return newmessages
	 
        def IterateSumProduct(self,mg,iteration_limit=10):
		observed = []
		for iter in range(iteration_limit):
		    print 'Sum Product Belief Propagation: iteration',iter
		    mg.messages = self.sumproduct_round(mg.messages)

		    ## Just here to decide to stop early
		    assignments = [self.beliefs2rotamer(resid,self.calc_beliefs(resid,mg.messages)) for resid in self.resids]
		    print assignments,eg.ComputeEnergy(IntVector(assignments))
		    if assignments in observed:
		        print 'Already found this solution. Terminating'
		        break
		    observed.append(assignments)
 		return mg

class MessageGraph:
	def __init__(self,bp):
		self.bp = bp
		self.messages = {}
		for i,j in bp.resedges:
			self.messages[(i,j)] = bp.normalize([1.0 for x in bp.res2rots[j]])

	def calc_rot_beliefs(self,beliefs):
	    rot_beliefs = {}
	    for i in beliefs.keys():
	        for ii,irot in enumerate(self.bp.res2rots[i]):
	            rot_beliefs[irot] = beliefs[i][ii]
	    return rot_beliefs
	 
	def calc_rot_rot_beliefs(self,res_res_beliefs):
	    rotrot_beliefs = {}
	    for i,j in res_res_beliefs.keys():
	        for ii,irot in enumerate(self.bp.res2rots[i]):
	            for jj,jrot in enumerate(self.bp.res2rots[j]):
	                rotrot_beliefs[(irot,jrot)] = res_res_beliefs[(i,j)][ii][jj]
	    return rotrot_beliefs
		
 
	def calc_single_beliefs(self):
	    """Calculate beliefs given the potential for each choice and the current messages"""
	    beliefs = dict((i,[]) for i in self.bp.resids)
	    for i in self.bp.resids:
	        incoming_messages = [self.messages[(j,i)] for j in self.bp.res2partners[i]]
	        incoming_product  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages)
	        for xii,xi in enumerate(self.bp.res2rots[i]):
	            beliefs[i].append( self.bp.P(xi)*incoming_product[xii] )
	        beliefs[i] /= sum(beliefs[i])
	    return beliefs
 
	def calc_two_node_beliefs(self):
	    """Calculate two-node beliefs given the potentials and current messages"""
	    resid_pairs = [(i,j) for i in self.bp.resids for j in self.bp.resids]
	    beliefs = dict( (pair,[]) for pair in resid_pairs)
	    for i,j in resid_pairs:
	        incoming_messages_i = [self.messages[(k,i)] for k in self.bp.res2partners[i] if k != j]
	        incoming_product_i  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages_i)
	        incoming_messages_j = [self.messages[(l,j)] for l in self.bp.res2partners[j] if l != i]
	        incoming_product_j  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages_j)
	        for xii,xi in enumerate(self.bp.res2rots[i]):
	            row = []
	            for xjj,xj in enumerate(self.bp.res2rots[j]):
	                row.append( self.bp.P(xi)*self.bp.P(xj)*self.bp.P(xi,xj)*incoming_product_i[xii]*incoming_product_j[xjj] )
	            beliefs[(i,j)].append( row[:] )
	        beliefs[(i,j)] /= sum(sum(numpy.array(beliefs[(i,j)])))
	    return beliefs
	 
	def CalcBeliefGraph(self):
	    B = BeliefGraph()
	    singles = self.calc_single_beliefs()
	    pairs = self.calc_two_node_beliefs()
	    singles_rot = self.calc_rot_beliefs(singles)
	    pairs_rot = self.calc_rot_rot_beliefs(pairs)
	    for x in singles_rot.keys(): B.graph[x] = singles_rot[x]
	    for x in pairs_rot.keys(): B.graph[x] = pairs_rot[x]
	    return B
 

class BeliefGraph:
	def __init__(self):
		self.graph = {}

	def MeanFieldApproximation(self,bp):
		epsilon=1e-100
		B = self.graph
		edgeE = {}
		for edge in bp.rotedges:
			try:
				edgeE[edge] = bp.eg.GetEdgeEnergy(*edge)
			except:
				edgeE[edge] = 0

		MFpair = sum([ (B[edge[0]]*B[edge[1]]* edgeE[edge]) for edge in bp.rotedges ])
		MFsingle = sum([ B[rotid]*bp.eg.GetVertexEnergy(rotid) for rotid in bp.rotids ])
		MFentropy = -sum([B[rotid]*ln(B[rotid]+epsilon) for rotid in bp.rotids ])
		return MFpair+MFsingle, MFentropy
 
	def BetheApproximation(self,bp):
		epsilon=1e-100
		B = self.graph
		edgeE = {}
		for edge in bp.rotedges:
			try:
				edgeE[edge] = bp.eg.GetEdgeEnergy(*edge)
			except:
				edgeE[edge] = 0

		Bethepair = sum([ (B[edge[0],edge[1]]* edgeE[edge]) for edge in bp.rotedges ])
		Bethesingles = sum([ B[rotid]*bp.eg.GetVertexEnergy(rotid) for rotid in bp.rotids ])
		BetheSpairs = -sum([ B[edge]*ln(B[edge]+epsilon) for edge in bp.rotedges ])
		BetheSsingles = sum([ (bp.res2q[resid]-1)*sum([ B[rotid]*ln(B[rotid]+epsilon) for rotid in bp.res2rots[resid] ]) for resid in bp.resids ])
		return Bethepair+Bethesingles, BetheSpairs+BetheSsingles



bp = BeliefPropagator(rl,n,eg)
mg = MessageGraph(bp)
mg = bp.IterateSumProduct(mg, iteration_limit=8)
bg = mg.CalcBeliefGraph()
U,S = bg.MeanFieldApproximation(bp)
print 'Mean Field',U,S,U-S
U,S = bg.BetheApproximation(bp)
print 'Bethe',U,S,U-S

