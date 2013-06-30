#!/usr/bin/env python
from CHOMP import *
import sys,os,random,time
from numpy import log as ln
from time import time

def P(i,j=-1):
    if j>=0:
        try:
            return numpy.exp(-eg.GetEdgeEnergy(i,j))
        except RuntimeError:
            #print 'Did not find GetEdgeEnergy(%d,%d)' % (i,j)
            return 1.0
    else:
        try:
            return numpy.exp(-eg.GetVertexEnergy(i))
        except RuntimeError:
            print 'Failed to GetVertexEnergy(%d)' % i
            sys.exit()


try:
    inputpdb = sys.argv[1]
except:
    print 'Provide pdb'
    sys.exit()

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

activeresids = [r.ID for r in s]
activeres = [s.GetResidue(resid) for resid in activeresids]

print 'generating rotamers...'
n = NeighborMap(s, data_dir.PairCutoffTable)
rl = RotamerLibrary(s)
rg = RotamerGenerator(n, dunbrack_lib)
for res in activeres:
    rl.GenerateRotamers(res,rg)
ri = RotamerIterationMap(rl)

size = 1
for dim in range(rl.NumDimensions()):
    r = rl.DimensionIDRange(dim)
    size *= (r.second - r.first + 1)
print 'problem size',size

print 'filling energy graph'
eg = EnergyGraph()
energy_function.FillEnergyGraph(rl, n, ri, eg)


def bp_setup_lookups(rl,n):
    """Construct useful lookups residue->rotamers, residue->interacting neighbors, and an edge list"""
    residues = [r for r in rl.GetSystems()[0]]
    res2rots = dict((r.ID,[x for x in range(rl.GetRotamerIDRange(r).first,rl.GetRotamerIDRange(r).second+1)]) for r in residues)
    res2neighbors = dict((r.ID,[x for x in n.GetNeighborIDs(r)]) for r in residues)
    rot2res = {}
    for rid in [r.ID for r in residues]:
        for rot in res2rots[rid]:
            rot2res[rot]=rid
    res2partners = dict( (r.ID,[]) for r in residues )
    resedges = []
    rotedges = []
    for r in res2neighbors.keys():
        for neighbor in res2neighbors[r]:
            if max_rot_rot_interaction(r,neighbor,res2rots)>0:
                if (r,neighbor) in resedges: continue
                resedges.append((r,neighbor))
                res2partners[r].extend([neighbor])

                if (neighbor,r) in resedges: continue
                resedges.append((neighbor,r))
                res2partners[neighbor].extend([r])

                for rot1 in res2rots[r]:
                    for rot2 in res2rots[neighbor]:
                        rotedges.append((rot1,rot2))
    return res2rots,res2partners,resedges,rot2res,rotedges

def max_rot_rot_interaction(resA,resB,res2rots):
    """Find the interaction between rotamers of resA and rotamers of resB with the largest amplitude"""
    A,B = res2rots[resA],res2rots[resB]
    max = -1
    for a in A:
        #print 'Checking rotamer',a
        partners = [x.first for x in eg.GetVertexEdges(a)]
        for p in partners:
            if p in B:
                #print 'against',p
                if abs(P(a,p)) > max: 
                    max = abs(P(a,p))
    return max            

def normalize(vec):
    tot = sum(vec)
    if tot>0:
        return list( numpy.array(vec)/tot )
    else:
        print 'Cannot normalize',vec
        sys.exit()
        return vec

def uniform_messages():
    """Define a uniform message for each edge for initialization purposes"""
    messages = {}
    for i,j in resedges:
        messages[(i,j)] = normalize([1.0 for x in res2rots[j]])
    return messages

def sumproduct_message(i,j,messages):
    """Calculate an updated message from node i to node j, given the graph potentials and the current messages"""
    new_message = []
    incoming_messages = [messages[(k,i)] for k in res2partners[i] if k != j]
    incoming_product  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages)
    for xj in res2rots[j]:
        sumproduct = sum([P(xi)*P(xi,xj)*incoming_product[xii]  for xii,xi in enumerate(res2rots[i])])
        new_message.append(sumproduct) 
    return normalize(new_message)

def sumproduct_round(messages):
    """Inefficient, random flood of edge updates"""
    newmessages = {}
    for i,j in messages.keys():
        newmessages[(i,j)] = sumproduct_message(i,j,messages)
    return newmessages
    
def calc_beliefs(i,messages):
    """Calculate the node beliefs given the potential for each choice and the current messages"""
    beliefs = []
    incoming_messages = [messages[(j,i)] for j in res2partners[i]]
    incoming_product  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages)
    for xii,xi in enumerate(res2rots[i]):
        beliefs.append( P(xi)*incoming_product[xii] ) 
    return normalize(beliefs)

def calc_single_beliefs(messages):
    """Calculate beliefs given the potential for each choice and the current messages"""
    beliefs = dict((i,[]) for i in resids) 
    for i in resids:
        incoming_messages = [messages[(j,i)] for j in res2partners[i]]
        incoming_product  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages)
        for xii,xi in enumerate(res2rots[i]):
            beliefs[i].append( P(xi)*incoming_product[xii] ) 
        beliefs[i] /= sum(beliefs[i])
    return beliefs

def calc_rot_beliefs(beliefs):
    rot_beliefs = {}
    for i in beliefs.keys():
        for ii,irot in enumerate(res2rots[i]):
            rot_beliefs[irot] = beliefs[i][ii]
    return rot_beliefs

def calc_two_node_beliefs(messages):
    """Calculate two-node beliefs given the potentials and current messages"""
    resid_pairs = [(i,j) for i in resids for j in resids]
    beliefs = dict( (pair,[]) for pair in resid_pairs) 
    for i,j in resid_pairs:
        incoming_messages_i = [messages[(k,i)] for k in res2partners[i] if k != j]
        incoming_product_i  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages_i)
        incoming_messages_j = [messages[(l,j)] for l in res2partners[j] if l != i]
        incoming_product_j  = reduce(lambda x,y: numpy.array(x)*numpy.array(y), incoming_messages_j)
        for xii,xi in enumerate(res2rots[i]):
            row = []
            for xjj,xj in enumerate(res2rots[j]):
                row.append( P(xi)*P(xj)*P(xi,xj)*incoming_product_i[xii]*incoming_product_j[xjj] ) 
            beliefs[(i,j)].append( row[:] ) 
        beliefs[(i,j)] /= sum(sum(numpy.array(beliefs[(i,j)])))
    return beliefs

def calc_all_beliefs(messages):
    B = {}
    singles = calc_single_beliefs(messages)
    pairs = calc_two_node_beliefs(messages)
    singles_rot = calc_rot_beliefs(singles)
    pairs_rot = calc_rot_rot_beliefs(pairs)
    for x in singles_rot.keys(): B[x] = singles_rot[x]
    for x in pairs_rot.keys(): B[x] = pairs_rot[x]
    return B

def calc_rot_rot_beliefs(res_res_beliefs):
    rotrot_beliefs = {}
    for i,j in res_res_beliefs.keys():
        for ii,irot in enumerate(res2rots[i]):
            for jj,jrot in enumerate(res2rots[j]):
                rotrot_beliefs[(irot,jrot)] = res_res_beliefs[(i,j)][ii][jj]
    return rotrot_beliefs

def sharpen_beliefs(beliefs):
    maxbelief = max(beliefs)
    nummax = beliefs.count(maxbelief)
    newbeliefs = []
    for belief in beliefs:
        if belief==maxbelief:
            newbeliefs.append(1.0/nummax)
        else:
            newbeliefs.append(0)
    return newbeliefs

def beliefs2rotamer(resid,beliefs):
    assignment = sharpen_beliefs(beliefs)
    rotids = res2rots[resid]
    ## Break ties in favor of the lower rotamer index. Dangerous!
    for index,rotid in enumerate(rotids):
        if assignment[index]>0:
            return rotid
    print 'failed to find an assignment for',resid,beliefs,assignment




################### Setup ##########################################
before = time()
res2rots,res2partners,resedges,rot2res,rotedges = bp_setup_lookups(rl,n)
rotids = eg.GetVertexIDs()
resids = [r.ID for r in rl.GetSystems()[0]] 
resedgeindices = [(resids.index(x[0]),resids.index(x[1])) for x in resedges]
res2q = dict( (resid,len(res2partners[resid])) for resid in resids)

messages = uniform_messages()


################### Finding the solution ###########################
observed = []
for iter in range(5):
    print 'Sum Product Belief Propagation: iteration',iter
    messages = sumproduct_round(messages)
    
    assignments = [beliefs2rotamer(resid,calc_beliefs(resid,messages)) for resid in resids]
    print assignments,eg.ComputeEnergy(IntVector(assignments))       
    if assignments in observed:
        print 'Already found this solution. Terminating'
        break
    observed.append(assignments)


################### Free energy approximations #####################
B = calc_all_beliefs(messages)

epsilon = 1e-100

edgeE = {}
for edge in rotedges:
    try:
        edgeE[edge] = eg.GetEdgeEnergy(*edge)
    except:
        edgeE[edge] = 0

MFpair = sum([ (B[edge[0]]*B[edge[1]]*edgeE[edge]) for edge in rotedges ])
MFsingle = sum([ (B[rotid] * eg.GetVertexEnergy(rotid)) for rotid in rotids ])
MFentropy = -sum([B[rotid]*ln(B[rotid]+epsilon) for rotid in rotids])

Bethepair = sum([ (B[edge[0],edge[1]]*edgeE[edge]) for edge in rotedges ])
Bethesingles = sum([ B[rotid]*eg.GetVertexEnergy(rotid) for rotid in rotids ])
Bethesingles2 = sum([ (res2q[resid]-1)*sum([ B[rotid]*eg.GetVertexEnergy(rotid) for rotid in res2rots[resid] ]) for resid in resids ])
BetheSpairs = -sum([ B[edge]*ln(B[edge]+epsilon) for edge in rotedges ])
BetheSsingles = sum([ (res2q[resid]-1)*sum([ B[rotid]*ln(B[rotid]+epsilon) for rotid in res2rots[resid] ]) for resid in resids ])
print MFpair,MFsingle,MFentropy
print 'Mean-field avg energy:',MFpair+MFsingle
print 'Mean-field entropy:',MFentropy
print 'Mean-field free energy:',MFpair+MFsingle-MFentropy

print 'Bethe avg energy:',Bethepair+Bethesingles
print 'Bethe entropy:',BetheSpairs+BetheSsingles
print 'Bethe free energy:',Bethepair+Bethesingles-(BetheSpairs+BetheSsingles)
print inputpdb, len(resids),time()-before
