import PyCHOMP, numpy
from PyCHOMP import *

#///////////////////////////////////////////////////////////////////////////
#//////                      CONTAINER ADAPTERS                      ///////
#///////////////////////////////////////////////////////////////////////////

# TODO: deal with IntVectorFloatPairVector

vectors, vecs = [], []

def get_el_type(item):
	split = item.split('Vector')
	
	# ending with "2" means it's a vector of vectors
	if split[1] == '2':
		return PyCHOMP.__dict__[item[:-1]]
	else:
		d = {
			'Int' 		: int,
			'Float' 	: float,
			'Bool'		: bool,
			'Str'		: str,
			'ConstSystem'	: PyCHOMP.System
		}
		try:
			return d[split[0]]
		except:
			return PyCHOMP.__dict__[split[0]]

for item, cls in PyCHOMP.__dict__.items():
	if item in ['EGVertexVector']:
		continue

	if item.endswith('Vector') or item.endswith('Vector2'):
		vectors.append(cls)
		cls.el_type = get_el_type(item)

	elif item.startswith('Vec') or item.startswith('IVec'):
		vecs.append(cls)
		if item[0] == 'I':
			cls.el_type = int
		else:
			cls.el_type = float

def vector_fill(self, args, old_init):
	old_init(self)
	for x in args: self.append(x)
	
def vec_fill(self, args, old_init):
	old_init(self, *args)

def wrapped_ctor(old_init, fill):
	def call(obj, *args):
		try:
			# first, try converting the full argument list
			cargs = [obj.el_type(x) for x in args]
			fill(obj, cargs, old_init)
		except:
			try:
				# maybe we can convert just the first argument
				cargs = [obj.el_type(x) for x in args[0]]
				fill(obj, cargs, old_init)
			except Exception,e:
				# no luck; just do the regular constructor
				old_init(obj, *args)
	return call

def iter_vec(v):
	for i in xrange(v.length()):
		yield v[i]
	raise StopIteration	

def str_vector(v):
	def format(sep, pre, start, end):
		ss = [pre + str(x) for x in v]
		return v.__class__.__name__ + ': ' + start + sep.join(ss) + end
	
	s = format(', ', '', '[', ']')
	if len(s) > 80:
		s = format(',\n', '\t', '[\n', '\n]')

	return s

for vector in vectors:
	vector.__init__ = wrapped_ctor(vector.__init__, vector_fill)
	vector.__str__ = str_vector
	vector.__repr__ = str_vector

for vec in vecs:
	vec.__init__ = wrapped_ctor(vec.__init__, vec_fill)
	vec.__iter__ = iter_vec
	vec.__len__ = vec.length
	vec.__str__ = str_vector
	vec.__repr__ = str_vector

#Overloading to allow easy addition, subtraction of coordinates
class Vec3(Vec3):
	def __sub__(self,other):
		return Vec3(self[0]-other[0],self[1]-other[1],self[2]-other[2])
	def __add__(self,other):
		return Vec3(self[0]+other[0],self[1]+other[1],self[2]+other[2])
	def __mul__(self,other):
		return Vec3(self[0]*other[0],self[1]*other[1],self[2]*other[2])
	def __div__(self,other):
		return Vec3(self[0]/other[0],self[1]/other[1],self[2]/other[2])
	def __pow__(self,other):
		return Vec3(self[0]**other,self[1]**other,self[2]**other)

#///////////////////////////////////////////////////////////////////////////
#//////                         METHOD ADAPTERS                      ///////
#///////////////////////////////////////////////////////////////////////////

def subst_arg(arg, container_type):
	container = container_type()
	for el in arg:
		container.append(el)
	return container

old_System_Join = System.Join
new_System_Join = lambda cls, sv: old_System_Join(subst_arg(sv, ConstSystemVector))
new_System_Join.__doc__ = old_System_Join.__doc__
System.Join = classmethod(new_System_Join)

old_NeighborMap_Create = NeighborMap.Create
def new_NeighborMap_Create(self, systems, cutoffs):
	if isinstance(systems, list) or isinstance(systems, SystemVector):
		return old_NeighborMap_Create(subst_arg(systems, ConstSystemVector), cutoffs)
	return old_NeighborMap_Create(systems, cutoffs)
new_NeighborMap_Create.__doc__ = old_NeighborMap_Create.__doc__
NeighborMap.Create = classmethod(new_NeighborMap_Create)

old_RotamerLibrary_Create = RotamerLibrary.Create
def new_RotamerLibrary_Create(self, systems):
	if isinstance(systems, list) or isinstance(systems, SystemVector):
		return old_RotamerLibrary_Create(subst_arg(systems, ConstSystemVector))
	return old_RotamerLibrary_Create(systems)
new_RotamerLibrary_Create.__doc__ = old_RotamerLibrary_Create.__doc__
RotamerLibrary.Create = classmethod(new_RotamerLibrary_Create)

old_Residue_SetCoords = Residue.SetAtomCoords
def new_Residue_SetCoords(self, i, coords):
	coords = Vec3(numpy.array(coords, dtype=numpy.float32))
	return old_Residue_SetCoords(self, i, coords)
new_Residue_SetCoords.__doc__ = old_Residue_SetCoords.__doc__
Residue.SetAtomCoords = new_Residue_SetCoords


#///////////////////////////////////////////////////////////////////////////
#//////                     CONVENIENCE FUNCTIONS                     //////
#///////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////
#//////                           AtomTypes                           //////
#///////////////////////////////////////////////////////////////////////////

def AtomTypesIndexer(self, args):
	item, typesys = args
	if isinstance(item, str):
		return self.ToType(item, typesys)
	else:
		return self.ToString(item, typesys)

AtomTypes.__getitem__ = AtomTypesIndexer	# a.ToType(s, sys) becomes a[s,sys] and 
						# a.ToString(i, sys) becomes a[i,sys]

#///////////////////////////////////////////////////////////////////////////
#//////                      AminoAcidDirectory                      ///////
#///////////////////////////////////////////////////////////////////////////

AminoAcidDirectory.__getitem__ = AminoAcidDirectory.Get # aad.Get('ALA') becomes aad['ALA']

#///////////////////////////////////////////////////////////////////////////
#//////                            System                            ///////
#///////////////////////////////////////////////////////////////////////////

def new_System_getitem(self, n):
	try:
		return self.GetResidueByStorageIndex(n)
	except:
		raise IndexError()

System.__getitem__ = new_System_getitem
System.__len__     = lambda self: self.size

#//////////////////////////////////////////////////////////////////////////
#/////                          CONSTRUCTORS                          /////
#//////////////////////////////////////////////////////////////////////////

# turn Create methods into constructors for convenience
create_fns = {}

for cls in PyCHOMP.__dict__.values():
	try: 	
		create_fns[cls] = cls.Create
		del cls.Create

		# make __new__ just call the Create function
		cls.__new__  = classmethod( lambda c, *args: create_fns[c](*args[1:]) )

		# the default __init__ for these classes throws an exception, so replace it
		cls.__init__ = lambda self, *args: None

	# ignore objects with no Create method
	except AttributeError: pass
