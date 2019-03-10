# INTERSECTION_CONES_TOP_DIRECTORY = TODO

# Currently, everything you could possibly want should be produced using
# raw_cones
# true_decomposition_matrices
# blocked (both of the above folders in blocks)
# computations (all of the computational data that is produced from the above)


class rec(object):
	# Duplicates the rec class in GAP.
	def __init__(self, **kargs):
		self._keys = set()
		self._keys.remove('_keys')
		for (key, value) in kargs.iteritems():
			setattr(self, key, value)
			#self._keys.append(key)
	def __str__(self):
		return "rec(" + ", ".join(" = ".join( (key, str( getattr(self, key)))) for key in self._keys) + ")"
	def __repr__(self):
		return "rec(" + ", ".join(" = ".join( (key, repr(getattr(self, key)))) for key in self._keys) + ")"
	def __setattr__(self, at, val):
		super(rec, self).__setattr__(at, val)
		self._keys.add(at)
	def __delattr__(self, at):
		super(rec, self).__delattr__(at)
		self._keys.remove(at)
f = open("settings.txt", "r")
s = f.read()
f.close()
s = s[6:].strip(" ;\n").replace(":=", "=")
settings = eval(s)


def get_top_filename(p, n=None, core=None, length=None):
	"""
	(p, n, core, length) -> "p_**n_**core_**_length_**.txt"
	with **'s replaced with representations of respective terms.
	Also does "p_**n_**.txt".
	"""
	if n is None: return p
	if core is None: return "".join(["p_", str(p), "n_", str(n), ".txt"])
	return "".join(
			["p_", 
				str(p),
				"n_",
				str(n),
				"core_",
				"_".join(str(i) for i in core),
				"_length_",
				(str(length) if length not in (-1, "infinity") else "infinity"),
				".txt"])
def get_pncorelength_from_filename(filename):
	"""
	Inverse of get_top_filename().
	"""
	npos = filename.find("n", 2)
	p = int(filename[2:npos])
	corepos = filename.find("core", npos)
	if corepos == -1:
		n = int(filename[npos+2:filename.find(".", npos)])
		return (p, n)
	n = int(filename[npos+2:corepos])
	lengthpos = filename.find("length", corepos)
	dotpos = filename.find(".", corepos)
	if lengthpos == -1:
		core = [int(i) for i in filename[corepos+5:dotpos-1].split('_') if i != '']
		length = "infinity"
	else:
		core = [int(i) for i in filename[corepos+5:lengthpos-1].split('_') if i != '']
		length = filename[lengthpos+7:dotpos]
		if length != "infinity": length = int(length)
	return (p, n, core, length)

def matrix_to_list(m):
	return list(tuple(row) for row in m)

def cut_zero_columns_and_sort(m):
	if type(m) == list: m = matrix(m)
	m = m.transpose()
	new_matrix = matrix([row for row in m if not all(i == 0 for i in row)])
	li = matrix_to_list(new_matrix)
	li.sort()
	li.reverse()
	return matrix_to_list(matrix(li).transpose())

import gc, time
def block_truncate_raw_cones(p, n=None):
	# print "block_truncate_raw_cones(%s, %s)" % (p, n)
	if n is None:
		filename = p
		(p, n) = get_pncorelength_from_filename(filename)
	else:
		filename = get_top_filename(p,n)
	if n >= 25:
		# print "\tSkipping."
		return
	principal_core = n % p
	if principal_core == 0:
		principal_core = []
	else:
		principal_core = [principal_core]
	# Look for principal block.
	if os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/" + get_top_filename(
	 p, n, principal_core, -1)):
		# print "\tSkipping."
		return
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/actually_raw_cones/" + filename, "r")
	s = f.read()
	f.close()
	s = s[6:].strip(" ;\n").replace(":=", "=")
	r = eval(s)
	s = None # gc
	gc.collect()
	for core in r.parts.p_cores:
		# print "\tCore %s" % core
		core_indices = [i for i in range(len(r.P_p)) if Partition(r.parts.partitions[i]).core(p) == core]
		lower = core_indices[0]
		upper = core_indices[-1]+1
		assert core_indices == range(lower, upper)
		parts = r.parts.partitions[lower:upper]
		P_p = cut_zero_columns_and_sort(r.P_p[lower:upper])
		gc.collect()
		P_es = [cut_zero_columns_and_sort(P_e[lower:upper]) for P_e in r.P_es]
		gc.collect()
		assert len(P_p) == len(parts)
		assert all(len(P_e) == len(parts) for P_e in P_es)
		for length in ["infinity", 5, 4, 3, 2]:
			# print "\t\tLength %s" % length
			if length != "infinity":
				cutoff = len(tuple(part for part in parts if len(part) <= length))
				# print "\t\t %s rows" % cutoff
				P_p = cut_zero_columns_and_sort(P_p[:cutoff])
				gc.collect()
				P_es = [cut_zero_columns_and_sort(P_e[:cutoff]) for P_e in P_es]
				gc.collect()
				parts = parts[:cutoff]
				gc.collect()
				assert len(P_p) == len(parts)
				assert all(len(P_e) == len(parts) for P_e in P_es)
			else:
				pass # print "\t\t %s rows" % len(parts)
			f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/" + get_top_filename(
	 		 p, n, core, length), "w")
			f.write(repr(rec(P_p = P_p, P_es=P_es, parts=parts)))
			f.close()
			gc.collect()
def block_truncate_all_raw_cones():
	for filename in sorted(os.listdir("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/actually_raw_cones/"),
	 key = get_pncorelength_from_filename):
		block_truncate_raw_cones(filename)

def generate_order_function(p):
	def order_function(mu, nu):
		ecore_mu = mu.core(p)
		ecore_nu = nu.core(p)
		c = cmp(ecore_mu, ecore_nu)
		if c != 0:
			return True
		elif len(mu) != len(nu):
			return cmp(len(mu), len(nu))
		else:
			return cmp(nu, mu)
	return order_function
ordering = generate_order_function(2) # Just for example. This should be changed when p is selected.

def block_truncate_dec_mats(p, n=None):
	# print "block_truncate_dec_mats(%s, %s)" % (p, n)
	if n is None:
		filename = p
		(p, n) = get_pncorelength_from_filename(filename)
	else:
		filename = get_top_filename(p,n)
	if n >= 25:
		# print "\tSkipping."
		return
	principal_core = n % p
	if principal_core == 0:
		principal_core = []
	else:
		principal_core = [principal_core]
	# Look for principal block.
	if os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_decomposition_matrices/" + get_top_filename(
	 p, n, principal_core, -1)):
		# print "\tSkipping."
		return
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/true_decomposition_matrices/" + filename, "r")
	s = f.read()
	f.close()
	s = s[6:].strip(" ;\n").replace(":=", "=")
	r = eval(s)
	s = None # gc
	gc.collect()
	cores = list(set(Partition(part).core(p) for part in r.parts))
	cores.sort()
	for core in cores:
		# print "\tCore %s" % core
		core_indices = [i for i in range(len(r.decmat)) if Partition(r.parts[i]).core(p) == core]
		lower = core_indices[0]
		upper = core_indices[-1]+1
		assert core_indices == range(lower, upper)
		parts = r.parts[lower:upper]
		decmat = cut_zero_columns_and_sort(r.decmat[lower:upper])
		gc.collect()
		assert len(decmat) == len(parts)
		for length in ["infinity", 5, 4, 3, 2]:
			# print "\t\tLength %s" % length
			if length != "infinity":
				cutoff = len(tuple(part for part in parts if len(part) <= length))
				# print "\t\t %s rows" % cutoff
				decmat = cut_zero_columns_and_sort(decmat[:cutoff])
				gc.collect()
				parts = parts[:cutoff]
				gc.collect()
				assert len(decmat) == len(parts)
			else:
				pass # print "\t\t %s rows" % len(parts)
			f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_decomposition_matrices/" + get_top_filename(
	 		 p, n, core, length), "w")
			f.write(repr(rec(decmat = decmat, parts=parts)))
			f.close()
			gc.collect()
def block_truncate_all_dec_mats():
	for filename in sorted(os.listdir("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/true_decomposition_matrices/"),
	 key = get_pncorelength_from_filename):
		block_truncate_dec_mats(filename)

# Converts an arbitrary matrix of entries into a human-readable table in string format.
def get_nice_string_rep_of_matrix(m):
	m = [[str(e) for e in row] for row in m]
	lens = [max(map(len, col)) for col in zip(*m)]
	fmt = ' '.join('{{:{}}}'.format(x) for x in lens)
	table = [fmt.format(*row) for row in m]
	return '\n'.join(table)

# Assumes everything is given as a list of lists or list of list of lists (as opposed to a matrix or list of matrices).
def nonreduced_dual_cone_gens_raw(P_p, higher_order_cones, parts=None):
	k = len(P_p)
	if k == 0:
		if parts is not None:
			return rec(P_p = P_p, p_part = [], parts =parts)
		return rec(P_p = P_p, p_part = [])
	l_p = len(P_p[0])
	P_p = matrix(P_p)
	p_part = list(matrix.identity(l_p)) # The positive part of the dual cone in adjustment matrix form.

	#e = p*p
	for P_e in higher_order_cones:
		l_e = len(P_e[0])
		P_e = matrix(P_e)
		Q_e = P_e.solve_left(matrix.identity(l_e))
		p_part.extend(Q_e*P_p)
	p_part.sort();
	i = 1;
	while i < len(p_part):
		if p_part[i] == p_part[i-1]:
			del p_part[i];
		else:
			i += 1;
	if parts is not None:
		return rec(P_p = matrix_to_list(P_p), p_part = matrix_to_list(p_part), parts=parts)
	return rec(P_p = matrix_to_list(P_p), p_part = p_part)

# Dualizes r.p_part to obtain r.intersection
def dualize(r):
	if len(r.p_part) == 0:
		r.intersection = None
		r.intersection_adjustment_matrix = []
	poly = Polyhedron(ieqs = ((0,)+tuple(gen) for gen in r.p_part))
	r.intersection = poly
	r.intersection_adjustment_matrix = matrix_to_list(
	 matrix(list(reversed(sorted(tuple(row.vector()) for row in r.intersection.Vrepresentation() if row.is_ray())))).transpose())

def get_blocked_raw_cones(p, n=None, core=None, length=None):
	if n is not None:
		filename = "../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/" + get_top_filename(p, n, core, length)
	else:
		filename = p
	f = open(filename, "r")
	s = f.read()
	f.close()
	return eval(s)

def is_p_regular(p, partition):
	last_row = None
	last_row_count = 0
	for row in partition:
		if row == last_row:
			last_row_count += 1
			if last_row_count >= p:
				return False
		else:
			last_row = row
			last_row_count = 1
	return True
def dualize_raw_cones(p, n=None, core=None, length=None):
	# print "dualize_raw_cones(%s, %s, %s, %s)" %(p, n, core, length)
	if os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/dual_raw_cones/" + get_top_filename(
	 p, n, core, length)):
		# print "\tSkipping."
		return
	if n is None:
		(p, n, core, length) = get_pncorelength_from_filename(p)
	#if p <= 5: # TODO this is just to skip some hard cases.
	#	# print "\tSkipping."
	#	return
	if n > 20:
		# print "\tSkipping."
		return
	try:
		r = get_blocked_raw_cones(p, n, core, length)
	except:
		# print "\tFile retrieval Error."
		return
	try:
		r = nonreduced_dual_cone_gens_raw(r.P_p, r.P_es, r.parts)
	except:
		# print "\tError occurred in processing."
		return
	filename = get_top_filename(p, n, core, length)
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/dual_raw_cones/" + filename, "w")
	f.write(repr(r))
	f.close()
	# f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/displayed_intersection_raw_cones/" + filename, "w")
	# f.write("Intersection of cones in adjustment matrix form.\n")
	# f.write(get_nice_string_rep_of_matrix([ ( (r.preg_parts[i],) + r.intersection_adjustment_matrix[i])
	#  for i in range(len(r.preg_parts))]))
	# f.write("\nHecke algebra decomposition matrix.\n")
	# f.write(get_nice_string_rep_of_matrix([ ( (r.parts[i],) + r.P_p[i])
	#  for i in range(len(r.parts))]))
	# f.close()
def dualize_all_raw_cones():
	for filename in sorted(os.listdir("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/"),
	 key = get_pncorelength_from_filename):
		# print(filename)
		dualize_raw_cones(filename)
		# print("Garbage collecting...")
		gc.collect()
def get_dualized_raw_cones(p, n=None, core=None, length=None):
	if n is not None:
		filename = "../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/dual_raw_cones/" + get_top_filename(
		 p, n, core, length)
	else:
		filename = "../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/dual_raw_cones/" + p
	f = open(filename, "r")
	s = f.read()
	f.close()
	return eval(s)
def intersect_raw_cones(p, n=None, core=None, length=None):
	# print "intersect_raw_cones(%s, %s, %s, %s)" %(p, n, core, length)
	if os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/intersection_raw_cones/" + get_top_filename(
	 p, n, core, length)):
		# print "\tSkipping."
		return
	if n is None:
		(p, n, core, length) = get_pncorelength_from_filename(p)
	#if p <= 5: # TODO this is just to skip some hard cases.
	#	print "\tSkipping."
	#	return
	if n > 20:
		# print "\tSkipping."
		return
	try:
		r = get_dualized_raw_cones(p, n, core, length)
	except:
		# print "\tFile retrieval Error."
		return
	try:
		dualize(r)
		del r.intersection
		r.preg_parts = [part for part in r.parts if is_p_regular(p, part)]
	except:
		# print "\tError occurred in processing."
		return
	assert len(r.preg_parts) == len(r.intersection_adjustment_matrix)
	filename = get_top_filename(p, n, core, length)
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/intersection_raw_cones/" + filename, "w")
	f.write(repr(r))
	f.close()
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/displayed_intersection_raw_cones/" + filename, "w")
	f.write("Intersection of cones in adjustment matrix form.\n")
	f.write(get_nice_string_rep_of_matrix([ ( (r.preg_parts[i],) + r.intersection_adjustment_matrix[i])
	 for i in range(len(r.preg_parts))]))
	f.write("\nHecke algebra decomposition matrix.\n")
	f.write(get_nice_string_rep_of_matrix([ ( (r.parts[i],) + r.P_p[i])
	 for i in range(len(r.parts))]))
	f.close()

def intersect_all_raw_cones():
	for filename in sorted(os.listdir("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/"),
	 key = get_pncorelength_from_filename):
		# print(filename)
		intersect_raw_cones(filename)
		# print("Garbage collecting...")
		gc.collect()

def get_part(partition):
	if isinstance(partition, Partition): return partition
	else: return Partition(partition)

def get_indexing_set(p, n, core, length, regular = True):
	if regular is True: regular = p
	core = get_part(core)
	ordering = generate_order_function(p)
	if regular:
		if length not in (-1, "infinity"):
			parts = [part for part in Partitions(n, regular=regular) if (len(part) <= length and part.core(p) == core)]
		else:
			parts = [part for part in Partitions(n, regular=regular) if (part.core(p) == core)]
	else:
		if length not in (-1, "infinity"):
			parts = [part for part in Partitions(n) if (len(part) <= length and part.core(p) == core)]
		else:
			parts = [part for part in Partitions(n) if (part.core(p) == core)]
	parts.sort(cmp=ordering)
	return parts

# Relevant partition tools
# corners_residue
# remove_cell
# (seems like corners and inside corners are the same thing)
# addable_cells
# add_cell
# addable_cells_residue


# Compute the matrix corresponding to induction f_i where i=residue.
def get_induct_matrix_to(p, n, core, length, residue):
	# It is easier to create rows because that's how our matrices work.
	# This corresponds to restricting partitions.
	index_upstairs = get_indexing_set(p, n, core, length, False)
	for part in index_upstairs:
		remove = part.corners_residue(residue, p)
		if len(remove)>0:
			core_downstairs = part.remove_cell(*remove[0]).core(p)
			break
	else:
		return None
	index_downstairs = get_indexing_set(p, n-1, core_downstairs, length, False)
	return rec(
		induct = [ [ 1 if part_downstairs in set(part_upstairs.remove_cell(*cell) for cell in part_upstairs.corners_residue(residue, p))
 		 else 0
 			for part_downstairs in index_downstairs]
 			for part_upstairs in index_upstairs],
		index_upstairs = index_upstairs,
		index_downstairs = index_downstairs,
		core_upstairs = core,
		core_downstairs = core_downstairs)

def get_blocked_dec_mat(p, n, core, length):
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_decomposition_matrices/" + 
	 get_top_filename(p, n, core, length), "r")
	s = f.read()
	f.close()
	return eval(s)

def get_intersected_cones(p, n, core, length):
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/intersection_raw_cones/" + 
	 get_top_filename(p, n, core, length), "r")
	s = f.read()
	f.close()
	return eval(s)

def induce_true_dec_mats_to(p, n=None, core=None, length=None):
	# print "induce_true_dec_mats_to(%s, %s, %s, %s)" %(p, n, core, length)
	if n is None:
		(p, n, core, length) = get_pncorelength_from_filename(p)
	if not os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/true_decomposition_matrices/" + 
	 get_top_filename(p, n-1)):
		# print "\tSkipping (don't have true dec mat)."
		return
	filename = get_top_filename(p, n, core, length)
	if os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/induced_cone_human/" + filename):
		# print "\tSkipping (already did)."
		return
	f_is = []
	for residue in range(p):
		f_is.append(get_induct_matrix_to(p, n, core, length, residue))
	if all(f_i is None for f_i in f_is):
		# print "\tSkipping (block upstairs is boooooooorring)."
		return
	induced_dec_matrices = []
	for f_i in f_is:
		if f_i is not None:
			r = get_blocked_dec_mat(p, n-1, f_i.core_downstairs, length)
			induced_dec_matrices.append(matrix(f_i.induct)*matrix(r.decmat))
		else:
			induced_dec_matrices.append(None)
	try:
		r = get_intersected_cones(p, n, core, length)
	except:
		# Might not have obtained intersection of cones.
		r = get_blocked_raw_cones(p, n, core, length)
		r.preg_parts = [part for part in r.parts if is_p_regular(p, part)]
	P_p = matrix(r.P_p)
	induced_adjustment_matrices = [matrix_to_list(P_p.solve_right(ind)) if ind is not None else None for ind in induced_dec_matrices]
	# detailed f_i breakdown in machine readable format
	# print "\tf_i_machine"
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/f_i_machine/" + filename, "w")
	f.write(repr(rec(induced_adjustment_matrices=induced_adjustment_matrices, P_p = r.P_p, parts = r.preg_parts)))
	f.close()
	# detailed f_i breakdown in machine readable format
	# print "\tf_i_human"
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/f_i_human/" + filename, "w")
	for i in range(p):
		if induced_adjustment_matrices[i] is not None:
			f.write("f_%s in adjustment matrix format:\n" %i)
			f.write(get_nice_string_rep_of_matrix([ ( (r.preg_parts[j],) + induced_adjustment_matrices[i][j])
			 for j in range(len(r.preg_parts))]))
			f.write("\n\n")
	f.close()
	# print "\tGetting cone."
	cone_generators = []
	for ind in induced_adjustment_matrices:
		if ind is not None:
			cone_generators.extend(matrix_to_list(matrix(ind).transpose()))
	cone_generators = [row for row in cone_generators if not all(i == 0 for i in row)]
	poly = Polyhedron(rays = cone_generators)
	cone = matrix_to_list(
	 matrix(list(reversed(sorted(tuple(row.vector()) for row in poly.Vrepresentation() if row.is_ray())))).transpose())
	# print "\t\tFor machine."
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/induced_cone_machine/" + filename, "w")
	f.write(repr(rec(induced_cone = cone, P_p = r.P_p, preg_parts = r.preg_parts, parts = r.parts)))
	f.close()
	# print "\t\tFor hooman."
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/induced_cone_human/" + filename, "w")
	f.write(get_nice_string_rep_of_matrix([ ( (r.preg_parts[j],) + cone[j])
	 for j in range(len(r.preg_parts))]))
	f.write("\n")
	f.close()
def induce_all_dec_mats():
	for filename in sorted(os.listdir("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/"),
	 key = get_pncorelength_from_filename):
		# print(filename)
		induce_true_dec_mats_to(filename)
		# print("Garbage collecting...")
		gc.collect()

def generate_contains_function(p_part):
	def f(row_vector):
		return all(sum(dual[i]*row_vector[i] for i in range(len(row_vector))) >= 0 for dual in p_part)
	return f
def generate_poss_adj(p, n=None, core=None, length=None):
	# print "generate_poss_adj(%s, %s, %s, %s)" %(p, n, core, length)
	if n is None:
		(p, n, core, length) = get_pncorelength_from_filename(p)
	filename = get_top_filename(p, n, core, length)
	if not os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/induced_cone_machine/" + filename):
		# print "\tSkipping (don't have induced cone)."
		return
	if not os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/dual_raw_cones/" + filename):
		# print "\tSkipping (don't have dual intersection of cones)."
		return
	if os.path.isfile("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/poss_adj_human/" + filename):
		# print "\tSkipping (already did)."
		return
	#if p <= 5: # TODO this is just to skip some hard cases.
	#	print "\tSkipping."
	#	return
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/induced_cone_machine/" + filename, "r")
	s = f.read()
	f.close()
	r_induced_cone = eval(s)
	r_dualized_intersection_cones = get_dualized_raw_cones(p, n, core, length)
	induced_cone = r_induced_cone.induced_cone
	dualized_cones = r_dualized_intersection_cones.p_part
	# assert r_induced_cone.preg_parts == r_intersection_cones.preg_parts
	parts = r_induced_cone.preg_parts
	assert len(parts)>0
	assert len(parts) == len(induced_cone) == len(dualized_cones[0])
	r_induced_cone = None
	r_dualized_intersection_cones = None
	s = None
	f = None
	gc.collect()
	induced_transposed = [tuple(induced_cone[j][i] for j in range(len(induced_cone))) for i in range(len(induced_cone[0]))]
	##old method
	##intersect_transposed = [tuple(intersection_cones[j][i] for j in range(len(intersection_cones)))
	## for i in range(len(intersection_cones[0]))]
	# I don't think I actually use this one.
	# induced_poly = Polyhedron(rays = induced_transposed)
	intersect_poly = rec(contains = generate_contains_function(dualized_cones))
	# This part is being changed.
	# intersect_poly = Polyhedron(rays = intersect_transposed)
	# first nonzero entries of each column of the induced cone.
	# Also note that the "columns" here became rows because I transposed, but I'm still going to call them columns.
	first_nonzero_entries = tuple(min(i for i in range(len(column)) if column[i]) for column in induced_transposed)
	# The induced cone gives direct upper bounds on the entries of the adjustment matrix.
	# I'll name this "by column" to mean I'm not yet taking a min over all duplicate columns.
	implied_upper_bounds_by_column = [
	 tuple(induced_transposed[i][j]//induced_transposed[i][first_nonzero_entries[i]] for j in range(len(induced_transposed[i])))
	 for i in range(len(induced_transposed))]
	square_dim = len(parts)
	upper_bounds = [
	 tuple(min(implied_upper_bounds_by_column[k][j] for k in range(len(induced_transposed)) if first_nonzero_entries[k] == i)
	 for j in range(square_dim))
	 for i in range(square_dim)]
	# BAM. Now we have our upper bounds, and we compute lower bounds entry by entry using intersect_poly.
	# # Old method: this was horribly wrong - the lower bounds generated are way too strong.
	# # This is a list of lists instead of a list of tuples because I need to mutate it.
	# lower_bounds = [ [ 1 if i == j else 0
	#  for j in range(square_dim)]
	#  for i in range(square_dim)]
	# for i in range(square_dim):
	# 	for j in range(i+1, square_dim):
	# 		column = list(upper_bounds[i])
	# 		# Now keep trying to decrease it until I fall outside the intersection of cones.
	# 		while intersect_poly.contains(column):
	# 			column[j] -= 1
	# 		lower_bounds[i][j] = column[j] + 1
	# New method: Use entries of either
	# (1) the rays of the intersection of cones, or
	# (2) the hecke algebra at p^2, p^3, etc. roots of unity (this will be harder).
	try:
		f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/intersection_raw_cones/" + filename, "r")
		s = f.read()
		f.close()
		r_intersection_cones = eval(s)
		s = None
	except:
		#print "\tSkipping; I'm not ready to try the case when I don't have the intersection of cones yet."
		#raw_input()
		#return
		# print "\tUsing default lower bounds given by Hecke at a pth root of unity."
		lower_bounds = [ [ 1 if i == j else 0
		 for j in range(square_dim)]
		 for i in range(square_dim)]
	else:
		intersection_cones = r_intersection_cones.intersection_adjustment_matrix
		intersect_transposed = [tuple(intersection_cones[j][i] for j in range(len(intersection_cones)))
		 for i in range(len(intersection_cones[0]))]
		first_nonzero_entries = tuple(min(i for i in range(len(column)) if column[i]) for column in intersect_transposed)
		implied_lower_bounds_by_column = [
		 tuple(intersect_transposed[i][j]//intersect_transposed[i][first_nonzero_entries[i]]
		 for j in range(len(intersect_transposed[i])))
		 for i in range(len(intersect_transposed))]
		square_dim = len(parts)
		lower_bounds = [
		 tuple(min(implied_lower_bounds_by_column[k][j] for k in range(len(intersect_transposed)) if first_nonzero_entries[k] == i)
		 for j in range(square_dim))
		 for i in range(square_dim)]
	try:
		# now I want to enumerate the implied possibilities for each column. This involves taking a Cartesian product of ranges.
		import itertools
		# This list will contain a tuple of all possible tuples in each position.
		possible_columns = []
		for i in range(square_dim):
			possible_columns.append(
			 tuple(column for column in itertools.product(*(range(lower_bounds[i][j], upper_bounds[i][j]+1)
			 for j in range(square_dim)))
			 if intersect_poly.contains(column)))
		# Now we take a Cartesian product over all the possible columns.
		possible_adjustment_matrices = []
		for adj_mat in itertools.product(*possible_columns):
			poly = Polyhedron(rays = adj_mat)
			# Now I need to check that the induced cone is generated by this one.
			for column in induced_transposed:
				if not poly.contains(column):
					break
			else:
				possible_adjustment_matrices.append( tuple(
				 tuple(adj_mat[j][i] for j in range(square_dim))
				 for i in range(square_dim)))
	except:
		# print "\tError occurred in processing."
		return
	# Whew! Finally, I need to record this in poss_adj_machine and poss_adj_human.
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/poss_adj_machine/" + filename, "w")
	f.write(repr(rec(possible_adjustment_matrices = possible_adjustment_matrices, parts = parts)))
	f.close()
	f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/poss_adj_human/" + filename, "w")
	for adj in possible_adjustment_matrices:
		f.write(get_nice_string_rep_of_matrix([ ( (parts[j],) + adj[j])
		 for j in range(square_dim)]))
		f.write("\n\n")
	f.close()
	# Now pray all that works.
	# Holy hell. All I had to do was fix a reference to first_nonzero_entry.
def generate_all_poss_adj():
	for filename in sorted(os.listdir("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/induced_cone_machine/"),
	 key = get_pncorelength_from_filename):
		# print(filename)
		generate_poss_adj(filename)
		# print("Garbage collecting...")
		gc.collect()

def clean_up_files():
	for filename in sorted(os.listdir("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/"),
	 key = get_pncorelength_from_filename):
		f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/" + filename, "r")
		s = f.read()
		f.close()
		r = eval(s)
		s = None
		if len(r.P_p) == 0 or len(r.P_p[0]) <= 1:
			# Now cleanup.
			for directory in ("blocked_decomposition_matrices", "blocked_raw_cones", "displayed_intersection_raw_cones",
			 "dual_raw_cones", "f_i_human", "f_i_machine", "induced_cone_human", "induced_cone_machine", "intersection_raw_cones",
			 "poss_adj_human", "poss_adj_machine"):
				toats = "../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/%s/%s" % (directory, filename)
				if os.path.isfile(toats):
					os.remove(toats)

def are_files_clean():
	count0 = 0
	count1 = 0
	countplus = 0
	for filename in sorted(os.listdir("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/"),
	 key = get_pncorelength_from_filename):
		f = open("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/blocked_raw_cones/" + filename, "r")
		s = f.read()
		f.close()
		r = eval(s)
		s = None
		if len(r.P_p) == 0:
			count0 += 1
		elif len(r.P_p[0]) <= 1:
			count1 += 1
		else:
			countplus += 1
