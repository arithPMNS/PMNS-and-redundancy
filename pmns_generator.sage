from multiprocessing import Process, Queue

########################################################################

def build_lattice_base(p, n, gmm):
	b = []
	l = [p] + [0]*(n-1)
	b.append(l)
	m = identity_matrix(n)
	for i in range(1,n):
		t = (-gmm.powermod(i, p))%p
		if t%2 == 1 :
			t += p
		l = [t] + list(m[i][1:])
		b.append(l)
	bb = matrix(b)
	return bb.LLL(delta=1, algorithm='NTL:LLL')

########################################################################

def roots_computation(ext_pol, queue):
	res = []
	try:
		res = ext_pol.roots(multiplicities=False)
	except ValueError:  # ValueError, msg:
		#~ print msg
		res = []
	finally:
		queue.put(res)

#~ tries to find a root of 'ext_pol' within 'roots_max_duration_checks' seconds
def find_roots_with_timeout(ext_pol, roots_max_duration_checks):
	
	res = []
	queue = Queue() # used to get the result
	proc = Process(target=roots_computation, args=(ext_pol, queue)) # creation of a process calling longfunction with the specified arguments
	proc.start() # lauching the processus on another thread
	try:
		res = queue.get(timeout=roots_max_duration_checks) # getting the resultat under 'max_duration' seconds or stop
		proc.join() # proper delete if the computation has take less than timeout seconds
	except Exception:  #Exception, msg:
		proc.terminate() # kill the process
		#~ print ("Timed out!")
	return res

########################################################################
	
def extPolys_sorting_criteria(ext_pol_coeffs):
	return compute_redExtPol_w(len(ext_pol_coeffs) - 1, ext_pol_coeffs)


def gen_extPol_set(n, lambda_max):
	
	extPol_set = []
	
	# ~ extPol_set.append([-1]+[0]*(n-1)+[1])  # very unlikely to be useful
	if Integer(n).is_power_of(2):
		extPol_set.append(tuple([1]+[0]*(n-1)+[1]))
		
	extPol_set.append(tuple([-1, -1]+[0]*(n-2)+[1]))
	extPol_set.append(tuple([1, -1]+[0]*(n-2)+[1]))
	extPol_set.append(tuple([-1, 1]+[0]*(n-2)+[1]))
	extPol_set.append(tuple([1, 1]+[0]*(n-2)+[1]))
	extPol_set.append(tuple([1]*(n+1)))
	
	tpol=[1]
	coeff = -1
	for i in range(n):
		tpol = [coeff] + tpol
		coeff *= -1
	extPol_set.append(tuple(tpol))
	
	if n%2 == 0 :
		
		zero = [0]*((n//2)-1)
		extPol_set.append(tuple([1] + zero + [1] + zero + [1]))
		extPol_set.append(tuple([1] + zero + [-1] + zero + [1]))
		
		tpol = [1]
		for i in range(0,n,2):
			tpol = [1,0] + tpol
		extPol_set.append(tuple(tpol))
		
		tpol = [1]
		coeff = -1
		for i in range(0,n,2):
			tpol = [coeff,0] + tpol
			coeff *= -1
		extPol_set.append(tuple(tpol))
		
	for lmbd in range(2, lambda_max+1):
		extPol_set.append(tuple([-lmbd]+[0]*(n-1)+[1]))
		extPol_set.append(tuple([lmbd]+[0]*(n-1)+[1]))
	
	extPol_set = list(set(extPol_set))
	extPol_set.sort(key=extPolys_sorting_criteria)
	# ~ print(extPol_set)
	return extPol_set

def compute_redExtPol_w(n, ext_polC):
	R.<x> = ZZ[]
	x = R.gen()
	ext_pol = R(list(ext_polC))
	c = n-1
	tmp = x**n
	l1 = list(tmp%ext_pol)
	l1 = [abs(k) for k in l1]
	V = c * vector(l1 + [0]*(n-len(l1)))
	for d in range(n-2):
		c -= 1
		tmp *= x
		l1 = list(tmp%ext_pol)
		l1 = [abs(k) for k in l1]
		V += c * vector(l1 + [0]*(n-len(l1)))
	V += vector(range(1, n+1))
	return max(V)

########################################################################

def compute_param_u(p, n, delta, rho, ext_polC, wE, G):
	
	k = ceil(log(p,2)/n)
	bt = 1 << k
	c1 = n * (bt - 1)
	
	f_add = (delta + 1)**2 
	c2 = wE * f_add * (rho - 1)
	
	m = max(c1, c2)
	iG = G.inverse()
	u = ceil(m * (rho-1) * iG.norm(1))
	
	return u


def check_gamma(gmm, p, n, delta, ext_polC, wE, phi_log2):
	
	G = build_lattice_base(p, n, gmm)
	
	rho = Integer(G.norm(1)) + 1
	
	u = compute_param_u(p, n, delta, rho, ext_polC, wE, G)
	
	phi_log2min = ceil(log(2*u, 2))
	
	if phi_log2min > phi_log2:
		return -1

	pmns_found = [p, n, gmm, ext_polC, phi_log2, delta] # G is re-computed later
	
	return pmns_found


########################################################################

#~ checks if 'extPol' produces some PMNS and returns them if so. 
def quick_check_ext_pol(ZX, F, p, n, redExtPol_coeffs, roots_max_duration_checks, nb_free_add, phi_log2):

	redExtPol_coeffs = list(redExtPol_coeffs)
	
	T.<y> = F[]
	redExtPolT = T(redExtPol_coeffs)
	
	wE = compute_redExtPol_w(n, redExtPol_coeffs)
	
	print("Starting: " + str(ZX(redExtPol_coeffs)))
	print("E = " + str(redExtPol_coeffs))
	print("w = " + str(wE))
	
	gmms = find_roots_with_timeout(redExtPolT, roots_max_duration_checks)
	
	if gmms == []:
		print("Done: NO ROOT FOUND!\n")
		return []
	
	pmns_listt = []
	for gmm in gmms :
		
		gmm = Integer(gmm)
		
		if (gmm == 0) or (gmm == 1) or (gmm == (p-1)):
			continue # not useful roots
		
		pmns = check_gamma(gmm, p, n, nb_free_add, redExtPol_coeffs, wE, phi_log2)
		if pmns != -1 :
			pmns_listt.append(pmns)
			
	print("Done: " + str(len(pmns_listt)) + " PMNS found.\n")
	return pmns_listt


def build_pmns_candidates_for_n(p, n, phi_log2, nb_free_add, lambda_max, roots_max_duration_checks):
	
	F = GF(p)
	ZX.<X> = ZZ[]
	
	extPol_set = gen_extPol_set(n, lambda_max)
	
	print("STARTING\nCandidate external polynomials generation done, nbre cdts: " + str(len(extPol_set)) + "\n") 
	
	pmns_cands = []	
	for extPolC in extPol_set:
		rep = quick_check_ext_pol(ZX, F, p, n, extPolC, roots_max_duration_checks, nb_free_add, phi_log2)
		if rep != [] :
			pmns_cands.append(rep)
	
	pmns_cands = flatten(pmns_cands, max_level=1)
	
	print("DONE. Total number of PMNS found: " + str(len(pmns_cands)))
		
	return pmns_cands


########################################################################











