load("pmns_arith_ops.sage")

R.<x> = ZZ[]

#data struct: [p, n, gmm, E, phi_log2, delta]

#pmns = [65311477604033896723481007574120302187674052674795197557479610069165116413291, 5, 30496124510329138387329131429910118297009679881561008645677630224367532768779, (-1, -1, 0, 0, 0, 1), 64, 0]
pmns = [291791, 2, 11810, (-2, 0, 1), 16, 0]  
 
pmns_q = pmns[0] 	 
pmns_n = pmns[1]
pmns_gmm = pmns[2]
pmns_E = R(pmns[3])
pmns_phi_log2 = pmns[4]
pmns_dlt = pmns[5]

########################################################################

Zq = GF(pmns_q)

pmns_convbase_log2 = ceil(log(pmns_q,2)/pmns_n)
pmns_beta = 1 << pmns_convbase_log2

pmns_phi = 1 << pmns_phi_log2

pmns_iM_dom.<t> = ZZ.quo(pmns_phi)[]

pmns_in_mont_domain = True	# if True, non-exact conversions will be put in Montgomery domain; also used for converion from PMNS

pmns_j = 1 # can be updated "freely"

########################################################################

pmns_Mmat = build_lattice_base(pmns_q, pmns_n, pmns_gmm)  #WARNING: 'build_lattice_base' should be the exact same function used for PMNS generation.

pmns_Mmat_norm1 = Integer(pmns_Mmat.norm(1))

pmns_rho = 1 + pmns_Mmat_norm1

pmns_iMmat = compute_neg_inv_ri_mat(pmns_Mmat, pmns_iM_dom)

pmns_u = compute_param_u(pmns_q, pmns_n, pmns_dlt, pmns_rho, pmns_E, pmns_Mmat, R)

pmns_transV = compute_translation_vect(pmns_n, pmns_u, pmns_Mmat)

(pmns_Omega, pmns_Delta) = compute_trans_to_dHs_params(pmns_q, pmns_n, pmns_phi, pmns_Mmat, pmns_iMmat, pmns_iM_dom, Zq, pmns_transV, pmns_E, R)

pmns_conv_polys = compute_conv_polys__in_D1(pmns_q, pmns_n, pmns_convbase_log2, pmns_phi, pmns_Mmat, pmns_iMmat, pmns_iM_dom, Zq, pmns_transV, pmns_E, R, pmns_in_mont_domain)

pmns_zero_reps = compute_zero_reps_in_Dj(pmns_j, pmns_n, pmns_Mmat)

print(pmns_zero_reps)

########################################################################

v = randint(0, pmns_q - 1)

v_rep = exact_conv_to_pmns__to_D1(v, pmns_n, pmns_q, pmns_phi, pmns_Mmat, pmns_iMmat, pmns_iM_dom, Zq, pmns_transV)

v_repH = compute_rep_in_H(v_rep, pmns_n, pmns_Omega, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, R)

v_repHH = compute_rep_in_HH(v_rep, pmns_n, pmns_Delta, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, pmns_phi_log2, R)

print(v_rep,"\n\n",v_repH,"\n\n",v_repHH)


v_reps = compute_rep_set(v, pmns_q, pmns_n, pmns_Omega, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, pmns_zero_reps, Zq, R)

print(v_reps)


########################################################################

a = randint(0, pmns_q - 1)
b = randint(0, pmns_q - 1)

a_rep = conv_to_pmns__to_D1(a, pmns_n, pmns_convbase_log2, pmns_phi, pmns_conv_polys, pmns_Mmat, pmns_iMmat, pmns_iM_dom, Zq, pmns_transV)
b_rep = conv_to_pmns__to_D1(b, pmns_n, pmns_convbase_log2, pmns_phi, pmns_conv_polys, pmns_Mmat, pmns_iMmat, pmns_iM_dom, Zq, pmns_transV)

c_rep = pmns_trans_mont_mult(a_rep, b_rep, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, R)

c = conv_from_pmns(c_rep, pmns_n, pmns_gmm, pmns_phi, Zq, pmns_in_mont_domain) 

print(c_rep)

print(c == ((a*b)%pmns_q)) # caution: 'c' should be multiplied by 'pmns_phi' if 'pmns_in_mont_domain=False', because it is a product and the factor '1/phi' from GMont-like


########################################################################
######################### equality check ###############################

check1 = equality_check(a_rep, b_rep, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_transV)

a_rep2 = a_rep + pmns_zero_reps[randint(0, len(pmns_zero_reps)-1)]

check2 = equality_check(a_rep, a_rep2, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_transV)

print(a, b)
print(a_rep,"\n\n",a_rep2,"\n\n",b_rep)
print(check1,"\n", check2)



########################################################################












