load("pmns_generator.sage")

roots_max_duration_checks = 60  # seconds

bSize = 256
p = random_prime(2**bSize, lbound=2**(bSize-1)) # or p = Integer('a_value')

phi_log2 = 64 
lambda_max = 3

delta = 0

n = floor(p.nbits()/phi_log2) + 1

pmns_list = build_pmns_candidates_for_n(p, n, phi_log2, delta, lambda_max, roots_max_duration_checks)

if len(pmns_list) != 0 : print(pmns_list[0])


for pmns in pmns_list:
	print(pmns)
	print()
	
#~ NOTE: data struct: [p, n, gmm, E, phi_log2, delta]

########################################################################

# ~ An example: 256-bit prime

[65130057217352493063258531782528987527216027708801919615305112039069287150861, 5, 
40031617528150199516543576638109793972900977382860490744805902362180553667328, 
(1, -1, 0, 0, 0, 1), 64, 1]

########################################################################

