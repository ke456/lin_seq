# generates a polynomial P = p_1^m p_2^(m+1) p_3^(m+2) ... p_t^(m+t-1)
# where deg(p_i) = d over field
def gen_test(d,m,t,p):
	M.<x> = PolynomialRing(GF(p))
	P = 1
	for i in range(t):
		rand = M.random_element(d)^(m+i)
		P = P*rand
	return P
	
def prepare_test(d,m,t,p):
	P = gen_test(d,m,t,p)
	P = P.monic()
	print(P.squarefree_decomposition())
	init = []
	for i in range(P.degree()):
		init.append(1)
	name = "t_"+str(d)+"_"+str(m)+"_"+str(t)
	f = open(name, "w");
	f.write(str(p)+"\n")
	for i in P.coefficients(sparse=false):
		f.write(str(i) + " ")
	f.write("\n")
	for i in init:
		f.write(str(i) + " ")
	f.write("\n")
	seq = [P.degree()*100,P.degree()*200,P.degree()*300]
	for i in seq:
		f.write(str(i) + "\n")
	return P
	
