import numpy as np
import math

def normalizacija(x, y):

	# Prvi korak je trazenje tezista tacaka
	sumaX = [0,0,0]
	sumaY = [0,0,0]
	for i in range(8):
		sumaX = [sumaX[j]+x[i][j] for j in range(3)]
		sumaY = [sumaY[j]+y[i][j] for j in range(3)]

	sumaX = [sumaX[i]/8 for i in range(3)]
	sumaY = [sumaY[i]/8 for i in range(3)]

	# Sada sve tacke transliramo oko koordinatnog pocetka,
	# tako da je njihovo teziste upravo u koordinatnom pocetku.
	xn = x[:]
	yn = y[:]

	for i in range(8):
		xn[i] = [xn[i][j]-sumaX[j] for j in range(3)]
		yn[i] = [yn[i][j]-sumaY[j] for j in range(3)]

	# Za svaku tacku racunamo udaljenost od koordinatnog pocetka
	# i trazimo prosek tih rasrojanja, zelimo da namestimo da bude koren iz 2 na kraju

	rastojanjaX = [math.sqrt(xn[i][0]*xn[i][0] + xn[i][1]*xn[i][1]) for i in range(8)]
	rastojanjaY = [math.sqrt(yn[i][0]*yn[i][0] + yn[i][1]*yn[i][1]) for i in range(8)]

	# formula je koeficijent = sqrt(2) / (sum(rastojanja) / 8)
	koeficijentX = math.sqrt(2) * 8 / np.sum(rastojanjaX)
	koeficijentY = math.sqrt(2) * 8 / np.sum(rastojanjaY)

	for i in range(8):
		xn[i] = [xn[i][j] * koeficijentX for j in range(3)]
		yn[i] = [yn[i][j] * koeficijentY for j in range(3)]
		xn[i][2] = 1
		yn[i][2] = 1

	return xn, yn


def fundamentalna_matrica(x, y):

	# Matrica koeficijenata uz nepoznate(f00, f01,..)
	mat = np.array([ [x[i][0]*y[i][0] for i in range(8)],  
				   [x[i][1]*y[i][0] for i in range(8)],
				   [x[i][2]*y[i][0] for i in range(8)],
				   [x[i][0]*y[i][1] for i in range(8)],
				   [x[i][1]*y[i][1] for i in range(8)],
				   [x[i][2]*y[i][1] for i in range(8)],
				   [x[i][0]*y[i][2] for i in range(8)],
				   [x[i][1]*y[i][2] for i in range(8)],
				   [x[i][2]*y[i][2] for i in range(8)]  ]
		)
	matT = mat.transpose()

	# Ovaj nacin ne radi
	#b = np.array([0 for i in range(8)])
	#F = np.linalg.solve(matT, b)

	#Odredjujemo fundamentalnu matricu preko svd-a
	U, D, V = np.linalg.svd(matT)

	# Poslednja kolona matrice V predstavlja nase resenje
	rez = V[:][8]
	F = np.array([	[rez[i]   for i in range(3)],
					[rez[i+3] for i in range(3)],
					[rez[i+6] for i in range(3)]  ]
		)
	print("\n\n")
	print("Determinanta pre: ")
	print(np.linalg.det(F))

	U1, D1, V1 = np.linalg.svd(F)

	D2 = np.array([ [D1[0], 	0, 0],
					[	 0, D1[1], 0],
					[	 0,		0, 0]   ]
		)

	F1 = np.matmul(U1, np.matmul(D2, V1))
	
	print("Determinanta pposle: ")
	print(np.linalg.det(F1))
	return F1

def epipolovi(F):
	U, D, V = np.linalg.svd(F)

	# poslednja kolona matrice U je e2,
	# dok je poslednja kolona matrice V e1
	# Razlika je samo u tome sto je V matrica transponovana

	e1 = [V[2][i] / V[2][2] for i in range(3)]
	e2 = [U[i][2] / U[2][2] for i in range(3)]

	return e1,e2


def kamere(F, e2):
	
	# ~Prva kamera: T1 = [E|0]~
	T1 = np.array([ [1, 0, 0, 0],
					[0, 1, 0, 0],
					[0, 0, 1, 0]   ]
		)

	# ~Druga kamera: T2 = [E2F|e2]~

	# Prvo racunamo matricu mnozenja sa epipolom e2
	E2 = np.array([	[  	  0,-e2[2], e2[1]],
					[ e2[2], 	0, -e2[0]],
					[-e2[1], e2[0],		0]   ]
		)
	E2F = np.matmul(E2, F)

	T2 = np.array([ [E2F[0][0], E2F[0][1], E2F[0][2], e2[0]],
					[E2F[1][0], E2F[1][1], E2F[1][2], e2[1]],
					[E2F[2][0], E2F[2][1], E2F[2][2], e2[2]]   ]
		)

	return T1,T2

def koordinate_u_prostoru(x, y, T1, T2):
	iteracija = len(x)
	X = [[0,0,0,0] for i in range(iteracija)]
	mat = [[0, 0, 0, 0] for i in range(4)]

	# Za svaku tacku racunamo X(a : b : c : d) na osnovu 4 jednacine sa 4 nepoznate
	for i in range(iteracija):
		#pomocna matrica koja sadrzi koeficijente koji stoje uz homogene koordinate tacke X u jednacinama
		#mat = np.array([  [( x[i][1]*T1[2][0] - T1[1][0]), ( x[i][1]*T1[2][1] - T1[1][1]), ( x[i][1]*T1[2][2] - T1[1][2]), ( x[i][1]*T1[2][3] - T1[1][3])],
		#				[( y[i][1]*T2[2][0] - T2[1][0]), ( y[i][1]*T2[2][1] - T2[1][1]), ( y[i][1]*T2[2][2] - T2[1][2]), ( y[i][1]*T2[2][3] - T2[1][3])],
	 	#				[(-x[i][0]*T1[2][0] + T1[0][0]), (-x[i][0]*T1[2][1] + T1[0][1]), (-x[i][0]*T1[2][2] + T1[0][0]), (-x[i][0]*T1[2][3] + T1[0][3])],
		#				[(-y[i][0]*T2[2][0] + T2[0][0]), (-y[i][0]*T2[2][1] + T2[0][1]), (-y[i][0]*T2[2][2] + T2[0][0]), (-y[i][0]*T2[2][3] + T2[0][3])]
		#])

		mat[0] = x[i][1] * T1[2] - T1[1]
		mat[1] = -x[i][0] * T1[2] + T1[0]
		mat[2] = y[i][1] * T2[2] - T2[1]
		mat[3] = -y[i][0] * T2[2] + T2[0]

		# Poslednja kolona matrice V je nase resenje, u ovom slucaju poslednja vrsta, 
		# posto je matrica V transponovana
		U, D, V = np.linalg.svd(mat)

		X[i] = [V[3][i] / V[3][3] for i in range(4)]
		X[i][2] = X[i][2] * 100
		
		print(str(i) + ": \t" + str(X[i]))

	return X

def nedostajuce_tacke(X, Y):
	#Za X nam fale x5, x13 i x16
	#X5
	#Xbesk = 37 X 26
	X37 = np.cross(X[2], X[6])
	X26 = np.cross(X[1], X[5])
	Xb5 = np.cross(X37, X26)

	#Ybesk = 32 X 14
	X32 = np.cross(X[2], X[1])
	X14 = np.cross(X[0], X[3])
	Yb5 = np.cross(X32, X14)

	#X5 = (1 * Xbesk) X (8 * Ybesk)
	m1Xb5 = np.cross(X[0], Xb5)
	m8Yb5 = np.cross(X[7], Yb5)
	X5 = np.cross(m1Xb5, m8Yb5)
	X[4] = [X5[i]/X5[2] for i in range(3)]

	#x13 i x16
	#Ybesk = (11 15) X (10 14)
	X1115 = np.cross(X[10], X[14])
	X1014 = np.cross(X[9], X[13])
	Yb = np.cross(X1115, X1014)

	#Zbesk = (11 12) X (10 9)
	X1112 = np.cross(X[10], X[11])
	X109 = np.cross(X[9], X[8])
	Zb = np.cross(X1112, X109)

	Yb12 = np.cross(Yb, X[11])
	Yb9 = np.cross(Yb, X[8])
	Zb15 = np.cross(Zb, X[14])
	Zb14 = np.cross(Zb, X[13])

	#X16 = (12 Ybesk) X (15 Zbesk	)
	X16 = np.cross(Yb12, Zb15)
	X[15] = [X16[i]/X16[2] for i in range(3)]

	#X13 = (9 Ybesk) X (14 Zbesk)
	X13 = np.cross(Yb9, Zb14)
	X[12] = [X13[i]/X13[2] for i in range(3)]

	#Za Y nam fale y5 i y16
	#y16
	#Ybesk = (11 15) X (10 14)
	Y1115 = np.cross(Y[10], Y[14])
	Y1014 = np.cross(Y[9], Y[13])
	Yb = np.cross(Y1115, Y1014)

	#Zbesk = (11 12) X (10 9)
	Y1112 = np.cross(Y[10], Y[11])
	Y109 = np.cross(Y[9], Y[8])
	Zb = np.cross(Y1112, Y109)

	Yb12 = np.cross(Yb, X[11])
	Zb15 = np.cross(Zb, X[14])

	#Y16 = (12 Ybesk) X (15 Zbesk)
	Y16 = np.cross(Yb12, Zb15)
	Y[15] = [Y16[i]/Y16[2] for i in range(3)]

	#Y5
	#Xbesk = (3 7) X (2 6)
	Y37 = np.cross(Y[2], Y[6])
	Y26 = np.cross(Y[1], Y[5])
	Xb5 = np.cross(Y37, Y26)

	#Ybesk = (3 2) X (1 4)
	Y32 = np.cross(Y[2], Y[1])
	Y14 = np.cross(Y[0], Y[3])
	Yb5 = np.cross(Y32, Y14)

	#Y5 = (1 Xbesk) X (8 Ybesk)
	Xb51 = np.cross(Y[0], Xb5)
	Yb58 = np.cross(Y[7], Yb5)
	Y5 = np.cross(Xb51, Yb58)
	Y[4] = [Y5[i]/Y5[2] for i in range(3)]

	return X, Y


def main():

	#Sve tacke vidljive na slici
	#Tacke koje se ne vide cemo izracunati
	Xsve = [  [812, 98, 1],	#x1
			  [916, 141, 1],	#x2
			  [759, 255, 1],	#x3
			  [655, 212, 1],	#x4
			  [0, 0, 1], 		#x5
			  [903, 422, 1],	#x6
			  [750, 547, 1],	#x7
			  [649, 497, 1],	#x8
			  [372, 443, 1],	#x9 
			  [691, 710, 1],	#x10
			  [689, 1037, 1],	#x11
			  [388, 765, 1],	#x12
			  [0, 0, 1],		#x13
			  [1163, 465, 1],	#x14
			  [1144, 783, 1],	#x15
			  [0, 0, 1]			#x16
	]

	Ysve = [  [903, 144, 1],	#y1
			  [962, 207, 1],	#y2
			  [741, 262, 1],	#y3
			  [680, 201, 1],	#y4
			  [0, 0, 1],		#y5
			  [933, 473, 1],	#y6
			  [718, 541, 1],	#y7
			  [660, 473, 1],	#y8
			  [464, 356, 1],	#y9
			  [571, 653, 1],	#y10
			  [558, 970, 1],	#y11
			  [457, 663, 1],	#y12
			  [999, 298, 1],	#y13
			  [1145, 580, 1],	#y14
			  [1114, 890, 1],	#y15
			  [0, 0, 1]			#y16
	]

	#Za zadatak korisitmo tacke 1, 2, 3, 4, 9, 10, 11, 12
	x = [	[812, 98, 1],	#x1
			[916, 141, 1],	#x2
			[759, 255, 1],	#x3
			[655, 212, 1],	#x4
			[372, 443, 1],	#x9 
			[691, 710, 1],	#x10
			[689, 1037, 1],	#x11
			[388, 765, 1]	#x12
	]

	y = [   [903, 144, 1],	#y1
			[962, 207, 1],	#y2
			[741, 262, 1],	#y3
			[680, 201, 1],	#y4
			[464, 356, 1],	#y9
			[571, 653, 1],	#y10
			[558, 970, 1],	#y11
			[457, 663, 1]	#y12
	]

	xn, yn = normalizacija(x, y)

	print("\n~Normalizovane tacke~\nNormalizovane X koordinate:")
	for i in range(8):
		print(xn[i])
	
	print("Normalizovane Y koordinate:")
	for i in range(8):
		print(yn[i])


	Fn = fundamentalna_matrica(xn, yn)
	F = fundamentalna_matrica(x,  y)
	print("\n~Fundamentalna matrica~")
	print(F)

	e1, e2 = epipolovi(F)
	print("\n~Epipolovi~")
	print(e1)
	print(e2)

	# Matrice kamera
	T1, T2 = kamere(F, e2)
	print("\n~Kamere~\nPrva kamera:")
	print(T1)
	print("Druga kamera:")
	print(T2)

	matrica = [[T2[i][j] for j in range(3)] for i in range(3)]

	#Racunamo Q i R matrice kamere
	Q, R = np.linalg.qr(np.linalg.pinv(matrica))
	K = np.linalg.inv(R)
	print("\nMatrica K:")
	print(K/K[2][2])
	print("------------------\nMatrica A:")
	print(Q.T)
	print("------------------")

	Xsve, Ysve = nedostajuce_tacke(Xsve, Ysve)
	print("\n~Dodate nedostajuce tacke~\nKoordinate svih 16 temena jedne slike:")
	for i in range(len(Xsve)):
		print(Xsve[i])
	print("\nKoordinate svih 16 temena druge slike:")
	for i in range(len(Ysve)):
		print(Ysve[i])

	print("\n~Krajnje tacke~")
	Rezultat = koordinate_u_prostoru(Xsve, Ysve, T1, T2)



if __name__ == '__main__':
	main()