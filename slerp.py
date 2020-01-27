import math as m
import numpy as np

def ispisMatrice(A):
	for i in range(len(A)):
		print("    "+str(A[i]))


#Proveravamo da li je determinanta jednaka 1 i matrica nije jedinicna
def proveraMatrice(A):
	E = [[1, 0, 0],
		 [0, 1, 0],
		 [0, 0, 1]]

	if np.linalg.det(A) != 1 :
		return True

	if (A == E) == True:
		return True

	return False


# Prihvata kao argument tri Ojlerova ugla, vraca matricu
def Euler2A(fi, teta, psi):

	#Odredjujemo sinuse i kosinuse
	sinX = m.sin(fi)
	cosX = m.cos(fi)
	sinY = m.sin(teta)
	cosY = m.cos(teta)
	sinZ = m.sin(psi)
	cosZ = m.cos(psi)

	#odredjujemo matrice
	Rx = np.array([[1, 0, 0],
		  [0, cosX, -sinX],
		  [0, sinX,  cosX]])

	Ry = np.array([[ cosY, 0, sinY],
		  [0, 1, 0],
		  [-sinY, 0, cosY]])

	Rz = np.array([[cosZ, -sinZ, 0],
		  [sinZ,  cosZ, 0],
		  [0, 0, 1]])


	#Mnozenje matrica

	#Rez1 = Rz * Ry	
	Rez1 = [[0,0,0],
			[0,0,0],
			[0,0,0]]
	
	for i in range(len(Rz)):
		for j in range(len(Ry[0])):
			for k in range(len(Ry)):
				Rez1[i][j] += Rz[i][k] * Ry[k][j]


	#Rez2 = (Rz * Ry) * Rx
	Rez2 = [[0,0,0],
			[0,0,0],
			[0,0,0]]

	for i in range(len(Rez1)):
		for j in range(len(Rx[0])):
			for k in range(len(Rx)):
				Rez2[i][j] += Rez1[i][k] * Rx[k][j]
	return Rez2


#Prihvata kao argument matricu, vraca vektor i ugao
def AxisAngle(A):
	E = [[1, 0, 0],
		 [0, 1, 0],
		 [0, 0, 1]]

	if proveraMatrice(A):
		print("Matrica za funkciju AxisAngle ne valja!")
		exit(1)

	Ap = A
	Ap[0][0] -=1
	Ap[1][1] -=1
	Ap[2][2] -=1

	p =np.array(np.cross(Ap[0], Ap[1]))
	pom = m.sqrt(p[0]**2 + p[1]**2 + p[2]**2)
	pN= [i/pom for i in p]

	u = np.array(Ap[1])

	#Posto se gore kopirao pokazivac na matricu, sad moramo vratiti vrednosti, jer nam treba
	#matrica A!
	Ap[0][0] +=1
	Ap[1][1] +=1
	Ap[2][2] +=1

	Up = [0,0,0]
	for i in range(3):
		for j in range(3):
			Up[i] += A[i][j] * u[j]


	#odredjivanje ugla
	brojilac = u[0] * Up[0] + u[1] * Up[1] + u[2] * Up[2]
	imenilac = m.sqrt(u[0]**2 + u[1]**2 + u[2]**2)  *  m.sqrt(Up[0]**2 + Up[1]**2 + Up[2]**2 )
	ugao = m.acos(brojilac/imenilac)

	#na kraju proveravamo mesoviti proizvod
	if 0 > np.linalg.det([u, Up, pN]):
		pN = [-i for i in pN]
	
	return pN, ugao


#Prihvata kao argument vektor i ugao, vraca matricu
def Rodrigez(p, ugao):
	E = [[1, 0, 0],
		 [0, 1, 0],
		 [0, 0, 1]]

	pp = np.array([p[0], p[1], p[2]])

	#inicijalizujemo tri matrice,ciji zbir, pomnozen sa koeficijentima trazimo

	#1
	#ppT = pp.dot(p)
	ppT = [ [pp[0]*p[0], pp[0]*p[1], pp[0]*p[2]],
			[pp[1]*p[0], pp[1]*p[1], pp[1]*p[2]],
			[pp[2]*p[0], pp[2]*p[1], pp[2]*p[2]]]

	#2
	EMinusPPT = [ [1-pp[0]*p[0],  -pp[0]*p[1],  -pp[0]*p[2]],
				  [ -pp[1]*p[0], 1-pp[1]*p[1],  -pp[1]*p[2]],
				  [ -pp[2]*p[0],  -pp[2]*p[1], 1-pp[2]*p[2]]]

	#3
	mat = np.array([[    0, -p[2],  p[1]],
		   			[ p[2],     0, -p[0]],
		   			[-p[1],  p[0],  0]])

	R = [[0,0,0],
		 [0,0,0],
		 [0,0,0]]

	for i in range(3):
		for j in range(3):
			R[i][j] = ppT[i][j] + EMinusPPT[i][j] * m.cos(ugao) + mat[i][j] * m.sin(ugao)

	#R = ppT + EMinusPPT*m.cos(ugao) + m.sin(ugao)*mat
	return R


#Prihvata kao argument matricu, vraca Ojlerove uglove
def A2Euler(A):
	if proveraMatrice(A):
		print("Matrica za funkciju A2Euler ne valja!")
		exit(1)

	if A[2][0] < 1:
		if (-1 < A[2][0]):		#jedinstveno resenje
			psi  = m.atan2(A[1][0], A[0][0])
			teta = m.asin(-A[2][0])
			fi   = m.atan2(A[2][1], A[2][2])
		else:					#nije jedinstveno, slucaj Ox3 = -Oz
			psi  = m.atan2(-A[0][1], A[1][1])
			teta = m.pi / 2
			fi   = 0
	else:						#nije jedinstveno, slucaj Ox3 = Oz
		psi  = m.atan2(-A[0][1], A[1][1])
		teta = (-m.pi) / 2
		fi   = 0

	return fi, teta, psi


#Prihvata kao argument vektor i ugao, vraca kvaternion
def AxisAngle2Q(p, ugao):
	if ugao == 0:
		return [0, 0, 0, 1]

	sinus = m.sin(ugao / 2)

	return [p[0]*sinus, p[1]*sinus, p[2]*sinus, m.cos(ugao / 2)]


#Prihvata kao argument kvaternion, vraca vektor i ugao
def Q2AxisAngle(q):
	pom = m.sqrt(q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2)
	qP = [i/pom for i in q]
	if 0 > qP[3]:
		qP = [-i for i in qP]

	ugao = 2 * m.acos(qP[3])
	if qP[3] == 1 or qP[3] == -1:
		return [1, 0, 0], ugao
	else:
		pom1 = m.sqrt(qP[0]**2 + qP[1]**2 + qP[2]**2)
		p = [qP[0]/pom1, qP[1]/pom1, qP[2]/pom1]
		return p, ugao


#DRUGI DEO DOMACEG ZADATKA

def lerp(q1, q2, tm, t):
	q1 = [(i * (1-t/tm)) for i in q1]
	q2 = [(i * (t/tm)) 	 for i in q2]

	Rez = [q1[i] + q2[i] for i in range(4)]

	pom = m.sqrt(sum(i ** 2 for i in Rez))
	return [i / pom for i in Rez]

def slerp(q1, q2, tm, t):
	if tm < t and t < 0:
		print("vrednost t mora da bude izmedju 0 i tm!")
		exit(1)

	cosP = sum([q1[i] * q2[i] for i in range(4)])

	if cosP < 0:		#idi po kracem luku sfere
		q1 = [-i for i in q1]
		cosP = -cosP

	if cosP > 0.95:		#kvaternioni q1 i q2 previse blizu
		return lerp(q1, q2, tm, t)

	fiP = m.acos(cosP)
	alfa = m.sin(fiP * (1 - t/tm) / m.sin(fiP))
	beta = m.sin(fiP * t/tm) / m.sin(fiP)
	
	return [q1[i]*alfa + q2[i]*beta for i in range(4)]




def main():

	print("Uglovi:")

	fi   = m.radians(25)
	teta = m.radians(50)
	psi  = m.radians(75)
	print("  " + str(fi) + ", " +str(teta) + ", " + str(psi))

	print("\n1) Euler2A:")
	A = Euler2A(fi, teta, psi)
	print("  A: ")
	ispisMatrice(A)
	#------------------------------------------------------

	print("\n2) AxisAngle:")
	p, ugao = AxisAngle(Euler2A(fi, teta, psi))
	print("  P:  " + str(p) + "\n  fi: " + str(ugao))
	#------------------------------------------------------

	print("\n3) Rodrigez:")
	Ap = Rodrigez(p, ugao)
	print("  A: ")
	ispisMatrice(Ap)
	#------------------------------------------------------

	print("\n4) A2Euler:")
	fiP, tetaP, psiP = A2Euler(A)
	print("  Uglovi: " + str(fiP) + ", " +str(tetaP) + ", " + str(psiP))
	#------------------------------------------------------

	print("\n5) AxisAngle2Q:")
	q = AxisAngle2Q(p, ugao)
	print("  q: " + str(q))

	#------------------------------------------------------

	print("\n6) Q2AxisAngle:")
	pQ, ugaoQ = Q2AxisAngle(q)
	print("  p:    " + str(pQ) + "\n  ugao: " + str(ugaoQ))




if __name__ == '__main__':
	main()