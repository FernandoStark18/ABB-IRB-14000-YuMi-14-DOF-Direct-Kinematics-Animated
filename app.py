# Espericueta Fonseca Fernando Simón
# Mecatrónica 6to 6
# Cinemática de robots
# Análisis cinemático por la convención Denavit-Hartenberg del robot ABB IRB 14000 YuMi con 14 GDL

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm

fig, ax = plt.subplots()
ax = plt.axes(projection = "3d")
#Función para rotar la función en el eje x
def matriz_rotacion_y(grados):
	rad = grados/180*np.pi
	#Matriz de rotación
	rotacion=np.array([[np.cos(rad),0,-np.sin(rad),0],
					   [0,1,0,0],
					   [np.sin(rad),0,np.cos(rad),0],
					   [0,0,0,1]])
	return rotacion
#Función para rotar la función en el eje y
def matriz_rotacion_x(grados):
	rad = grados/180*np.pi
	#Matriz de rotacuón
	rotacion=np.array([[1,0,0,0],
					   [0,np.cos(rad),-np.sin(rad),0],
					   [0,np.sin(rad),np.cos(rad),0],
					   [0,0,0,1]])
	return rotacion
def matriz_rotacion_z(grados):
	rad = grados/180*np.pi
	#Matriz de rotacuón
	rotacion=np.array([[np.cos(rad),-np.sin(rad),0,0],
					   [np.sin(rad),np.cos(rad),0,0],
					   [0,0,1,0],
					   [0,0,0,1]])
	return rotacion

def matriz_traslacion_x(x):
	traslacion = np.array([[1,0,0,x],
						   [0,1,0,0],
						   [0,0,1,0],	   
						   [0,0,0,1]])
	return traslacion

def matriz_traslacion_y(y):
	traslacion = np.array([[1,0,0,0],
						   [0,1,0,y],
						   [0,0,1,0],	   
						   [0,0,0,1]])
	return traslacion

def matriz_traslacion_z(z):
	traslacion = np.array([[1,0,0,0],
						   [0,1,0,0],
						   [0,0,1,z],	   
						   [0,0,0,1]])
	return traslacion

def configuracion_grafica():
	plt.title("ABB IRB 14000 YuMi 14 GDL")
	ax.set_xlim(-70,60)
	ax.set_ylim(-70,60)
	ax.set_zlim(-70,60)
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")
	ax.view_init(elev=25,azim=-90)

def sistema_coordenadas(a,b,c,a_1,b_1,c_1):
	x = [a,a_1]
	y = [b,b_1]
	z = [c,c_1]
	ax.plot3D(x,[b,b],[c,c],color='red')
	ax.plot3D([a,a],y,[c,c],color='blue')
	ax.plot3D([b,a],[b,b],z,color='green')

def sistema_coordenadas_movil(matriz_rotacion):
	r_11 = matriz_rotacion[0,0]
	r_12 = matriz_rotacion[1,0]
	r_13 = matriz_rotacion[2,0]
	r_21 = matriz_rotacion[0,1]
	r_22 = matriz_rotacion[1,1]
	r_23 = matriz_rotacion[2,1]
	r_31 = matriz_rotacion[0,2]
	r_32 = matriz_rotacion[1,2]
	r_33 = matriz_rotacion[2,2]

	dx = matriz_rotacion[0,3]
	dy = matriz_rotacion[1,3]
	dz = matriz_rotacion[2,3]


	ax.plot3D([dx,dx+r_11],[dy,dy+r_22],[dz,dz+r_13], color='m')
	ax.plot3D([dx,dx+r_21],[dy,dy+r_22],[dz,dz+r_23], color='c')
	ax.plot3D([dx,dx+r_31],[dy,dy+r_32],[dz,dz+r_33], color='y')
	plt.draw()

def DH(theta_i, di, ai, alpha_i):
	MT = matriz_rotacion_z(theta_i)@matriz_traslacion_z(di)@matriz_traslacion_x(ai)@matriz_rotacion_x(alpha_i)
	#MT = matriz_rotacion_x(theta_i)@matriz_traslacion_x(di)@matriz_traslacion_y(ai)@matriz_rotacion_y(alpha_i)
	return MT

def Robot_RR(theta_1, d1, a1, alpha_1, theta_2, d2, a2, alpha_2, theta_3, d3, a3, alpha_3, theta_4, d4, a4, alpha_4, theta_5, d5, a5, alpha_5, theta_6, d6, a6, alpha_6, theta_7, d7, a7, alpha_7,theta_12, d12, a12, alpha_12, theta_22, d22, a22, alpha_22, theta_32, d32, a32, alpha_32, theta_42, d42, a42, alpha_42, theta_52, d52, a52, alpha_52, theta_62, d62, a62, alpha_62, theta_72, d72, a72, alpha_72):
	A0 = np.eye(4)
	_0A1 = DH(theta_1, d1, a1, alpha_1)
	_1A2 = DH(theta_2, d2, a2, alpha_2)
	_2A3 = DH(theta_3, d3, a3, alpha_3)
	_3A4 = DH(theta_4, d4, a4, alpha_4)
	_4A5 = DH(theta_5, d5, a5, alpha_5)
	_5A6 = DH(theta_6, d6, a6, alpha_6)
	_6A7 = DH(theta_7, d7, a7, alpha_7)
	_0A1 = _0A1@A0
	_0A2 = _0A1@_1A2
	_0A3 = _0A1@_1A2@_2A3
	_0A4 = _0A1@_1A2@_2A3@_3A4
	_0A5 = _0A1@_1A2@_2A3@_3A4@_4A5
	_0A6 = _0A1@_1A2@_2A3@_3A4@_4A5@_5A6
	_0A7 = _0A1@_1A2@_2A3@_3A4@_4A5@_5A6@_6A7

	sistema_coordenadas_movil(A0)
	sistema_coordenadas_movil(_0A1)
	sistema_coordenadas_movil(_0A2)
	sistema_coordenadas_movil(_0A3)
	sistema_coordenadas_movil(_0A4)
	sistema_coordenadas_movil(_0A5)
	sistema_coordenadas_movil(_0A6)
	sistema_coordenadas_movil(_0A7)

	A02 = np.array([[1,0,0,-15],
					[0,1,0,0],
			    	[0,0,1,0],	   
				    [0,0,0,1]])
	_0A12 = DH(theta_12, d12, a12, alpha_12)
	_1A22 = DH(theta_22, d22, a22, alpha_22)
	_2A32 = DH(theta_32, d32, a32, alpha_32)
	_3A42 = DH(theta_42, d42, a42, alpha_42)
	_4A52 = DH(theta_52, d52, a52, alpha_52)
	_5A62 = DH(theta_62, d62, a62, alpha_62)
	_6A72 = DH(theta_72, d72, a72, alpha_72)
	_0A12 = A02@_0A12
	_0A22 = _0A12@_1A22
	_0A32 = _0A12@_1A22@_2A32
	_0A42 = _0A12@_1A22@_2A32@_3A42
	_0A52 = _0A12@_1A22@_2A32@_3A42@_4A52
	_0A62 = _0A12@_1A22@_2A32@_3A42@_4A52@_5A62
	_0A72 = _0A12@_1A22@_2A32@_3A42@_4A52@_5A62@_6A72

	sistema_coordenadas_movil(A02)
	sistema_coordenadas_movil(_0A12)
	sistema_coordenadas_movil(_0A22)
	sistema_coordenadas_movil(_0A32)
	sistema_coordenadas_movil(_0A42)
	sistema_coordenadas_movil(_0A52)
	sistema_coordenadas_movil(_0A62)
	sistema_coordenadas_movil(_0A72)

	ax.plot3D([A0[0,3],_0A1[0,3]],[A0[1,3],_0A1[1,3]],[A0[2,3],_0A1[2,3]], color = 'red')
	ax.plot3D([_0A1[0,3],_0A2[0,3]],[_0A1[1,3],_0A2[1,3]],[_0A1[2,3],_0A2[2,3]], color = 'green')
	ax.plot3D([_0A2[0,3],_0A3[0,3]],[_0A2[1,3],_0A3[1,3]],[_0A2[2,3],_0A3[2,3]], color = 'red')
	ax.plot3D([_0A3[0,3],_0A4[0,3]],[_0A3[1,3],_0A4[1,3]],[_0A3[2,3],_0A4[2,3]], color = 'green')
	ax.plot3D([_0A4[0,3],_0A5[0,3]],[_0A4[1,3],_0A5[1,3]],[_0A4[2,3],_0A5[2,3]], color = 'red')
	ax.plot3D([_0A5[0,3],_0A6[0,3]],[_0A5[1,3],_0A6[1,3]],[_0A5[2,3],_0A6[2,3]], color = 'green')
	ax.plot3D([_0A6[0,3],_0A7[0,3]],[_0A6[1,3],_0A7[1,3]],[_0A6[2,3],_0A7[2,3]], color = 'red')

	ax.plot3D([A02[0,3],_0A12[0,3]],[A02[1,3],_0A12[1,3]],[A02[2,3],_0A12[2,3]], color = 'red')
	ax.plot3D([_0A12[0,3],_0A22[0,3]],[_0A12[1,3],_0A22[1,3]],[_0A12[2,3],_0A22[2,3]], color = 'green')
	ax.plot3D([_0A22[0,3],_0A32[0,3]],[_0A22[1,3],_0A32[1,3]],[_0A22[2,3],_0A32[2,3]], color = 'red')
	ax.plot3D([_0A32[0,3],_0A42[0,3]],[_0A32[1,3],_0A42[1,3]],[_0A32[2,3],_0A42[2,3]], color = 'green')
	ax.plot3D([_0A42[0,3],_0A52[0,3]],[_0A42[1,3],_0A52[1,3]],[_0A42[2,3],_0A52[2,3]], color = 'red')
	ax.plot3D([_0A52[0,3],_0A62[0,3]],[_0A52[1,3],_0A62[1,3]],[_0A52[2,3],_0A62[2,3]], color = 'green')
	ax.plot3D([_0A62[0,3],_0A72[0,3]],[_0A62[1,3],_0A72[1,3]],[_0A62[2,3],_0A72[2,3]], color = 'red')


def Robot_RR_animado(theta_1, d1, a1, alpha_1, theta_2, d2, a2, alpha_2, theta_3, d3, a3, alpha_3, theta_4, d4, a4, alpha_4, theta_5, d5, a5, alpha_5, theta_6, d6, a6, alpha_6, theta_7, d7, a7, alpha_7,theta_12, d12, a12, alpha_12, theta_22, d22, a22, alpha_22, theta_32, d32, a32, alpha_32, theta_42, d42, a42, alpha_42, theta_52, d52, a52, alpha_52, theta_62, d62, a62, alpha_62, theta_72, d72, a72, alpha_72):
	
	#configuracion_grafica()
	#Robot_RR(theta_1, d1, a1, alpha_1, theta_2, d2, a2, alpha_2, theta_3, d3, a3, alpha_3, theta_4, d4, a4, alpha_4, theta_5, d5, a5, alpha_5, theta_6, d6, a6, alpha_6, theta_7, d7, a7, alpha_7,-theta_12, d12, -a12, alpha_12, -theta_22, d22, -a22, alpha_22, -theta_32, d32, -a32, alpha_32, -theta_42, d42, -a42, alpha_42, -theta_52, d52, -a52, alpha_52, -theta_62, d62, -a62, alpha_62, -theta_72, d72, -a72, alpha_72)
	#plt.draw()

	if theta_1 > 0:
		min = 1
		max = theta_1+1
		step = 1
	else: 
		min = -1
		max = theta_1-1
		step = -1

	for i in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, 0, d2, a2, alpha_2, 0, d3, a3, alpha_3, 0, d4, a4, alpha_4, 0, d5, a5, alpha_5, 0, d6, a6, alpha_6, 0, d7, a7, alpha_7,0, d12, -a12, alpha_12, 0, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_2 > 0:
		min = 1
		max = theta_2+1
		step = 1
	else: 
		min = -1
		max = theta_2-1
		step = -1

	for j in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, 0, d3, a3, alpha_3, 0, d4, a4, alpha_4, 0, d5, a5, alpha_5, 0, d6, a6, alpha_6, 0, d7, a7, alpha_7,0, d12, -a12, alpha_12, 0, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_3 > 0:
		min = 1
		max = theta_3+1
		step = 1
	else: 
		min = -1
		max = theta_3-1
		step = -1

	for k in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, 0, d4, a4, alpha_4, 0, d5, a5, alpha_5, 0, d6, a6, alpha_6, 0, d7, a7, alpha_7,0, d12, -a12, alpha_12, 0, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_4 > 0:
		min = 1
		max = theta_4+1
		step = 1
	else: 
		min = -1
		max = theta_4-1
		step = -1

	for l in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, 0, d5, a5, alpha_5, 0, d6, a6, alpha_6, 0, d7, a7, alpha_7,0, d12, -a12, alpha_12, 0, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_5 > 0:
		min = 1
		max = theta_5+1
		step = 1
	else: 
		min = -1
		max = theta_5-1
		step = -1

	for m in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, 0, d6, a6, alpha_6, 0, d7, a7, alpha_7,0, d12, -a12, alpha_12, 0, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_6 > 0:
		min = 1
		max = theta_6+1
		step = 1
	else: 
		min = -1
		max = theta_6-1
		step = -1

	for n in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, 0, d7, a7, alpha_7,0, d12, -a12, alpha_12, 0, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_7 > 0:
		min = 1
		max = theta_7+1
		step = 1
	else: 
		min = -1
		max = theta_7-1
		step = -1

	for o in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, o, d7, a7, alpha_7,0, d12, -a12, alpha_12, 0, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)
		
	if theta_12 > 0:
		min = 1
		max = theta_12+1
		step = 1
	else: 
		min = -1
		max = theta_12-1
		step = -1

	for i2 in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, o, d7, a7, alpha_7,-i2, d12, -a12, alpha_12, 0, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_22 > 0:
		min = 1
		max = theta_22+1
		step = 1
	else: 
		min = -1
		max = theta_22-1
		step = -1

	for j2 in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, o, d7, a7, alpha_7,-i2, d12, -a12, alpha_12, -j2, d22, -a22, alpha_22, 0, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_32 > 0:
		min = 1
		max = theta_32+1
		step = 1
	else: 
		min = -1
		max = theta_32-1
		step = -1

	for k2 in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, o, d7, a7, alpha_7,-i2, d12, -a12, alpha_12, -j2, d22, -a22, alpha_22, -k2, d32, -a32, alpha_32, 0, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_42 > 0:
		min = 1
		max = theta_42+1
		step = 1
	else: 
		min = -1
		max = theta_42-1
		step = -1

	for l2 in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, o, d7, a7, alpha_7,-i2, d12, -a12, alpha_12, -j2, d22, -a22, alpha_22, -k2, d32, -a32, alpha_32, -l2, d42, -a42, alpha_42, 0, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_52 > 0:
		min = 1
		max = theta_52+1
		step = 1
	else: 
		min = -1
		max = theta_52-1
		step = -1

	for m2 in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, o, d7, a7, alpha_7,-i2, d12, -a12, alpha_12, -j2, d22, -a22, alpha_22, -k2, d32, -a32, alpha_32, -l2, d42, -a42, alpha_42, -m2, d52, -a52, alpha_52, 0, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_6 > 0:
		min = 1
		max = theta_62+1
		step = 1
	else: 
		min = -1
		max = theta_62-1
		step = -1

	for n2 in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, o, d7, a7, alpha_7,-i2, d12, -a12, alpha_12, -j2, d22, -a22, alpha_22, -k2, d32, -a32, alpha_32, -l2, d42, -a42, alpha_42, -m2, d52, -a52, alpha_52, -n2, d62, -a62, alpha_62, 0, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

	if theta_72 > 0:
		min = 1
		max = theta_72+1
		step = 1
	else: 
		min = -1
		max = theta_72-1
		step = -1

	for o2 in range(min,max,step):
		ax.cla()
		configuracion_grafica()
		Robot_RR(i, d1, a1, alpha_1, j, d2, a2, alpha_2, k, d3, a3, alpha_3, l, d4, a4, alpha_4, m, d5, a5, alpha_5, n, d6, a6, alpha_6, o, d7, a7, alpha_7,-i2, d12, -a12, alpha_12, -j2, d22, -a22, alpha_22, -k2, d32, -a32, alpha_32, -l2, d42, -a42, alpha_42, -m2, d52, -a52, alpha_52, -n2, d62, -a62, alpha_62, -o2, d72, -a72, alpha_72)
		plt.draw()
		plt.pause(1e-3)

Robot_RR_animado(1,16.6,0,0  ,1,0,3,90  ,1,25.15,-3,-30  ,1,0,4.05,90  ,1,26.5,-4.05,-90  ,1,0,2.7,90  ,1,3.6,-2.7,-90  ,1,16.6,0,0  ,1,0,3,90  ,1,25.15,-3,-30  ,1,0,4.05,90  ,1,26.5,-4.05,-90  ,1,0,2.7,90  ,1,3.6,-2.7,-90)

plt.show()
print("Programado por: Fernando Espericueta")