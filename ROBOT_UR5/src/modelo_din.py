#!/usr/bin/env python3
import rbdl
import numpy as np

def MatM(q):
    for i in range(ndof):
        tau_m = np.zeros(ndof)
        rbdl.InverseDynamics(modelo, q, zeros, e[i], tau_m)
        M[i] = tau_m - g
    return M

def MatC(q,dq):
    tau_c = np.zeros(ndof)
    rbdl.InverseDynamics(modelo, q, dq, zeros, tau_c)
    return tau_c - g

def Matg(q):
    rbdl.InverseDynamics(modelo, q, zeros, zeros, g)
    return g

if __name__ == '__main__':

  # Lectura del modelo del robot a partir de URDF (parsing)
  modelo = rbdl.loadModel('lab_ws/src/ProyectoRobotica/ROBOT_URDF/urdf/ROBOT_URDF.urdf')
  # Grados de libertad
  ndof = modelo.q_size

  # Configuracion articular
  q = np.array([0.5, 0.2, 0.3, 0.8, 0.5, 0.6])
  # Velocidad articular
  dq = np.array([0.8, 0.7, 0.8, 0.6, 0.9, 1.0])
  # Aceleracion articular
  ddq = np.array([0.2, 0.5, 0.4, 0.3, 1.0, 0.5])
  
  # Armys
  zeros = np.zeros(ndof)          # Vector de ceros
  tau   = np.zeros(ndof)          # Para torque
  g     = np.zeros(ndof)          # Para la gravedad
  c     = np.zeros(ndof)          # Para el vector de Coriolis+centrifuga
  M     = np.zeros([ndof, ndof])  # Para la matriz de inercia
  e     = np.eye(ndof)            # Vector identidad
  
  # Torque dada la configuracion del robot
  rbdl.InverseDynamics(modelo, q, dq, ddq, tau)
  g=Matg(q)
  C=MatC(q,dq)
  M=MatM(q); 

  # Parte 2: Calcular M y los efectos no lineales b usando las funciones
  # CompositeRigidBodyAlgorithm y NonlinearEffects. Almacenar los resultados
  # en los arreglos llamados M2 y b2
  b2 = np.zeros(ndof)          # Para efectos no lineales
  M2 = np.zeros([ndof, ndof])  # Para matriz de inercia
  
  rbdl.CompositeRigidBodyAlgorithm(modelo, q, M2)

  rbdl.NonlinearEffects(modelo, q, dq, b2)
  
  # Parte 3: Verificacion de la expresion de la dinamica
  print("\n Verificación:")
  print("tau")
  print(tau)
  tau_verif = M @ ddq + C + g
  print("tau_verif")
  print(tau_verif)