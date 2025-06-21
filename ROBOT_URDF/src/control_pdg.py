#!/usr/bin/env python3

import rospy
import matplotlib.pyplot as plt
from sensor_msgs.msg import JointState
from markers import *
from robot_functions import *
from roslib import packages
import rbdl

# ***************************************************************************************************
if __name__ == '__main__':

  rospy.init_node("control_pdg")
  pub = rospy.Publisher('joint_states', JointState, queue_size=10)
  bmarker_actual  = BallMarker(color['RED'])
  bmarker_deseado = BallMarker(color['GREEN'])

  # Archivos donde se almacenara los datos
  fqact = open("/home/user/qactual.txt", "w")
  fqdes = open("/home/user/qdeseado.txt", "w")
  fxact = open("/home/user/xactual.txt", "w")
  fxdes = open("/home/user/xdeseado.txt", "w")
  
  # Nombres de las articulaciones
  jnames = ['joint_1', 'joint_2', 'joint_3','joint_4', 'joint_5', 'joint_6', 'joint_7', 'joint_8']
 
  # Objeto (mensaje) de tipo JointState
  jstate = JointState()
  # Valores del mensaje
  jstate.header.stamp = rospy.Time.now()
  jstate.name = jnames
  
  # =============================================================
  # Configuracion articular inicial (en radianes)
  q = np.array([0., 0., 0., 0., 0., 0., 0.,0.])
  # Velocidad inicial
  dq = np.array([0., 0., 0., 0., 0., 0., 0.,0.])
  # Configuracion articular deseada
  xd = np.array([1.3 ,-1.12 , 1.60110448])
  #q0 = np.array([pi/3,-pi/6,pi/5,-pi/12,pi/12, pi/15,0,0])
  q0=np.array([0,0,0,0,0,0,0,0])
  qdes = ikine(xd, q0)
  # =============================================================
  
  # Posicion resultante de la configuracion articular deseada
  xdes = fkine(qdes)[0:3,3]
  # Copiar la configuracion articular en el mensaje a ser publicado
  jstate.position = q
  pub.publish(jstate)
  
  # Modelo RBDL
  modelo = rbdl.loadModel('lab_ws/src/ProyectoRobotica/ROBOT_URDF/urdf/ROBOT_URDF.urdf')
  ndof   = modelo.q_size     # Grados de libertad
  
  # Frecuencia del envio (en Hz)
  freq = 20
  dt = 1.0/freq
  rate = rospy.Rate(freq)
  #valores maximas y minimas de las articulaciones
  joint2_min=-1.57;  joint2_max=1.57
  joint3_min=-3.142; joint3_max=0.785
  joint4_min=-1.57;  joint4_max=1.57
  joint5_min=0;      joint5_max=0.4
  joint6_min=-4.04;  joint6_max=0.785
  joint7_min=0;      joint7_max=0.05
  joint8_min=-0.05;  joint8_max=0
  
  # Simulador dinamico del robot
  robot = Robot(q, dq, ndof, dt)

  # Se definen las ganancias del controlador
  Kp = np.eye(ndof)*np.array([60, 60, 55, 20, 45, 40, 5, 5])
  Kd = np.eye(ndof)*np.array([30, 30, 30, 30, 39, 38, 3, 3])

  # Inicializar acumulador para el error integral
  t,i=0,0
  MM=[]; CC=[]; GG=[]
  while not rospy.is_shutdown():
    if (t>=15): break

    # Calcular el error actual
    ea = qdes - q

    # Leer valores del simulador
    q  = robot.read_joint_positions()
    dq = robot.read_joint_velocities()

    q[1] = max(min(q[1], joint2_max), joint2_min)
    q[2] = max(min(q[2], joint3_max), joint3_min)
    q[3] = max(min(q[3], joint4_max), joint4_min)
    q[4] = max(min(q[4], joint5_max), joint5_min)
    q[5] = max(min(q[5], joint6_max), joint6_min)
    q[6] = max(min(q[6], joint7_max), joint7_min)
    q[7] = max(min(q[7], joint8_max), joint8_min)     

    # Posicion actual del efector final
    x = fkine(q)[0:3,3]
    # Tiempo actual (necesario como indicador para ROS)
    jstate.header.stamp = rospy.Time.now()

    # Almacenamiento de datos
    fxact.write(str(t)+' '+str(x[0])+' '+str(x[1])+' '+str(x[2])+'\n')
    fxdes.write(str(t)+' '+str(xdes[0])+' '+str(xdes[1])+' '+str(xdes[2])+'\n')
    fqact.write(str(t)+' '+str(q[0])+' '+str(q[1])+' '+ str(q[2])+' '+ str(q[3])+' '+str(q[4])+' '+str(q[5])+'\n ')
    fqdes.write(str(t)+' '+str(qdes[0])+' '+str(qdes[1])+' '+ str(qdes[2])+' '+ str(qdes[3])+' '+str(qdes[4])+' '+str(qdes[5])+'\n ')

    # ----------------------------
    # Control dinamico (COMPLETAR)
    # ----------------------------
    zeros = np.zeros(ndof)
    g = np.zeros(ndof)
    rbdl.InverseDynamics(modelo, q, zeros, zeros, g)
    
    q[1] = max(min(q[1], joint2_max), joint2_min)
    q[2] = max(min(q[2], joint3_max), joint3_min)
    q[3] = max(min(q[3], joint4_max), joint4_min)
    q[4] = max(min(q[4], joint5_max), joint5_min)
    q[5] = max(min(q[5], joint6_max), joint6_min)
    q[6] = max(min(q[6], joint7_max), joint7_min)
    q[7] = max(min(q[7], joint8_max), joint8_min)
    
    u = g + Kp@(qdes-q) - Kd@dq # Reemplazar por la ley de control

    # Simulacion del robot
    robot.send_command(u)

    # Publicacion del mensaje
    jstate.position = q
    pub.publish(jstate)
    bmarker_deseado.xyz(xdes)
    bmarker_actual.xyz(x)
    t+=dt
    # Esperar hasta la siguiente  iteracion
    rate.sleep()
    e = q-qdes    

    # ***********************************************

    # Lectura del modelo del robot a partir de URDF (parsing)
    modelo = rbdl.loadModel('lab_ws/src/ProyectoRobotica/ROBOT_URDF/urdf/ROBOT_URDF.urdf')
    
    # Armys
    tau   = np.zeros(ndof)          # Para torque
    c     = np.zeros(ndof)          # Para el vector de Coriolis+centrifuga
    M     = np.zeros([ndof, ndof])  # Para la matriz de inercia
    eye     = np.eye(ndof)               # Vector identidad
    
    # Calcular vector de gravedad, vector de Coriolis/centrifuga,
    # y matriz M usando InverseDynamics
    tau_c = np.zeros(ndof)
    rbdl.InverseDynamics(modelo, q, dq, zeros, tau_c)
    c = tau_c - g
    
    for i in range(ndof):
        tau_m = np.zeros(ndof)
        rbdl.InverseDynamics(modelo, q, zeros, eye[i], tau_m)
        M[i] = tau_m - g

    GG.append(np.linalg.norm(g))
    CC.append(np.linalg.norm(c))
    MM.append(np.linalg.norm(M))

    if (np.linalg.norm(e)<0.01): break

  fqact.close()
  fqdes.close()
  fxact.close()
  fxdes.close()

  

  # Cargar los datos desde los archivos de registro
  xcurrent = np.loadtxt("/home/user/xactual.txt")
  xdesired = np.loadtxt("/home/user/xdeseado.txt")
  q = np.loadtxt("/home/user/qactual.txt")
  
  # Calcular el número de muestras y el vector de tiempo
  num_samples = xcurrent.shape[0]
  dt = 1.0 / freq  # Intervalo de tiempo, basado en la frecuencia de control
  time_steps = np.linspace(0, num_samples * dt, num_samples)
  
  # ****************************************************************************************************

  # Crear subplots
  fig, axs = plt.subplots(3, 1, figsize=(10, 8))
  fig.suptitle("Posición actual y deseada en los ejes x, y, z")
  
  axs[0].plot(time_steps, MM, color="blue")
  axs[0].set_ylabel("Modulo de M")
  
  axs[1].plot(time_steps, CC, color="blue")
  axs[1].set_ylabel("Modulo de C")
  
  axs[2].plot(time_steps, GG, color="blue")
  axs[2].set_ylabel("Modulo de g")
  
  # Configurar el eje x como el tiempo en cada subplot
  for ax in axs:
      ax.set_xlabel("Tiempo")

  # Mostrar el gráfico
  plt.tight_layout(rect=[0, 0.03, 1, 0.95])
  plt.show()

  # ****************************************************************************************************

  # Crear subplots
  fig, axs = plt.subplots(3, 1, figsize=(10, 8))
  fig.suptitle("Posición actual y deseada en los ejes x, y, z")
  
  # Coordenada X
  axs[0].plot(time_steps, xcurrent[:, 1], label="Posición actual X", color="blue")
  axs[0].plot(time_steps, xdesired[:, 1], label="Posición deseada X", color="red", linestyle="--")
  axs[0].set_ylabel("Posición X")
  axs[0].legend()
  
  # Coordenada Y
  axs[1].plot(time_steps, xcurrent[:, 2], label="Posición actual Y", color="blue")
  axs[1].plot(time_steps, xdesired[:, 2], label="Posición deseada Y", color="red", linestyle="--")
  axs[1].set_ylabel("Posición Y")
  axs[1].legend()
  
  # Coordenada Z
  axs[2].plot(time_steps, xcurrent[:, 3], label="Posición actual Z", color="blue")
  axs[2].plot(time_steps, xdesired[:, 3], label="Posición deseada Z", color="red", linestyle="--")
  axs[2].set_ylabel("Posición Z")
  axs[2].legend()
  
  # Configurar el eje x como el tiempo en cada subplot
  for ax in axs:
      ax.set_xlabel("Tiempo")

  # Mostrar el gráfico
  plt.tight_layout(rect=[0, 0.03, 1, 0.95])
  plt.show()