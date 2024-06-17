from qiskit_aer import AerSimulator
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import UnitaryGate

import numpy as np
from scipy.linalg import sqrtm

# setting parameters
x0 = np.array([0.5,0.5,0.5,0.5]) #initial distribution
output_file_name = 'results/qsout.txt'
shot_num = 2e5
ntmax = 10
dt = 0.05

C = 0.25
D = 0.75
tau01 = 1
tau10 = 5
tau12 = 1
tau02 = 2
tau23 = 3
tau20 = 5
tau30 = 5

def Gamma(x):
  if(x>0):
    return 1-np.exp(-x)
  else:
    return 0

R01 = Gamma(C)*Gamma(D)/tau01
R02 = Gamma(C)*(1-Gamma(D))/tau02
R10 = Gamma(D)/tau10
R12 = Gamma(C)*(1-Gamma(D))/tau12
R20 = (1-Gamma(C))/tau20
R23 = 1/tau23
R30 = 1/tau30

P01 = R01*dt
P02 = R02*dt
P10 = R10*dt
P12 = R12*dt
P20 = R20*dt
P23 = R23*dt
P30 = R30*dt

M = np.array([[1-P01-P02, P01, P02, 0],
                [P10, 1-P10-P12, P12, 0],
                [P20, 0, 1-P20-P23, P23],
                [P30, 0, 0, 1-P30]], dtype=np.complex128)

def gram_schmidt(vectors):
  basis = []
  for v in vectors:
    w = v - np.sum( np.dot(v,b)*b  for b in basis )
    if (w > 1e-10).any():  
      basis.append(w/np.linalg.norm(w)) 
  return np.array(basis).T

M = M.T/(np.linalg.norm(M,ord=2)+0.01)
B = (M+M.T)/2
C = (M-M.T)/(2j)

F1 = (B + 1j*sqrtm(np.eye(4)-B@B.T))
F2 = (B - 1j*sqrtm(np.eye(4)-B@B.T))
F3 = (1j*C - sqrtm(np.eye(4)-C@np.conjugate(C.T)))
F4 = (1j*C + sqrtm(np.eye(4)-C@np.conjugate(C.T)))

F1_gate = UnitaryGate(F1,label="F1").control(2,ctrl_state=0b00)
F2_gate = UnitaryGate(F2,label="F2").control(2,ctrl_state=0b01)
F3_gate = UnitaryGate(F3,label="F3").control(2,ctrl_state=0b10)
F4_gate = UnitaryGate(F4,label="F4").control(2,ctrl_state=0b11)

def create_circ(x0):
  q = QuantumRegister(4)
  c = ClassicalRegister(4)
  qc = QuantumCircuit(q,c)
  init_mat = gram_schmidt(np.vstack((np.array(x0.reshape((1,4))),np.random.rand(3,4))))
  qc.unitary(init_mat,[2,3],label='Init')
  # qc.initialize(np.kron(x0,[1,0,0,0]))
  qc.h(0)
  qc.h(1)
  qc.append(F1_gate,[0,1,2,3])
  qc.append(F2_gate,[0,1,2,3])
  qc.append(F3_gate,[0,1,2,3])
  qc.append(F4_gate,[0,1,2,3])
  qc.h(0)
  qc.h(1)
  qc.measure(q,c)
  return qc

# main calculation
with open(output_file_name, 'w') as f:
  for t in range(ntmax):
    qc = create_circ(x0)
    backend = AerSimulator()
    qc_compiled = transpile(qc, backend)
    job_sim = backend.run(qc_compiled, shots=shot_num)
    result_sim = job_sim.result()
    counts = result_sim.get_counts(qc_compiled)
    
    result = []
    
    for i in range(4):
      index = format(i,"02b")[0]+format(i,"02b")[1]+'00'
      try:
        result.append(np.sqrt(counts[index]))
      except:
        result.append(0)
    prob2 = np.array(result)
    prob2 = prob2/np.linalg.norm(prob2)
    prob = prob2/np.sum(prob2)
    line = '\t'.join(['{}'.format(val) for val in prob])
    x0 = prob2
    
    f.write(line + '\n')
    print("t = {:5.2f}h :{}".format((t+1)*dt, line))