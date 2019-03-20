import numpy as np

# assume x is in polar form, with a radius and an angle
def P2R(x):
    radius = x[0]
    angle = x[1]
    return radius * np.exp(1j*angle)

# assume x is in cartesian form, with a real and imaginary component
def R2P(x):
    return np.abs(x), np.angle(x)

R = [100]*4
va = np.array([1, 0])
vb = np.array([1, np.deg2rad(90)])
vc = np.array([1, np.deg2rad(180)])
vd = np.array([1, np.deg2rad(270)])
v_gens = [P2R(va), P2R(vb), P2R(vc), P2R(vd)]
# print(v_gens)


# C_net0 = 3.4710153285650266e-19F
# C_net1 = 3.5368579719566134e-19F
# C_net2 = 3.353951273969931e-19F
# C_net3 = 3.2749458120962368e-19F
# 1.9533199170382535e-18 -8.108209265564725e-19 -1.8728167290854811e-19 -6.081157847167302e-19
# -8.108209265564725e-19 2.1734322659530066e-18 -8.163561133609602e-19 -1.9256942883991256e-19
# -1.8728167290854811e-19 -8.163561133609602e-19 2.1730140796244824e-18 -8.33981165957981e-19
# -6.081157847167302e-19 -1.9256942883991256e-19 -8.33981165957981e-19 1.9621609607242476e-18
f = 1E15
omega = 2*np.pi*f

Z = np.zeros([4,4])
Z = -np.array([[1.9533199170382535e-18, -8.108209265564725e-19, -1.8728167290854811e-19, -6.081157847167302e-19],
[-8.108209265564725e-19, 2.1734322659530066e-18, -8.163561133609602e-19, -1.9256942883991256e-19],
[-1.8728167290854811e-19, -8.163561133609602e-19, 2.1730140796244824e-18, -8.33981165957981e-19],
[-6.081157847167302e-19, -1.9256942883991256e-19, -8.33981165957981e-19, 1.9621609607242476e-18]])
for i in range(len(Z)):
    Z[i][i] = -sum(Z[i])
#Z is now the mutual capacitance matrix
Z = np.reciprocal(omega*Z)
# print(Z)
#Z is now the impedance

# A = np.zeros([4,4])
# A = np.reciprocal(-Z)
A = -Z
for i in range(len(Z)):
    A[i][i] = sum(Z[i]) + (1/R[i])
# print(A)

v_nodes = np.linalg.solve(A,np.divide(v_gens,R))
print(v_nodes)
