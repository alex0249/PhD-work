import numpy as np
import matplotlib.pyplot as plt

N1 = 2        # Number of elements
N0 = N1 + 1       # Number of vertices
h = 1.0 / N1      # Element length
i = np.arange(1, N0 + 1)  # i = 1, 2, ..., N0
x = (i - 1) * h           # Xi = (i - 1)h global nodes

# Each element connects two neighbouring vertices
elements = np.column_stack((
    np.arange(N1),
    np.arange(1, N0),
))


A = np.zeros((N0, N0))  # Global matrix: N0 × N0
B = np.zeros((N0, 1))  # Global vector: N0 × 1
# Same parent shape functions for every element
def phi(xi):
    return np.array([
        (1 - xi) / 2,
        (1 + xi) / 2
    ])

# Derivatives with respect to xi: constant for linear elements
dphi_dxi = np.array([-0.5, 0.5])

for e in range(N1):
    dx_dxi = h/2 #Jacobian
    J_e = h/2 # det
    J_e_inv = 2/h # det inverse
    dphi_dx = dphi_dxi * J_e_inv # derivative with respect to x

    x_left = x[e]
    x_right = x[e + 1]
    # Local stiffness matrix for element e
    A_e = (1 / h) * np.array([
        [ 1, -1],
        [-1,  1]
    ])

    B_e = (h / 6) * np.array([
        [2 * x_left + x_right],
        [x_left + 2 * x_right]
    ])

    # Assemble local matrices into global matrices
    A[e,     e]     += A_e[0, 0]
    A[e,     e + 1] += A_e[0, 1]
    A[e + 1, e]     += A_e[1, 0]
    A[e + 1, e + 1] += A_e[1, 1]

    B[e,     0] += B_e[0, 0]
    B[e + 1, 0] += B_e[1, 0]

# Prescribe u = 0 at the first and last nodes
def dirichlet_bc(A, B):

    for k in [0, N0 - 1]:
        T_k = 0.0

        B -= A[:, k:k+1] * T_k  # Subtract column k × prescribed value
        A[k, :] = 0             # Zero row k
        A[:, k] = 0             # Zero column k
        A[k, k] = 1             # Set diagonal entry to 1
        B[k, 0] = T_k           # Set prescribed value
    return A, B

# Case 1: Dirichlet at both endpoints
A_dir, B_dir = dirichlet_bc(A.copy(), B.copy())
U = np.linalg.solve(A_dir, B_dir)
y= np.linspace(0, 1, 100)
u = y*(1-y**2)/6 #exact solution

plt.plot(x, U, 'o-', label='FEM solution (N1 = {})'.format(N1))
plt.plot(y, u, 'x--', label='Exact solution')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title( "Dirichlet BC")
plt.legend()
plt.show()


def neumann_dirichlet_bc(A, B, g=0.0):
    # Neumann at x=1: T'(1) = g
    B[-1, 0] += g

    # Dirichlet at x=0: T(0) = 0
    T_left = 0.0
    B -= A[:, 0:1] * T_left

    A[0, :] = 0.0
    A[:, 0] = 0.0
    A[0, 0] = 1.0
    B[0, 0] = T_left

    return A, B

A_both, B_both = neumann_dirichlet_bc(A, B)
print(B_both)

U_both = np.linalg.solve(A_both, B_both)
y= np.linspace(0, 1, 100)
u = y*(3-y**2)/6 #exact solution

plt.plot(x, U_both, 'o-', label='FEM solution (N1 = {})'.format(N1))
plt.plot(y, u, 'x--', label='Exact solution')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title( "Dirichlet + Neumann BCs")
plt.legend()
plt.show()