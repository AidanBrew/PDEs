import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def numerical_solution(D, R, f, x, t):

    x = np.asarray(x, dtype=float)
    t = np.asarray(t, dtype=float)

    N = len(x) # Number of nodes
    h = np.diff(x) # Node spacings

    # Control volume lengths
    V = np.zeros(N) 
    V[0] = h[0]/2
    V[1:N-1] = (h[0:N-2] + h[1:N-1])/2
    V[N-1] = h[N-2]/2

    u0 = f(x)

    # Solves G
    solution = solve_ivp(
        fun=lambda time, u: Gfunc(u, D, R, N, h, V),
        t_span=(t[0], t[-1]),
        y0=u0,
        t_eval=t,
        method="BDF",
        rtol=1e-6,
        atol=1e-8
    )

    return solution.y

def Gfunc(u, D, R, N, h, V):

    G = np.zeros(N)

    # First node 
    DE = 0.5 * (D(u[0]) + D(u[1]))
    EAvg = (u[1]-u[0])/h[0]

    G[0] = (1 / V[0]) * (-DE * EAvg) + R(u[0])

    # Interior nodes
    for i in range(1, N-1):
        
        DE = 0.5 * (D(u[i]) + D(u[i+1]))
        DW = 0.5 * (D(u[i-1]) + D(u[i]))

        EAvg = (u[i+1]-u[i])/h[i]
        WAvg = (u[i] - u[i-1])/h[i-1]
        
        G[i] = (1 / V[i]) * (DW * WAvg - DE * EAvg) + R(u[i])
    
    # Last node
    DW = 0.5 * (D(u[N-2]) + D(u[N-1]))
    WAvg = (u[N-1] - u[N-2])/h[N-2]

    G[N-1] = (1 / V[N-1]) * (DW * WAvg) + R(u[N-1])

    return G

def main():
    
    m = 4 # Diffusion exponent
    K = 1 # Carrying capacity
    
    # Non linear diffusivity
    alpha = 1     
    def D(u): 
        return -alpha * (u/K)**m
    
    # Non linear reaction
    beta = 1
    def R(u):
        return beta*u*(1-u/K)
    
    L = 40 # Length of domain (mm)
    T = 36 # End time (hours)
    N = 201 # Number of nodes 
    M = 360 # Number of time steps

    def f(x):
        return 1.0 * (x<=L/4) + 1.0 * (x >= 3*L/4)
    
    x = np.linspace(0, L, N)

    atol = 1e-10 # Absolute tolerance for Newtons method
    rtol = 1e-10 # Relative tolerance for Newtons method
    maxiters = 100 # Maximum number of iterations

    # Animation of numerical solution
    t = np.arange(0, T + 0.5, 0.5)

    u1 = numerical_solution(D, R, f, x, t)

    fig, ax = plt.subplots()

    line, = ax.plot(x, u1[:, 0], "b-", linewidth=2)

    ax.set_ylim(0, 1)
    ax.set_xlabel("x", fontsize=14)
    ax.set_ylabel("u(x,t)", fontsize=14)

    time_text = ax.text(32, 0.05, f"t = {t[0]:.2f}", fontsize=14)

    def update(n):
        line.set_ydata(u1[:, n])
        time_text.set_text(f"t = {t[n]:.2f}")
        return line, time_text

    ani = FuncAnimation(
        fig,
        update,
        frames=len(t),
        interval=100,
        blit=True
    )

    plt.show()

    plt.figure()

    for time_value in [0, 6, 12, 24, 36]:
        j = np.argmin(np.abs(t - time_value))
        plt.plot(x, u1[:, j], linewidth=2, label=f"t = {t[j]:.1f}")

    plt.ylim(0, 1)
    plt.xlabel("x")
    plt.ylabel("u(x,t)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
