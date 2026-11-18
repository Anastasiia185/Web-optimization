import numpy as np
import matplotlib.pyplot as plt

def F(x):
    return 2*x[0]**2 + x[1]**2 - 2*x[0]*x[1] + 2*x[0] + x[1]

def grad_F(x):
    return np.array([4*x[0] - 2*x[1] + 2, -2*x[0] + 2*x[1] + 1])

def g1(x): return 1.0 - x[0]
def grad_g1(x): return np.array([-1.0, 0.0])

def g2(x): return x[1]
def grad_g2(x): return np.array([0.0, 1.0])

def golden_section(phi, a, b, eps=1e-6):
    tau = (np.sqrt(5) - 1) / 2
    x1, x2 = b - tau * (b - a), a + tau * (b - a)
    f1, f2 = phi(x1), phi(x2)
    while (b - a) > eps:
        if f1 <= f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - tau * (b - a); f1 = phi(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + tau * (b - a); f2 = phi(x2)
    return (a + b) / 2

def svenn(phi, x0=0.0, delta=0.01):
    x, f0 = x0, phi(x0)
    x_prev, f_prev = x, f0
    x += delta; f_curr = phi(x)
    if f_curr > f0:
        delta = -delta; x, f_curr = x0 + delta, phi(x0 + delta)
    while f_curr < f_prev:
        delta *= 2.0; x_prev, f_prev = x, f_curr
        x += delta; f_curr = phi(x)
    return min(x_prev, x), max(x_prev, x)

def line_search(F, x, d):
    phi = lambda h: F(x + h*d)
    a, b = svenn(phi, 0.0, 0.01)
    return golden_section(phi, a, b, 1e-6)

def fletcher_reeves(F_func, grad_func, x0, eps=0.01, n_dim=2):
    x = np.array(x0, float)
    g_val = grad_func(x)
    d = -g_val
    k = 0
    while np.linalg.norm(g_val) > eps and k < 100:
        h = line_search(F_func, x, d)
        x = x + h * d
        g_new = grad_func(x)
        beta = (g_new @ g_new) / (g_val @ g_val)
        if (k + 1) % n_dim == 0: beta = 0.0
        d = -g_new + beta * d
        g_val = g_new
        k += 1
    return x, k

def sphere_descent(F_func, grad_func, x0, eps=0.01):
    x = np.array(x0, float)
    h = 1.0
    k = 0
    while True:
        g_val = grad_func(x)
        if np.linalg.norm(g_val) <= eps or k > 1000: break
        while True:
            x_new = x - h * (g_val / np.linalg.norm(g_val))
            if F_func(x_new) < F_func(x):
                break
            h /= 2.0
            if h < 1e-10: return x, k
        x = x_new
        k += 1
        h = min(h * 1.2, 1.0)
    return x, k

if __name__ == "__main__":
    print("="*85)
    
    x0_barrier = np.array([5.0, -2.0]) 
    mu = 10.0
    eps = 0.01
    x = x0_barrier.copy()
    
    print(f"{'k':>2} | {'mu':>8} | {'x1':>8} | {'x2':>8} | {'F(x)':>8} | {'g1(x)':>8} | {'g2(x)':>8} | {'lam1':>8} | {'lam2':>8}")
    print("-" * 85)
    
    for k in range(10):
        def Phi(xx):
            g1v, g2v = g1(xx), g2(xx)
            if g1v >= 0 or g2v >= 0: return float('inf')
            return F(xx) - mu * (np.log(-g1v) + np.log(-g2v))
            
        def grad_Phi(xx):
            g1v, g2v = g1(xx), g2(xx)
            return grad_F(xx) - mu / g1v * grad_g1(xx) - mu / g2v * grad_g2(xx)
            
        x, iters = fletcher_reeves(Phi, grad_Phi, x, eps)
        
        g1x, g2x = g1(x), g2(x)
        lam1 = -mu / g1x if g1x < 0 else 0
        lam2 = -mu / g2x if g2x < 0 else 0
        
        print(f"{k:>2} | {mu:>8.4f} | {x[0]:>8.4f} | {x[1]:>8.4f} | {F(x):>8.4f} | {g1x:>8.4f} | {g2x:>8.4f} | {lam1:>8.4f} | {lam2:>8.4f}")
        
        if mu * 2 <= eps: break
        mu /= 10.0

    print("Аналітично точно: x* = (1; 0), F* = 4, lambda = (5; 1)\n")

    print("="*85)
    
    def g_circle(x): return 16.0 - x[0]**2 - x[1]**2
    def grad_g_circle(x): return np.array([-2.0*x[0], -2.0*x[1]])
    
    x0_circle = np.array([5.0, 1.0]) 
    mu = 10.0
    x = x0_circle.copy()
    
    print(f"{'k':>2} | {'mu':>8} | {'x1':>8} | {'x2':>8} | {'F(x)':>8} | {'g(x)':>8} | {'||x||':>8} | {'lam':>8}")
    print("-" * 75)
    
    for k in range(10):
        def Phi_c(xx):
            gv = g_circle(xx)
            if gv >= 0: return float('inf')
            return F(xx) - mu * np.log(-gv)
            
        def grad_Phi_c(xx):
            gv = g_circle(xx)
            return grad_F(xx) - mu / gv * grad_g_circle(xx)
            
        x, iters = fletcher_reeves(Phi_c, grad_Phi_c, x, eps)
        
        gx = g_circle(x)
        norm_x = np.linalg.norm(x)
        lam = -mu / gx if gx < 0 else 0
        
        print(f"{k:>2} | {mu:>8.4f} | {x[0]:>8.4f} | {x[1]:>8.4f} | {F(x):>8.4f} | {gx:>8.4f} | {norm_x:>8.4f} | {lam:>8.4f}")
        
        if mu * 1 <= eps: break
        mu /= 10.0

    print("="*85)
    
    r_test = 100.0
    def Phi_ext(xx): 
        gv = max(0.0, g1(xx))
        return F(xx) + r_test * gv**2
    def grad_Phi_ext(xx):
        gv = max(0.0, g1(xx))
        return grad_F(xx) + 2 * r_test * gv * grad_g1(xx)
        
    x_start = np.array([5.0, 2.0])
    
    x_fr, it_fr = fletcher_reeves(Phi_ext, grad_Phi_ext, x_start, eps)
    print(f"Флетчер-Рівс:             ітерацій = {it_fr}, x = ({x_fr[0]:.4f}; {x_fr[1]:.4f})")
    
    x_sp, it_sp = sphere_descent(Phi_ext, grad_Phi_ext, x_start, eps)
    print(f"Сферична метрика (ЛР 8):  ітерацій = {it_sp}, x = ({x_sp[0]:.4f}; {x_sp[1]:.4f})")

    x1_vals = np.linspace(-2, 6, 100)
    x2_vals = np.linspace(-4, 3, 100)
    X1, X2 = np.meshgrid(x1_vals, x2_vals)
    Z = 2*X1**2 + X2**2 - 2*X1*X2 + 2*X1 + X2
    plt.figure(figsize=(10, 6))
    plt.contour(X1, X2, Z, levels=30, cmap="viridis", alpha=0.5)
    plt.axvline(x=1, color='g', lw=2, label="g1(x)=0 (x1=1)")
    plt.axhline(y=0, color='b', lw=2, label="g2(x)=0 (x2=0)")
    plt.axvspan(1, 6, color='green', alpha=0.1)
    plt.axhspan(-4, 0, color='blue', alpha=0.1)
    plt.plot(1, 0, "rX", ms=10, label="Умовний мінімум (1, 0)")
    plt.xlabel("x1"); plt.ylabel("x2")
    plt.title("Два обмеження (кут)")
    plt.legend(); plt.grid(alpha=0.3)
    plt.savefig("lab 11/part2.png", dpi=150, bbox_inches="tight")
    plt.show()