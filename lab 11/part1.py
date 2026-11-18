import numpy as np
import matplotlib.pyplot as plt

def F(x):
    return 2*x[0]**2 + x[1]**2 - 2*x[0]*x[1] + 2*x[0] + x[1]

def grad_F(x):
    return np.array([4*x[0] - 2*x[1] + 2, -2*x[0] + 2*x[1] + 1])

def g(x):
    return 1.0 - x[0]

def grad_g(x):
    return np.array([-1.0, 0.0])

def simple_gradient_descent(F_func, grad_func, x0, eps=0.01):
    x = np.array(x0, float)
    h = 1.0
    k = 0
    while True:
        g_val = grad_func(x)
        if np.linalg.norm(g_val) <= eps or k > 1000:
            break
        while True:
            x_new = x - h * g_val
            if F_func(x_new) < F_func(x):
                break
            h /= 2.0
            if h < 1e-10:
                return x, k
        x = x_new
        k += 1
        h = min(h * 1.2, 1.0) 
    return x, k

if __name__ == "__main__":
    x0 = np.array([5.0, 2.0])
    r0, c, eps = 1.0, 10.0, 0.01
    x = x0.copy()
    r = r0
    
    history_x = [x.copy()]
    history_g = []
    history_r = []
    
    print(f"{'k':>2} | {'r':>6} | {'x(r)':>14} | {'ітер':>4} | {'x1':>8} | {'x2':>8} | {'F(x)':>8} | {'g(x)':>8} | {'r*P':>8} | {'2r*g':>8}")
    print("-" * 105)
    
    for k in range(10):
        def Phi(xx): 
            g_val = max(0.0, g(xx))
            return F(xx) + r * g_val**2
        def grad_Phi(xx):
            g_val = max(0.0, g(xx))
            return grad_F(xx) + 2 * r * g_val * grad_g(xx)
            
        x, iters = simple_gradient_descent(Phi, grad_Phi, x, eps)
        
        gx = g(x)
        Px = max(0.0, gx)**2
        rP = r * Px
        lam = 2 * r * max(0.0, gx)
        
        history_x.append(x.copy())
        history_g.append(gx)
        history_r.append(r)
        
        print(f"{k:>2} | {r:>6.0f} | ({x[0]:>6.4f}; {x[1]:>6.4f}) | {iters:>4} | {x[0]:>8.4f} | {x[1]:>8.4f} | {F(x):>8.4f} | {gx:>8.4f} | {rP:>8.4f} | {lam:>8.4f}")
        
        if rP <= eps: break
        r *= c

    x_hat = x - g(x) * grad_g(x) / (grad_g(x) @ grad_g(x))
    print(f"\nВідповідь x* ≈ ({x[0]:.4f}, {x[1]:.4f}), порушення {g(x):.4f}")
    print(f"Після проєкції на межу: ({x_hat[0]:.4f}, {x_hat[1]:.4f}), g = {g(x_hat):.4f}, F = {F(x_hat):.4f}")
    print(f"Аналітично точно: (1.0000, 0.5000), F = 3.7500, lambda = 5")

    x1_vals = np.linspace(-3, 6, 100)
    x2_vals = np.linspace(-4, 4, 100)
    X1, X2 = np.meshgrid(x1_vals, x2_vals)
    Z = 2*X1**2 + X2**2 - 2*X1*X2 + 2*X1 + X2
        
    plt.figure(figsize=(10, 6))
    plt.contour(X1, X2, Z, levels=30, cmap="viridis", alpha=0.5)
    plt.axvline(x=1, color='g', lw=2, label="Межа g(x)=0 (x1=1)")
    plt.axvspan(1, 6, color='green', alpha=0.1, label="Допустима область Ω")
    hist = np.array(history_x)
    plt.plot(hist[:, 0], hist[:, 1], "ro-", label="Траєкторія x(r)")
    plt.plot(x0[0], x0[1], "k*", ms=12, label="Старт (5, 2)")
    plt.plot(-1.5, -2.0, "bx", ms=10, label="Безумовний мінімум")
    plt.plot(1.0, 0.5, "rX", ms=10, label="Умовний мінімум")
    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.title("Зовнішній штраф: Траєкторія x(r)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig("lab 11/trajectory.png", dpi=150, bbox_inches="tight")
    plt.figure(figsize=(8, 5))
    plt.semilogx(history_r, history_g, "bo-")
    plt.xlabel("r (логарифмічна шкала)")
    plt.ylabel("Порушення g(x(r))")
    plt.title("Зменшення порушення зі зростанням r")
    plt.grid(alpha=0.3)
    plt.savefig("lab 11/g_vs_r.png", dpi=150, bbox_inches="tight")
    plt.show()