import numpy as np
import matplotlib.pyplot as plt

def F(x):
    x1, x2 = x[0], x[1]
    return 2*(x1**2) + (x2**2) - 2*x1*x2 + 2*x1 + x2

def grad_F(x):
    x1, x2 = x[0], x[1]
    return np.array([4*x1 - 2*x2 + 2, -2*x1 + 2*x2 + 1])

def golden_section(phi, a, b, eps=1e-6):
    tau = (np.sqrt(5) - 1) / 2
    x1 = b - tau * (b - a)
    x2 = a + tau * (b - a)
    f1, f2 = phi(x1), phi(x2)
    while (b - a) > eps:
        if f1 <= f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - tau * (b - a)
            f1 = phi(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + tau * (b - a)
            f2 = phi(x2)
    return (a + b) / 2

def svenn(phi, x0=0.0, delta=0.01):
    x = x0
    f0 = phi(x)
    x_prev, f_prev = x, f0
    x += delta
    f_curr = phi(x)
    if f_curr > f0:
        delta = -delta
        x, f_curr = x0 + delta, phi(x0 + delta)
    while f_curr < f_prev:
        delta *= 2.0
        x_prev, f_prev = x, f_curr
        x += delta
        f_curr = phi(x)
    return min(x_prev, x), max(x_prev, x)

def line_search(F, x, d):
    phi = lambda h: F(x + h*d)
    a, b = svenn(phi, 0.0, 0.01)
    return golden_section(phi, a, b, 1e-6)

def fletcher_reeves(F, grad, x0, eps=0.01, n_dim=2):
    x = np.array(x0, float)
    g = grad(x)
    d = -g
    k = 0
    history_x = [x.copy()]
    
    print(f"\nФлітчер-Рівс")
    print(f"{'k':>2} | {'x1':>9} | {'x2':>9} | {'F(x)':>10} | {'||g||':>8} | {'d1':>8} | {'d2':>8} | {'h':>8} | {'beta':>8}")
    print("-" * 90)
    
    while np.linalg.norm(g) > eps and k < 100:
        h = line_search(F, x, d)
        x = x + h * d
        history_x.append(x.copy())
        g_new = grad(x)
        beta = (g_new @ g_new) / (g @ g)
        if (k + 1) % n_dim == 0: beta = 0.0
        print(f"{k:>2} | {x[0]:>9.4f} | {x[1]:>9.4f} | {F(x):>10.4f} | {np.linalg.norm(g_new):>8.4f} | {d[0]:>8.4f} | {d[1]:>8.4f} | {h:>8.4f} | {beta:>8.4f}")
        d = -g_new + beta * d
        g = g_new
        k += 1
    return x, k, np.array(history_x)

def anti_gradient(F, grad, x0, eps=0.01):
    x = np.array(x0, float)
    g = grad(x)
    k = 0
    hist = [x.copy()]
    print(f"\nАнтиградієнт + точний пошук")
    print(f"{'k':>2} | {'x1':>9} | {'x2':>9} | {'F(x)':>10} | {'||g||':>8} | {'d1':>8} | {'d2':>8} | {'h':>8}")
    print("-" * 80)
    while np.linalg.norm(g) > eps and k < 100:
        d = -g
        h = line_search(F, x, d)
        print(f"{k:>2} | {x[0]:>9.4f} | {x[1]:>9.4f} | {F(x):>10.4f} | {np.linalg.norm(g):>8.4f} | {d[0]:>8.4f} | {d[1]:>8.4f} | {h:>8.4f}")
        x = x + h * d
        hist.append(x.copy())
        g = grad(x)
        k += 1
    return x, k, np.array(hist)

def sphere_metric(F, grad, x0, h0=1.0, eps=0.01):
    x = np.array(x0, float)
    Fx = F(x)
    h = h0
    k = 0
    hist = [x.copy()]
    
    #print(f"\nСферична метрика")
    #print(f"{'k':>2} | {'x1':>9} | {'x2':>9} | {'F(x)':>10} | {'||g||':>8} | {'h':>8}")
    #print("-" * 60)
    
    while True:
        g = grad(x)
        if np.linalg.norm(g) <= eps:
            break
        d = -g / np.linalg.norm(g) 
        
        while True:
            xb = x + h * d
            Fb = F(xb)
            if Fb < Fx:
                break
            h /= 2.0
            if h < 1e-10: return x, k, np.array(hist)
        
        x, Fx = xb, Fb
        hist.append(x.copy())
        k += 1
    return x, k, np.array(hist)

if __name__ == "__main__":
    x0 = np.array([5.0, 2.0])
    eps = 0.01
    
    x_cg, n_cg, hist_cg = fletcher_reeves(F, grad_F, x0, eps)
    x_ag, n_ag, hist_ag = anti_gradient(F, grad_F, x0, eps)
    x_sp, n_sp, hist_sp = sphere_metric(F, grad_F, x0, h0=1.0, eps=eps)
    
    print(f"\n[Довідково] Сферична метрика (з ЛР8) зробила {n_sp} ітерацій.")
    
    x1_vals = np.linspace(-3, 6, 100)
    x2_vals = np.linspace(-4, 4, 100)
    X1, X2 = np.meshgrid(x1_vals, x2_vals)
    Z = 2*X1**2 + X2**2 - 2*X1*X2 + 2*X1 + X2
        
    plt.figure(figsize=(10, 6))
    plt.contour(X1, X2, Z, levels=30, cmap="viridis", alpha=0.4)
    
    plt.plot(hist_sp[:, 0], hist_sp[:, 1], "r-.", label=f"Сферична метрика ({n_sp} іт.)", lw=1.5, alpha=0.7)
    plt.plot(hist_ag[:, 0], hist_ag[:, 1], "b--o", label=f"Антиградієнт + точний пошук ({n_ag} іт.)", lw=1.5, alpha=0.8)
    plt.plot(hist_cg[:, 0], hist_cg[:, 1], "go-", label=f"Флетчер-Рівс ({n_cg} іт.)", lw=2)
    
    plt.plot(x0[0], x0[1], "k*", ms=12, label="Старт (5, 2)")
    plt.plot(-1.5, -2.0, "rX", ms=10, label="Точний мінімум")
    
    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.title("Порівняння 3 траєкторій")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig("lab 9/Графік_траєкторій.png", dpi=150, bbox_inches="tight")
    plt.show()