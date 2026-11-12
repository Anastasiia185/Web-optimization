import numpy as np
import matplotlib.pyplot as plt

def F(x):
    x1, x2 = x[0], x[1]
    return 2*(x1**2) + (x2**2) - 2*x1*x2 + 2*x1 + x2

nF = 0
def Feval(x):
    global nF
    nF += 1
    return F(x)

def explore(y, Fy, delta):
    global nF
    z = y.copy()
    Fz = Fy
    for i in range(len(z)):
        t_plus = z.copy()
        t_plus[i] += delta
        F_plus = Feval(t_plus)
        if F_plus < Fz:
            z, Fz = t_plus, F_plus
            continue
        
        t_minus = z.copy()
        t_minus[i] -= delta
        F_minus = Feval(t_minus)
        if F_minus < Fz:
            z, Fz = t_minus, F_minus
            
    return z, Fz

def hooke_jeeves(x0, delta0, eps):
    global nF
    
    base = np.array(x0, dtype=float)
    Fb = Feval(base)
    y = base.copy()
    Fy = Fb
    delta = delta0
    k = 0
    from_pattern = False
    
    history_base = [base.copy()]
    
    print(f"{'k':>2} | {'Дія':<25} | {'Звідки':>12} | {'F':>6} | {'Delta':>6} | {'Куди':>12} | {'F':>6} | {'проб':>4} | {'N_F':>4}")
    print("-" * 110)
    
    while True:
        old_nF = nF
        z, Fz = explore(y, Fy, delta)
        probs = nF - old_nF
        
        if Fz < Fb:
            new_y = z + (z - base)
            F_new_y = Feval(new_y)
            
            print(f"{k:>2} | {'нова база -> зразок':<25} | ({base[0]:>4.1f};{base[1]:>4.1f}) | {Fb:>6.2f} | {delta:>6.2f} | ({new_y[0]:>4.1f};{new_y[1]:>4.1f}) | {F_new_y:>6.2f} | {probs:>4} | {nF:>4}")
            
            base, Fb = z, Fz
            y, Fy = new_y, F_new_y
            from_pattern = True
            k += 1
            history_base.append(base.copy())
            
        elif from_pattern:
            print(f"{k:>2} | {'повернення в базу':<25} | ({y[0]:>4.1f};{y[1]:>4.1f}) | {Fy:>6.2f} | {delta:>6.2f} | {'--':>12} | {'--':>6} | {probs:>4} | {nF:>4}")
            y, Fy = base.copy(), Fb
            from_pattern = False
            
        else:
            delta /= 2.0
            print(f"{k:>2} | {'дроблення Delta':<25} | ({base[0]:>4.1f};{base[1]:>4.1f}) | {Fb:>6.2f} | {delta:>6.2f} | {'--':>12} | {'--':>6} | {probs:>4} | {nF:>4}")
            if delta < eps:
                print("-" * 110)
                print(f"СТОП: Delta = {delta:.4f} < eps = {eps}")
                break
            y, Fy = base.copy(), Fb
            
    return base, Fb, np.array(history_base), nF

if __name__ == "__main__":
    x0 = np.array([5.0, 2.0])
    eps = 0.01
    
    # Запуск 1: Delta0 = 1.0 
    print("="*110)
    print(" ЗАПУСК 1: Delta0 = 1.0")
    print("="*110)
    x_opt, F_opt, hist1, nF1 = hooke_jeeves(x0, delta0=1.0, eps=eps)
    print(f"\nРезультат: x* = ({x_opt[0]:.4f}, {x_opt[1]:.4f}), F* = {F_opt:.6f}")
    print(f"Всього обчислень F: {nF1}")

    # Запуск 2: Delta0 = 0.7
    print("\n" + "="*110)
    print(" ЗАПУСК 2: Delta0 = 0.7")
    print("="*110)
    x_opt2, F_opt2, hist2, nF2 = hooke_jeeves(x0, delta0=0.7, eps=eps)
    print(f"\nРезультат: x* = ({x_opt2[0]:.4f}, {x_opt2[1]:.4f}), F* = {F_opt2:.6f}")
    print(f"Всього обчислень F: {nF2}")

    x1_vals = np.linspace(-3, 6, 100)
    x2_vals = np.linspace(-4, 3, 100)
    X1, X2 = np.meshgrid(x1_vals, x2_vals)
    Z = 2*X1**2 + X2**2 - 2*X1*X2 + 2*X1 + X2
        
    plt.figure(figsize=(10, 6))
    plt.contour(X1, X2, Z, levels=30, cmap="viridis", alpha=0.5)
    plt.plot(hist1[:, 0], hist1[:, 1], "ro-", label=f"Хук-Дживс (Delta0=1.0) - {len(hist1)-1} баз")
    plt.plot(hist2[:, 0], hist2[:, 1], "bs--", label=f"Хук-Дживс (Delta0=0.7) - {len(hist2)-1} баз", alpha=0.7)
    plt.plot(x0[0], x0[1], "k*", ms=12, label="Старт (5, 2)")
    plt.plot(-1.5, -2.0, "rX", ms=10, label="Точний мінімум")
    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.title("Траєкторія методу Хука-Дживса")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig("lab 10/trajectories.png", dpi=150, bbox_inches="tight")
    plt.show()