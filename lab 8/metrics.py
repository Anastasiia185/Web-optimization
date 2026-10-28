import numpy as np
import matplotlib.pyplot as plt

def F(x):
    x1, x2 = x[0], x[1]
    return 2*(x1**2) + (x2**2) - 2*x1*x2 + 2*x1 + x2

def grad_F(x):
    x1, x2 = x[0], x[1]
    df_dx1 = 4*x1 - 2*x2 + 2
    df_dx2 = 2*x2 - 2*x1 + 1
    return np.array([df_dx1, df_dx2])

def check_gradients(x, h=1e-5):
    print("\nПеревірка градієнта скінченними різницями")
    g = grad_F(x)
    for i in range(2):
        e = np.zeros(2)
        e[i] = 1.0
        num_deriv = (F(x + h*e) - F(x - h*e)) / (2*h)
        print(f"dF/dx{i+1}: аналітично = {g[i]:.6f}, чисельно = {num_deriv:.6f}, різниця = {abs(g[i]-num_deriv):.2e}")

def get_direction(g, metric):
    if metric == "sphere":
        return -g / np.linalg.norm(g)
    elif metric == "cube":
        return -np.sign(g)
    elif metric == "octa":
        j = np.argmax(np.abs(g))
        d = np.zeros_like(g)
        d[j] = -np.sign(g[j])
        return d

def gradient_descent(x0, h0, eps, metric):
    x = np.array(x0, dtype=float)
    Fx = F(x)
    h = h0
    nF = 1
    ng = 0
    
    history_x = [x.copy()]
    
    for k in range(5000):
        g = grad_F(x)
        ng += 1
        
        if np.linalg.norm(g) <= eps:
            break
            
        d = get_direction(g, metric)
        
        while True:
            xb = x + h * d
            Fb = F(xb)
            nF += 1
            if Fb < Fx:
                break
            h /= 2.0
            if h < 1e-10:
                raise RuntimeError("Крок здрібнів до нуля")
                
        x, Fx = xb, Fb
        history_x.append(x.copy())
        
    return x, Fx, k+1, nF, ng, np.array(history_x)

if __name__ == "__main__":
    x0 = np.array([5.0, 2.0])
    eps = 0.01
    exact_min = np.array([-1.5, -2.0])
    
    print(f"Старт: x0 = {x0}, eps = {eps}")
    check_gradients(x0)
    
    metrics = ["sphere", "cube", "octa"]
    names = {"sphere": "Сферична", "cube": "Кубічна", "octa": "Октаедрична"}
    colors = {"sphere": "blue", "cube": "orange", "octa": "green"}
    
    for h0 in [1.0, 0.7]:
        print(f" ЗАПУСК ДЛЯ h0 = {h0}")
        
        print(f"\n{'Метрика':<15} | {'n':>3} | {'NF':>4} | {'N_grad':>6} | {'x*':>18} | {'||x*-x_точн||':>14}")
        print("-" * 75)
        
        results_hist = {}
        for metric in metrics:
            x_star, F_star, n, nF, ng, hist = gradient_descent(x0, h0, eps, metric)
            results_hist[metric] = hist
            err = np.linalg.norm(x_star - exact_min)
            x_str = f"({x_star[0]:.4f}, {x_star[1]:.4f})"
            print(f"{names[metric]:<15} | {n:>3} | {nF:>4} | {ng:>6} | {x_str:>18} | {err:>14.2e}")
            
            print(f"\nДетальна таблиця для {names[metric]} (h0={h0}):")
            print(f"{'k':>3} | {'x1':>8} | {'x2':>8} | {'F(x)':>10} | {'g1':>8} | {'g2':>8} | {'||g||':>8}")
            print("-" * 70)
            for i in range(min(6, len(hist))):
                x_c = hist[i]
                g_c = grad_F(x_c)
                print(f"{i:>3} | {x_c[0]:>8.4f} | {x_c[1]:>8.4f} | {F(x_c):>10.4f} | {g_c[0]:>8.4f} | {g_c[1]:>8.4f} | {np.linalg.norm(g_c):>8.4f}")
            if len(hist) > 6:
                print("...")
                x_c = hist[-1]
                g_c = grad_F(x_c)
                print(f"{len(hist)-1:>3} | {x_c[0]:>8.4f} | {x_c[1]:>8.4f} | {F(x_c):>10.4f} | {g_c[0]:>8.4f} | {g_c[1]:>8.4f} | {np.linalg.norm(g_c):>8.4f}")
            print()
            
        x1_vals = np.linspace(-3, 6, 100)
        x2_vals = np.linspace(-4, 3, 100)
        X1, X2 = np.meshgrid(x1_vals, x2_vals)
        Z = 2*X1**2 + X2**2 - 2*X1*X2 + 2*X1 + X2
            
        plt.figure(figsize=(10, 6))
        plt.contour(X1, X2, Z, levels=30, cmap="viridis", alpha=0.5)
        
        for metric in metrics:
            hist = results_hist[metric]
            plt.plot(hist[:, 0], hist[:, 1], marker="o", ms=4, label=names[metric], color=colors[metric])
            
        plt.plot(x0[0], x0[1], "k*", ms=12, label="Старт")
        plt.plot(exact_min[0], exact_min[1], "rX", ms=10, label="Точний мінімум")
        
        plt.xlabel("x1")
        plt.ylabel("x2")
        plt.title(f"Траєкторії методів (h0 = {h0})")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.savefig(f"lab 8/Графік_траєкторій_h{h0}.png", dpi=150, bbox_inches="tight")
        plt.close()