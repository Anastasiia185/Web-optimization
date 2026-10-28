import math
import matplotlib.pyplot as plt

def f(x):
    return math.exp(x - 1.0) - 1.0 - x - (x**2)/2.0 - (x**3)/6.0

def df(x):
    return math.exp(x - 1.0) - 1.0 - x - (x**2)/2.0

def secant_method(a0, b0, eps=0.001):
    print(f"Старт: x0 = {a0:.4f}, x1 = {b0:.4f}, eps = {eps}")
    
    xm = a0  
    xk = b0  
    k = 0
    
    pm = df(xm)
    pk = df(xk)
    
    print(f"f'(x0) = {pm:.6f}, f'(x1) = {pk:.6f}")
    print(f"\n{'k':>2} | {'x_k-1':>9} | {'f\"(x_k-1)':>10} | {'x_k':>9} | {'f\"(x_k)':>9} | {'x_k+1':>9} | {'f\"(x_k+1)':>10} | {'|x_k+1 - x*|':>12} | Дія")
    print("-" * 105)
    
    x_prev_vals = []
    
    for k in range(1, 50):
        xn = (pk * xm - pm * xk) / (pk - pm)
        
        pn = df(xn)
        
        x_prev_vals.append((xm, pm))
        x_prev_vals.append((xk, pk))
        
        exact_min = 3.19944
        err = abs(xn - exact_min)

        if abs(pn) <= eps:
            print(f"{k:>2} | {xm:>9.6f} | {pm:>10.6f} | {xk:>9.6f} | {pk:>9.6f} | {xn:>9.6f} | {pn:>10.6f} | {err:>12.1e} | |f'| <= eps -> СТОП")
            break
        if pn > 0:
            print(f"{k:>2} | {xm:>9.6f} | {pm:>10.6f} | {xk:>9.6f} | {pk:>9.6f} | {xn:>9.6f} | {pn:>10.6f} | {err:>12.1e} | x_k <- x_k+1")
            xk = xn
            pk = pn
        else:
            print(f"{k:>2} | {xm:>9.6f} | {pm:>10.6f} | {xk:>9.6f} | {pk:>9.6f} | {xn:>9.6f} | {pn:>10.6f} | {err:>12.1e} | x_k-1 <- x_k+1")
            xm = xn
            pm = pn
            
            
    print(f"\nРезультат:")
    print(f"x* = {xn:.6f}")
    print(f"f(x*) = {f(xn):.6f}")
    print(f"n (ітерацій) = {k}")
    print(f"N_f' = {k + 2} (2 на старті + по 1 на ітерацію)")

    plt.figure(figsize=(10, 6))
    x_plot = [2.5 + i*0.01 for i in range(100)]
    y_plot = [df(x) for x in x_plot]
    plt.plot(x_plot, y_plot, label="f'(x)", color="blue")
    plt.axhline(0, color="black", lw=1)
    
    colors = ["orange", "green"]
    for i in range(2):
        if i * 2 + 1 < len(x_prev_vals):
            p1 = x_prev_vals[i*2]
            p2 = x_prev_vals[i*2 + 1]
            slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
            x_line = [2.6, 3.3]
            y_line = [p1[1] + slope * (x - p1[0]) for x in x_line]
            plt.plot(x_line, y_line, linestyle="--", color=colors[i], label=f"Січна {i+1}")
            plt.scatter([p1[0], p2[0]], [p1[1], p2[1]], color=colors[i], zorder=5)

    plt.scatter([3.19944], [0], color="red", zorder=6, label="Корінь (x*)")
    plt.xlim(2.6, 3.35)
    plt.ylim(-2, 5)
    plt.xlabel("x")
    plt.ylabel("f'(x)")
    plt.title("Метод січних: графік f'(x) та перші дві січні")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig("lab 7/Графік_січних.png", dpi=150, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    secant_method(2.70, 3.20, 0.001)