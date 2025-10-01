import cmath

def f(x):
    return x**3 - 13*x - 12

def muller(f, x0, x1, x2, tol=1e-6, max_iter=100):
    for iteration in range(max_iter):
        h0 = x1 - x0
        h1 = x2 - x1
        d0 = (f(x1) - f(x0)) / h0
        d1 = (f(x2) - f(x1)) / h1
        a = (d1 - d0) / (h1 + h0)
        b = a * h1 + d1
        c = f(x2)

        # Compute discriminant
        rad = cmath.sqrt(b*b - 4*a*c)

        # Denominator with larger magnitude
        if abs(b + rad) > abs(b - rad):
            denom = b + rad
        else:
            denom = b - rad

        dx = -2 * c / denom
        x3 = x2 + dx

        # Convergence check
        if abs(dx) < tol:
            return x3, iteration + 1

        # Update guesses
        x0, x1, x2 = x1, x2, x3

    raise ValueError("Muller's method did not converge")

# Initial guesses
x0, x1, x2 = 4.5, 5.5, 5.0
root, iters = muller(f, x0, x1, x2)

print(f"Root found: {root}")
print(f"Iterations: {iters}")
