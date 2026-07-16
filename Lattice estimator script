from estimator import *

# secp256k1 elliptic curve
q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
l = ceil(log(q, 2))
E = EllipticCurve(GF(q), [0, 7])
print(f"{E}\n")
r = E.cardinality()
assert r.is_prime()
assert l == ceil(log(r, 2))

lam = 0
n = 0
while lam < 128:
    n = n + 1
    mpr = n*(l + 1)
    m = ceil(1.5*mpr)
    c = RR(m/mpr)    
    d = m - mpr
    print(f"n={n}, m={m}, mpr={mpr}, c={c}, d={d}")
    
    params = SIS.Parameters(n=n, m=m+1, length_bound=1, q=q, norm=oo)
    est = SIS.estimate(params)
    rop = est['lattice']['rop']
    lam = log(rop, 2)
    print()
    
print(f"n*(m - 1)={n*(m - 1)}, m/2={RR(m/2)}, n*(m/2 - 1)={RR(n*(m/2 - 1))}, l={l}, 2*(l + l/2)={2*(l + l/2)}") 
