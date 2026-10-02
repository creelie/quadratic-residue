"""Third-order coefficients of the expansions in Theorems 1.3 and 1.4 (Section 5.3).
E(x) = C x (log x)^(-1/2) (1 + kappa/L + kappa3/L^2 + ...), and similarly for sum rho(n) with
kappa2, kappa2_3. These are used only to interpret the residual columns of Tables 2 and 3.
Run: python3 code/third_order.py   (about 5 minutes)"""
import mpmath as mp, math
mp.mp.dps=30
chi4=[0,1,0,-1]
def R(s): return (1-mp.power(2,-s))*mp.zeta(s)/mp.dirichlet(s,chi4)
def logP3(s):
    tot,k=mp.mpf(0),1
    while True:
        t=mp.log(R(s*2**(k-1)))/2**k; tot+=t
        if abs(t)<mp.mpf(10)**-28: return tot
        k+=1
def s1z(s):
    h=s-1
    if abs(h)<mp.mpf('1e-3'):
        return 1+sum((-1)**n*mp.stieltjes(n)*h**(n+1)/mp.factorial(n) for n in range(12))
    return h*mp.zeta(s)
# E: Z(s) = sqrt(s1z * (1-2^-s) P(2s)/L) / s
def ZE(s): return mp.sqrt(s1z(s)*(1-mp.power(2,-s))*mp.e**logP3(2*s)/mp.dirichlet(s,chi4))/s
Z0=ZE(mp.mpf(1)); Z1=mp.diff(ZE,1); Z2=mp.diff(ZE,1,2)/2
l0=Z0/mp.gamma(0.5); l1=Z1/mp.gamma(-0.5); l2=Z2/mp.gamma(-1.5)
print('E lambdas',l0,l1,l2)
# A(y) = y (log y)^{-1/2}(l0+l1/Ly+l2/Ly^2); y=x/2, Ly = L - log2 =: L(1-a/L), a=log2
a=mp.log(2)
# (Ly)^{-1/2} = L^{-1/2}(1 + a/(2L) + 3a^2/(8L^2)); 1/Ly = 1/L(1+a/L); 1/Ly^2 = 1/L^2
c0=l0; c1=l1+l0*a/2; c2=l2+l1*a+l1*a/2+l0*3*a*a/8
print('E: C=',c0/2,' kappa=',c1/c0,' kappa3=',c2/c0)

# rho: G(s)=prod_p (1-((1+1/p)/2)p^-s)(1-p^-s)^-1/2 ; log G derivatives via prime zeta expansions
# log E_p(s) = log(1 - u/2 - u/(2p)) - (1/2)log(1-u), u=p^-s. Expand in powers of p^-k: write
# log(1-(u/2)(1+1/p)) = -sum_{j>=1} (1/j)(1/2)^j u^j (1+1/p)^j ; u^j (1+1/p)^j = sum_i C(j,i) p^{-js - i}
def logG(s):
    tot=mp.mpf(0)
    # term: -sum_j (1/j) 2^-j sum_i C(j,i) P(js+i)  +  (1/2) sum_j (1/j) P(js)
    tot=-mp.primezeta(s+1)/2
    for j in range(2,90):
        t=mp.mpf(0)
        for i in range(0,j+1):
            t+=mp.binomial(j,i)*mp.primezeta(j*s+i)
        term=-(mp.power(2,-j)/j)*t + mp.primezeta(j*s)/(2*j)
        tot+=term
        if j>3 and abs(term)<mp.mpf(10)**-22: break
    return tot
mp.mp.dps=25
G1=mp.e**logG(1)
d1=mp.diff(logG,1); d2=mp.diff(logG,1,2)
print('G(1)=',G1,' dlogG=',d1)
def logZr(s): return mp.log(s1z(s))/2 + logG(s) - mp.log(s)
Z0=G1; Zd1=mp.diff(logZr,1); Zd2=mp.diff(logZr,1,2)
Z1=Z0*Zd1; Z2=Z0*(Zd2+Zd1**2)/2
l0=Z0/mp.gamma(0.5); l1=Z1/mp.gamma(-0.5); l2=Z2/mp.gamma(-1.5)
print('rho lambdas',l0,l1,l2)
print('C2=',l0/2,' kappa2=',(l1-l0/4)/l0,' kappa2_3=',(l2-3*l1/4-3*l0/16)/l0)
