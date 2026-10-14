"""Exact audit of C6, its algebraic maximum, and the boundary certificate.
Uses rational polynomial arithmetic. No floating comparison is accepted.
"""
from pathlib import Path
from math import isqrt
import hashlib,json
import sympy as s
if not __debug__:
    raise RuntimeError('Run without -O: proof assertions must remain enabled')
X,Y,U,V,t,e,x=s.symbols('X Y U V t e x')
ps=[s.S(3),U,U*U-2*V]
for k in range(3,7):ps.append(s.expand(U*ps[-1]-V*ps[-2]+ps[-3]))
formula=U**6-6*U**4*V+9*U**2*V**2+6*U**3-2*V**3-12*U*V
assert s.expand(ps[6]-3-formula)==0
A=t**6*e**2+(6*t**3-6*t**5)*e+9*t**4-12*t*t-2*t**3/e
H=s.expand(A*A.subs(e,1/e)/(t*t))
T=H.coeff(e,0)+sum(2*H.coeff(e,j)*s.chebyshevt(j,x) for j in range(1,4))
Q=-Y*(4*X-Y-4)*(4*X*X*Y-8*X*Y*Y+4*X*Y+12*X+Y**3+14*Y*Y-47*Y+36)
assert s.expand(T-Q.subs({X:t*x,Y:t*t}))==0
F=25*Y**3-50*Y*Y-5*Y+12
qx,qy=s.diff(Q,X),s.diff(Q,Y)
expected=110592*Y**4*(Y-4)**5*(Y-1)**2*F
assert s.expand(s.resultant(qx,qy,X)-expected)==0
relation=2*X-5*Y*Y+11*Y-4
gb=s.groebner([qx,qy,F],X,Y,order='lex',domain=s.QQ)
want=s.groebner([relation,F],X,Y,order='lex',domain=s.QQ)
assert gb==want
atX=(5*Y*Y-11*Y+4)/2
crit=36*(545*Y*Y-1420*Y+704)/125
for expr in [qx.subs(X,atX),qy.subs(X,atX),Q.subs(X,atX)-crit]:
    assert s.rem(s.together(expr),F,Y)==0
assert s.expand(27-Q.subs(Y,1)-(2*X-1)**2*(4*X+7))==0
gp=t*(3-t)*(1+t)*(2-t)**2
gm=t*(1-t)*(t+2)**2*(t+3)
assert s.expand(Q.subs({X:t,Y:t*t})-gp*gp)==0
assert s.expand(Q.subs({X:-t,Y:t*t})-gm*gm)==0
rows=[]
for name,g,a,b,minimum in [
    ('plus',gp,s.S(0),s.S(1),s.Rational(9,10)),
    ('minus',gm,s.S(0),s.Rational(1,2),s.Rational(17,32)),
    ('minus',gm,s.Rational(1,2),s.Rational(3,4),s.Rational(61,1280)),
    ('minus',gm,s.Rational(3,4),s.S(1),s.Rational(699,1024))]:
    pol=s.Poly(s.expand(6-g.subs(t,a+(b-a)*t)),t)
    bc=[sum(pol.nth(k)*s.binomial(j,k)/s.binomial(5,k) for k in range(j+1)) for j in range(6)]
    assert min(bc)==minimum and minimum>0
    recovered=sum(bc[j]*s.binomial(5,j)*t**j*(1-t)**(5-j) for j in range(6))
    assert s.expand(recovered-pol.as_expr())==0
    rows.append({'branch':name,'interval':[str(a),str(b)],'coefficients':list(map(str,bc)),'minimum':str(minimum)})
lo,hi=s.Rational(503,1000),s.Rational(63,125)
assert F.subs(Y,lo)>0>F.subs(Y,hi)
assert s.expand(s.diff(F,Y)-(-25*Y-5)-75*Y*(Y-1))==0
assert F.subs(Y,0)>0>F.subs(Y,1)
assert s.diff(atX,Y).subs(Y,hi)<0
assert -s.Rational(1,7)<atX.subs(Y,hi)<atX.subs(Y,lo)<0
assert s.Rational(1,49)<lo<hi<1
assert s.diff(crit,Y).subs(Y,hi)<0
assert crit.subs(Y,hi)-36==s.Rational(197856,390625)>0
for step in range(100):
    mid=(lo+hi)/2
    if F.subs(Y,mid)>0:lo=mid
    else:hi=mid
assert F.subs(Y,lo)>0>F.subs(Y,hi)
vlo,vhi=crit.subs(Y,hi),crit.subs(Y,lo)
den=10**12
low_int=isqrt(int(s.floor(vlo*den*den)))
high_int=isqrt(int(s.floor(vhi*den*den)))+1
clo,chi=s.Rational(low_int,den),s.Rational(high_int,den)
assert clo*clo<=vlo<vhi<=chi*chi
out={'all_passed':True,'arithmetic':'exact rational polynomial arithmetic',
     'stationary_cubic':str(F),'stationary_x':str(atX),
     'squared_optimal_constant':str(crit),
     'certified_y_interval':[str(lo),str(hi)],
     'certified_C6_interval':[str(clo),str(chi)],
     'C6_decimal_interval_for_reading':[str(s.N(clo,16)),str(s.N(chi,16))],
     'bernstein_boundary_certificates':rows,
     'equality_traces':'abs(s)^2=y_star; Re(s^3)=y_star*x_star; exactly six',
     'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
dest=Path(__file__).with_name('sixth_power_constant_verification.json')
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'all_passed':True,'C6_interval':out['C6_decimal_interval_for_reading'],'output':str(dest)}))
