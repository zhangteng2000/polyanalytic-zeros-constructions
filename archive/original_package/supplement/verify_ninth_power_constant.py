"""Exact polynomial audit for the ninth-order unit-triple power constant.
All acceptance comparisons are rational; decimal output is for reading.
"""
from pathlib import Path
from math import isqrt
import sympy as s
import hashlib,json
if not __debug__:raise RuntimeError('Run without -O: proof assertions must remain enabled')
X,Y,R,x,v,U,V=s.symbols('X Y R x v U V')
a=Y*X-3*Y+3
b=s.expand(a*a+Y*Y*(Y-X*X))
Q=s.cancel(b*(b*b+9*b-24*a**3+18*a*b)/Y)
assert s.Poly(Q,X,Y)
p=[s.S(3),U,U*U-2*V]
for k in range(3,10):p.append(s.expand(U*p[-1]-V*p[-2]+p[-3]))
assert s.expand(p[9]-3-(p[3]**3-3*p[3]*(V**3-3*U*V+3)))==0
t=s.symbols('t')
laurent=s.expand((p[9]-3).subs({U:R*t,V:R/t}))
A=0
for term in s.Add.make_args(laurent):
    power=term.as_powers_dict().get(t,0)
    assert power%3==0
    A+=term/t**power*v**(power/3)
H=s.expand(A*A.subs(v,1/v)/(R*R))
realH=H.coeff(v,0)
for j in range(1,7):
    assert s.expand(H.coeff(v,j)-H.coeff(v,-j))==0
    realH+=2*H.coeff(v,j)*s.chebyshevt(j,x)
assert s.expand(realH-Q.subs({X:R*x,Y:R*R}))==0
F=512*Y**9-7616*Y**8+43848*Y**7-123203*Y**6+172055*Y**5-97935*Y**4-9396*Y**3+23040*Y**2-189*Y-972
G=Y**3-9*Y**2+18*Y-9
N=15207936*Y**8-197849152*Y**7+915234552*Y**6-1698080469*Y**5+594170401*Y**4+1588963305*Y**3-1428234222*Y**2+164181132*Y+22342419
relation=22654044*X+N
qx,qy=s.diff(Q,X),s.diff(Q,Y)
assert s.expand(s.resultant(qx,qy,X)-429981696*Y**8*(Y-3)**5*G**5*F)==0
assert s.groebner([qx,qy,F],X,Y,order='lex',domain=s.QQ)==s.groebner([relation,F],X,Y,order='lex',domain=s.QQ)
xp=-N/s.Integer(22654044)
assert s.rem(qx.subs(X,xp),F,Y)==0
assert s.rem(qy.subs(X,xp),F,Y)==0
upper=s.Rational(49,100)
assert s.diff(G,Y,2).subs(Y,upper)<0
assert s.diff(G,Y).subs(Y,upper)>0 and G.subs(Y,upper)<0
sturm=s.sturm(F,Y)
def variations(at):
    signs=[s.sign(h.subs(Y,at)) for h in sturm]
    signs=[z for z in signs if z]
    return sum(a!=b for a,b in zip(signs,signs[1:]))
vv=[variations(s.S(0)),variations(upper)]
assert vv[0]-vv[1]==1
def bernstein(g,lo,hi):
    q=s.Poly(s.expand(g.subs(R,lo+(hi-lo)*R)),R)
    degree=8
    cs=[s.factor(sum(q.nth(k)*s.binomial(j,k)/s.binomial(degree,k) for k in range(j+1))) for j in range(degree+1)]
    assert s.expand(sum(cs[j]*s.binomial(degree,j)*R**j*(1-R)**(degree-j) for j in range(degree+1))-q.as_expr())==0
    return cs
gp=R*(3-R)*(R**3-3*R**2+3)**2
gm=R*(3+R)*(3-3*R**2-R**3)**2
assert s.expand(Q.subs({X:R,Y:R*R})-gp*gp)==0
assert s.expand(Q.subs({X:-R,Y:R*R})-gm*gm)==0
boundary=[]
for name,g in [('plus',gp),('minus',gm)]:
    pending=[(s.S(0),s.Rational(7,10))]
    while pending:
        lo,hi=pending.pop(0);cs=bernstein(s.Rational(60,7)-g,lo,hi)
        if min(cs)<=0:
            assert hi-lo>s.Rational(1,10000)
            mid=(lo+hi)/2;pending[:0]=[(lo,mid),(mid,hi)];continue
        boundary.append({'branch':name,'interval':[str(lo),str(hi)],'coefficients':list(map(str,cs)),'minimum':str(min(cs))})
lo,hi=s.S(0),upper
assert F.subs(Y,lo)<0<F.subs(Y,hi)
for k in range(130):
    mid=(lo+hi)/2
    if F.subs(Y,mid)<0:lo=mid
    else:hi=mid
def interval(poly,variables,bounds):
    # Bounds used here may straddle zero in X; multiplication is general.
    def times(q,r):
        vals=[c*d for c in q for d in r];return min(vals),max(vals)
    def power(q,n):
        ans=(s.S(1),s.S(1))
        for k in range(n):ans=times(ans,q)
        return ans
    low=high=s.S(0)
    for powers,co in s.Poly(poly,*variables).terms():
        q=(co,co)
        for box,n in zip(bounds,powers):q=times(q,power(box,n))
        low+=q[0];high+=q[1]
    return low,high
xl,xh=interval(xp,[Y],[(lo,hi)])
assert -s.Rational(1,10)<xl<xh<0
assert xh*xh<xl*xl<lo<hi<upper
vl,vh=interval(Q,[X,Y],[(xl,xh),(lo,hi)])
assert vl>76>s.Rational(3600,49)
den=10**12
clo=s.Rational(isqrt(int(s.floor(vl*den*den))),den)
chi=s.Rational(isqrt(int(s.floor(vh*den*den)))+1,den)
assert clo*clo<=vl<vh<=chi*chi and chi-clo==s.Rational(1,den)
shortden=10**9
yls=s.floor(lo*shortden)/shortden;yhs=s.ceiling(hi*shortden)/shortden
assert F.subs(Y,yls)<0<F.subs(Y,yhs)
xls,xhs=interval(xp,[Y],[(yls,yhs)])
xls=s.floor(xls*10**7)/10**7;xhs=s.ceiling(xhs*10**7)/10**7
shortvl,shortvh=interval(Q,[X,Y],[(xls,xhs),(yls,yhs)])
assert -s.Rational(1,10)<xls<xhs<0 and xls*xls<yls and shortvl>76
out={'all_passed':True,'arithmetic':'exact rational polynomial arithmetic',
     'stationary_polynomial':str(F),'x_numerator':str(N),'x_denominator':22654044,
     'sturm_variations_0_49over100':vv,'boundary_certificates':boundary,
     'certified_y_interval':[str(lo),str(hi)],'certified_x_interval':[str(xl),str(xh)],
     'certified_C9_interval':[str(clo),str(chi)],
     'short_stationary_box':{'Y':[str(yls),str(yhs)],'X':[str(xls),str(xhs)],'Q_lower_exceeds':76},
     'C9_decimal_for_reading':[str(s.N(clo,16)),str(s.N(chi,16))],
     'equality_traces':'abs(s)^2=y_star and Re(s^3)=y_star*x_star; exactly six',
     'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
dest=Path(__file__).with_name('ninth_power_constant_verification.json')
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'all_passed':True,'C9_interval':out['C9_decimal_for_reading'],'sturm_variations':vv,'stationary_box':out['short_stationary_box'],'boundary_minima':[(r['branch'],r['interval'],r['minimum']) for r in boundary],'output':str(dest)}))
