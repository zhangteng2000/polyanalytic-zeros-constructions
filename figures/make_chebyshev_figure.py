"""Generate the scientific figure and its embedded TikZ source."""
from pathlib import Path
import csv
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
n=5
epsilon=1/n**2
def phi(t):
    return math.sinh(t)/math.cosh(n*t)
def bisect(fun,a,b):
    fa=fun(a)
    assert fa*fun(b)<0
    for _ in range(100):
        c=(a+b)/2
        if fa*fun(c)<=0:b=c
        else:a=c;fa=fun(a)
    return (a+b)/2
peak=bisect(lambda t:1/math.tanh(t)-n*math.tanh(n*t),.001,1)
roots=[(math.cos((2*k-1)*math.pi/(2*n)),0.0,"real",k)
       for k in range(1,n+1)]
times={}
for j in range(1,n):
    theta=j*math.pi/n
    fun=lambda t:2*math.sin(theta)*phi(t)-epsilon
    t1=bisect(fun,0,peak)
    t2=bisect(fun,peak,2)
    times[j]=(t1,t2)
    side=1 if j%2==0 else -1
    for t in (t1,t2):
        roots.append((math.cos(theta)*math.cosh(t),
                      side*math.sin(theta)*math.sinh(t),"nonreal",j))
with (HERE/"chebyshev_n5_roots.csv").open("w",newline="",encoding="utf-8") as stream:
    writer=csv.writer(stream)
    writer.writerow(["real_part","imaginary_part","type","index"])
    writer.writerows(roots)

# Matplotlib publication figure, independent of the LaTeX implementation.
plt.rcParams.update({"font.family":"serif","font.size":10,
                     "axes.spines.top":False,"axes.spines.right":False})
fig,(ax,bx)=plt.subplots(1,2,figsize=(9.0,3.5),layout="constrained")
ts=np.linspace(0,1.18,150)
for j in range(1,n):
    theta=j*math.pi/n
    side=1 if j%2==0 else -1
    x=np.cos(theta)*np.cosh(ts)
    y=np.sin(theta)*np.sinh(ts)
    ax.plot(x,side*y,color="#285a84",lw=1)
    ax.plot(x,-side*y,color="0.78",lw=.8,ls="--")
ax.axhline(0,color=".5",lw=.7)
ax.axvline(0,color=".5",lw=.7)
ax.scatter([x for x,y,kind,j in roots if kind=="real"],
           [y for x,y,kind,j in roots if kind=="real"],s=22,c="black",zorder=3)
ax.scatter([x for x,y,kind,j in roots if kind=="nonreal"],
           [y for x,y,kind,j in roots if kind=="nonreal"],
           s=32,facecolor="white",edgecolor="#ae452e",lw=1.2,zorder=4)
ax.set(xlim=(-1.4,1.4),ylim=(-1.25,1.25),
       xlabel=r"$\mathrm{Re}\,z$",ylabel=r"$\mathrm{Im}\,z$")
ax.set_aspect("equal")
ax.set_title(r"$n=5,\quad \varepsilon=1/25:\quad 13$ zeros")
t=np.linspace(0,1.2,300)
A=2*math.sin(math.pi/n)*np.sinh(t)/np.cosh(n*t)
bx.plot(t,A,color="#285a84")
bx.axhline(epsilon,color="#ae452e",ls="--",lw=1)
bx.scatter(times[1],[epsilon,epsilon],s=28,facecolor="white",edgecolor="#ae452e",zorder=4)
bx.set(xlim=(0,1.2),ylim=(0,.18),xlabel=r"$t$",
       ylabel=r"$2\sin(\pi/5)\,\phi_5(t)$")
bx.text(.74,epsilon+.006,r"$\varepsilon=1/25$",color="#ae452e")
bx.set_title("Two roots on a selected branch")
fig.savefig(HERE/"chebyshev_n5.pdf")
fig.savefig(HERE/"chebyshev_n5.png",dpi=200)
plt.close(fig)

# Inline source keeps the manuscript a standalone TeX file.
lines=[r"\begin{figure}[htbp]",r"\centering",
       r"\begin{tikzpicture}[font=\small]",
       r"\definecolor{curveblue}{RGB}{40,90,132}",
       r"\definecolor{rootred}{RGB}{174,69,46}",
       r"\begin{scope}[x=1.8cm,y=1.8cm]",
       r"\clip (-1.52,-1.48) rectangle (1.64,1.55);",
       r"\draw[gray!65,->] (-1.43,0)--(1.43,0) node[below] {\(\Rea z\)};",
       r"\draw[gray!65,->] (0,-1.25)--(0,1.28) node[right] {\(\Ima z\)};",
       r"\node at (0,1.47) {\(n=5,\ \eps=1/25\): thirteen zeros};"]
for j in range(1,n):
    degree=180*j/n
    side=1 if j%2==0 else -1
    for sign,color in ((side,"curveblue"),(-side,"gray!45,densely dashed")):
        lines.append(
            r"\draw["+color+r",thin,domain=0:1.18,samples=90,variable=\t] plot "
            +"({cos("+str(degree)+r")*(exp(\t)+exp(-\t))/2},"
            +"{"+str(sign)+"*sin("+str(degree)+r")*(exp(\t)-exp(-\t))/2});")
for x,y,kind,j in roots:
    style="black" if kind=="real" else "draw=rootred,fill=white,line width=.7pt"
    lines.append(r"\fill["+style+f"] ({x:.9f},{y:.9f}) circle (1.65pt);")
for x in (-1,1):
    lines.append(rf"\node[below] at ({x},-.025) {{\({x}\)}};")
lines.extend([r"\end{scope}",
              r"\begin{scope}[shift={(4.4cm,-2.02cm)},x=4.4cm,y=22cm]",
              r"\draw[gray!65,->] (0,0)--(1.21,0) node[below] {\(t\)};",
              r"\draw[gray!65,->] (0,0)--(0,.177);",
              r"\node[anchor=south west] at (0,.18) {\(2\sin(\pi/5)\,\phi_5(t)\)};",
              r"\draw[curveblue,domain=0:1.15,samples=140,variable=\t] plot "
              r"(\t,{2*sin(36)*(exp(\t)-exp(-\t))/(exp(5*\t)+exp(-5*\t))});",
              r"\draw[rootred,dashed] (0,.04)--(1.15,.04);",
              r"\node[rootred,above] at (.92,.04) {\(\eps=1/25\)};"])
for t in times[1]:
    lines.append(rf"\fill[draw=rootred,fill=white,line width=.7pt] ({t:.9f},.04) circle (1.65pt);")
lines.extend([r"\node[below] at (0,0) {\(0\)};",
              r"\node[below] at (1,0) {\(1\)};",
              r"\end{scope}",r"\end{tikzpicture}",
              r"\caption{The degree-five example. Filled dots are the five real zeros;",
              r"open circles are the eight nonreal zeros. Solid curves are the selected",
              r"branches of \(\Ima T_5=0\); their conjugates are dashed.",
              r"The right panel shows the two intersections that occur on each",
              r"selected branch, here for \(j=1\) and \(j=4\).}",
              r"\label{fig:chebyshev}",r"\end{figure}"])
(HERE/"chebyshev_n5.tikz").write_text("\n".join(lines)+"\n",encoding="utf-8")
print("Generated figure, inline TikZ, and 13 root coordinates.")
