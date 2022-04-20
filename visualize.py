import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ArrowStyle


def readFile(time):
    xs = []
    thisx = []
    ys = []
    thisy = []
    vxs = []
    thisvx = []
    vys = []
    thisvy = []
    ps = []
    thisp = []
    curx = -1
    with open("data "+time,"r") as f:
        while line := f.readline():
            dat = line[:-2].split(",")
            datF = [float(d) for d in dat]
            if datF[0] != curx and curx != -1:
                xs.append(thisx)
                ys.append(thisy)
                vxs.append(thisvx)
                vys.append(thisvy)
                ps.append(thisp)
                thisx = []
                thisy = []
                thisvx = []
                thisvy = []
                thisp = []
            thisx.append(datF[0])
            thisy.append(datF[1])
            thisvx.append(datF[2])
            thisvy.append(datF[3])
            thisp.append(datF[4])
            curx = datF[0]
    return np.array(xs), np.array(ys), np.array(vxs), np.array(vys), np.array(ps)

xs, ys, vxs, vys, ps = readFile("5.6800e+00")
xReg = np.linspace(xs[0,0],xs[-1,-1],len(xs))
yReg = np.linspace(ys[0,0],ys[-1,-1],len(ys[0]))
xReg2 = np.linspace(xs[0,0],xs[-1,-1],len(xs)*4)
yReg2 = np.linspace(ys[0,0],ys[-1,-1],len(xs)*4)
print(xs)
xs,ys = np.meshgrid(xReg,yReg)
print(xs)
vs = np.sqrt(vxs*vxs + vys*vys)
fig, (ax1) = plt.subplots(nrows=1)
#ax1.contour(xs, ys, vs.T, levels=40, linewidths=0.5, colors='k')
cntr1 = ax1.contourf(xs, ys, vs.T, levels=40)
fig.colorbar(cntr1,ax=ax1)
lw = np.power(vs.T/np.max(vs),0.3)*0.5
seeds = np.zeros((len(xReg2),2))
seeds[:,0] = xReg2
seeds[:,1] = yReg2
print(seeds)
ax1.streamplot(xs,ys,vxs.T,vys.T,color='r',linewidth=lw, density=(4,2), arrowstyle="-",maxlength=106)
plt.show()

#print(np.meshgrid(np.arange(7),np.arange(10)))
