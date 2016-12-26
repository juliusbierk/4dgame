import numpy as np
import matplotlib.pyplot as plt

col1 = np.array([0.6,0.1,0.8])
col2 = np.array([0.3,0.6,0.1])
size = 75
for _ in range(100):
    x0 = 0
    z0 = 0
    w0 = 0
    x = [w0]
    z = [z0]
    w = [w0]
    col = [col1]
    for i in range(500):
        r = np.random.random()
        if r>0.5:
            x0 += 2*(np.random.random()>0.5) - 1
        elif r>0.2:
            z0 += 2*(np.random.random()>0.5) - 1
        else:
            w0 += 2*(np.random.random()>0.5) - 1

        if (x0, z0, w0) in zip(x, z, w):
            x0 = x[-1]
            z0 = z[-1]
            w0 = x[-1]
        else:
            x.append(x0)
            z.append(z0)
            w.append(w0)
            col.append((col1*(size-len(x)-1) + col2*len(x))/size)
        if len(x)>=size:
            break
    if len(x)>=size:
        break

print len(x)

X, Z, W = np.meshgrid(np.arange(min(x)-1, max(x)+2),
            np.arange(min(z)-1, max(z)+2),
            np.arange(min(w)-1, max(w)+2))

X = X.ravel(); Z = Z.ravel(); W = W.ravel()
grid = zip(X, Z, W)
path = zip(x, z, w)
for i in reversed(range(len(grid))):
    if grid[i] in path:
        del grid[i]

x = np.array(x)
z = np.array(z)
w = np.array(w)
for i in reversed(range(len(grid))):
    dist = np.min((grid[i][0] - x)**2 + (grid[i][1] - z)**2 +
                  (grid[i][2] - w)**2)
    if dist>1:
        del grid[i]
X, Z, W = zip(*grid)

s = ''
for i in range(len(x)):
    s += '%f %f %f %f\n'%(10*x[i], 0, 10*z[i], 10*w[i])
    s += '5 0 0 0\n'
    s += '0 0 5 0\n'
    s += '0 0 0 5\n'
    s += '%f %f %f\n'%tuple(col[i])

def variate_color(col, x):
    if x>0:
        return col*(1-x) + x*np.array([1,1,1])
    else:
        return col*(1+x) - x*np.array([0,0,0])

# x, z direction walls
count = 0
for j in range(len(x)):
    for i in range(len(X)):
        if abs(X[i]-x[j])+abs(Z[i]-z[j])+abs(W[i]-w[j]) == 0:
            raise AssertionError("")
        if abs(X[i]-x[j])+abs(Z[i]-z[j])+abs(W[i]-w[j]) == 1:
            posx = X[i] + 0.5*(x[j]-X[i])
            posz = Z[i] + 0.5*(z[j]-Z[i])
            posw = W[i] + 0.5*(w[j]-W[i])
            s += '%f %f %f %f\n'%(10*posx, -5, 10*posz, 10*posw)
            s += '0 5 0 0\n'
            if abs(x[j]-X[i]) == 0:
                s += '5 0 0 0\n'
            else:
                thiscol = '%f %f %f\n'%tuple(variate_color(col[j], 0.1))
            if abs(z[j]-Z[i]) == 0:
                s += '0 0 5 0\n'
            else:
                thiscol = '%f %f %f\n'%tuple(variate_color(col[j], -0.1))
            if abs(w[j]-W[i]) == 0:
                s += '0 0 0 5\n'
            else:
                thiscol = '%f %f %f\n'%tuple(variate_color(col[j], 0.2))
            s += thiscol

with open("Game4D/map3d.txt", "w") as f:
    f.write(s)

plt.scatter(X, Z)
plt.scatter(x, z, c='r')
plt.show()
