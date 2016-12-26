import numpy as np
import matplotlib.pyplot as plt

x = []
z = []
x0 = 0
z0 = 0
for i in range(1000):
    if np.random.random()>0.5:
        x0 += 2*(np.random.random()>0.5) - 1
    else:
        z0 += 2*(np.random.random()>0.5) - 1

    if (x0, z0) in zip(x, z):
        x0 = x[-1]
        z0 = z[-1]
    else:
        x.append(x0)
        z.append(z0)
    if len(x)>=50:
        break

X, Z = np.meshgrid(np.arange(min(x)-1, max(x)+2),
            np.arange(min(z)-1, max(z)+2))
X = X.ravel(); Z = Z.ravel()
grid = zip(X, Z)
path = zip(x, z)
for i in reversed(range(len(grid))):
    if grid[i] in path:
        del grid[i]

x = np.array(x)
z = np.array(z)
for i in reversed(range(len(grid))):
    dist = np.min((grid[i][0] - x)**2 + (grid[i][1] - z)**2)
    if dist>1:
        del grid[i]
X, Z = zip(*grid)

print len(x)

s = ''
for i in range(len(x)):
    s += '%f %f %f %f\n'%(10*x[i], 0, 10*z[i], 0)
    s += '5 0 0 0\n'
    s += '0 0 5 0\n'
    s += '0 0 0 5\n'
    s += '%f %f %f\n'%tuple(np.random.random(3))

# w direction walls
for i in range(len(x)):
    s += '%f %f %f %f\n'%(10*x[i], -5, 10*z[i], 5)
    s += '5 0 0 0\n'
    s += '0 5 0 0\n'
    s += '0 0 5 0\n'
    s += '%f %f %f\n'%tuple(np.random.random(3))

    s += '%f %f %f %f\n'%(10*x[i], -5, 10*z[i], -5)
    s += '5 0 0 0\n'
    s += '0 5 0 0\n'
    s += '0 0 5 0\n'
    s += '%f %f %f\n'%tuple(np.random.random(3))

# x, z direction walls
count = 0
for j in range(len(x)):
    for i in range(len(X)):
        if abs(X[i]-x[j])+abs(Z[i]-z[j]) == 0:
            raise AssertionError("")
        if abs(X[i]-x[j])+abs(Z[i]-z[j]) == 1:
            posx = X[i] + 0.5*(x[j]-X[i])
            posy = Z[i] + 0.5*(z[j]-Z[i])
            s += '%f %f %f %f\n'%(10*posx, -5, 10*posy, 0)
            s += '0 5 0 0\n'
            s += '0 0 0 5\n'
            if abs((x[j]-X[i])) == 0:
                s += '5 0 0 0\n'
            else:
                s += '0 0 5 0\n'
            s += '%f %f %f\n'%tuple(np.random.random(3))

with open("Game4D/map2d.txt", "w") as f:
    f.write(s)

plt.scatter(X, Z)
plt.scatter(x, z, c='r')
plt.show()
