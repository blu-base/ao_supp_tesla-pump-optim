#!/usr/bin/env python3
import numpy as np
import pandas as pd

#import chart_studio.plotly as plt

import matplotlib.cm as cm

import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
from mpl_toolkits.mplot3d import Axes3D

from scipy.spatial import Delaunay
from collections import defaultdict


def alpha_shape_3D(pos, alpha):
    """
    Compute the alpha shape (concave hull) of a set of 3D points.
    Parameters:
        pos - np.array of shape (n,3) points.
        alpha - alpha value.
    return
        outer surface vertex indices, edge indices, and triangle indices
    """

    tetra = Delaunay(pos)
    # Find radius of the circumsphere.
    # By definition, radius of the sphere fitting inside the tetrahedral needs 
    # to be smaller than alpha value
    # http://mathworld.wolfram.com/Circumsphere.html
    tetrapos = np.take(pos,tetra.vertices,axis=0)
    normsq = np.sum(tetrapos**2,axis=2)[:,:,None]
    ones = np.ones((tetrapos.shape[0],tetrapos.shape[1],1))
    a = np.linalg.det(np.concatenate((tetrapos,ones),axis=2))
    Dx = np.linalg.det(np.concatenate((normsq,tetrapos[:,:,[1,2]],ones),axis=2))
    Dy = -np.linalg.det(np.concatenate((normsq,tetrapos[:,:,[0,2]],ones),axis=2))
    Dz = np.linalg.det(np.concatenate((normsq,tetrapos[:,:,[0,1]],ones),axis=2))
    c = np.linalg.det(np.concatenate((normsq,tetrapos),axis=2))
    r = np.sqrt(Dx**2+Dy**2+Dz**2-4*a*c)/(2*np.abs(a))

    # Find tetrahedrals
    tetras = tetra.vertices[r<alpha,:]
    # triangles
    TriComb = np.array([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
    Triangles = tetras[:,TriComb].reshape(-1,3)
    Triangles = np.sort(Triangles,axis=1)
    # Remove triangles that occurs twice, because they are within shapes
    TrianglesDict = defaultdict(int)
    for tri in Triangles:TrianglesDict[tuple(tri)] += 1
    Triangles=np.array([tri for tri in TrianglesDict if TrianglesDict[tri] ==1])
    #edges
    EdgeComb=np.array([(0, 1), (0, 2), (1, 2)])
    Edges=Triangles[:,EdgeComb].reshape(-1,2)
    Edges=np.sort(Edges,axis=1)
    Edges=np.unique(Edges,axis=0)

    Vertices = np.unique(Edges)
    return Vertices,Edges,Triangles


def normalVectors(points, faces):
    norm = np.zeros(points.shape)
    tris = points[faces] 

    n = np.cross( tris[::,1] - tris[::,0], tris[::,2] - tris[::,0])
    
    return n


gen = pd.read_csv('genAll.csv')

pts = gen[[ 'ns', 'specificDiaBalje', 'efficiency']].to_numpy()

vertices,edges,triangles = alpha_shape_3D(pts, 200)

norms = normalVectors(pts,triangles)
coveringTris = []
for i in range(0, len(triangles)):
    if norms[i][2] < 0:
        coveringTris.append(triangles[i].tolist())


coveringTris = np.array(coveringTris)
coveringVerts = []

for triangle in coveringTris:
   
    missing = set(triangle.tolist())-set(coveringVerts)
    coveringVerts.extend(missing)

coveringVerts.sort()
coveringVerts = np.array(coveringVerts)
coveringPts = pts[coveringVerts]


X = coveringPts[:,0]
Y = coveringPts[:,1]
Z = coveringPts[:,2]

triangulation = Triangulation(pts[:,0],pts[:,1], triangles=triangles)
levels = np.arange(0.4,0.6, 0.05)

cmap = cm.get_cmap(name='gnuplot2_r', lut=None)

fig,ax = plt.subplots()
CS = ax.tricontour(triangulation, pts[:,2], levels=levels, colors='grey', linewidths=[.5,.5,1.])

strs = ['0.4', '0.45', '0.5', '0.55', '0.6']
fmt = {}
for l,s in zip(CS.levels, strs):
    fmt[l] = s
    
ax.clabel(CS, CS.levels[::], inline=True, fmt=fmt)

#ax.set_xlim(.1,.4)
#ax.set_ylim(10,40)


#ax.set_xscale('log')
#ax.set_yscale('log')

ax.set_ylabel('specific diameter')
ax.set_xlabel('specific speed ns')

ax.set_title('Cordier diagram')

ax.xaxis.grid(color='grey', alpha=0.5, zorder=-2, linestyle='dashed')
ax.yaxis.grid(color='grey', alpha=0.5, zorder=-2, linestyle='dashed')
plt.grid(True)
cbar = ax.scatter(pts[:,0], pts[:,1], c=pts[:,2], zorder=-1, s=1, cmap=cmap, marker='o')

cbarobj = fig.colorbar(cbar, ax=ax)
cbarobj.set_label('hydr. efficiency ')
#plt.show()

#plt.rcParams['svg.fonttype'] = 'svgfont'

#plt.savefig('balje.svg')

plt.savefig('balje.pdf')
