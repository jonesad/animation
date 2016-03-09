# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:54:56 2016

Animation test
make some scripts for running and testing animations

@author: Adam
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def regPolyVert(nSides, fCircRad=1.0, npaTransVect=[0,0], fOrient=0, coord='cart'):
    '''
        Return a numpy array with the cartesian coordinates of the vertices of
        the specified polygon.    
    '''
    fDelta=2.*np.pi/float(nSides)
    npaOut = np.zeros([nSides + 1,2])
#    first calculate cartesian coordinates to handle translations.
    for nVert in range(nSides):
        npaOut[nVert, 0] = fCircRad*np.cos(float(nVert)*fDelta+ fOrient) + npaTransVect[0]
        npaOut[nVert, 1] = fCircRad*np.sin(float(nVert)*fDelta+ fOrient) + npaTransVect[1]
    npaOut[-1,:]=npaOut[0,:] #get back to the beginning with no rounding error
    if coord == 'cart':
        return npaOut[:, 0], npaOut[:, 1]
    if coord == 'polar':
        npaTheta = np.arctan2(npaOut[:,1], npaOut[:,0])
        npaTheta
        npaR = np.sqrt(np.power(npaOut[:,1], 2), np.power(npaOut[:,0], 2))
        return npaR, npaTheta

def regPolyPol(nSides, fCR=1., fOrient=0.):
    fDelta=2.*np.pi/float(nSides)
    npaOut = np.zeros([nSides + 1,2])
#    first calculate cartesian coordinates to handle translations.
    for nVert in range(nSides):
        npaOut[nVert, 0] = fCR
        npaOut[nVert, 1] = float(nVert)*fDelta+ fOrient
    npaOut[-1,1]=2.*np.pi #get back to the beginning with no rounding error and distinct theta
    return npaOut[:, 0], npaOut[:, 1]
    
def linterpSeq(npaX, npaY, fRes=4):
    '''
        Return a list of points on the line interpolating the points given.
        Returns res*(len(x)-1) points
    '''
    if len(x) != len(y):
        print 'Error: len(x) != len(y)'
        return None
    if fRes == 0:
        fRes = 2*len(x)
    from scipy.interpolate import interp1d
    lXOut=[]
    lYOut=[]
    for nVert in range(len(x) - 1):
        f=interp1d(x[nVert:nVert + 2], y[nVert:nVert + 2])
        newx = np.linspace(x[nVert], x[nVert + 1], fRes)
        newy = f(newx)
        lXOut.extend(newx)
        lYOut.extend(newy)
    return lXOut, lYOut

def lsClose(npaX1, npaX2, tol=10.**(-4)):
    '''
        Combine two lists discarding any two elements that are within tol of 
        eachothers value. Returns teh sorted result.    
    '''
    npaXnew = np.append(npaX1,npaX2)
    lDiscard = []
    for nIdx in range(npaXnew.size):
        for nJIdx in range(nIdx + 1, npaXnew.size):
            if abs(npaXnew[nIdx] - npaXnew[nJIdx]) <= tol:
                lDiscard.append(nJIdx)
    lDiscard = list(set(lDiscard))
    npaXnew = npaXnew.take([elem for elem in range(npaXnew.size) if elem not in lDiscard])
    npaXnew.sort(0)    
    return npaXnew
    
def matchUp(npaX1, npaY1, npaX2, npaY2):
    '''
        Take a pair of sequences that define 2 segments cartesian coordinates 
        and interpolate the loop (in a parameter) such that:
        1) The resulting interpolations have even spacing in the parameter.
        2) The resulting interpolation variable is in the interval [0, 1].
        3) The resulting interpolations have the same number of points.
        4) All the original points are in the interpolation.
    '''
    if len(npaX1) != len(npaY1) or len(npaX2) != len(npaY2):
        print 'Error in match up: len(x) != len(y)'
        return None
    nPoints = lcm(len(npaX1) - 1, len(npaX2) - 1) + 1
    npaT1 = np.linspace(0, 1, len(npaX1))
    npaT2 = np.linspace(0, 1, len(npaX2))
    npaTnew = np.linspace(0, 1, nPoints)
    from scipy import interp
    lY1new = []
    lX1new = []
    lY2new = []
    lX2new = []
    for nVert in range(len(npaX1) - 1):
        if nVert < len(npaX1) - 2:
            tempt = [elem for elem in npaTnew if (elem < npaT1[nVert + 1] and elem >= npaT1[nVert])]
        else: #get the rest
            tempt = [elem for elem in npaTnew if elem >= npaT1[nVert]]
        tempx = interp(tempt, npaT1[nVert:nVert + 2], npaX1[nVert:nVert + 2])
        tempy = interp(tempt, npaT1[nVert:nVert + 2], npaY1[nVert:nVert + 2])
        lY1new.extend(tempy)
        lX1new.extend(tempx)
    for nVert in range(len(npaX2) - 1):
        if nVert < len(npaX2) - 2:
            tempt = [elem for elem in npaTnew if (elem < npaT2[nVert + 1] and elem >= npaT2[nVert])]
        else: #get the rest
            tempt = [elem for elem in npaTnew if elem >= npaT2[nVert]]
        tempx = interp(tempt, npaT2[nVert:nVert + 2], npaX2[nVert:nVert + 2])
        tempy = interp(tempt, npaT2[nVert:nVert + 2], npaY2[nVert:nVert + 2])
        lY2new.extend(tempy)
        lX2new.extend(tempx)
    return     lX1new, lY1new, lX2new, lY2new, npaTnew



def gcd(*numbers):
    """Return the greatest common divisor of the given integers"""
    from fractions import gcd
    return reduce(gcd, numbers)

# Least common multiple is not in standard libraries? It's in gmpy, but this is simple enough:

def lcm(*numbers):
    """Return lowest common multiple."""    
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    return reduce(lcm, numbers, 1)

def intermFFT(npaR1, npaR2, nSteps):
    '''
        Take a set of data and do a fast fourier transform, then linearly 
        interpolate the the fft to generate intermediate steps. Use inverse fft 
        to untransform intermediate steps.
    '''
    if len(npaR1) != len(npaR2):
        print 'Error: len(npaR1) != len(npaR2)'
        return None
    npaR1T = np.fft.rfft(npaR1)
    npaR2T = np.fft.rfft(npaR2)
    npaRreal, npaS = np.array(intermSimp(np.real(npaR1T), np.real(npaR2T), nSteps))
    npaRimag, npaS = np.array(intermSimp(np.imag(npaR1T), np.imag(npaR2T), nSteps))
    npaOut = np.zeros([len(npaR1), nSteps])
    for nColIdx in range(npaRreal.shape[1]):
        temp = npaRreal[:, nColIdx] + 1.0j*npaRimag[:, nColIdx]
        npaOut[:,nColIdx] = np.fft.irfft(temp, n=len(npaR1))
    return npaOut, npaS

def intermSimp(npaX1, npaX2, nSteps):
    '''
        Take two sequences and smoothly transform the the first  sequence into 
        the second sequence by interpolating the coordinates in the 
        intermediate steps.
    '''
    if len(npaX1) != len(npaX2):
        print 'Error in intemSimpCart: sequences must be same length'
        return None
    elif nSteps < 2:
        print 'Error in intermSimpCart: must draw at lest 2 frames'
        return None
    npaXs = np.zeros([len(npaX1), nSteps], dtype=float)
    npaS = np.linspace(0,1, nSteps)
    from scipy.interpolate import interp1d
    for nIdx in range(len(npaX1)):
        fX= interp1d([0, 1], [npaX1[nIdx], npaX2[nIdx]])
        npaXs[nIdx, :] = fX(npaS)
    return npaXs, npaS

def toPol(npaX, npaY, bInv = False):
    '''
        Convert ordered pairs in the arrays x and y and convert into polar
        coordinates. bInv = True will do the inverse operation with x acting as 
        theta and y acting as r. 
    '''
    if not bInv:
        npaR = np.sqrt(np.power(npaX, 2) + np.power(npaY, 2))
        npaTheta = np.arctan2(npaY, npaX)
        return npaR, npaTheta
    else:
        npaXnew = np.multiply(npaY, np.cos(npaX))
        npaYnew = np.multiply(npaY, np.sin(npaX))
        return npaXnew, npaYnew

def loopCycle(npaTheta1, npaTheta2):
    '''
        Take two ordered pairs that represent loops and then return the cyclic 
        permutation of the loops that minimizes the difference between theta1
        and theta2.
    '''
    if len(npaTheta1) != len(npaTheta2):
        print 'Error in loopCycle: both arguments must be same length.'
        return None
    if npaTheta1[-1] != npaTheta1[-1]:
        print 'Error in loopCycle: First sequence is not a loop.'
        return None
    if npaTheta2[-1] != npaTheta2[-1]:
        print 'Error in loopCycle: Second sequence is not a loop.'
        return None
    nSize = len(npaTheta1) - 1
    fLowNorm = np.linalg.norm(npaTheta1[:-1]) + np.linalg.norm(npaTheta2[:-1])
    lnBest = []
    for i in range(nSize):
        lnCyc = [j % nSize for j in range(i, i + nSize)]
        temp = np.linalg.norm(npaTheta1[lnCyc]-npaTheta2[:-1])
        if temp < fLowNorm:
            fLowNorm = temp
            lnBest = [elem for elem in lnCyc]
    return lnBest

def reOrderLoop(lnPermute, npaData):
    '''
        Take a permutation and reorder and restore a loop based on the indices
        passed.
    '''
    newData = [npaData[i] for i in lnPermute]
    newData.append(newData[0])
    return newData
    
def centroid(x, y):
    '''
        Find the centroid of a finite set of points and return it. I will 
        probably replace this with a more sophisticated algorithm that is less
        sensitive to intepolating the boundary. 
    '''
    x0 = float(sum(x)) / float(len(x))
    y0 = float(sum(y)) / float(len(y))
    return [x0,y0]

def plotFrame(nFrame, npaXData, npaYData, lines):
    for line, xdat, ydat in zip(lines, npaXData, npaYData):
        line.set_data([xdat[:, nFrame], ydat[:, nFrame]])
    return lines
    
    
##########################Test code

x1,y1 = regPolyVert(3)
x2,y2 = regPolyVert(4,fOrient = np.pi/4.)
x1new, y1new, x2new, y2new, t = matchUp(x1, y1, x2, y2)

r1, theta1 = toPol(np.array(x1new), np.array(y1new))
r2, theta2 = toPol(np.array(x2new), np.array(y2new))

perm = loopCycle(theta1, theta2)
r1 = reOrderLoop(perm, r1)
theta1 = reOrderLoop(perm, theta1)
x1new = reOrderLoop(perm, x1new)
y1new = reOrderLoop(perm, y1new)

nframes = 120
npaThetas, npaS = intermSimp(theta1, theta2, nframes)
npaRs, npaS = intermSimp(r1, r2, nframes)
npaRfft, npaS = intermFFT(r1, r2, nframes)

npaXpol, npaYpol = toPol(npaThetas, npaRs, bInv=True)
npaXfft, npaYfft = toPol(npaThetas, npaRfft, bInv=True)
npaXs, npaS = intermSimp(x1new, x2new, nframes)
npaYs, npaS = intermSimp(y1new, y2new, nframes)

fig, ax = plt.subplots(3, 1)
line1 = ax[0].plot(npaXpol[:, 0],npaYpol[:, 0], marker='o')
line2 = ax[1].plot(npaXs[:, 0],npaYs[:, 0], marker='o')
line3 = ax[2].plot(npaXfft[:,0], npaYfft[:,0], marker='o')

lines = [line1[0], line2[0], line3[0]]
npaXl = [npaXpol, npaXs, npaXfft]
npaYl = [npaYpol, npaYs, npaYfft]

for i in range(len(ax)):
    ax[i].set_aspect(aspect='equal')
line_ani = animation.FuncAnimation(fig, plotFrame, nframes, fargs=(npaXl, npaYl, lines), interval=16, blit=True)
plt.show()

