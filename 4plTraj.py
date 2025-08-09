#!/usr/bin/env python3
import numpy as np  # for array and matrix operations
import matplotlib.pyplot as plt  # for graphing


def windowCoord(W1,W2,W3,W4):
    W = np.zeros((4,3))

    W[0,0] = W1[0]
    W[0,1] = W2[1]
    W[0,2] = W4[2]

    W[1,0] = W2[0]
    W[1,1] = W1[1]
    W[1,2] = W2[2]

    W[2,0] = W3[0]
    W[2,1] = W4[1]
    W[2,2] = W3[2]

    W[3,0] = W4[0]
    W[3,1] = W3[1]
    W[3,2] = W1[2]

    return W

class trajFeasibility:
    def __init__(self, bounds: dict, eps: float, i, f):
        self.openBounds(bounds)
        self.eps = eps
        self.i = i
        self.f = f
    
    def openBounds(self, bounds: dict):
        self.tf = bounds['tf']
        self.Vbnd = bounds['Vbnd']
        self.Abnd = bounds['Abnd']
        self.Jbnd = bounds['Jbnd']

    def velMax(self, B):
        return np.abs(((B+1)**((1/B)+1))/(4*B*self.Vbnd*((B-1)**((1/B)-1))))
    
    def accMax(self, B):
        a = 4*(B-1)*(B+1)
        b = 2*B*np.sqrt(3*(B-1)*(B+1))
        c = 2*(B+1)*(B+2)

        k1 = (a+b)/c
        k2 = (a-b)/c

        k = np.array([k1,k2])
        f = np.zeros((2))

        for i in range(len(k)):
            f[i] = np.abs((B*(k[i]**(1-(2/B)))*((B-1)-k[i]*(B+1)))/(self.Abnd*(1+k[i])**3))

        return max(np.abs(f))
    
    def jerkMax(self, B):
        k1 = 24*(B**3)
        k2 = 36*(B**2)*(B-1)
        k3 = 2*B*(B-1)*(7*B-11)
        k4 = (B-1)*(B-2)*(B-3)

        a = k1 - k2 + k3 - k4
        b = -(k2 - 2*k3 - 3*k4)
        c = k3 - 3*k4
        d = -k4

        # Step 1: Reduce the cubic equation to depressed form: t^3 + pt + q = 0
        if a != 1:
            b, c, d = b / a, c / a, d / a
        
        p = c - (b ** 2) / 3
        q = (2 * b ** 3) / 27 - (b * c) / 3 + d
        
        # Step 2: Calculate discriminant
        discriminant = (q / 2) ** 2 + (p / 3) ** 3
        
        # Step 3: Use Cardano's formula to find roots
        if discriminant > 0:
            # One real root, two complex conjugates
            u = (-q / 2 + np.sqrt(discriminant)) ** (1 / 3)
            v = (-q / 2 - np.sqrt(discriminant)) ** (1 / 3)
            
            root1 = u + v
            root2 = -(u + v) / 2 + (u - v) * np.sqrt(3) * 1j / 2
            root3 = -(u + v) / 2 - (u - v) * np.sqrt(3) * 1j / 2
            root = np.array([root1])

        
        elif discriminant == 0:
            # All roots real, at least two equal
            u = (-q / 2) ** (1 / 3)
            
            root1 = 2 * u
            root2 = -u
            root3 = -u
            root = np.array([root1, root2])

        else:
            # All roots real and distinct
            theta = np.arccos(-q / (2 * (-p / 3) ** (3 / 2)))
            root1 = 2 * (-p / 3) ** 0.5 * np.cos(theta / 3)
            root2 = 2 * (-p / 3) ** 0.5 * np.cos((theta + 2 * np.pi) / 3)
            root3 = 2 * (-p / 3) ** 0.5 * np.cos((theta + 4 * np.pi) / 3)
            root = np.array([root1,root2,root3])
        
        # Step 4: Convert roots back to original variable x
        root -= b / 3

        kj = []

        for i in range(len(root)):
            kj.append(B*(root[i]**(1-(3/B)))*((6*B*(root[i]**2)/((1+root[i])**2)) - 6*B*(B-1)*root[i]/(1+root[i]) + (B-1)*(B-2))/((1+root[i])**2)/self.Jbnd)

        return max(np.abs(np.array(kj)))
    
    def upperBoundC(self, B):
        return self.tf*((self.eps/(np.abs(self.i - self.f) - self.eps))**(1/B))
    
    def lowerBoundC(self, B):
        c1 = self.velMax(B)*np.abs(self.i-self.f)
        c2 = (np.abs(self.i-self.f)*self.accMax(B))**(1/2)
        c3 = (np.abs(self.i-self.f)*self.jerkMax(B))**(1/3)

        return max(c1,c2,c3)
    
    def algo1(self):
        sol = []
        for B in range(4,11):
            uC = self.upperBoundC(B)
            lC = self.lowerBoundC(B)
            if uC > lC:
                sol.append(np.array([B,uC,lC]))
        
        sol = np.array(sol)

        return sol
    
class windowCheck:
    def __init__(self, W, i, f, r, SdynX, SdynY, SdynZ):
        self.W = W
        self.i = i
        self.f = f
        self.SdynX = SdynX
        self.SdynY = SdynY
        self.SdynZ = SdynZ
        self.r = r

    def Sx(self):
        SxWnd = []
        Sx = []
        for posX in self.SdynX:
            for Cx in np.linspace(posX[2], posX[1], num=int((posX[1] - posX[2]) / 0.1) + 1):
                for posY in self.SdynY:
                    for posZ in self.SdynZ:
                        lCxyWnd = Cx*(((self.W[3,0] + self.r - self.i[0])/(self.f[0] - self.W[3,0] - self.r))**(1/posX[0]))*(((self.f[1] - self.W[2,1] + self.r)/(self.W[2,1] - self.r - self.i[1]))**(1/posY[0]))
                        uCxyWnd = Cx*(((self.W[0,0] + self.r - self.i[0])/(self.f[0] - self.W[0,0] - self.r))**(1/posX[0]))*(((self.W[1,1] + self.r - self.i[1])/(self.f[1] - self.W[1,1] - self.r))**(1/posY[0]))
                        lCxzWnd = Cx*(((self.W[3,0] + self.r - self.i[0])/(self.f[0] - self.W[3,0] - self.r))**(1/posX[0]))*(((self.f[2] - self.W[2,2] + self.r)/(self.W[2,2] - self.r - self.i[2]))**(1/posZ[0]))
                        uCxzWnd = Cx*(((self.W[0,0] + self.r - self.i[0])/(self.f[0] - self.W[0,0] - self.r))**(1/posX[0]))*(((self.W[1,2] + self.r - self.i[2])/(self.f[2] - self.W[1,2] - self.r))**(1/posZ[0]))
                        
                        if (lCxyWnd < uCxyWnd) and (lCxzWnd < uCxzWnd):
                            lCxyfes = max(lCxyWnd,posY[2])
                            uCxyfes = min(uCxyWnd,posY[1])
                            lCxzfes = max(lCxzWnd,posZ[2])
                            uCxzfes = min(uCxzWnd,posZ[1])

                            if (lCxyfes < uCxyfes) and (lCxzfes < uCxzfes):
                                a = np.array([posX[0],Cx,posY[0],lCxyfes,uCxyfes,posZ[0],lCxzfes,uCxzfes])
                                SxWnd.append(a)
        
        for q in SxWnd:
            for Cy in np.linspace(q[3], q[4], num=int((q[4] - q[3]) / 0.1) + 1):
                for Cz in np.linspace(q[6], q[7], num=int((q[7] - q[6]) / 0.1) + 1):
                    if Cz == 0 or Cy == 0 or q[1] == 0:
                        continue
                    sol = np.array([q[0],q[1],q[1],Cy,q[5],Cz])
                    Sx.append(sol)
        
        Sx = np.array(Sx)

        if Sx.size == 0:
            Sx = np.empty((1,6))

        return Sx
    
    def Sy(self):
        SyWnd = []
        Sy = []
        for posY in self.SdynY:
            for Cy in np.linspace(posY[2], posY[1], num=int((posY[1] - posY[2]) / 0.1) + 1):
                for posZ in self.SdynZ:
                    for posX in self.SdynX:
                        lCxyWnd = Cy*(((self.W[3,1] + self.r - self.i[1])/(self.f[1] - self.W[3,1] - self.r))**(1/posY[0]))*(((self.f[0] - self.W[2,0] + self.r)/(self.W[2,0] - self.r - self.i[0]))**(1/posX[0]))
                        uCxyWnd = Cy*(((self.W[0,1] + self.r - self.i[1])/(self.f[1] - self.W[0,1] - self.r))**(1/posY[0]))*(((self.W[1,0] + self.r - self.i[0])/(self.f[0] - self.W[1,0] - self.r))**(1/posX[0]))
                        lCyzWnd = Cy*(((self.W[3,1] + self.r - self.i[1])/(self.f[1] - self.W[3,1] - self.r))**(1/posY[0]))*(((self.f[2] - self.W[2,2] + self.r)/(self.W[2,2] - self.r - self.i[2]))**(1/posZ[0]))
                        uCyzWnd = Cy*(((self.W[0,1] + self.r - self.i[1])/(self.f[1] - self.W[0,1] - self.r))**(1/posY[0]))*(((self.W[1,2] + self.r - self.i[2])/(self.f[2] - self.W[1,2] - self.r))**(1/posZ[0]))
                        
                        if (lCxyWnd < uCxyWnd) and (lCyzWnd < uCyzWnd):
                            lCxyfes = max(lCxyWnd,posX[2])
                            uCxyfes = min(uCxyWnd,posX[1])
                            lCyzfes = max(lCyzWnd,posZ[2])
                            uCyzfes = min(uCyzWnd,posZ[1])

                            if (lCxyfes < uCxyfes) and (lCyzfes < uCyzfes):
                                a = np.array([posY[0],Cy,posX[0],lCxyfes,uCxyfes,posZ[0],lCyzfes,uCyzfes])
                                SyWnd.append(a)
        
        for q in SyWnd:
            for Cx in np.linspace(q[3], q[4], num=int((q[4] - q[3]) / 0.1) + 1):
                for Cz in np.linspace(q[6], q[7], num=int((q[7] - q[6]) / 0.1) + 1):
                    if Cx == 0 or Cz == 0 or q[1] == 0:
                        continue
                    sol = np.array([q[2],Cx,q[0],q[1],q[5],Cz])
                    Sy.append(sol)
        
        Sy = np.array(Sy)

        if Sy.size == 0:
            Sy = np.empty((1,6))
        
        return Sy
    
    def Sz(self) -> np.ndarray:
        SzWnd = []
        Sz = []
        for posZ in self.SdynZ:
            for Cz in np.linspace(posZ[2], posZ[1], num=int((posZ[1] - posZ[2]) / 1) + 1):
                for posY in self.SdynY:
                    for posX in self.SdynX:
                        lCzyWnd = Cz*(((self.W[3,2] + self.r - self.i[2])/(self.f[2] - self.W[3,2] - self.r))**(1/posZ[0]))*(((self.f[1] - self.W[2,1] + self.r)/(self.W[2,1] - self.r - self.i[1]))**(1/posY[0]))
                        uCzyWnd = Cz*(((self.W[0,2] + self.r - self.i[2])/(self.f[2] - self.W[0,2] - self.r))**(1/posZ[0]))*(((self.W[1,1] + self.r - self.i[1])/(self.f[1] - self.W[1,1] - self.r))**(1/posY[0]))
                        lCxzWnd = Cz*(((self.W[3,2] + self.r - self.i[2])/(self.f[2] - self.W[3,2] - self.r))**(1/posZ[0]))*(((self.f[0] - self.W[2,0] + self.r)/(self.W[2,0] - self.r - self.i[0]))**(1/posX[0]))
                        uCxzWnd = Cz*(((self.W[0,2] + self.r - self.i[2])/(self.f[2] - self.W[0,2] - self.r))**(1/posZ[0]))*(((self.W[1,0] + self.r - self.i[0])/(self.f[0] - self.W[1,0] - self.r))**(1/posX[0]))
                        
                        if (lCzyWnd < uCzyWnd) and (lCxzWnd < uCxzWnd):
                            lCzyfes = max(lCzyWnd,posY[2])
                            uCzyfes = min(uCzyWnd,posY[1])
                            lCxzfes = max(lCxzWnd,posZ[2])
                            uCxzfes = min(uCxzWnd,posZ[1])

                            if (lCzyfes < uCzyfes) and (lCxzfes < uCxzfes):
                                a = np.array([posZ[0],Cz,posY[0],lCzyfes,uCzyfes,posX[0],lCxzfes,uCxzfes])
                                SzWnd.append(a)
        
        for q in SzWnd:
            for Cy in np.linspace(q[3], q[4], num=int((q[4] - q[3]) / 0.1) + 1):
                for Cx in np.linspace(q[6], q[7], num=int((q[7] - q[6]) / 0.1) + 1):
                    if Cx == 0 or Cy == 0 or q[1] == 0:
                        continue
                    sol = np.array([q[5],Cx,q[2],Cy,q[0],q[1]])
                    Sz.append(sol)
        
        Sz = np.array(Sz)

        if Sz.size == 0:
            return np.empty((1,6))

        else: return Sz


def plot3D(S,W,i,f,tf):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    Wx = W[:,0]
    Wy = W[:,1]
    Wz = W[:,2]

    t = np.linspace(0,tf,50)

    for s in S:
        x = f[0] + (i[0] - f[0])/(1 + (t/s[1])**s[0])
        y = f[1] + (i[1] - f[1])/(1 + (t/s[3])**s[2])
        z = f[2] + (i[2] - f[2])/(1 + (t/s[5])**s[4])
        ax.plot(x,y,z,color='g')

    ax.plot(Wx,Wy,Wz,color='r')
    ax.plot(i[0],i[1],i[2],color='b')
    ax.plot(f[0],f[1],f[2],color='b')

    plt.show()


def main():
    bounds = {
        'tf': 10,
        'Vbnd': 5,
        'Abnd': 10,
        'Jbnd': 20,
    }

    i = [0,0,0]
    f = [5,3,3]
    
    r = 0.045

    eps = 0.01

    # Case 1
    W1 = [2.5,0.05,0.05]
    W2 = [2.5,2.95,0.05]
    W3 = [2.5,2.95,2.95]
    W4 = [2.5,0.05,2.95]

    # Case 2
    # W1 = [2.51,1.43,1.43]
    # W2 = [2.51,1.56,1.43]
    # W3 = [2.48,1.56,1.56]
    # W4 = [2.48,1.43,1.56]

    # Case 3a
    # W1 = [4.20,2.25,2.41]
    # W2 = [4.28,2.68,2.58]
    # W3 = [3.79,2.77,2.58]
    # W4 = [3.71,2.31,2.41]

    # Case 3b
    # W1 = [1.25,1.87,0.78]
    # W2 = [1.25,2.12,1.12]
    # W3 = [0.75,2.12,1.21]
    # W4 = [0.75,1.87,0.78]

    # Case 3c
    # W1 = [4.6,0.9,1]
    # W2 = [4.6,1.1,1]
    # W3 = [4.4,1.1,1]
    # W4 = [4.4,0.9,1]

    W = windowCoord(W1,W2,W3,W4)
    W_ = np.array([W1,W2,W3,W4,W1])

    feasibleX = trajFeasibility(bounds, eps, i[0], f[0])
    feasibleY = trajFeasibility(bounds, eps, i[1], f[1])
    feasibleZ = trajFeasibility(bounds, eps, i[2], f[2])

    SdynX = feasibleX.algo1()
    SdynY = feasibleY.algo1()
    SdynZ = feasibleZ.algo1()

    # print(SdynX)
    # print(SdynY)
    # print(SdynZ)

    wndSx = windowCheck(W,i,f,r,SdynX,SdynY,SdynZ)
    Sx = wndSx.Sx()
    wndSy = windowCheck(W,i,f,r,SdynX,SdynY,SdynZ)
    Sy = wndSy.Sy()
    wndSz = windowCheck(W,i,f,r,SdynX,SdynY,SdynZ)
    Sz = wndSz.Sz()

    S = np.concatenate((Sx,Sy,Sz),axis=0)
    print(S)
    print(S.shape)

    plot3D(Sz,W_,i,f,bounds['tf'])


if __name__ == "__main__":
    main()