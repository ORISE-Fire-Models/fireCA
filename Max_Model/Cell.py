import numpy as np
from scipy.optimize import minimize

class Cell:
    def __init__(self, ignPt, x, y, aos, res):
        # directions are in order sw, w, nw, n, ne, e, se, s
        self.dirs = np.array([[0, 0], [0, 0.5], [0, 1], [0.5, 1], [1, 1], [1, 0.5], [1, 0], [0.5, 0]])
        self.x = x
        self.y = y
        self.findIgnPt(ignPt)
        #self.RoS = np.array([np.nan] * 8)
        self.aos = aos
        self.phi = [self.dist(d,self._ignPt)*res for d in self.dirs]
        

    def __str__(self):
        return f"<{self.x},{self.y}>"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (self.x==other.x and self.y==other.y)

    def __lt__(self, other):
        return (self.x + self.y) < (other.x + other.y)

    def dist(self, a, b):
        return np.sqrt((b[0]-a[0])**2+(b[1]-a[1])**2)
  
    # returns a list of times it would take to the fire to get to each edge
    def getTravelTimes(self):
        tt = np.array([])
        for i in range(len(self.phi)):
            if self.phi[i] > 0 and self.RoS[i] > self.phi[i]:
                tt = np.append(tt, self.phi[i]/self.RoS[i])
        return tt
    
    # updates phi by subtracting the RoS* step size from the distance to the edge of the cell
    # step size units is hour
    # returns the edges that the fire has reached
    def updatePhi(self, s):
        burnt = []
        if any(np.isnan(self.RoS)):
            #print("RoS has NaN values")
            return np.array(burnt)
        else:
            for p in range(len(self.phi)):
                if self.phi[p] > 0:
                    self.phi[p] = self.phi[p] - self.RoS[p]*s
                    if self.phi[p] <= 0:
                        burnt.append(self.dirs[p])
            return np.array(burnt)
    
    # TODO: may want to adjust so negative thetas are made to be positive, idk
    # currently output of arctan2 is [-pi, pi]
    def calcThetas(self):
        ix = self._ignPt[0]
        iy = self._ignPt[1]
        thetas = np.array([np.arctan2(e[1]-iy, e[0]-ix) for e in self.dirs])
        for t in range(len(thetas)):
            if thetas[t] < 0:
                thetas[t] = thetas[t] + (2*np.pi)
        return thetas
    
    # TODO: check that this is doing the projections properly
    def calcProjections(self, v):
        u = np.abs(v[0])
        uv = v[1]
        projs = [u*np.cos(uv-t) for t in self.thetas]
        # set any negative values to 0
        for i in range(len(projs)):
            if projs[i] < 0:
                projs[i] = 0
        return projs

    # calculate f
    def f(self, d):
        s = 0
        for i in range(len(d)-1):
            area = d[i]*d[i+1]*np.abs(np.sin(np.abs(self.thetas[i]-self.thetas[i+1])))
            #if np.isinf(area):
                #print('\n area is inf @',self.x, self.y, i, d[i],d[i+1], self.thetas[i], self.thetas[i+1], np.abs(np.sin(np.abs(self.thetas[i]-self.thetas[i+1]))))
            s = s + area
        # possibly only do this part if it's a center ignition
        if self.isCenter: #should always happen
            s = s + d[-1]*d[0]*np.sin(np.abs(self.thetas[0]-self.thetas[-1])) #should be thetas[0], not thetas[1]
        return 0.5*s-self.aos

    def partialDiff(self, d1, d3, t1, t2, t3):
        return d1*np.sin(np.abs(t2-t1))+d3*np.sin(np.abs(t3-t2))
    
    # calculate the gradient of f
    # TODO: might want to check this later, but seems fiine
    def df(self, d):
        gradf = [None]*(len(d))
        if self.isCenter:
            m = 0
        else:
            m = 1
        for i in range(m, len(d)-1):
            gradf[i] = self.partialDiff(d[i-1], d[i+1], self.thetas[i-1], self.thetas[i], self.thetas[i+1])
        
        # possibly only do this part if it's a center ignition
        if self.isCenter:
            gradf[-1] = (d[-2]*np.sin(np.abs(self.thetas[-1]-self.thetas[-2]))+
                                d[0]*np.sin(np.abs(self.thetas[0]-self.thetas[-1])))
        else:
            gradf[0] = self.partialDiff(0, d[1], 0, self.thetas[0], self.thetas[1])
            gradf[-1] = self.partialDiff(d[-2], 0, self.thetas[-2], self.thetas[-1], 0)

        #gradf = np.array(gradf) #+ _f/np.abs(_f)
        return np.array(gradf)

    def fStep(self, a, xs, dxs):
        out = self.f(xs-(a*dxs))**2
        #print(out)
        return out

    def getSlope(self, x1, x2, y1, y2):
        if (x2-x1) == 0:
            return np.nan
        return (y2-y1)/(x2-x1)

    def getYIntercept(self, x, y, m):
        return y - (m*x)
    
    def findIntersectionDistance(self, d1, t1, d2, t2, t):
        t1 = t1-t
        t2 = t2-t
        x1 = d1*np.cos(t1)
        y1 = d1*np.sin(t1)
        x2 = d2*np.cos(t2)
        y2 = d2*np.sin(t2)
        m = self.getSlope(x1, x2, y1, y2)
        if m == np.nan:
            return x1
        if m == 0:
            return 0
        b = self.getYIntercept(x1, y1, m)
        return b/(-m)

    def convexify(self, d, t):
        for i in range(len(d)-1):
            minDist = self.findIntersectionDistance(d[i+1], t[i+1], d[i-1], t[i-1], t[i])
            if d[i] < minDist:
                d[i] = minDist
                
        if self.isCenter:
            minDist = self.findIntersectionDistance(d[0], t[0], d[-2], t[-2], t[-1])
            if d[-1] < minDist:
                d[-1] = minDist

        return d
        
    # v is combined wind slope vector
    # v[0] is the length and v[1] is the angle
    
    # TODO: check the units for the combined wind slope vector
    #       currently assuming the combined wind slope vector has units m/hr
    
    # aos (area of spread) is the rate of spread calculated by Zack in acres/day, we'll assume here that it's already been converted to acres/hr
    # should this be converted to m^2/hr? 
    
    # here we will calculate RoS based on the fire spread after 1 hour
    # in the future maybe we can make this time variable based on aos or just change it to find a good time that works for 'most' cells

    # step is the step size, maybe will remove this and add some kind of variable stepping
    def calcRoS(self, v, step):
        self.thetas = self.calcThetas()
        d = self.calcProjections(v)
        g = self.f(d)
        df = self.df(d)*g
        i = 0
        minError = np.linalg.norm(df)
        minD = d
        gIsZero = False
        while(np.linalg.norm(df)>0.0001 and i<500):
            a = minimize(self.fStep, step, args=(d, df), options={'maxiter':100}) #gets a guess at best stepsize for gradient descent
            step = a.x
            d = d-a.x*df
            #d = d - step*df

            # constraints
            # all RoS >= 0
            d = self.replaceZeros(d)
            # the RoS shape must be convex
            d = self.convexify(d, self.thetas)

            g = self.f(d)
                
            df = self.df(d)*g
            if np.linalg.norm(df) < minError:
                minError = np.linalg.norm(df)
                minD = d
            if g == 0:
                self.RoS = d
                gIsZero = True
                break
            i = i+1
        if i == 500 and not gIsZero:
            self.RoS = minD
        elif not gIsZero:
            self.RoS = d
    
    def replaceZeros(self, lis):
        for i in range(len(lis)):
            if lis[i] < 0:
                lis[i] = 0
        return lis

    def clearDirs(self):
        x = self._ignPt[0]
        y = self._ignPt[1]
        if not self.isCenter:
            mask = [np.array_equal(a1, [x, y]) for a1 in self.dirs]
            #self.dirs = self.dirs[mask]
            if self.isDiagonal:
                mask2 = [np.array_equal(a1, [x, 0.5]) for a1 in self.dirs]
                mask3 = [np.array_equal(a1, [0.5, y]) for a1 in self.dirs]
                mask = [a or b or c for a, b, c in zip(mask, mask2, mask3)]
            mask = [not val for val in mask]
            self.dirs = self.dirs[mask]
    
    def findIgnPt(self,ignPt):
        self.isCenter = False
        self.isDiagonal = False
        
        match ignPt:
            case 'sw':
                self._ignPt = [0, 0]
                self.isDiagonal = True
            case 'w':
                self._ignPt = [0, 0.5]
            case 'nw':
                self._ignPt = [0, 1]
                self.isDiagonal = True
            case 'n':
                self._ignPt = [0.5, 1]
            case 'ne':
                self._ignPt = [1, 1]
                self.isDiagonal = True
            case 'e':
                self._ignPt = [1, 0.5]
            case 'se':
                self._ignPt = [1, 0]
                self.isDiagonal = True
            case 's':
                self._ignPt = [0.5, 0]
            case 'c':
                self.isCenter = True
                self._ignPt = [0.5, 0.5]
            case _:
                self._ignPt = np.nan
                
        self.clearDirs()
                
    def get_ignPt(self): 
        return self._ignPt
       

    def set_ignPt(self, a): 
        self.findIgnPt(a)
  

    def del_ignPt(self): 
        del self._ignPt 
     
    ignPt = property(get_ignPt, set_ignPt, del_ignPt) 