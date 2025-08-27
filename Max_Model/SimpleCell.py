from Cell import Cell
import numpy as np

class SimpleCell(Cell):
    def __init__(self, ignPt, x, y, aos, res):
        super().__init__(ignPt, x, y, aos, res)

    def calcMaximum(self, d):
        m = [i + 10 for i in d]
        area = self.f(m)
        while area < self.aos:
            m = [i + 10 for i in m]
            area = self.f(m)
        return m

    def calcMinimum(self, d):
        m = [i - 10 for i in d]
        m = self.replaceZeros(m)
        area = self.f(m)
        while area > self.aos:
            m = [i - 10 for i in m]
            m = self.replaceZeros(m)
            area = self.f(m)
        return m

    # calculate f
    def f(self, d):
        s = 0
        for i in range(len(d)-1):
            area = d[i]*d[i+1]*np.abs(np.sin(np.abs(self.thetas[i]-self.thetas[i+1])))
            s = s + area
        if self.isCenter:
            s = s + d[-1]*d[0]*np.sin(np.abs(self.thetas[1]-self.thetas[-1]))
        return 0.5*s
    
    def calcRoS(self, v, step):
        self.thetas = self.calcThetas()
        #print(self.thetas)
        #print(self.dirs)
        d = self.calcProjections(v)
        #print(d)
        area = self.f(d)
        error = np.abs(area - self.aos)
        self.RoS = d
        if area < self.aos:
            minimum = d
            maximum = self.calcMaximum(d)
        elif area > self.aos:
            minimum = self.calcMinimum(d)
            maximum = d
        num_iters = 0
        while error > 0.0001 and num_iters < 500:
            d = np.add(maximum,minimum)
            d = [i/2 for i in d]
            area = self.f(d)
            if area > self.aos:
                maximum = d
            elif area < self.aos:
                minimum = d
            error = np.abs(area - self.aos)
            num_iters = num_iters + 1
            self.RoS = d
            
            
            
    
