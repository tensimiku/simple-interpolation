import matrix
class Point:
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

class Interpolation:
    def __init__(self, points):
        self.points = points
        self.plen = len(points)
        self.divtable = []

    def __lag__(self, x, index):
        xu = 1
        xl = 1
        for i in range(0, self.plen):
            if index == i:
                continue
            xu *= x - self.points[i].x
        for i in range(0, self.plen):
            if index == i:
                continue
            xl *= self.points[index].x - self.points[i].x
        return xu/xl

    def lagrange(self, x):
        px = 0
        for i in range(0, self.plen):
            px += self.points[i].y*self.__lag__(x, i)
        return px

    def __divdif__(self, n):
        if self.divtable:
            return self.divtable[n]
        rtable = []
        dtable = []
        for i in range(1, self.plen):
            diff = (self.points[i].y-self.points[i-1].y)/(self.points[i].x - self.points[i-1].x)
            dtable.append(diff)
        rtable.append(dtable)
        for i in range(2, self.plen):
            dtable = []
            for x in range(i, self.plen):
                yval = rtable[i-2]
                ydiff = (yval[x-i+1] - yval[x-i]) / (self.points[i].x - self.points[0].x)
                dtable.append(ydiff)
            rtable.append(dtable)
        for val in rtable:
            self.divtable.append(val[0])
        return self.divtable[n]

    def newton(self, x):
        plen = self.plen
        px = (x - self.points[plen-2].x) * self.__divdif__(plen-2)
        for i in range(plen-3, -1, -1):
            px += self.__divdif__(i)
            px *= x - self.points[i].x
        px += self.points[0].y
        return px

    def __sindex__(self, x):
        index = 0
        for i in range(0, self.plen-1):
            xi1 = self.points[i+1].x
            xi = self.points[i].x
            if x >= xi and x <= xi1:
                index = i
                break
        return index

    def splinesquare(self, x):
        z = [0] #init value
        for i in range(0, self.plen-1):
            pi1 = self.points[i+1]
            pi = self.points[i]
            _z = -z[i] + 2*((pi1.y - pi.y)/(pi1.x-pi.x))
            z.append(_z)

        i = self.__sindex__(x)
        pi1 = self.points[i+1]
        pi = self.points[i]
        a = ((z[i+1]-z[i])/(pi1.x - pi.x))/2
        #print(str(a)+"(x-"+str(pi.x)+")^2+"+str(z[i])+"(x-"+str(pi.x)+")+"+str(pi.y))
        return a*((x-pi.x)**2)+z[i]*(x-pi.x) + pi.y

    def splinecubic(self, x):
        plen = self.plen-1
        p1 = self.points[1]
        p0 = self.points[0]
        z = [0] * self.plen
        h = [0] * plen
        u = [0] * plen
        c = [0] * plen
        v = [0] * plen
        h[0] = p1.x-p0.x
        c[0] = (p1.y-p0.y)/h[0]
        v = [0] * self.plen
        for i in range(1, plen):
            pi1 = self.points[i+1]
            pi = self.points[i]
            pi_1 = self.points[i-1]
            h[i] = (pi1.x - pi.x)
            c[i] = ((pi1.y - pi.y)/h[i])
            u[i] = 2*(h[i-1]+h[i])
            v[i] = 6*(c[i]-c[i-1])
        mat = [[0]*(plen-1) for i in range(0, plen-1)]
        mat[0][0] = u[1]
        if plen > 2:
            mat[0][1] = h[1]
            mat[plen-2][plen-3] = h[plen-2]
            mat[plen-2][plen-2] = u[plen-1]
        for i in range(1, plen-2):
            mat[i][i-1] = h[i]
            mat[i][i] = u[i+1]
            mat[i][i+1] = h[i+1]

        vecc = [[vec] for vec in v[1:plen]]
        solve = matrix.Matrix(mat, vecc)
        zmat = solve.gaussjordan()
        for  i in range(0, len(zmat)):
            z[i+1] = zmat[i][0]

        i = self.__sindex__(x)
        pi1 = self.points[i+1]
        pi = self.points[i]
        weight = []
        weight.append(z[i+1]/(6*h[i]))
        weight.append(z[i]/(6*h[i]))
        weight.append(pi1.y/h[i]-h[i]*z[i+1]/6)
        weight.append(pi.y/h[i]-h[i]*z[i]/6)
        xi1 = pi1.x - x
        xi = x - pi.x
        return weight[0]*xi**3+weight[1]*xi1**3+weight[2]*xi+weight[3]*xi1

    def lsapprox(self, x):
        mat = [[0]*2 for i in range(0,2)]
        vecc = [[0] for i in range(0,2)]
        for i in range(0, self.plen):
            mat[0][0] += self.points[i].x**2
            mat[0][1] += self.points[i].x
            vecc[0][0] += self.points[i].x * self.points[i].y
            vecc[1][0] += self.points[i].y
        mat[1][0] = mat[0][1]
        mat[1][1] = self.plen
        solve = matrix.Matrix(mat, vecc)
        ab = [ x[0] for x in solve.gaussjordanwithpivot()]
        return ab[0]*x + ab[1]

if __name__ == "__main__":
    tuples = [(0,2.5), (1, 0.5), (2, 1.5), (2.5, 1.5), (3, 1.5), (3.5, 1.25), (4, 0)]
    dataset = [Point(tuples[_i][0], tuples[_i][1]) for _i in range(0,len(tuples))]
    tet = Interpolation(dataset)
    varx = 1.5
    print("x="+str(varx))
    print("lagrange:")
    print(tet.lagrange(varx))
    print("newton:")
    print(tet.newton(varx))
    print("square spline:")
    print(tet.splinesquare(varx))
    print("cubic spline:")
    print(tet.splinecubic(varx))
    print("least square approx:")
    print(tet.lsapprox(varx))
    varx = 3.2
    print("x="+str(varx))
    print("lagrange:")
    print(tet.lagrange(varx))
    print("newton:")
    print(tet.newton(varx))
    print("square spline:")
    print(tet.splinesquare(varx))
    print("cubic spline:")
    print(tet.splinecubic(varx))
    print("least square approx:")
    print(tet.lsapprox(varx))