from tridiag import TRIDIAG_SOLVER

import math
import os

# a = [0, -1, -9, -1, 9]
# b = [-6, 13, -15, -7, -18]
# c = [5, 6, -4, 1, 0]
# d = [51, 100, -12, 47, -90]

# print(TRIDIAG_SOLVER(a,b,c,d).solve())

# Вариант 7
class PARAB_SOLVER:
    def __init__(self, x_steps = 10, scheme_type: int = 1, approx_type: int = 1):
        """
        Схема scheme_type: 1 - явная, 2 - неявная, 3 - Кранка-Николсона.
        
        Аппроксимация approx_type: 1 - 1п2т, 2 - 2п3т, 3 - 2п2т.
        """
        if scheme_type not in [1,2,3] and approx_type not in [1,2,3]:
            raise ValueError("Неверно указаны параметры решателя!")
        
        self._scheme = scheme_type
        self._approx = approx_type
        self._n = x_steps

        self._xd = math.pi / x_steps
        self._td = 0.25*self._xd**2

        self._start_cond = lambda x: math.sin(x)
        self._anal_sol = lambda x, t: math.exp(-0.5*t)*math.sin(x)

    def _write_res(self, t:float):
        PATH = os.path.split(os.path.realpath(__file__))[0] + '/t' + str(round(t,5)) + ".txt"
        with open(PATH, "w") as f:
            for i in range (self._n+1):
                cur_x = i*self._xd
                f.write(str(cur_x) + ' ' + str(self._u[i]) + ' ' + str(self._anal_sol(cur_x, t)) + '\n')

    def solve(self):
        self._u = [self._start_cond(i*self._xd) for i in range(0, self._n+1)]
        if self._scheme == 1:
            for i in range(1, self._n):
                self._u[i] = self._u[i] + self._td*(self._u[i+1]-2*self._u[i]+self._u[i-1])/self._xd**2 
                + 0.5*self._td*math.exp(-0.5*self._td)*math.cos(self._xd)
            self._u[0] = self._u[1] - self._xd*math.exp(-0.5*self._td)
            self._u[self._n] = self._u[self._n-1] - self._xd*math.exp(-0.5*self._td)
            self._write_res(self._td)
        

if __name__ == "__main__":
    solver = PARAB_SOLVER()
    solver.solve()
    print("DONE!")