from tridiag import TRIDIAG_SOLVER
from visual import visualise

import math
import os
import shutil
from copy import deepcopy
import sys

data_folder = "results"
DATA_PATH = os.path.join(os.path.split(os.path.realpath(__file__))[0], data_folder)
if os.path.exists(DATA_PATH):
        shutil.rmtree(DATA_PATH)
os.makedirs(DATA_PATH, exist_ok=True)

# Вариант 7
class PARAB_SOLVER:
    def __init__(self, x_steps = 20, t_steps = 1000, scheme_type: int = 1, approx_type: int = 1):
        """
        Схема scheme_type: 1 - явная, 2 - неявная, 3 - Кранка-Николсона.
        
        Аппроксимация approx_type: 1 - 1п2т, 2 - 2п3т, 3 - 2п2т.
        """
        if scheme_type not in [1,2,3] and approx_type not in [1,2,3]:
            raise ValueError("Неверно указаны параметры решателя!")
        
        self._scheme = scheme_type
        self._approx = approx_type
        self._n = x_steps

        self._xd = math.pi / self._n
        self._td = 0.1*self._xd**2
        self._t_steps = t_steps

        self._start_cond = lambda x: math.sin(x)
        self._anal_sol = lambda x, t: math.exp(-0.5*t)*math.sin(x)

    def _write_res(self, u:list[float], t:float, num: int):
        
        PATH = DATA_PATH + '/' + str(num) + ".txt"
        with open(PATH, "w") as f:
            f.write(str(round(t,5)) + '\n')
            for i in range (self._n+1):
                cur_x = i*self._xd
                f.write(str(cur_x) + ' ' + str(u[i]) + ' ' + str(self._anal_sol(cur_x, t)) + '\n')

    def solve(self):

        u = [self._start_cond(i*self._xd) for i in range(0, self._n+1)]

        u_prev = deepcopy(u)
        
        if self._scheme == 1:
            
            for j in range (1, self._t_steps+1):
                cur_t = j*self._td
                for i in range(1, self._n):
                    cur_x = self._xd*i

                    u[i] = u_prev[i] + self._td*(u_prev[i+1]-2*u_prev[i]+u_prev[i-1])/(self._xd**2) \
                        + 0.5*self._td*math.exp(-0.5*(cur_t-self._td))*math.sin(cur_x)

                if self._approx == 1:
                    u[0] = u[1] - self._xd*math.exp(-0.5*cur_t)
                    u[self._n] = u[self._n-1] - self._xd*math.exp(-0.5*cur_t)

                elif self._approx == 2:
                    u[0] = 1/3*(4*u[1] - u[2] - 2*self._xd*math.exp(-0.5*cur_t))
                    u[self._n] = 1/3*(-u[self._n-2] + 4*u[self._n-1] - 2*self._xd*math.exp(-0.5*cur_t))

                elif self._approx == 3:
                    u[0] = u[1] - self._xd*math.exp(-0.5*cur_t) + self._xd**2/2*u_prev[0]/self._td
                    u[0] /= 1 + (self._xd**2)/(2*self._td)

                    u[self._n] = u[self._n-1] - self._xd*math.exp(-0.5*cur_t) \
                        + self._xd**2/2*(u_prev[self._n]/self._td + 0.5*math.exp(-0.5*cur_t)*math.sin(self._xd*self._n))
                    u[self._n] /= 1 + (self._xd**2)/(2*self._td)

                u_prev = deepcopy(u)

                self._write_res(u, cur_t, j)

        else:
            param = 1
            if self._scheme == 3:
                param = 0.5

            for j in range (1, self._t_steps+1):
                cur_t = j*self._td
                a = [0]*(self._n+1); b = [0]*(self._n+1); c = [0]*(self._n+1); d = [0]*(self._n+1)
                    
                for i in range(1, self._n):
                    cur_x = self._xd*i

                    a[i] = -self._td*param
                    b[i] = self._xd**2 + 2*self._td*param
                    c[i] = -self._td*param
                    d[i] = u[i]*self._xd**2 + self._td*(1-param)*(u[i+1]-2*u[i]+u[i-1]) \
                        + 0.5*self._xd**2*self._td*math.sin(cur_x)*(param*math.exp(-0.5*cur_t)+(1-param)*math.exp(-0.5*(cur_t-self._td)))

                if self._approx == 1:
                    b[0] = -1
                    c[0] = 1
                    d[0] = self._xd*math.exp(-0.5*cur_t)

                    a[self._n] = -1
                    b[self._n] = 1
                    d[self._n] = -self._xd*math.exp(-0.5*cur_t)

                elif self._approx == 2:
                    b[0] = -3 - 1/self._td*a[1]
                    c[0] = 4 - 1/self._td*b[1]
                    d[0] = 2*self._xd*math.exp(-0.5*cur_t) - 1/self._td*d[1]

                    a[self._n] = -4 + 1/self._td*b[self._n-1]
                    b[self._n] = 3 + 1/self._td*c[self._n-1]
                    d[self._n] = -2*self._xd*math.exp(-0.5*cur_t) + 1/self._td*d[self._n-1]

                elif self._approx == 3:
                    b[0] = 1 + (self._xd**2)/(2*self._td)
                    c[0] = -1
                    d[0] = -self._td*(u_prev[i+1]-2*u_prev[i]+u_prev[i-1])/(self._xd**2) \
                        + 0.5*self._td*math.exp(-0.5*(cur_t-self._td))*math.sin(cur_x)

                    a[self._n] = -1
                    b[self._n] = 1 + (self._xd**2)/(2*self._td)
                    d[self._n] = - self._xd*math.exp(-0.5*cur_t) \
                        + self._xd**2/2*(u_prev[self._n]/self._td + 0.5*math.exp(-0.5*cur_t)*math.sin(self._xd*self._n))
                
                progon = TRIDIAG_SOLVER(a,b,c,d)
                
                try:
                    u = progon.solve()
                except Exception as e:
                    print(e)
                    sys.exit()

                u_prev = deepcopy(u)

                self._write_res(u, cur_t, j)

if __name__ == "__main__":
    solver = PARAB_SOLVER(scheme_type=3, approx_type=1)
    solver.solve()
    visualise(path=DATA_PATH)