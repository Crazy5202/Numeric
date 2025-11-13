import os

os.environ["QT_LOGGING_RULES"] = "qt.qpa.wayland=false;qt.qpa.socketnotifier=false"

import math
from copy import deepcopy
import sys

from visual import visualise

data_folder = "results"
DATA_PATH = os.path.join(os.path.split(os.path.realpath(__file__))[0], data_folder)
os.makedirs(DATA_PATH, exist_ok=True)

# Вариант 7
class ELLIP_SOLVER:
    def __init__(self, saving_path, x_steps = 20, y_steps = 20):
        """
        Папка сохранения результатов saving_path (не должно быть других .txt).
        
        Args:

            x_steps (int): Разбиение по x-координате.

            y_steps (int): Разбиение по y-координате.
        """
        
        self._path = saving_path
        
        self._nx = x_steps
        self._xd = math.pi / 2 / self._nx

        self._ny = y_steps
        self._yd = math.pi / 2 / self._ny

        self._start_left = lambda y: math.cos(y)
        self._start_right = lambda y: 0*y
        self._start_bottom = lambda x: math.cos(x)
        self._start_top = lambda x: 0*x

        self._true_sol = lambda x, y: math.cos(x)*math.cos(y)

    def _cleanup_dir(self):
        """
        Очистить папку результатов от "*.txt".
        """
        txt_files = [f for f in os.listdir(self._path) if f.endswith('.txt')]
    
        for file in txt_files:
            file_path = os.path.join(self._path, file)
            os.remove(file_path)
            # print(f"Removed: {file_path}")

    def _write_res(self, u:list[list[float]], ts: list[list[float]], iter:int, pogr: float):
        """
        Записать время, вычисленное и точное значение для каждой точки в файл.
        """
        
        ind_path = self._path + '/' + str(iter) + ".txt"
        with open(ind_path, "w") as f:
            f.write(str(iter) + '\n')
            for j in range (self._ny+1):
                for i in range (self._nx+1):
                    f.write(str(j*self._yd) + ' ' + str(i*self._xd) + ' ' +  str(u[j][i]) + ' ' + str(ts[j][i]) + '\n')

        pogr_path = self._path + "/p.txt"
        with open(pogr_path, "a") as f:
            f.write(str(iter) + ' ' + str(pogr) + '\n')

    def _pogr_step(self, u: list[float], ts: list[float]):
        """
        Рассчитать погрешность для данного времени.
        """
        pogr = 0
        for j in range (self._ny+1):
            for i in range (self._nx+1):
                pogr = max(pogr, abs(u[j][i]-ts[j][i]))
        return pogr
        
    def _post_solution(self, u:list[float], iter:float):
        """
        Сделать действия после шага решения.
        """
        cur_true = []
        for j in range(self._ny+1):
            cur_true.append([self._true_sol(self._xd*i, self._yd*j) for i in range(self._nx+1)])
                
        pogr = self._pogr_step(u, cur_true)
        self._write_res(u, cur_true, iter, pogr)

    def solve(self, scheme_type: int = 1, eps: float = 0.1, interpol: bool = True, w: int = 1.5):
        """
        Вычислить и сохранить решение.

        Args:

            scheme_type (str): Тип схемы. 1 - МПИ, 2 - Зейдель, 3 - Верхняя релаксация.
            
            w (float): параметр верхней релаксации (1;2).

            interpol (bool): Производить ли интерполяцию по границам.
        """

        if scheme_type not in [1,2,3] or w <= 1 or w >= 2:
            raise ValueError("Неверно указаны параметры решателя!")
        
        self._cleanup_dir()

        u = []

        # Нижняя граница
        u.append([self._start_bottom(self._xd*i) for i in range(self._nx+1)])

        # Всё кроме нижней и верхней границы
        for _ in range(1, self._ny):
            u.append([0]*(self._nx+1))

        # Верхняя граница
        u.append([self._start_top(self._xd*i) for i in range(self._nx+1)])

        # Левая и правая граница
        for j in range(self._ny+1):
            u[j][0] = self._start_left(self._yd*j)
            u[j][self._nx] = self._start_right(self._yd*j)

        # Интерпроляция
        if (interpol == True):
            for j in range(1, self._ny):
                for i in range(1, self._nx):
                    f1 = (u[0][0] * (self._nx+1-i) + u[0][self._nx] * i) / (self._nx + 1)
                    f2 = (u[self._ny][0] * (self._nx+1-i) + u[self._ny][self._nx] * i) / (self._nx + 1)
                    u[j][i] = (f1 * (self._ny+1-j) + f2 * j) / (self._ny + 1)

        u_prev = deepcopy(u)

        cur_iter = 1
        end_crit = 1e9
        
        if scheme_type == 1:

            while (end_crit > eps):
            
                for j in range(1, self._ny):

                    for i in range(1, self._nx):

                        u[j][i] = (u_prev[j][i+1] + u_prev[j][i-1]) / (self._xd**2) + (u_prev[j+1][i] + u_prev[j-1][i]) / (self._yd**2)

                        u[j][i] /= (2 / self._xd**2) + (2 / self._yd**2) - 2

                self._post_solution(u, cur_iter)

                end_crit = self._pogr_step(u_prev, u)

                u_prev = deepcopy(u)

                cur_iter += 1

            print(f"Закончили вычисления на итерации {cur_iter-1} с значением критерия остановки {end_crit}")

        # else:
        #     param = 1
        #     if scheme_type == 3:
        #         param = 0.5

        #     for j in range (1, self._t_steps+1):
        #         cur_t = j*self._td
        #         a = [0]*(self._n+1); b = [0]*(self._n+1); c = [0]*(self._n+1); d = [0]*(self._n+1)
                    
        #         for i in range(1, self._n):
        #             cur_x = self._xd*i

        #             a[i] = -self._td*param
        #             b[i] = self._xd**2 + 2*self._td*param
        #             c[i] = -self._td*param
        #             d[i] = u[i]*(self._xd**2) + self._td*(1-param)*(u[i+1]-2*u[i]+u[i-1]) \
        #                 + 0.5*(self._xd**2)*self._td*math.sin(cur_x)*(param*math.exp(-0.5*cur_t)+(1-param)*math.exp(-0.5*(cur_t-self._td)))

        #         if approx_type == 1:
        #             b[0] = -1
        #             c[0] = 1
        #             d[0] = self._xd*math.exp(-0.5*cur_t)

        #             a[self._n] = -1
        #             b[self._n] = 1
        #             d[self._n] = -self._xd*math.exp(-0.5*cur_t)

        #         elif approx_type == 2:
        #             b[0] = -3 - a[1] / self._td / param
        #             c[0] = 4 - b[1] / self._td / param
        #             d[0] = 2*self._xd*math.exp(-0.5*cur_t) - d[1] / self._td / param

        #             a[self._n] = -4 + b[self._n-1] / self._td / param
        #             b[self._n] = 3 + c[self._n-1] / self._td / param
        #             d[self._n] = -2*self._xd*math.exp(-0.5*cur_t) + d[self._n-1] / self._td / param

        #         elif approx_type == 3:
        #             b[0] = 1 + (self._xd**2)/(2*self._td)
        #             c[0] = -1
        #             d[0] = -self._xd*math.exp(-0.5*cur_t)

        #             a[self._n] = -1
        #             b[self._n] = 1 + (self._xd**2)/(2*self._td)
        #             d[self._n] = - self._xd*math.exp(-0.5*cur_t) \
        #                 + self._xd**2/2*(u_prev[self._n]/self._td)
                
        #         # progon = TRIDIAG_SOLVER(a,b,c,d)
                
        #         # # try:
        #         # #     u = progon.solve()
        #         # # except Exception as e:
        #         # #     print(e)
        #         # #     sys.exit(1)

        #         # u = progon.solve()

        #         u_prev = deepcopy(u)

        #         self._post_solution(u, cur_t, j)

if __name__ == "__main__":
    solver = ELLIP_SOLVER(saving_path=DATA_PATH, x_steps = 10, y_steps = 10)

    for i in range (1,2):
        solver.solve(i, eps = 0.01)
        visualise(path=DATA_PATH)