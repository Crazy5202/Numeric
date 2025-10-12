from tridiag import TRIDIAG_SOLVER

a = [0, -1, -9, -1, 9]
b = [-6, 13, -15, -7, -18]
c = [5, 6, -4, 1, 0]
d = [51, 100, -12, 47, -90]

print(TRIDIAG_SOLVER(a,b,c,d).solve())