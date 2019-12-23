from sympy import *
import csv
from math import *
import random

x,y,z = symbols('x y z')

def get_file_name(func):
    if func == poli:
        file_name 
    elif func == T:
        pass
    elif func == U:
        pass

def poli(n, x):
    t = n
    k = 1/((2**n)*factorial(n))*(x**2-1)**n
    for i in range(t):
        k = diff(k,x)
    return k

def T(n, x):
    v,l = x,1
    if(n==0):
        return l
    for m in range(1,n):
        t = v
        v = 2*x*t-l
        l = t
    return v

def U(n, x):
    v,l = 2*x,1
    if(n==0):
        return l
    for m in range(1,n):
        t = v
        v = 2*x*t-l
        l = t
    return v

def integr(n,m):
    integral = integrate(n*m,(x,-1,1)) 
    integral = integrate(integral,(y,-1,1))
    integral = integrate(integral,(z,-1,1))
    return integral

def sys_pol(n,funk = poli):
    sys = []
    for i in range(n):
        sys.append(funk(i))
    return sys

def sys_diff(k):
    tam = []
    for i in k:
        tam.append(rotor(i))
    return tam

def matrix(x):
    matrix = []
    t = sys_diff(x)
    for i in range(len(x)):
        k = []
        for j in range(len(x)):
            k.append(integr(t[j],x[i])/integr(x[i],x[i]))
        matrix.append(k)
    return matrix

def poli_vector(n, func = None):
    vect = []
    n += 1
    for i in range(0, n):
        for j in range (0, n):
            for k in range(0, n):
                if func is None:
                    function = x**(i)*y**(j)*z**(k)
                else:
                    function = func(i, x)*func(j, y)*func(k, z)
                vect.append(function)
    return vect

def poli_vector_ijk(poli_vector):
    vect = list()
    n = len(poli_vector)
    for i in range(0, 3):
        for j in range(0, n):
            mini_vect = [x*0,x*0,x*0]
            mini_vect[i] = (poli_vector[j])
            vect.append(mini_vect)
    return vect

def rot(vect):
    rotor = list()
    for i in range(0, len(vect)):
        roti = diff(vect[i][2], y) - diff(vect[i][1], z)
        rotj = diff(vect[i][0], z) - diff(vect[i][2], x)
        rotk = diff(vect[i][1], x) - diff(vect[i][0], y)
        rotor.append([(roti), (rotj), (rotk)])
    return rotor

def evaluate_als(Fn):
    rotor = rot(Fn)
    N = len(Fn)
    L_max = N
    S_max = N//3
    a = []
    for j in range(S_max):
        aj = []
        for i in range(L_max):
            aji = (integr(rotor[i][0], Fn[j][0])/integr(Fn[j][0], Fn[j][0])) 
            aj.append(aji)
        a.append(aj)
        print(aj)
    print("-----------")
    S_min= N//3
    S_max = 2*N//3
    L_max = N
    for j in range(S_min, S_max):
        aj = []
        for i in range(L_max):
            aji = (integr(rotor[i][1], Fn[j][1])/integr(Fn[j][1], Fn[j][1])) 
            aj.append(aji)
        a.append(aj)
        print(aj)
    print("----------")
    S_min= 2*N//3
    S_max = N
    L_max = N
    for j in range(S_min, S_max):
        aj = []
        for i in range(L_max):
            aji = (integr(rotor[i][2], Fn[j][2])/integr(Fn[j][2], Fn[j][2])) 
            aj.append(aji)
        a.append(aj)
        print(aj)
    return a

def get_y(a, x):
    M = len(a[0])
    N = len(a)
    y = []
    for i in range(N):
        yi = 0
        for j in range(M):
            yi += float(a[i][j])*x[j] 
        y.append(yi)
    return y

def make_x(N):
    random.seed()
    N = ((N+1)**3)*3
    x = [0 for i in range(N)]
    for i in range(0, N):
        if i%3 == 0: 
            x[i] = random.randint(0, 10)
    return x

def make_x_ijk(base, x):
    xijk = []
    for i in range(len(x)):
        xijk_1 = []
        for j in range(len(base[i])):
            xijk_1.append(x[i] * base[i][j])
        xijk.append(xijk_1)
    return xijk

def simplify(xijk):
    xi = 0
    xj = 0
    xk = 0
    for i in range(len(xijk)):
       xi += xijk[i][0]
       xj += xijk[i][1]
       xk += xijk[i][2]
    return [xi, xj, xk]

def save_csv(a, file_name):
    N = len(a)
    with open(file_name, mode='w+') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter = ',', quotechar='"',
                quoting=csv.QUOTE_MINIMAL)
        for i, row in enumerate(a):
            if i%5 == 0:
                print('{}-s of {} rows writed'.format(i, N))
            csv_writer.writerow(row)
        
def read_csv(file_name):
    table = []
    with open(file_name, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = ',')
        for row in csv_reader:
            table.append(row)
        return table

def test(N, func, file_name='', save_matrix=True, details=False):
    N = N
    eps = 10**(-10)
    base_vector = poli_vector(N, func = func)
    base_vector_ijk = poli_vector_ijk(base_vector)
    
    x_ = make_x(N)
    #The process of evaluating rotor directly 
    rotor = rot(base_vector_ijk)
    fake_r = [0.0*x,0.0*x,0.0*x] 
    for i,j in enumerate(x_):
        if j != 0:
            fake_r[0] += rotor[i][0]*float(j)
            fake_r[1] += rotor[i][1]*float(j)
            fake_r[2] += rotor[i][2]*float(j)
    
    if details:
        for i in range(len(base_vector_ijk)):
            print("{}th basic vector and rotor, ".format(i), base_vector_ijk[i],"  ",rotor[i])

    #The process of evaluationg rotor via matrix
    if save_matrix == True:
        a = evaluate_als(base_vector_ijk)
        save_csv(a, file_name)
    a = read_csv(file_name)
    y_ = get_y(a, x_) 
    #Striping 
    for j,i in enumerate(y_):
        if abs(i) < eps:y_[j]= 0

    #Printing results
    simplify_y = simplify(make_x_ijk(base_vector_ijk, y_))
    if details: print("Vector Y in Legandr Basis:", y_) 
    if details: print("Vector Y in xyz Basis:", simplify_y) 
    if fake_r == simplify_y:
        print('Test passed successfully') 
    else:
        print('Test failed')
    print('Direct rotor', fake_r)
    print('Matrix rotor', simplify_y)

def main():
    test(2, poli, file_name='rotor_matrix_legandre.csv', save_matrix=False,
            details=False)
    test(1, T, file_name='rotor_T_test.csv', save_matrix=False, 
            details=False)
    test(1, U, file_name='rotor_U_test.csv', save_matrix=False, 
            details=False)

if __name__ == '__main__':
    main()
