#!/usr/bin/python
import os

xl = -3.99
xr = 3.99
yl = -1.99
yr = 1.99
zl = -1.99
zr = 1.99

#size = [20]
size = [20, 40, 80, 160]

t_start = 0
t_finish = 1

#t_step = [0.1, 0.05, 0.025, 0.01, 0.008, 0.005, 0.0025, 0.001, 0.0008, 0.0005, 0.00025, 0.0001]
t_step = [0.0001]

surface = 'E'

equation = 2
beta = 0

accuracy = 2


def main(method, mib_method):
    data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    if not os.path.exists(data):
        os.makedirs(data)
    if (len(t_step) != 1):
        for i in range(len(t_step)):
            file_name = os.path.join(data, "data" + str(i + 1) + ".txt")
            f = open(file_name, "w+")
            write_txt(f, method, mib_method)
            f.write("nx = %s \n" % size[0])
            f.write("ny = %s \n" % size[0])
            f.write("nz = %s \n" % size[0])
            f.write("t_step = %s \n" % t_step[i])
            f.close()
    elif (len(size) != 1):
        for i in range(len(size)):
            file_name = os.path.join(data, "data" + str(i + 1) + ".txt")
            f = open(file_name, "w+")
            write_txt(f, method, mib_method)
            f.write("nx = %s \n" % size[i])
            f.write("ny = %s \n" % size[i])
            f.write("nz = %s \n" % size[i])
            f.write("t_step = %s \n" % t_step[0])
            f.close()
    else:
        print('Error!')
        os._exit()


def write_txt(f, method, mib_method):
    f.write("xl = %s \n" % xl)
    f.write("xr = %s \n" % xr)
    f.write("yl = %s \n" % yl)
    f.write("yr = %s \n" % yr)
    f.write("zl = %s \n" % zl)
    f.write("zr = %s \n" % zr)

    f.write("t_start = %s \n" % t_start)
    f.write("t_finish = %s \n" % t_finish)

    f.write("surface = %s \n" % surface)
    f.write("method = %s \n" % method)
    f.write("equation = %s \n" % equation)
    f.write("beta = %s \n" % beta)
    f.write("mib_method = %s \n" % mib_method)
    f.write("accuracy = %s \n" % accuracy)


if __name__ == '__main__':
    #method = ['A', 'T', 'I', 'C']
    #mib_method = [1, 2]
    method = ['A']
    mib_method = [1]
    for m in method:
        for mib_m in mib_method:
            main(m, mib_m)
