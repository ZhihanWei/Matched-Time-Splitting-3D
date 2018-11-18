#!/usr/bin/python
import os
import argparse

xl = -3.99
xr = 3.99
yl = -1.99
yr = 1.99
zl = -1.99
zr = 1.99

t_start = 0
t_finish = 1

surface = 'E'

equation = 2
beta = 0

accuracy = 2


def main(method, mib_method, t_step, size):
    data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    if not os.path.exists(data):
        os.makedirs(data)
    if (len(t_step) != 1):
        for i in range(len(t_step)):
            file_name = os.path.join(data, "data" + str(i + 1) + ".txt")
            f = open(file_name, "w+")
            write_txt(f, method, mib_method)
            f.write("t_step = %s \n \n" % t_step[i])
            f.write("nx = %s \n" % size[0])
            f.write("ny = %s \n" % size[0])
            f.write("nz = %s \n" % size[0])
            f.close()
    elif (len(size) != 1):
        for i in range(len(size)):
            file_name = os.path.join(data, "data" + str(i + 1) + ".txt")
            f = open(file_name, "w+")
            write_txt(f, method, mib_method)
            f.write("t_step = %s \n \n" % t_step[0])
            f.write("nx = %s \n" % size[i])
            f.write("ny = %s \n" % size[i])
            f.write("nz = %s \n" % size[i])
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
    f.write("zr = %s \n \n" % zr)

    f.write("surface = %s \n" % surface)
    f.write("method = %s \n" % method)
    f.write("equation = %s \n" % equation)
    f.write("beta = %s \n" % beta)
    f.write("mib_method = %s \n" % mib_method)
    f.write("accuracy = %s \n \n" % accuracy)

    f.write("t_start = %s \n" % t_start)
    f.write("t_finish = %s \n" % t_finish)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = "create input data")
    parser.add_argument('-t', '--temporal_method', type = str, required = True, help = 'Temporal method, A: ADI; T: TS; I: LODIE; C: LODCN')
    parser.add_argument('-s', '--spatial_method', type = str, required = True, help = 'Spatial method, 1: MIBL1; 2: MIBL2')
    parser.add_argument('-f', '--focus', type = str, required = True, help = 'Focus on accuracy or speed, a: accuracy; s: speed') 
    args = parser.parse_args()

    if (args.focus == 'a'):
        t_step = [0.0001]
        size = [20, 40, 80, 160]
    elif (args.focus == 's'):
        t_step = [0.1, 0.05, 0.025, 0.01, 0.008, 0.005, 0.0025, 0.001, 0.0008, 0.0005, 0.00025, 0.0001]
        size = [160]
    else:
        print('Error!')
        os._exit()

    main(args.temporal_method, args.spatial_method, t_step, size)
