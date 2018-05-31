#!/usr/bin/python
import os
import shutil
import subprocess

if __name__ == '__main__':
    root_dir = os.path.dirname(os.path.abspath(__file__))
    focus = ['a', 's']
    for f in focus:
        if f == 'a':
            folder = os.path.join(root_dir, 'Accuracy')
        elif f == 's':
            folder = os.path.join(root_dir, 'Time')
        else:
            print('Error!')
            os._Exit()

        if not os.path.exists(folder):
            os.makedirs(folder)

        file_dir = os.path.join(root_dir, 'get_results.sh')
        shutil.copy(file_dir, folder)

        temporal_method = ['A', 'I', 'C', 'T']
        spacial_method = ['1', '2']
        for t in temporal_method:
            for s in spacial_method:
                if t == 'A' and s == '1':
                    subfolder = os.path.join(folder, 'ADI-L1')
                elif t == 'A' and s == '2':
                    subfolder = os.path.join(folder, 'ADI-L2')
                elif t == 'I' and s == '1':
                    subfolder = os.path.join(folder, 'LODIE-L1')
                elif t == 'I' and s == '2':
                    subfolder = os.path.join(folder, 'LODIE-L2')
                elif t == 'C' and s == '1':
                    subfolder = os.path.join(folder, 'LODCN-L1')
                elif t == 'C' and s == '2':
                    subfolder = os.path.join(folder, 'LODCN-L2')
                elif t == 'T' and s == '1':
                    subfolder = os.path.join(folder, 'TS-L1')
                elif t == 'T' and s == '2':
                    subfolder = os.path.join(folder, 'TS-L2')
                else:
                    print('Error!')
                    os._Exit()

                source_folder = os.path.join(root_dir, 'original')
                shutil.copytree(source_folder, subfolder)
                executable = os.path.join(subfolder, 'set_paras.py')
                subprocess.call(['python', executable, '-t', t, '-s', s, '-f', f])
