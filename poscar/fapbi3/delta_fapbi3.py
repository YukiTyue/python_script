#!/usr/bin/env Python3
# coding=utf-8
"""

@author: CuanRan(Meiyukichan@163.com)
@file: delta_fapbi3.py
@time: 2020/12/28 20:16

"""

import os
import copy
import math


def poscar_input(file):  # 一个函数是可以返还多个值的！！！
    # lattice_vector:
    lattice_vector = []
    # scale_factor = []   # poscar里面的缩放系数，比如1，因为读取行line虽然是列表，但是里面的元素是字符串，所以需要定义为列表
    poscar_tempt = []
    scale_and_lat = []
    other_lines = []
    coordinate_atom = []
    coord_num_atom = []
    # coordinate_tempt = []  如果不对这个空列表做操作，则不需要声明，可以后面直接赋值
    if os.path.exists(file):  # os.path.exists(path) 判断系统路径是否存在
        poscar = open(file, 'r')  # 使用 open() 方法一定要保证关闭文件对象，即调用 close() 方法。open() 函数常用形式是接收两个参数：文件名(file)和模式(mode)。
    else:
        raise IOError('POSCAR does not exist!')  # IOError主要是指要打开的文件不存在的错误提示。
        # 一旦执行了raise语句，raise后面的语句将不能执行。
        # 所以循环可以写在if判断完成之后，增加效率！！！
    for line in poscar.readlines():  # readlines()方法读取整个文件所有行，保存在一个列表(list)变量中，每行作为一个元素，但读取大文件会比较占内存
        # 依次读取每行
        poscar_tempt.append(line)
    poscar.close()
    scale_factor = poscar_tempt[1].split()   # 因为缩放系数在poscar第二行，所以是poscar_tempt[1]，因为只有一个字符串，所以split的时候不需要加空格
    lattice = poscar_tempt[2:5]   # 晶格矢量在POSCAR的第3-5行，所以为[2:5],但是这时的lattice是一行的内容，包括空格和换行符
    for axis in lattice:
        lattice_vector.append(axis.split())  # 如果不用split，将会保留换行符\n，这里相当于取了整个第二行为一个字符串,然后对这个字符串切片
        #                                         但是直接用split(' ')也不行，因为不知道两个数之间有多少空格，而且还会多出\n
    scale_and_lat.append(scale_factor)
    scale_and_lat.append(lattice_vector)

    # other_lines:
    line1 = poscar_tempt[0]  # POSCAR第一行
    line6 = poscar_tempt[5].split()  # POSCAR第六行，元素！
    line7 = poscar_tempt[6].split()  # POSCAR第七行，元素个数！
    if poscar_tempt[7].strip().upper().startswith('S'):
        # 如果POSCAR第七行去掉首尾（空格和尾部的换行符）后的字符串的第一个字母大写后为‘S’则为真！
        # 7 Selective dynamics
        # 8 Cartesian 的情况：
        line8 = poscar_tempt[7]
        line9 = poscar_tempt[8]
    else:
        line8 = ''
        line9 = poscar_tempt[7]
    other_lines.append(line1)
    other_lines.append(line6)
    other_lines.append(line7)
    other_lines.append(line8)
    other_lines.append(line9)
    if len(other_lines[3]) < 2:
        del other_lines[3]       # 后面用len(other_lines)就可以重新确定index

    # coordinates_of_atoms
    atom_composition = poscar_tempt[6]
    atom_composition = atom_composition.split()
    atom_number = 0   # 总的原子数
    for a_i in range(len(atom_composition)):
        atom_number += int(atom_composition[a_i])
    if poscar_tempt[7].strip().upper().startswith('S'):
        coordinate_tempt = poscar_tempt[9:9 + atom_number]
    else:
        coordinate_tempt = poscar_tempt[8:8 + atom_number]   # 是赋值而不是引用，因为没有alist = blist
    for atom in coordinate_tempt:
        coordinate_atom.append(atom.split())
    coord_num_atom.append(atom_number)
    coord_num_atom.append(coordinate_atom)

    return scale_and_lat, other_lines, coord_num_atom


def cord_elem_reset(dict_plus, dict_moto, cord_num_moto, cord_num_plus):
    # 这是一个index的字典
    dict_plus_index = {}
    atom_index = 0
    keys_plus = list(dict_plus.keys())
    for plus_elements in keys_plus:
        dict_plus_index[plus_elements] = atom_index
        atom_index += int(dict_plus[plus_elements])

    for elements in keys_plus:
        # 更新后的pbi3中元素的后移index
        dict_moto_taildex = {}
        atom_index = 0
        keys_moto = list(dict_moto.keys())
        for moto_elements in keys_moto:
            atom_index += int(dict_moto[moto_elements])
            dict_moto_taildex[moto_elements] = atom_index
        if elements in keys_moto:
            dict_moto[elements] = str(int(dict_moto[elements]) + int(dict_plus[elements]))
            for atom_elments in range(dict_plus_index[elements], dict_plus_index[elements] + int(dict_plus[elements])):
                cord_num_moto[1].insert(dict_moto_taildex[elements], cord_num_plus[1][atom_elments])
        else:
            dict_moto[elements] = dict_plus[elements]
            for atom_elments in range(dict_plus_index[elements], dict_plus_index[elements] + int(dict_plus[elements])):
                cord_num_moto[1].append(cord_num_plus[1][atom_elments])

    return cord_num_moto, dict_moto
            

def poscar_output(scal_lattice, other_line, cordin_num):
    # 写入文件，说白了就是把这个文件做成一个字符串的列表，然后循环依次写入。
    # writing file!
    poscar = [other_line[0], '   ' + scal_lattice[0][0] + '\n']    # 相当于声明列表list=[a,b]
    for vector in range(3):
        line_tempt = ''
        for cord_index in range(3):
            line_tempt += '    '
            line_tempt += scal_lattice[1][vector][cord_index].rjust(13)
        line_tempt += '\n'
        poscar.append(line_tempt)
        # poscar.append('    ' + '   '.join(scal_lattice[1][vector]) + '\n')
    poscar.append('   ' + '    '.join(other_line[1]) + '\n')
    poscar.append('    ' + '    '.join(other_line[2]) + '\n')
    poscar.append(other_line[3])
    for atom_num in range(cordin_num[0]):
        cord_tempt = ''
        for cord_index in range(3):
            cord_tempt += '     '
            cord_tempt += cordin_num[1][atom_num][cord_index].rjust(12)
        cord_tempt += '\n'
        poscar.append(cord_tempt)
        # poscar.append('    ' + '   '.join(cordin_num[1][atom_num]) + '\n')

    return poscar


def rotation_trans_2d(origin, ropoints, roangle):
    roangle = roangle * math.pi / 180
    rotation_matrix = [[math.cos(roangle), -math.sin(roangle)], [math.sin(roangle), math.cos(roangle)]]
    for cordinate in range(len(ropoints)):
        reduce_x_0 = float(ropoints[cordinate][0]) - float(origin[0])
        reduce_y_0 = float(ropoints[cordinate][1]) - float(origin[1])
        reduce_x_1 = reduce_x_0 * rotation_matrix[0][0] + reduce_y_0 * rotation_matrix[0][1]
        reduce_y_1 = reduce_x_0 * rotation_matrix[1][0] + reduce_y_0 * rotation_matrix[1][1]
        rotation_x = reduce_x_1 + float(origin[0])
        rotation_y = reduce_y_1 + float(origin[1])
        ropoints[cordinate][0] = '{:.9f}'.format(rotation_x)
        ropoints[cordinate][1] = '{:.9f}'.format(rotation_y)

    return ropoints


if __name__ == '__main__':
    poscar_c = 'D_FA2'
    poscar_fa = 'A-FA2'
    poscar_pbi3 = 'D-FA-BX3'
    scal_lat_c, other_lin_c, cord_num_c = poscar_input(poscar_c)
    scal_lat_fa, other_lin_fa, cord_num_fa = poscar_input(poscar_fa)
    scal_lat_pbi3, other_lin_pbi3, cord_num_pbi3 = poscar_input(poscar_pbi3)
    scal_lat_fapbi3, other_lin_fapbi3, cord_num_fapbi3 = scal_lat_pbi3, other_lin_pbi3, cord_num_pbi3
    cord_num_1 = copy.deepcopy(cord_num_fa)
    cord_num_2 = copy.deepcopy(cord_num_fa)
    for i in range(cord_num_fa[0]):
        tempt1_x = float(cord_num_fa[1][i][0]) - float(cord_num_fa[1][0][0]) + float(cord_num_c[1][0][0])
        tempt1_y = float(cord_num_fa[1][i][1]) - float(cord_num_fa[1][0][1]) + float(cord_num_c[1][0][1])
        tempt2_x = float(cord_num_fa[1][i][0]) - float(cord_num_fa[1][0][0]) + float(cord_num_c[1][1][0])
        tempt2_y = float(cord_num_fa[1][i][1]) - float(cord_num_fa[1][0][1]) + float(cord_num_c[1][1][1])
        cord_num_1[1][i][0] = '{:.9f}'.format(tempt1_x)   # 输出类型为str
        cord_num_1[1][i][1] = '{:.9f}'.format(tempt1_y)
        cord_num_1[1][i][2] = cord_num_c[1][0][2]
        cord_num_2[1][i][0] = '{:.9f}'.format(tempt2_x)
        cord_num_2[1][i][1] = '{:.9f}'.format(tempt2_y)
        cord_num_2[1][i][2] = cord_num_c[1][1][2]
    cord_num_pbi3[0] += cord_num_1[0] + cord_num_2[0]

    cord_num_1[1] = rotation_trans_2d(cord_num_1[1][0], cord_num_1[1], -120)
    cord_num_2[1] = rotation_trans_2d(cord_num_2[1][0], cord_num_2[1], 60)

    dict_fa1 = dict(zip(other_lin_fa[1], other_lin_fa[2]))
    dict_fa2 = dict(zip(other_lin_fa[1], other_lin_fa[2]))
    dict_pbi3 = dict(zip(other_lin_pbi3[1], other_lin_pbi3[2]))

    cord_num_pbi3, dict_pbi3 = cord_elem_reset(dict_fa1, dict_pbi3, cord_num_pbi3, cord_num_1)
    cord_num_pbi3, dict_pbi3 = cord_elem_reset(dict_fa2, dict_pbi3, cord_num_pbi3, cord_num_2)

    elements_pbi3 = list(dict_pbi3.keys())
    elements_comp = list(dict_pbi3.values())
    other_lin_pbi3[1] = elements_pbi3
    other_lin_pbi3[2] = elements_comp
    other_lin_pbi3[0] = 'Delta-FAPbI3  generated by delta_fapbi3.py\n'

    poscar_fapbi3 = poscar_output(scal_lat_fapbi3, other_lin_fapbi3, cord_num_fapbi3)
    poscar_fin_fapbi3 = open('POSCAR_DFA', 'w+')
    for lines in poscar_fapbi3:
        poscar_fin_fapbi3.write(lines)
    poscar_fin_fapbi3.close()

    print(
        '''
             へ　　　　　／|
        　　/＼7　　　 ∠＿/
        　 /　│　　 ／　／
        　│　Z ＿,＜　／　　 /`ヽ
        　│　　　　　ヽ　　 /　　〉
        　 Y　　　　　`　 /　　/
        　ｲ●　､　●　　⊂⊃〈　　/
        　()　 へ　　　　|　＼〈
        　　>ｰ ､_　 ィ　 │ ／／
        　 / へ　　 /　ﾉ＜| ＼＼
        　 ヽ_ﾉ　　(_／　 │／／
        　　7　　　　　　　|／
        　　＞―r￣￣`ｰ―＿_|
        '''
    )  # '''可以在print里面换行输出!
    print('Pika chu!')
