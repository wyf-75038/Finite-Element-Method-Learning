import numpy as np
from scipy import linalg


def input_2d(file_name):  # 输入

    def evaluate(var, line_start, line_last):  # 把某几行数据赋给矩阵var
        j = 0
        for i in range(line_start, line_start + line_last):
            var[j] = lines[i].split(',')
            j += 1
        return var

    with open(file_name) as FInput:
        lines = FInput.readlines()

    line_len = len(lines)

    for i in range(0, line_len):
        if lines[i] == 'A\n':
            a = float(lines[i + 1])
        elif lines[i] == 'E\n':
            e = float(lines[i + 1])
        elif lines[i] == 'NodeNum\n':
            node_num = int(lines[i + 1])
        elif lines[i] == 'Node\n':
            node_line = i + 1  # 记录节点数据开始的行数
        elif lines[i] == 'StickNum\n':
            stick_num = int(lines[i + 1])
        elif lines[i] == 'Stick\n':
            stick_line = i + 1  # 记录杆件数据开始的行数
        elif lines[i] == 'ForceNode\n':
            force_node = [int(n) for n in lines[i + 1].split(',')]
        elif lines[i] == 'Force\n':
            force_line = i + 1  # 记录力数据开始的行数
        elif lines[i] == 'XCon\n':
            con_x = [int(n) for n in lines[i + 1].split(',')]
        elif lines[i] == 'YCon\n':
            con_y = [int(n) for n in lines[i + 1].split(',')]

    node = np.zeros((node_num, 2))
    stick = np.zeros((stick_num, 2))
    node = evaluate(node, node_line, node_num)
    stick = evaluate(stick, stick_line, stick_num)

    force_node_num = len(force_node)
    force_input = np.zeros((force_node_num, 2))
    force = np.zeros((node_num, 2))
    force_input = evaluate(force_input, force_line, force_node_num)
    force_node_index = [i - 1 for i in force_node]
    force[force_node_index] = force_input

    con = np.zeros((node_num, 2))
    for i in con_x:
        con[i - 1, 0] = 1
    for i in con_y:
        con[i - 1, 1] = 1

    stick = stick.astype(int)
    con = con.astype(int)
    return [a, e, node, node_num, stick, stick_num, force, con]


def stick_len_and_angle():  # 求杆件长度和角度
    stick_l = np.zeros(StickNum)
    c = np.zeros(StickNum)
    s = np.zeros(StickNum)

    for i in range(0, StickNum):
        x_sta = Node[Stick[i, 0] - 1, 0]
        y_sta = Node[Stick[i, 0] - 1, 1]
        x_end = Node[Stick[i, 1] - 1, 0]
        y_end = Node[Stick[i, 1] - 1, 1]
        stick_l[i] = ((x_end - x_sta) ** 2 +
                      (y_end - y_sta) ** 2) ** 0.5
        c[i] = (x_end - x_sta) / stick_l[i]
        s[i] = (y_end - y_sta) / stick_l[i]
    return [stick_l, c, s]


def k_element():
    k = [0] * StickNum
    for i in range(0, StickNum):
        k[i] = np.array([[C[i] ** 2, C[i] * S[i], -C[i] ** 2, -C[i] * S[i]],
                         [C[i] * S[i], S[i] ** 2, -C[i] * S[i], -S[i] ** 2],
                         [-C[i] ** 2, -C[i] * S[i], C[i] ** 2, C[i] * S[i]],
                         [-C[i] * S[i], -S[i] ** 2, C[i] * S[i], S[i] ** 2]])
        k[i] = A * E / L[i] * k[i]
    return k


def k_total(k_element):
    k_total = np.zeros((NodeNum * 2, NodeNum * 2))

    for i in range(0, StickNum):  # 计算总刚
        ind1 = Stick[i, 0]
        ind2 = Stick[i, 1]
        k_index = [2 * ind1 - 2, 2 * ind1 - 1, 2 * ind2 - 2, 2 * ind2 - 1]

        jj = 0
        for j in k_index:
            kk = 0
            for k in k_index:
                k_total[j, k] += k_element[i][jj, kk]
                kk += 1
            jj += 1
    return k_total


def u_solve():
    index_no_con = np.argwhere(Con == 0).ravel()

    f_cut = F[index_no_con]
    k_total_cut = KTotal[index_no_con]
    k_total_cut = k_total_cut[:, index_no_con]

    u_result = linalg.solve(k_total_cut, f_cut)

    U[index_no_con] = u_result
    urs = U.reshape(NodeNum, 2)
    return urs


def f_solve():
    f = KTotal @ U
    f_rs = f.reshape(NodeNum, 2)
    return f_rs


def sigma_solve():
    sigma = np.zeros(StickNum)

    for i in range(0, StickNum):
        ind1 = Stick[i, 0]
        ind2 = Stick[i, 1]
        u_rs_cut = np.array([[U[2 * ind1 - 2]], [U[2 * ind1 - 1]],
                             [U[2 * ind2 - 2]], [U[2 * ind2 - 1]]])
        cap = E / L[i] * np.array([-C[i], -S[i], C[i], S[i]])
        sigma[i] = cap @ u_rs_cut
    return sigma


def output_2d(file_name, dig):
    with open(file_name, "w") as FOutput:
        FOutput.write('2D Bracing System Results\n')
        FOutput.write('U/mm\n')
        FOutput.write(str((UResult * 1000).round(dig)))
        FOutput.write('\n')
        FOutput.write('F/kN\n')
        FOutput.write(str((FResult / 1000).round(dig)))
        FOutput.write('\n')
        FOutput.write('sigma/MPa\n')
        FOutput.write(str((Sigma / 10 ** 6).round(dig)))
        FOutput.write('\n')
    return 1


FileNameInput = input('请输入输入文件名\n')
Dig = int(input('请输入计算结果保留小数位数\n'))
FileNameOutput = 'Result_of_' + FileNameInput

[A, E, Node, NodeNum, Stick, StickNum, Force, Con] = input_2d(FileNameInput)

[L, C, S] = stick_len_and_angle()

KElement = k_element()

KTotal = k_total(KElement)

Con = Con.ravel()

U = np.zeros(NodeNum * 2)
F = Force.ravel()

UResult = u_solve()

FResult = f_solve()

Sigma = sigma_solve()

output_2d(FileNameOutput, Dig)

print('Done Successfully at', FileNameOutput)
