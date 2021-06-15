import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import cm


def input_2d(file_name):  # 输入

    def evaluate(var, line_start, line_last):  # 把某几行数据赋给矩阵var
        j = 0
        for i in range(line_start, line_start + line_last):
            var[j] = lines[i].split(',')  # 半角逗号分隔各元素
            j += 1
        return var

    with open(file_name) as FInput:
        lines = FInput.readlines()

    line_num = len(lines)

    for i in range(0, line_num):
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
            con_x = [int(n) for n in lines[i + 1].split(',')]  # x方向约束信息
        elif lines[i] == 'YCon\n':
            con_y = [int(n) for n in lines[i + 1].split(',')]  # y方向约束信息

    node = np.zeros((node_num, 2))
    stick = np.zeros((stick_num, 2))
    node = evaluate(node, node_line, node_num)
    stick = evaluate(stick, stick_line, stick_num)

    force_node_num = len(force_node)  # 有外载荷的节点数
    force_input = np.zeros((force_node_num, 2))
    force = np.zeros((node_num, 2))
    force_input = evaluate(force_input, force_line, force_node_num)
    force_node_index = [i - 1 for i in force_node]  # 每个元素减一
    force[force_node_index] = force_input

    con = np.zeros((node_num, 2))
    for i in con_x:
        con[i - 1, 0] = 1
    for i in con_y:
        con[i - 1, 1] = 1

    stick = stick.astype(int)  # 转化为整型
    con = con.astype(int)  # 转化为整型
    con = con.ravel()
    return [a, e, node, node_num, stick, stick_num, force_node, force, con]


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


def k_element():  # 单元刚度
    k = [0] * StickNum  # StickNum个元素的全零列表
    for i in range(0, StickNum):
        k[i] = np.array([[C[i] ** 2, C[i] * S[i], -C[i] ** 2, -C[i] * S[i]],
                         [C[i] * S[i], S[i] ** 2, -C[i] * S[i], -S[i] ** 2],
                         [-C[i] ** 2, -C[i] * S[i], C[i] ** 2, C[i] * S[i]],
                         [-C[i] * S[i], -S[i] ** 2, C[i] * S[i], S[i] ** 2]])
        k[i] = A * E / L[i] * k[i]
    return k


def k_total(k_element):  # 总刚度
    k_total = np.zeros((NodeNum * 2, NodeNum * 2))

    for i in range(0, StickNum):
        ind1 = Stick[i, 0]
        ind2 = Stick[i, 1]

        k_index = [2 * ind1 - 2, 2 * ind1 - 1, 2 * ind2 - 2, 2 * ind2 - 1]
        # 某个单元刚度矩阵每行每列在总刚度矩阵中的位置索引

        jj = 0
        for j in k_index:
            kk = 0
            for k in k_index:
                k_total[j, k] += k_element[i][jj, kk]
                kk += 1
            jj += 1
    return k_total


def u_solve():  # 求解位移
    index_no_con = np.argwhere(Con == 0).ravel()  # 找到没有约束的节点的索引

    f_input_vector = Force.ravel()  # 展开为向量格式
    f_cut = f_input_vector[index_no_con]  # 缩减力向量
    k_total_cut = KTotal[index_no_con]  # 缩减总刚度矩阵向量的行
    k_total_cut = k_total_cut[:, index_no_con]  # 缩减总刚度矩阵向量的列

    u_result = linalg.solve(k_total_cut, f_cut)  # 求解线性方程组

    u_vector = np.zeros(NodeNum * 2)
    u_vector[index_no_con] = u_result  # 将结果填写在没有约束的位置
    u_matrix = u_vector.reshape(NodeNum, 2)  # 将u_vector整形为2列的矩阵
    return u_matrix, u_vector


def f_solve():  # 求解力（包括载荷和支反力）
    f_vector = KTotal @ UVector  # 点乘
    f_matrix = f_vector.reshape(NodeNum, 2)  # 将f整形为2列的矩阵
    return f_matrix


def sigma_solve():  # 求解应力
    sigma = np.zeros(StickNum)

    for i in range(0, StickNum):
        ind1 = Stick[i, 0]
        ind2 = Stick[i, 1]
        u_cut = np.array([[UVector[2 * ind1 - 2]], [UVector[2 * ind1 - 1]],
                          [UVector[2 * ind2 - 2]], [UVector[2 * ind2 - 1]]])
        c_ap = E / L[i] * np.array([-C[i], -S[i], C[i], S[i]])  # [C']
        sigma[i] = c_ap @ u_cut  # 点乘
    return sigma


def txt_output_2d(txt_name, dig):
    def add_seq(matrix):
        line = matrix.shape[0]
        seq = np.array([i for i in range(1, line + 1)])
        seq = seq.reshape(line, 1)
        seq_add_matrix = np.block([[seq, matrix]])
        return seq_add_matrix

    u_output = add_seq(UMatrix * 1000).round(dig)
    f_output = add_seq(FMatrix / 1000).round(dig)

    sigma_line = Sigma.shape[0]
    sigma_t = Sigma.reshape(sigma_line, 1)

    sigma_output = add_seq(sigma_t / 10 ** 6).round(dig)

    with open(txt_name, 'w') as FOutput:
        FOutput.write('2D Truss System Results\n')
        FOutput.write('U/mm\n')
        FOutput.write(str(u_output))
        FOutput.write('\n')
        FOutput.write('F/kN\n')
        FOutput.write(str(f_output))
        FOutput.write('\n')
        FOutput.write('sigma/MPa\n')
        FOutput.write(str(sigma_output))
        FOutput.write('\n')


def fig_output_2d(fig_name, u_factor):
    def stick_draw():  # 绘制杆件原始位置
        for i in range(StickNum):
            x = [Node[Stick[i, 0] - 1, 0], Node[Stick[i, 1] - 1, 0]]
            y = [Node[Stick[i, 0] - 1, 1], Node[Stick[i, 1] - 1, 1]]
            plt.plot(x, y, color='gray', linestyle=':', dash_capstyle='round', linewidth=2, zorder=0)

    def node_draw(circle_size):  # 绘制节点原始位置
        node_x = Node[:, 0]
        node_y = Node[:, 1]
        plt.scatter(node_x, node_y, circle_size, facecolors='w', edgecolors='gray', zorder=1)

    def value_max(matrix, num):  # 求最大值
        value = np.zeros(num)
        for i in range(num):
            value[i] = (matrix[i, 0] ** 2 + matrix[i, 1] ** 2) ** 0.5
        value_max = np.max(value)
        return value_max

    def f_draw(f_factor):  # 绘制外载荷
        for i in ForceNode:
            plt.quiver(Node[i - 1, 0] - Force[i - 1, 0] / f_factor, Node[i - 1, 1] - Force[i - 1, 1] / f_factor,
                       Force[i - 1, 0] / f_factor, Force[i - 1, 1] / f_factor, angles='xy', scale_units='xy', scale=1,
                       color='gray', zorder=2)

    def stick_result_draw(node_result, cmap, cmap_par):  # 绘制杆件计算结果
        for i in range(StickNum):
            x = [node_result[Stick[i, 0] - 1, 0], node_result[Stick[i, 1] - 1, 0]]
            y = [node_result[Stick[i, 0] - 1, 1], node_result[Stick[i, 1] - 1, 1]]
            plt.plot(x, y, c=cmap(cmap_par[i]), linewidth=3, zorder=3)
        cb = plt.colorbar()
        cb.set_label('$\it{\sigma}$(normalized)')

    def normalize(var, num):  # 归一化
        var_normalized = np.zeros(num)
        for i in range(num):
            var_normalized[i] = (var[i] - np.min(var)) / (np.max(var) - np.min(var))
        return var_normalized

    def node_result_draw(node_result):  # 绘制节点计算结果
        node_result_x = node_result[:, 0]
        node_result_y = node_result[:, 1]
        u_value = np.zeros(NodeNum)
        for i in range(NodeNum):
            u_value[i] = 1000 * (UMatrix[i, 0] ** 2 + UMatrix[i, 1] ** 2) ** 0.5
        plt.scatter(node_result_x, node_result_y, circle_size, c=u_value, cmap='jet', edgecolors='k', zorder=4)

    circle_size = 60
    stick_draw()
    node_draw(circle_size)
    node_distance_max = value_max(Node, NodeNum)
    f_value_max = value_max(Force, NodeNum)

    f_factor = 5 * f_value_max / node_distance_max  # 确定外载荷缩小系数

    f_draw(f_factor)

    node_result = Node + u_factor * UMatrix

    sigma_abs = abs(Sigma)
    sigma_normalized = normalize(sigma_abs, StickNum)
    cmap_stick = cm.get_cmap('viridis')
    stick_result_draw(node_result, cmap_stick, sigma_normalized)

    node_result_draw(node_result)
    cb = plt.colorbar()
    cb.set_label('$\it{u}$/mm')
    plt.axis('equal')
    fig_title = 'Visualization Result, Amplification Factor of u = ' + str(u_factor)
    plt.title(fig_title)
    plt.savefig(fig_name, dpi=600)
    plt.show()


def whether_default(var, default, var_type):
    if len(var.strip(' ')) == 0 or var.lower() == 'default':  # “default”不区分大小写
        var = var_type(default)
    else:
        var = var_type(var)
    return var


DigDefault = 3
FacDefault = 50

FileNameInput = input('请输入输入文件名：\n')
print('以下内容不输入或输入空格或输入“default”，则使用默认值。')
print('请输入计算结果保留小数位数，默认值为', DigDefault, '：')
Dig = input()
print('请输入位移放大系数，默认值为', FacDefault, '：')
Fac = input()

if FileNameInput[-4:] != '.txt':  # 如果最后4个字符不是“.txt”
    FileNameInput = FileNameInput + '.txt'  # 防止用户忘记输入后缀名

Dig = whether_default(Dig, DigDefault, int)
Fac = whether_default(Fac, FacDefault, float)

TxtNameOutput = 'Result_of_' + FileNameInput
FigNameOutput = 'Result_of_' + FileNameInput.strip('.txt') + '.png'

[A, E, Node, NodeNum, Stick, StickNum, ForceNode, Force, Con] = input_2d(FileNameInput)

[L, C, S] = stick_len_and_angle()

KElement = k_element()

KTotal = k_total(KElement)

[UMatrix, UVector] = u_solve()

FMatrix = f_solve()

Sigma = sigma_solve()

txt_output_2d(TxtNameOutput, Dig)
fig_output_2d(FigNameOutput, Fac)

print('计算完成，输出文件请见', TxtNameOutput, ',可视化结果请见', FigNameOutput, '。')
