import numpy as np

# 生成示例数据
x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
y = np.array([1, 2, 3, 7, 6, 5, 8, 9, 10, 12])

# 将数据写入 input.txt 文件
with open('input2.txt', 'w') as f:
    for i in range(len(x)):
        f.write(f"{x[i]} {y[i]}\n")