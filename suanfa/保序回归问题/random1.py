import numpy as np

# 随机生成1000列的自变量数据，包含1000行，值在0到10000之间
num_samples = 1000
num_features = 1000
X = np.random.randint(0, 10001, size=(num_samples, num_features))

# 随机生成因变量数据，值在0到10000之间
y = np.random.randint(0, 10001, size=num_samples)

# 将数据写入 input.txt 文件
with open('input1.txt', 'w') as f:
    for i in range(X.shape[0]):
        row = list(X[i]) + [y[i]]
        f.write(' '.join(map(str, row)) + '\n')
