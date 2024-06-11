# 生成示例数据
x = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
y = [1.2, 1.8, 2.4, 2.9, 3.0, 3.8, 4.1, 4.5, 5.2, 5.4]

# 将数据写入 input.txt 文件
with open('input3.txt', 'w') as f:
    for i in range(len(x)):
        f.write(f"{x[i]} {y[i]}\n")
