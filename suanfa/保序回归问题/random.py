import random

# 生成1000个随机数
random_numbers = [random.randint(0, 100) for _ in range(1000)]

# 将随机数写入input.txt文件
with open('input.txt', 'w') as file:
    for number in random_numbers:
        file.write(f"{number} ")
