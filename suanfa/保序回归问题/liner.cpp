#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <limits>

// 线性规划保序回归的简单实现
std::vector<double> isotonic_regression_lp(const std::vector<double>& y) {
    int n = y.size();
    std::vector<double> ghat(y);
    std::vector<double> cost(n, 0.0);
    std::vector<int> index(n, 0);

    // 初始化
    for (int i = 0; i < n; ++i) {
        cost[i] = y[i];
        index[i] = i;
    }

    // 贪心算法
    for (int i = 1; i < n; ++i) {
        if (cost[i - 1] > cost[i]) {
            int j = i;
            double sum = 0.0;
            int count = 0;
            while (j >= 0 && cost[j] > cost[i]) {
                sum += cost[j];
                count++;
                j--;
            }
            double mean = sum / count;
            for (int k = j + 1; k <= i; ++k) {
                cost[k] = mean;
            }
        }
    }

    // 将结果复制到ghat
    for (int i = 0; i < n; ++i) {
        ghat[i] = cost[i];
    }

    return ghat;
}

int main() {
    // 记录开始时间
    auto start = std::chrono::high_resolution_clock::now();

    // 读取输入文件
    std::ifstream infile("input.txt");
    std::vector<double> y;
    double value;
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            while (iss >> value) {
                y.push_back(value);
            }
        }
        infile.close();
    } else {
        std::cerr << "无法打开输入文件" << std::endl;
        return 1;
    }

    // 应用线性规划保序回归算法
    std::vector<double> ghat;
    try {
        ghat = isotonic_regression_lp(y);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // 将结果写入输出文件
    std::ofstream outfile("output.txt");
    if (outfile.is_open()) {
        for (double val : ghat) {
            outfile << val << " ";
        }
        outfile.close();
    } else {
        std::cerr << "无法打开输出文件" << std::endl;
        return 1;
    }

    // 记录结束时间
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // 输出运行时间
    std::cout << "程序运行时间: " << duration.count() << " 秒" << std::endl;

    return 0;
}
