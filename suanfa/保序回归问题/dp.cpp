#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>

std::vector<double> isotonic_regression_dp(const std::vector<double>& y) {
    int n = y.size();
    std::vector<double> cost(n, 0.0);
    std::vector<int> index(n, 0);
    std::vector<double> ghat(n, 0.0);
    
    for (int i = 0; i < n; ++i) {
        cost[i] = y[i];
        index[i] = i;
        
        for (int j = 0; j < i; ++j) {
            if (cost[j] <= cost[i]) {
                if (cost[j] + y[i] < cost[i]) {
                    cost[i] = cost[j] + y[i];
                    index[i] = j;
                }
            }
        }
    }
    
    int i = n - 1;
    while (i >= 0) {
        int j = index[i];
        for (int k = j; k <= i; ++k) {
            ghat[k] = cost[i] / (i - j + 1);
        }
        i = j - 1;
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
    
    // 应用动态规划保序回归算法
    std::vector<double> ghat = isotonic_regression_dp(y);
    
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
