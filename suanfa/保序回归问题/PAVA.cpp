#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>

std::vector<double> pava(const std::vector<double>& y) {
    int n = y.size();
    std::vector<double> ghat(y);
    std::vector<int> stack(n, 0);
    int stack_size = 0;
    
    for (int i = 0; i < n; ++i) {
        stack[stack_size] = i;
        stack_size++;
        
        while (stack_size > 1 && ghat[stack[stack_size - 1]] < ghat[stack[stack_size - 2]]) {
            stack_size--;
            int k = stack[stack_size - 1];
            int j = stack[stack_size];
            
            double sum = 0.0;
            for (int m = k; m <= j; ++m) {
                sum += ghat[m];
            }
            
            double mean = sum / (j - k + 1);
            for (int m = k; m <= j; ++m) {
                ghat[m] = mean;
            }
        }
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
    
    // 应用PAVA算法
    std::vector<double> ghat = pava(y);
    
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
