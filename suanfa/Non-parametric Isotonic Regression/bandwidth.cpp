#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>

// 高斯核函数
double gaussian_kernel(double u) {
    return std::exp(-0.5 * u * u) / std::sqrt(2 * M_PI);
}

// 核回归
std::function<double(double)> kernel_regression(const std::vector<double>& x, const std::vector<double>& y, double bandwidth, std::function<double(double)> kernel = gaussian_kernel) {
    int n = x.size();
    
    return [=](double x0) {
        std::vector<double> weights(n);
        for (int j = 0; j < n; ++j) {
            weights[j] = kernel((x0 - x[j]) / bandwidth);
        }
        double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        double weighted_sum = std::inner_product(weights.begin(), weights.end(), y.begin(), 0.0);
        return weighted_sum / weight_sum;
    };
}

// Silverman's Rule of Thumb for bandwidth selection
double silverman_bandwidth(const std::vector<double>& x) {
    int n = x.size();
    double sigma = std::sqrt(std::accumulate(x.begin(), x.end(), 0.0, [=](double a, double b) { return a + std::pow(b - std::accumulate(x.begin(), x.end(), 0.0) / n, 2); }) / n);
    return 1.06 * sigma * std::pow(n, -1.0 / 5.0);
}

// Scott's Rule for bandwidth selection
double scott_bandwidth(const std::vector<double>& x) {
    int n = x.size();
    double sigma = std::sqrt(std::accumulate(x.begin(), x.end(), 0.0, [=](double a, double b) { return a + std::pow(b - std::accumulate(x.begin(), x.end(), 0.0) / n, 2); }) / n);
    return 3.49 * sigma * std::pow(n, -1.0 / 3.0);
}

// 计算均方误差
double mean_squared_error(const std::vector<double>& y_true, const std::vector<double>& y_pred) {
    double mse = 0.0;
    int n = y_true.size();
    for (int i = 0; i < n; ++i) {
        mse += std::pow(y_true[i] - y_pred[i], 2);
    }
    return mse / n;
}

int main() {
    // 记录开始时间
    auto start = std::chrono::high_resolution_clock::now();

    // 读取输入文件
    std::ifstream infile("input3.txt");
    std::vector<double> x, y;
    double value_x, value_y;
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (iss >> value_x >> value_y) {
                x.push_back(value_x);
                y.push_back(value_y);
            }
        }
        infile.close();
    } else {
        std::cerr << "无法打开输入文件" << std::endl;
        return 1;
    }

    // 使用Silverman's Rule和Scott's Rule选择带宽
    double silverman_h = silverman_bandwidth(x);
    double scott_h = scott_bandwidth(x);
    std::cout << "根据Silverman's Rule选择的带宽: " << silverman_h << std::endl;
    std::cout << "根据Scott's Rule选择的带宽: " << scott_h << std::endl;

    // 应用核回归算法
    auto kr_silverman = kernel_regression(x, y, silverman_h);
    auto kr_scott = kernel_regression(x, y, scott_h);

    // 生成新数据点
    std::vector<double> x_new(100);
    double step = (x.back() - x.front()) / (x_new.size() - 1);
    std::generate(x_new.begin(), x_new.end(), [n = 0, &step]() mutable { return n++ * step; });

    std::vector<double> y_silverman(x_new.size());
    std::vector<double> y_scott(x_new.size());
    for (size_t i = 0; i < x_new.size(); ++i) {
        y_silverman[i] = kr_silverman(x_new[i]);
        y_scott[i] = kr_scott(x_new[i]);
    }

    // 计算均方误差
    std::vector<double> y_pred_silverman(x.size());
    std::vector<double> y_pred_scott(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        y_pred_silverman[i] = kr_silverman(x[i]);
        y_pred_scott[i] = kr_scott(x[i]);
    }
    double mse_silverman = mean_squared_error(y, y_pred_silverman);
    double mse_scott = mean_squared_error(y,    y_pred_scott);
    std::cout << "Silverman's Rule 带宽的均方误差: " << mse_silverman << std::endl;
    std::cout << "Scott's Rule 带宽的均方误差: " << mse_scott << std::endl;

    // 将结果写入输出文件
    std::ofstream outfile("output3.txt");
    if (outfile.is_open()) {
        outfile << "x_new y_silverman y_scott\n";
        for (size_t i = 0; i < y_silverman.size(); ++i) {
            outfile << x_new[i] << " " << y_silverman[i] << " " << y_scott[i] << "\n";
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

