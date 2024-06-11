#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <fstream>
#include "third_party/eigen/Eigen/Dense"

using namespace Eigen;

// 计算两个向量之间的皮尔逊相关系数
double pearsonCorrelation(const VectorXd &x, const VectorXd &y) {
    double mean_x = x.mean();
    double mean_y = y.mean();
    double numerator = ((x.array() - mean_x) * (y.array() - mean_y)).sum();
    double denominator = std::sqrt(((x.array() - mean_x).square().sum()) * ((y.array() - mean_y).square().sum()));
    return numerator / denominator;
}

// SIS方法筛选出与响应变量最相关的特征
std::vector<int> sisMethod(const MatrixXd &X, const VectorXd &y, int num_features) {
    std::vector<std::pair<double, int>> correlations;
    for (int j = 0; j < X.cols(); ++j) {
        double corr = pearsonCorrelation(X.col(j), y);
        correlations.push_back(std::make_pair(std::abs(corr), j));
    }
    std::sort(correlations.rbegin(), correlations.rend());
    std::vector<int> selected_indices;
    for (int i = 0; i < num_features; ++i) {
        selected_indices.push_back(correlations[i].second);
    }
    return selected_indices;
}

int main() {
    // 生成高维数据集
    int n_samples = 100;
    int n_features = 500;
    MatrixXd X = MatrixXd::Random(n_samples, n_features);
    VectorXd y = VectorXd::Random(n_samples);

    // 选择前10个最相关的特征
    int num_features = 10;
    std::vector<int> selected_indices = sisMethod(X, y, num_features);

    // 构建筛选后的特征矩阵
    MatrixXd X_selected(n_samples, num_features);
    for (int i = 0; i < num_features; ++i) {
        X_selected.col(i) = X.col(selected_indices[i]);
    }

    // 拟合线性回归模型
    VectorXd beta = (X_selected.transpose() * X_selected).ldlt().solve(X_selected.transpose() * y);

    // 将结果输出到txt文件
    std::ofstream outputFile("output4.txt");

    // 输出高维数据集
    outputFile << "High-dimensional dataset (X):\n";
    outputFile << X << "\n\n";

    // 输出响应变量
    outputFile << "Response variable (y):\n";
    outputFile << y << "\n\n";

    // 输出选择的特征索引
    outputFile << "Selected feature indices:\n";
    for (int idx : selected_indices) {
        outputFile << idx << " ";
    }
    outputFile << "\n\n";

    // 输出线性回归系数
    outputFile << "Linear regression coefficients:\n";
    outputFile << beta.transpose() << "\n";

    outputFile.close();

    std::cout << "Results have been written to output.txt" << std::endl;

    return 0;
}
