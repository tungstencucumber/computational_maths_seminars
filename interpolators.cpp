#include <iostream>
#include <vector>

class LagrangeInterpolator {
private:
  std::vector<double> nodes; // набор узлов интерполяции
  std::vector<double> coefficients; // вектор коэффициентов в формуле интерполяционного многочлена Лагранжа, k-й элемент которого равен f(x_k) / П_{i = 0}^N(x_k - x_i)
public:
  LagrangeInterpolator(std::vector<double> _nodes, std::vector<double> _nodeFunctionValues) {
      this->nodes = std::move(_nodes);
      // если точек для построения интерполянта недостаточно, выбрасываем исключение
      if (nodes.size() != _nodeFunctionValues.size())
        throw -1;
      std::vector<double> coefficients = _nodeFunctionValues;
      for (int k = 0; k < nodes.size(); ++k) {
        double t = 1.;
        for (int i = 0; i < nodes.size(); ++i) {
          if (i == k) continue;
          t /= (nodes[k] - nodes[i]);
        }
        coefficients[k] *= t;
      }
      this->coefficients = std::move(coefficients);
  }
  double interpolate(double x) {
    if (x < nodes[0] || x > nodes[nodes.size() - 1])
      throw -1;
    double result = 0.;
    for (int k = 0; k < nodes.size(); ++k) {
      double t = 1.;
      for (int i = 0; i < nodes.size(); ++i) {
        if (i == k) continue;
        t *= (x - nodes[i]);
      }
      result += coefficients[k] * t;
    }
    return result;
  }
};

class NewtonInterpolator {
private:
  std::vector<double> nodes; // набор узлов интерполяции
  std::vector<double> coefficients; // вектор разделённых разностей, используемых в конечной формуле интерполяции Ньютона
public:
  NewtonInterpolator(std::vector<double> _nodes, std::vector<double> _nodeFunctionValues) {
    this->nodes = std::move(_nodes);
    // если точек для построения интерполянта недостаточно, выбрасываем исключение
    if (nodes.size() != _nodeFunctionValues.size())
      throw -1;
    std::vector<double> _coefficients;
    _coefficients.resize(nodes.size());
    _coefficients[0] = _nodeFunctionValues[0];
    // считаем разделённые разности
    std::vector<double> oldDifferences = _nodeFunctionValues;
    for (int k = 1; k < nodes.size(); ++k) {
      std::vector<double> newDifferences;
      newDifferences.resize(nodes.size() - k);
      for (int i = 1; i < newDifferences.size(); ++i) {
        newDifferences[i - 1] = (oldDifferences[i] - oldDifferences[i - 1])/(nodes[i - 1 + k] - nodes[i-1]);
      }
      _coefficients[k] = newDifferences[0];
      oldDifferences.resize(newDifferences.size());
      oldDifferences = newDifferences;
    }
    this->coefficients = std::move(_coefficients);
  }

  double interpolate(double x) {
    if (x < nodes[0] || x > nodes[nodes.size() - 1])
      throw -1;
    double result = 0.;
    for (int k = 0; k < nodes.size(); ++k) {
      double t = 1.;
      for (int i = 0; i < k; ++i) {
        t *= (x - nodes[i]);
      }
      result += coefficients[k] * t;
    }
    return result;
  }
};

int main() {
  std::vector<double> nodes = {0., 1., 2.};
  std::vector<double> function = {2., 4., 6.};
  NewtonInterpolator myInterpolator(nodes, function);
  LagrangeInterpolator dimasInterpolator(nodes, function);
  std::vector<double> referencePoints = {0.05, 0.1, 0.3, 0.77, 1., 1.22};
  for (int j = 0; j < referencePoints.size(); ++j)
    std::cout << myInterpolator.interpolate(referencePoints[j]) -
                dimasInterpolator.interpolate(referencePoints[j]) << std::endl;
  return 0;
}
