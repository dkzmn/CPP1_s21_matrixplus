#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 0 || cols < 0) {
    throw std::invalid_argument("Invalid size");
  }
  allocMemory();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  allocMemory();
  copyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix::~S21Matrix() { clear_matrix(); }

void S21Matrix::allocMemory() {
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::copyMatrix(const S21Matrix& other) {
  int min_rows = other.getRows();
  if (rows_ < min_rows) min_rows = rows_;
  int min_cols = other.getCols();
  if (cols_ < min_cols) min_cols = cols_;
  for (int i = 0; i < min_rows; i++)
    for (int j = 0; j < min_cols; j++) matrix_[i][j] = other(i, j);
}

void S21Matrix::clear_matrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
  }
  delete[] matrix_;
  rows_ = 0;
  cols_ = 0;
  matrix_ = nullptr;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  bool res = true;
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    res = false;
  } else {
    for (int i = 0; i < rows_ && res != false; i++)
      for (int j = 0; j < cols_ && res != false; j++)
        if (fabs(matrix_[i][j] - other(i, j)) > ACC) res = false;
  }
  return res;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.getRows() || cols_ != other.getCols()) {
    throw std::length_error("Different dimensions of matrices");
  }
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] += other(i, j);
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.getRows() || cols_ != other.getCols()) {
    throw std::length_error("Different dimensions of matrices");
  }
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] -= other(i, j);
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] *= num;
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.getRows()) {
    throw std::length_error(
        "The number of columns of the first matrix is ​​not equal to the "
        "number of rows of the second matrix");
  }
  S21Matrix result = S21Matrix(rows_, other.getCols());
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < other.getCols(); j++)
      for (int k = 0; k < cols_; k++)
        result.matrix_[i][j] += matrix_[i][k] * other(k, j);
  *this = result;
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix result = S21Matrix(cols_, rows_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) result.matrix_[j][i] = matrix_[i][j];
  return result;
}

S21Matrix S21Matrix::Minor(int col, int row) const {
  S21Matrix result = S21Matrix(cols_ - 1, rows_ - 1);
  int newi, newj;
  for (int i = 0; i < rows_ - 1; i++) {
    if (i >= row)
      newi = i + 1;
    else
      newi = i;
    for (int j = 0; j < cols_ - 1; j++) {
      if (j >= col)
        newj = j + 1;
      else
        newj = j;
      result.matrix_[i][j] = matrix_[newi][newj];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::length_error("Matrix is ​​not square");
  }
  double det_temp = 0;
  S21Matrix result = S21Matrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor = Minor(j, i);
      det_temp = minor.Determinant();
      minor.~S21Matrix();
      if ((i + j) % 2 == 1) det_temp = -det_temp;
      result.matrix_[i][j] = det_temp;
    }
  }
  return result;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::length_error("Matrix is ​​not square");
  }
  double det = 0, det_temp = 0;
  if (cols_ == 1) {
    det = matrix_[0][0];
  } else if (cols_ == 2) {
    det = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    for (int i = 0; i < cols_; i++) {
      S21Matrix minor = Minor(i, 0);
      det_temp = minor.Determinant();
      minor.~S21Matrix();
      if (i % 2 == 1) det_temp = -det_temp;
      det = det + matrix_[0][i] * det_temp;
    }
  }
  return det;
}

S21Matrix S21Matrix::InverseMatrix() const {
  double det = Determinant();
  if (fabs(det) < ACC) {
    throw std::length_error("Matrix determinant is 0");
  }
  S21Matrix result;
  if (cols_ == 1) {
    S21Matrix result1 = S21Matrix(1, 1);
    result1.matrix_[0][0] = 1 / matrix_[0][0];
    result = result1;
  } else {
    S21Matrix calc_compl = CalcComplements();
    S21Matrix result2 = calc_compl.Transpose();
    result2.MulNumber(1 / det);
    result = result2;
  }
  return result;
}

void S21Matrix::setRows(int rows) {
  if (rows < 0) {
    throw std::length_error("Rows count must be >= 0");
  }
  S21Matrix result = S21Matrix(rows, cols_);
  result.copyMatrix(*this);
  *this = result;
}

int S21Matrix::getRows() const { return rows_; }

void S21Matrix::setCols(int cols) {
  if (cols < 0) {
    throw std::length_error("Cols count must be >= 0");
  }
  S21Matrix result = S21Matrix(rows_, cols);
  result.copyMatrix(*this);
  *this = result;
}

int S21Matrix::getCols() const { return cols_; }

double S21Matrix::operator()(int i, int j) const {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::out_of_range("Incorrect index!");
  }
  return matrix_[i][j];
}

double& S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::out_of_range("Incorrect index!");
  }
  return matrix_[i][j];
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) {
    throw std::invalid_argument("Same matrix");
  }
  clear_matrix();
  rows_ = other.getRows();
  cols_ = other.getCols();
  allocMemory();
  copyMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}

S21Matrix S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}
