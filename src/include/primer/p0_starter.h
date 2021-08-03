//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// p0_starter.h
//
// Identification: src/include/primer/p0_starter.h
//
// Copyright (c) 2015-2020, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include <memory>

namespace bustub {

/*
 * The base class defining a Matrix
 */
template <typename T>
class Matrix {
 protected:
  // TODO(P0): Add implementation
  Matrix(int r, int c) {
    // 为矩阵分配内存空间
    linear = new T[r * c];

    // 设置矩阵行数
    this->rows = r;

    // 设置矩阵列数
    this->cols = c;
  }

  // # of rows in the matrix
  int rows;
  // # of Columns in the matrix
  int cols;
  // Flattened array containing the elements of the matrix
  // TODO(P0) : Allocate the array in the constructor. Don't forget to free up
  // the array in the destructor.
  T *linear;

 public:
  // Return the # of rows in the matrix
  virtual int GetRows() = 0;

  // Return the # of columns in the matrix
  virtual int GetColumns() = 0;

  // Return the (i,j)th  matrix element
  virtual T GetElem(int i, int j) = 0;

  // Sets the (i,j)th  matrix element to val
  virtual void SetElem(int i, int j, T val) = 0;

  // Sets the matrix elements based on the array arr
  virtual void MatImport(T *arr) = 0;

  // TODO(P0): Add implementation
  virtual ~Matrix() {
    // 释放矩阵所占空间
    delete[] linear;
  }
};

template <typename T>
class RowMatrix : public Matrix<T> {
 public:
  // TODO(P0): Add implementation
  RowMatrix(int r, int c) : Matrix<T>(r, c) {
    // 为矩阵分配内存空间
    data_ = new T*[r];
    for (int i = 0; i < r; i++) {
      data_[i] = new T[c];
    }
  }

  // TODO(P0): Add implementation
  int GetRows() override {
    // 返回矩阵行数
    return this->rows;
  }

  // TODO(P0): Add implementation
  int GetColumns() override {
    // 返回矩阵列数
    return this->cols;
  }

  // TODO(P0): Add implementation
  T GetElem(int i, int j) override {
    // 获取矩阵指定位置的元素
    return data_[i][j];
  }

  // TODO(P0): Add implementation
  void SetElem(int i, int j, T val) override {
    // 为矩阵的指定位置设置元素值
    data_[i][j] = val;
  }

  // TODO(P0): Add implementation
  void MatImport(T *arr) override {
    // 从一维数组向矩阵导入元素
    for (int i = 0; i < this->rows; i++) {
      for (int j = 0; j < this->cols; j++) {
        data_[i][j] = arr[this->cols * i + j];
      }
    }
  }

  // TODO(P0): Add implementation
  ~RowMatrix() override {
    // 释放矩阵所占空间
    for (int i = 0; i < this->rows; i++) {
      delete[] data_[i];
    }

    delete[] data_;
  }

 private:
  // 2D array containing the elements of the matrix in row-major format
  // TODO(P0): Allocate the array of row pointers in the constructor. Use these pointers
  // to point to corresponding elements of the 'linear' array.
  // Don't forget to free up the array in the destructor.
  T **data_;
};

template <typename T>
class RowMatrixOperations {
 public:
  // Compute (mat1 + mat2) and return the result.
  // Return nullptr if dimensions mismatch for input matrices.
  static std::unique_ptr<RowMatrix<T>> AddMatrices(std::unique_ptr<RowMatrix<T>> mat1,
                                                   std::unique_ptr<RowMatrix<T>> mat2) {
    // TODO(P0): Add code
    if (mat1 == nullptr || mat2 == nullptr) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }

    if (mat1->GetRows() != mat2->GetRows() || mat1->GetColumns() != mat2->GetColumns()) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }

    int rows = mat1->GetRows(), cols = mat1->GetColumns();
    auto result = std::unique_ptr<RowMatrix<T>>(new RowMatrix<T>(rows, cols));
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        T val1 = mat1->GetElem(i, j);
        T val2 = mat2->GetElem(i, j);
        result->SetElem(i, j, val1 + val2);
      }
    }

    return result;
  }

  // Compute matrix multiplication (mat1 * mat2) and return the result.
  // Return nullptr if dimensions mismatch for input matrices.
  static std::unique_ptr<RowMatrix<T>> MultiplyMatrices(std::unique_ptr<RowMatrix<T>> mat1,
                                                        std::unique_ptr<RowMatrix<T>> mat2) {
    if (mat1 == nullptr || mat2 == nullptr) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }

    if (mat1->GetColumns() != mat2->GetRows()) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }

    int M = mat1->GetRows(), K = mat2->GetRows(), N = mat2->GetColumns();

    auto result = std::unique_ptr<RowMatrix<T>>(new RowMatrix<T>(M, N));

    // 初始化返回矩阵
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        result->SetElem(i, j, 0);
      }
    }

    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < K; k++) {
          result->SetElem(i, j, result->GetElem(i, j) + mat1->GetElem(i, k) * mat2->GetElem(k, j));
        }
      }
    }

    return result;
  }

  // Simplified GEMM (general matrix multiply) operation
  // Compute (matA * matB + matC). Return nullptr if dimensions mismatch for input matrices
  static std::unique_ptr<RowMatrix<T>> GemmMatrices(std::unique_ptr<RowMatrix<T>> matA,
                                                    std::unique_ptr<RowMatrix<T>> matB,
                                                    std::unique_ptr<RowMatrix<T>> matC) {
    // TODO(P0): Add code
    if (matA == nullptr || matB == nullptr || matC == nullptr) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }

    int rowsA = matA->GetRows(), colsA = matA->GetColumns(), rowsB = matB->GetRows(), colsB = matB->GetColumns(), rowsC = matC->GetRows(), colsC = matC->GetColumns();
    if (colsA != rowsB || rowsC != rowsA || colsC != colsB) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }

    auto result = std::unique_ptr<RowMatrix<T>>(new RowMatrix<T>(rowsC, colsC));

    // 初始化返回矩阵
    for (int i = 0; i < rowsC; i++) {
      for (int j = 0; j < colsC; j++) {
        result->SetElem(i, j, 0);
      }
    }

    for (int i = 0; i < rowsA; i++) {
      for (int j = 0; j < colsB; j++) {
        for (int k = 0; k < colsA; k++) {
          matC->SetElem(i, j, matC->GetElem(i, j) + matA->GetElem(i, k) * matB->GetElem(k, j));
        }
      }
    }
    return result;
  }
};
}  // namespace bustub