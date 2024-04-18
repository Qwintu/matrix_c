#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error_code = 0;
  if (rows < 1 || columns < 1 || result == NULL) return 1;
  result->matrix = NULL;
  result->rows = rows;
  result->columns = columns;

  size_t matrix_size = rows * columns;
  if (!(result->matrix = (calloc(rows + matrix_size, sizeof(double))))) {
    error_code = 1;
  } else {
    double *rows_ptr = (double *)(result->matrix + rows);
    for (int i = 0; i < rows; i++) {
      (result->matrix[i]) = (rows_ptr + columns * i);
    }
  }
  return error_code;
}

void s21_remove_matrix(matrix_t *A) {
  free(A->matrix);
  A->matrix = NULL;
  A->rows = 0;
  A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  // #define SUCCESS 1
  // #define FAILURE 0

  if (check_matrix(A) || check_matrix(B)) return 0;

  double eps = 0.000001;
  int error_code = 1;

  if ((A->rows != B->rows) || (A->columns != B->columns)) return 0;

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= eps) {
        error_code = 0;
        break;
      }
    }
    if (error_code == 0) break;
  }
  return error_code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (check_matrix(A)) return 1;
  int error_code = 0;
  if (s21_create_matrix(A->rows, A->columns, result)) {
    error_code = 1;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  };
  return error_code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (check_matrix(A) || check_matrix(B)) return 1;
  if ((A->rows != B->rows) || (A->columns != B->columns)) return 2;

  int error_code = 0;
  if (s21_create_matrix(A->rows, A->columns, result)) {
    error_code = 1;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  };
  return error_code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (check_matrix(A) || check_matrix(B)) return 1;
  if ((A->rows != B->rows) || (A->columns != B->columns)) return 2;

  int error_code = 0;
  if (s21_create_matrix(A->rows, A->columns, result)) {
    error_code = 1;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  };
  return error_code;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (check_matrix(A) || check_matrix(B)) return 1;
  if ((A->columns != B->rows)) return 2;

  int error_code = 0;
  int res_rows = A->rows;
  int res_columns = B->columns;
  if (s21_create_matrix(res_rows, res_columns, result)) {
    error_code = 1;
  } else {
    for (int i = 0; i < res_rows; i++) {
      for (int j = 0; j < res_columns; j++) {
        double value_of_calc = 0;
        for (int n = 0; n < A->columns; n++) {
          value_of_calc += A->matrix[i][n] * B->matrix[n][j];
        }
        result->matrix[i][j] = value_of_calc;
      }
    }
  };
  return error_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (check_matrix(A)) return 1;
  int error_code = 0;

  if (s21_create_matrix(A->columns, A->rows, result)) {
    error_code = 1;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  };
  return error_code;
}

int s21_determinant(matrix_t *A, double *result) {
  if (check_matrix(A)) return 1;
  if (A->columns != A->rows) return 2;
  double det = 0;
  int error_code = 0;

  if (A->rows <= 3 && A->rows > 0) {
    determinant_3_2_1(A, &det);
  } else if (A->rows > 3) {
    for (int j = 0; j < A->columns; j++) {
      double minor_det = 0;
      matrix_t minor = {0};
      s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
      create_minor(A, 0, j, &minor);
      s21_determinant(&minor, &minor_det);
      det += pow((-1), (0 + j)) * A->matrix[0][j] * minor_det;
      s21_remove_matrix(&minor);
    }
  } else {
    error_code = 1;
  }

  *result = det;
  return error_code;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (check_matrix(A)) return 1;
  if (A->columns != A->rows) return 2;
  int error_code = 0;

  if (s21_create_matrix(A->rows, A->columns, result)) {
    error_code = 1;
  } else {
    if (A->rows == 1) {
      result->matrix[0][0] = A->matrix[0][0];
    } else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          double minor_det = 0;
          matrix_t minor = {0};
          s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
          create_minor(A, i, j, &minor);
          s21_determinant(&minor, &minor_det);
          minor_det *= pow((-1), (i + j));
          result->matrix[i][j] = minor_det;
          s21_remove_matrix(&minor);
        }
      }
    }
  }
  return error_code;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (check_matrix(A)) return 1;
  if (A->columns != A->rows) return 2;

  int error_code = 0;
  double det = 0;
  s21_determinant(A, &det);
  if (fabs(det - 0.0) < 1e-6) return 2;

  matrix_t complements_matrix = {0};
  error_code = s21_calc_complements(A, &complements_matrix);

  matrix_t transp_compl_matrix = {0};
  error_code += s21_transpose(&complements_matrix, &transp_compl_matrix);

  error_code += s21_mult_number(&transp_compl_matrix, 1 / det, result);
  s21_remove_matrix(&complements_matrix);
  s21_remove_matrix(&transp_compl_matrix);
  if (error_code > 1) error_code = 1;
  return error_code;
}

// -------

int check_matrix(matrix_t *A) {
  int error_code = 0;  // 0 - ok
  if (A == NULL || A->matrix == NULL || A->rows < 1 || A->columns < 1)
    error_code = 1;
  return error_code;
}

int create_minor(matrix_t *A, int i, int j, matrix_t *minor) {
  int minor_row = 0;
  int minor_columns = 0;

  for (int k = 0; k < A->rows; k++) {
    if (k == i) continue;
    minor_columns = 0;
    for (int m = 0; m < A->columns; m++) {
      if (m != j) {
        minor->matrix[minor_row][minor_columns] = A->matrix[k][m];
        minor_columns++;
      }
    }
    minor_row++;
  }
  return 0;
}

int determinant_3_2_1(matrix_t *A, double *result) {
  int error_code = 0;
  if (A->rows == 1) {
    *result = A->matrix[0][0];
  } else if (A->rows == 2) {
    // определитель матрицы 2*2
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[1][0] * A->matrix[0][1];
  } else if (A->rows == 3) {
    // определитель матрицы 3*3 методом треугольника
    *result = A->matrix[0][0] * A->matrix[1][1] * A->matrix[2][2] +
              +A->matrix[1][0] * A->matrix[2][1] * A->matrix[0][2] +
              +A->matrix[0][1] * A->matrix[1][2] * A->matrix[2][0] +
              -A->matrix[0][2] * A->matrix[1][1] * A->matrix[2][0] +
              -A->matrix[1][0] * A->matrix[0][1] * A->matrix[2][2] -
              A->matrix[2][1] * A->matrix[1][2] * A->matrix[0][0];
  } else {
    error_code = 1;
  }
  return error_code;
}

// -----------

void fill_matrix(int rows, int columns, double start_from, double step,
                 matrix_t *result) {
  double mi = start_from;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      (result->matrix)[i][j] = mi;
      mi += step;
    }
  }
}

void change_matrix_value(matrix_t *A, int i, int j, double new_value) {
  A->matrix[i][j] = new_value;
}

void fill_matrix_by_massive(int rows, int columns, double *values,
                            matrix_t *result) {
  int mass_ind = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      (result->matrix)[i][j] = values[mass_ind];
      mass_ind++;
      ;
    }
  }
}