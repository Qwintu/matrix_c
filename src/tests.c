#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "s21_matrix.h"

START_TEST(create_matrix) {
  matrix_t result;
  int rows = 2, columns = 10;
  int exit_code = s21_create_matrix(rows, columns, &result);
  ck_assert_int_eq(0, exit_code);
  s21_remove_matrix(&result);

  //
  rows = 2;
  columns = 2;
  exit_code = s21_create_matrix(rows, columns, &result);
  ck_assert_int_eq(0, exit_code);
  s21_remove_matrix(&result);

  //

  rows = 0;
  columns = 2;
  exit_code = s21_create_matrix(rows, columns, &result);
  ck_assert_int_eq(1, exit_code);
  s21_remove_matrix(&result);

  //

  rows = 0;
  columns = 0;
  exit_code = s21_create_matrix(rows, columns, &result);
  ck_assert_int_eq(1, exit_code);
  s21_remove_matrix(&result);
}
END_TEST

//  --------------------------

START_TEST(equal_matrix) {
  matrix_t A, B;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 5, &B);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  //

  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  double mass[9] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
  fill_matrix_by_massive(3, 3, mass, &A);
  fill_matrix_by_massive(3, 3, mass, &B);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ck_assert_double_eq(A.matrix[i][j], B.matrix[i][j]);
    }
  }

  ck_assert_int_eq(s21_eq_matrix(&A, &B), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  //

  s21_create_matrix(10, 15, &A);
  s21_create_matrix(10, 15, &B);
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 15; j++) {
      B.matrix[i][j] = A.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 15; j++) {
      ck_assert_double_eq(A.matrix[i][j], B.matrix[i][j]);
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&A, &B), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  //

  s21_create_matrix(10, 15, &A);
  s21_create_matrix(10, 15, &B);
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 15; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
      B.matrix[i][j] = A.matrix[i][j];
    }
  }
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 15; j++) {
      ck_assert_double_eq(A.matrix[i][j], B.matrix[i][j]);
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&A, &B), 1);  // ok
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  //

  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  double massA[9] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
  double massB[9] = {2, 5, 7, 6, 44, 4, 5, -2, -3};
  fill_matrix_by_massive(3, 3, massA, &A);
  fill_matrix_by_massive(3, 3, massB, &B);

  ck_assert_int_eq(s21_eq_matrix(&A, &B), 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  //

  s21_create_matrix(1, 2, &A);
  s21_create_matrix(2, 5, &B);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  //

  s21_create_matrix(1, 2, &A);
  B.matrix = NULL;
  ck_assert_int_eq(s21_eq_matrix(&A, &B), 0);
  s21_remove_matrix(&A);

  //

  s21_create_matrix(1, 2, &A);
  s21_create_matrix(0, 5, &B);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  //

  s21_create_matrix(0, 0, &A);
  s21_create_matrix(0, 0, &B);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

//  --------------------------

START_TEST(sum_matrix) {
  matrix_t A, B, result;
  int rows = 5, columns = 5;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(rows, columns, &B);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
      B.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), 0);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      ck_assert_double_eq_tol(A.matrix[i][j] + B.matrix[i][j],
                              result.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  //

  rows = 5;
  columns = 1;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(rows, columns, &B);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
      B.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), 0);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      ck_assert_double_eq_tol(A.matrix[i][j] + B.matrix[i][j],
                              result.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  //

  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 5, &B);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  //

  A.matrix = NULL, B.matrix = NULL;
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

//  --------------------------

START_TEST(sub_matrix) {
  matrix_t A, B, result;
  A.matrix = NULL, B.matrix = NULL;
  int rows = 5, columns = 5;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(rows, columns, &B);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
      B.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), 0);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      ck_assert_double_eq_tol(A.matrix[i][j] - B.matrix[i][j],
                              result.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  //

  rows = 8;
  columns = 5;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(rows, columns, &B);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
      B.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), 0);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      ck_assert_double_eq_tol(A.matrix[i][j] - B.matrix[i][j],
                              result.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  //

  rows = 5;
  columns = 1;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(rows, columns, &B);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
      B.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), 0);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      ck_assert_double_eq_tol(A.matrix[i][j] - B.matrix[i][j],
                              result.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  //

  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 5, &B);
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  //

  A.matrix = NULL, B.matrix = NULL;
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

//  --------------------------

START_TEST(mult_number_matrix) {
  matrix_t A, result;
  int rows = 5, columns = 5;
  double number = (double)(rand()) / RAND_MAX;
  s21_create_matrix(rows, columns, &A);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_mult_number(&A, number, &result), 0);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      ck_assert_double_eq_tol(A.matrix[i][j] * number, result.matrix[i][j],
                              1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  //

  rows = 5;
  columns = 7;
  number = (double)(rand()) / RAND_MAX;
  s21_create_matrix(rows, columns, &A);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_mult_number(&A, number, &result), 0);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      ck_assert_double_eq_tol(A.matrix[i][j] * number, result.matrix[i][j],
                              1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  //

  rows = 5;
  columns = 1;
  number = (double)(rand()) / RAND_MAX;
  s21_create_matrix(rows, columns, &A);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_mult_number(&A, number, &result), 0);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      ck_assert_double_eq_tol(A.matrix[i][j] * number, result.matrix[i][j],
                              1e-8);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  //

  A.matrix = NULL;
  number = (double)(rand()) / RAND_MAX;
  ck_assert_int_eq(s21_mult_number(&A, number, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  //

  number = (double)(rand()) / RAND_MAX;
  s21_create_matrix(0, 0, &A);
  ck_assert_int_eq(s21_mult_number(&A, number, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

//  --------------------------

START_TEST(mult_matrix) {
  matrix_t A, B, control, result;
  int RowsA = 3, ColumnsA = 4;
  int RowsB = 4, ColumnsB = 5;
  int RowsCont = 3, ColumnsCont = 5;
  s21_create_matrix(RowsA, ColumnsA, &A);
  s21_create_matrix(RowsB, ColumnsB, &B);
  s21_create_matrix(RowsCont, ColumnsCont, &control);
  double massA[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  double massB[20] = {3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                      13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
  double massC[15] = {130, 140, 150, 160, 170, 298, 324, 350,
                      376, 402, 466, 508, 550, 592, 634};
  fill_matrix_by_massive(3, 4, massA, &A);
  fill_matrix_by_massive(4, 5, massB, &B);
  fill_matrix_by_massive(3, 5, massC, &control);

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), 0);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&control);
  s21_remove_matrix(&result);

  //

  RowsA = 3, ColumnsA = 2;
  RowsB = 2, ColumnsB = 3;
  RowsCont = 3, ColumnsCont = 3;
  s21_create_matrix(RowsA, ColumnsA, &A);
  s21_create_matrix(RowsB, ColumnsB, &B);
  s21_create_matrix(RowsCont, ColumnsCont, &control);
  double massA1[6] = {
      1, 4, 2, 5, 3, 6,
  };
  double massB1[6] = {1, -1, 1, 2, 3, 4};
  double massC1[9] = {9, 11, 17, 12, 13, 22, 15, 15, 27};
  fill_matrix_by_massive(3, 2, massA1, &A);
  fill_matrix_by_massive(2, 3, massB1, &B);
  fill_matrix_by_massive(3, 3, massC1, &control);

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), 0);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&control);
  s21_remove_matrix(&result);

  //

  RowsA = 4;
  ColumnsA = 10;
  RowsB = 13;
  ColumnsB = 27;
  s21_create_matrix(RowsA, ColumnsA, &A);
  s21_create_matrix(RowsB, ColumnsB, &B);
  for (int i = 0; i < RowsA; i++) {
    for (int j = 0; j < ColumnsA; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  for (int i = 0; i < RowsB; i++) {
    for (int j = 0; j < ColumnsB; j++) {
      B.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);

  //

  A.matrix = NULL;
  B.matrix = NULL;
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

//  --------------------------

START_TEST(transpose_matrix) {
  matrix_t A, control, result;
  int rows = 3, columns = 5;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(columns, rows, &control);
  fill_matrix(3, 5, 1, 1, &A);
  double massC1[15] = {1, 6, 11, 2, 7, 12, 3, 8, 13, 4, 9, 14, 5, 10, 15};
  fill_matrix_by_massive(5, 3, massC1, &control);

  ck_assert_int_eq(s21_transpose(&A, &result), 0);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&control);

  //

  rows = 15;
  columns = 1;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(columns, rows, &control);
  fill_matrix(15, 1, 1, 1, &A);
  double massC2[15] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  fill_matrix_by_massive(1, 15, massC2, &control);

  ck_assert_int_eq(s21_transpose(&A, &result), 0);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&control);

  //

  A.matrix = NULL;
  ck_assert_int_eq(s21_transpose(&A, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

////////////////////////////////////////////////////////////////

START_TEST(calc_complements_matrix) {
  matrix_t A, control, result;
  int rows = 5, columns = 5;

  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(columns, rows, &control);
  fill_matrix(5, 5, 1, 1, &A);
  change_matrix_value(&A, 0, 0, 32);
  change_matrix_value(&A, 2, 1, 13);
  change_matrix_value(&A, 3, 2, 17);
  double massC3[25] = {15,  0,    0,   -60,  45,    -20, -310, 155, -75,
                       219, 0,    465, 0,    -1395, 930, 0,    0,   -465,
                       930, -465, 5,   -155, 310,   135, -264};
  fill_matrix_by_massive(5, 5, massC3, &control);
  ck_assert_int_eq(s21_calc_complements(&A, &result), 0);
  // print_matrix(result);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&control);

  //

  rows = 3;
  columns = 3;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(columns, rows, &control);
  double massC4[9] = {1, 2, 3, 0, 4, 2, 5, 2, 1};
  fill_matrix_by_massive(3, 3, massC4, &A);
  double massC5[9] = {0, 10, -20, 4, -14, 8, -8, -2, 4};
  fill_matrix_by_massive(3, 3, massC5, &control);
  ck_assert_int_eq(s21_calc_complements(&A, &result), 0);
  // print_matrix(control);
  // print_matrix(result);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&control);

  //

  rows = 1;
  columns = 1;
  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(columns, rows, &control);
  double massC6[1] = {14};
  fill_matrix_by_massive(1, 1, massC6, &A);
  double massC7[1] = {14};
  fill_matrix_by_massive(1, 1, massC7, &control);
  ck_assert_int_eq(s21_calc_complements(&A, &result), 0);
  // print_matrix(control);
  // print_matrix(result);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&control);

  //

  A.matrix = NULL;
  result.matrix = NULL;
  ck_assert_int_eq(s21_calc_complements(&A, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

////////////////////////////////////////////////////////////////

START_TEST(determinant_matrix) {
  matrix_t A;
  int rows = 5, columns = 5;
  double result = 0.0;

  s21_create_matrix(rows, columns, &A);
  fill_matrix(5, 5, 1, 1, &A);
  ck_assert_int_eq(s21_determinant(&A, &result), 0);
  ck_assert_int_eq(result, 0);
  s21_remove_matrix(&A);

  //

  s21_create_matrix(rows, columns, &A);
  fill_matrix(5, 5, 1, 1, &A);
  change_matrix_value(&A, 0, 0, -32);
  change_matrix_value(&A, 4, 4, 33);
  change_matrix_value(&A, 2, 2, -22);
  ck_assert_double_eq(s21_determinant(&A, &result), 0);
  ck_assert_int_eq(result, -184800);
  s21_remove_matrix(&A);

  //

  rows = 1;
  columns = 1;
  result = 0.0;
  s21_create_matrix(rows, columns, &A);
  fill_matrix(1, 1, 1, 1, &A);
  ck_assert_int_eq(s21_determinant(&A, &result), 0);
  ck_assert_int_eq(result, 1);
  s21_remove_matrix(&A);
  //

  rows = 2;
  columns = 2;
  result = 0.0;
  s21_create_matrix(rows, columns, &A);
  fill_matrix(2, 2, 1, 1, &A);
  ck_assert_int_eq(s21_determinant(&A, &result), 0);
  ck_assert_int_eq(result, -2);
  s21_remove_matrix(&A);

  //

  rows = 10;
  columns = 7;
  s21_create_matrix(rows, columns, &A);
  result = 0.0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_determinant(&A, &result), 2);
  s21_remove_matrix(&A);

  //

  A.matrix = NULL;
  result = 0.0;
  ck_assert_int_eq(s21_determinant(&A, &result), 1);
  s21_remove_matrix(&A);
}
END_TEST

////////////////////////////////////////////////////////////////

START_TEST(inverse_matrix) {
  matrix_t A, control, result;
  int rows = 3, columns = 3;

  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(columns, rows, &control);
  double massC1[9] = {3, 2, 5, 4, 5, 6, 7, 7, 7};
  fill_matrix_by_massive(3, 3, massC1, &A);
  double massC2[9] = {0.25,      -0.75, 0.464286, -0.5, 0.5,
                      -0.071429, 0.25,  0.25,     -0.25};
  fill_matrix_by_massive(3, 3, massC2, &control);
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), 0);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&control);

  //

  s21_create_matrix(rows, columns, &A);
  s21_create_matrix(columns, rows, &control);
  double massC3[9] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
  fill_matrix_by_massive(3, 3, massC3, &A);
  double massC4[9] = {1, -1, 1, -38, 41, -34, 27, -29, 24};
  fill_matrix_by_massive(3, 3, massC4, &control);
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), 0);
  ck_assert_int_eq(s21_eq_matrix(&result, &control), 1);  // ok

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&control);

  //

  rows = 10;
  columns = 7;
  s21_create_matrix(rows, columns, &A);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      A.matrix[i][j] = (double)(rand()) / RAND_MAX;
    }
  }
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  //

  rows = -10;
  columns = 7;
  s21_create_matrix(rows, columns, &A);
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  //

  A.matrix = NULL;
  result.matrix = NULL;
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);

  //
}
END_TEST

////////////////////////////////////////////////////////////////

Suite *matrix_suite(void) {
  Suite *s = suite_create("Matrix_testcase");
  TCase *tc_matrix = tcase_create("Core");

  suite_add_tcase(s, tc_matrix);
  tcase_add_test(tc_matrix, create_matrix);
  tcase_add_test(tc_matrix, equal_matrix);
  tcase_add_test(tc_matrix, sum_matrix);
  tcase_add_test(tc_matrix, sub_matrix);
  tcase_add_test(tc_matrix, mult_number_matrix);
  tcase_add_test(tc_matrix, mult_matrix);
  tcase_add_test(tc_matrix, transpose_matrix);
  tcase_add_test(tc_matrix, calc_complements_matrix);
  tcase_add_test(tc_matrix, determinant_matrix);
  tcase_add_test(tc_matrix, inverse_matrix);

  return s;
}

int main() {
  int number_failed;
  Suite *s;
  SRunner *sr;

  s = matrix_suite();
  sr = srunner_create(s);

  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}