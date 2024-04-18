CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c11 -g 

OS:=$(shell uname -s)
ifeq ($(OS), Darwin)
LDFLAGS = -lcheck -lm
OPEN_CM = open
endif
ifeq ($(OS),Linux)
LDFLAGS=-lcheck_pic -lrt -lpthread -lsubunit -lm
OPEN_CM=xdg-open
endif

all: clean s21_matrix.a gcov_report

test: clean s21_matrix.a
	$(CC) $(CFLAGS) tests.c s21_matrix.a -o s21_matrix_test $(LDFLAGS)
	./s21_matrix_test

s21_matrix.a:
	$(CC) $(CFLAGS) -c "s21_matrix.c" "s21_matrix.h"
	ar rc s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a

gcov_report: clean s21_matrix.a
	$(CC) $(CFLAGS) --coverage tests.c s21_matrix.c s21_matrix.a $(LDFLAGS) -o gcov_tests
	./gcov_tests
	lcov -t "s21_matrix" -o s21_matrix.info -c -d . --rc lcov_branch_coverage=1
	genhtml -o report s21_matrix.info --rc lcov_branch_coverage=1
	open report/index.html
	
clean:
	rm -rf vgcheck ./*.o ./*.a ./*.out ./*.dSYM gcov_test *.html *.css ./GcovReport ./test ./*.gcda ./*.gch ./*.gcno ./report 
	rm -rf ./*.info ./gcov_tests ./gcov_tests ./res_leak.txt ./s21_matrix_test

style:
	clang-format -style=google -n *.c ./bin/*.c
	clang-format -style=google -n *.h ./bin/*.h

make_style:
	clang-format -style=google -i ./bin/*.c *.c bin/*.c
	clang-format -style=google -i ./bin/*.h *.h bin/*.h

#memcheck: clean
#	$(CC) $(CFLAGS) s21_matrix.c -o vgcheck $(LDFLAGS)
#	valgrind -s --leak-check=full --show-leak-kinds=all --tool=memcheck --trace-children=yes --quiet --track-origins=yes --verbose --log-file=res_leak.txt ./vgcheck
