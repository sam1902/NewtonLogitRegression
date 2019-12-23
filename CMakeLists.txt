cmake_minimum_required(VERSION 3.15)
project(RegLog)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm")

add_executable(RegLog
        Base_apprentissage_Notes_Admis.dat
        Base_apprentissage_Notes_NonAdmis.dat
        main.c gauss.c gauss.h vector_utils.c vector_utils.h functions_init.c functions_init.h marcros.h)
