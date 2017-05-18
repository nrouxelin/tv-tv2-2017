###
#Compiler
###

CXX = g++
FLAGS = -g -Wall
LIBS =  -L/usr/X11R6/lib -lm -lpthread -lX11

###
#Files
###

#Directories
OBJ_DIR = obj
BIN_DIR = bin
SRC_DIR = src

#binaries
EXEC = main_denoiser main_inpainter

#source files
SOURCES = $(wildcard $(SRC_DIR)/[!m]*.cpp)
#Objects
OBJ1 = $(SOURCES:.cpp=.o)#extension
OBJ = $(patsubst $(SRC_DIR)/%, $(OBJ_DIR)/%, $(OBJ1))#directory


###
#Rules
###

#Making objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(FLAGS) $(LIBS) -c $< -o $@

#Making the main program
all: $(OBJ)
	$(CXX) $(FLAGS) $(OBJ) $(LIBS) -o $(BIN_DIR)/main

denoiser: obj/main_denoiser.o $(OBJ)
	$(CXX) $(FLAGS) $< $(OBJ) $(LIBS) -o $(BIN_DIR)/denoiser

inpainter: obj/main_inpainter.o $(OBJ)
	$(CXX) $(FLAGS) $< $(OBJ) $(LIBS) -o $(BIN_DIR)/inpainter

h1: obj/main_h1.o
		$(CXX) $(FLAGS) $<  $(LIBS) -o $(BIN_DIR)/h1
#Cleaning rule
clean:
	rm $(BIN_DIR)/* $(OBJ_DIR)/*
