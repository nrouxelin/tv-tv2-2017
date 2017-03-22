###
#Compiler
###

CC = g++
FLAGS = -g -Wall
LIBS =  -L/usr/X11R6/lib -lm -lpthread -lX11

###
#Files
###

#Directories
OBJ_DIR = obj
BIN_DIR = bin
SRC_DIR = src

#C files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
#Objects
OBJ1 = $(SOURCES:.cpp=.o)#extension
OBJ = $(patsubst $(SRC_DIR)/%, $(OBJ_DIR)/%, $(OBJ1))#directory


###
#Rules
###

#Making objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(FLAGS) $(LIBS) -c $< -o $@

#Making the main program
all: $(OBJ)
	$(CC) $(FLAGS) $(OBJ) $(LIBS) -o $(BIN_DIR)/main


#Cleaning rule
clean:
	rm $(BIN_DIR)/* $(OBJ_DIR)/* $(FOBJ_DIR)/*
