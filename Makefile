CC = cc
VPATH = src include
FLAGS = 

all: gendiplo mongrail2

debug: FLAGS += -g -Wall
debug: all

gendiplo: gendiplo.c
	$(CC) $(FLAGS) -o gendiplo $< -I include -lm

mongrail2: mongrail2.c
	$(CC) $(FLAGS) -o mongrail2 $< -lm `pkg-config --cflags --libs glib-2.0`


clean:
	$(RM) mongrail2 gendiplo
