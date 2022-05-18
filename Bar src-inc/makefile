CC := g++

CFLAGS := -std=c++17 -Wall -Wextra -pedantic -ggdb -Wno-unused-variable -Wno-unused-parameter -Wunknown-pragmas

OPTFLAGS := -O3 -ffast-math

CLI := cli
TEST := test

all : $(CLI) $(TEST)

test : $(TEST)
cli : $(CLI)

$(CLI) : unit.o IPlugin.o Cbar.o main.o
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ $^

$(TEST) : test.o
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ $^

%.o : %.cpp
	clear
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $<
	
clean :
	rm *.o $(CLI) $(TEST)
	clear

run :
	make cli
	./cli
	
