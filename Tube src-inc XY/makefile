CC := g++

CFLAGS := -std=c++17 -Wall -Wextra -pedantic -ggdb -Wno-unused-variable -Wno-unused-parameter -Wunknown-pragmas -O3 -ffast-math

CLI := cli

all : $(CLI)

$(CLI) : unit.o IPlugin.o Ctube.o main.o
	$(CC) $(CFLAGS) -o $@ $^

%.o : %.cpp
	clear
	$(CC) $(CFLAGS) -o $@ -c $<
	
clean :
	rm *.o $(CLI)
	clear
	
