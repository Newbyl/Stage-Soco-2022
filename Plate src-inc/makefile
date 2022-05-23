CC := g++

CFLAGS := -std=c++17 -Wall -Wextra -pedantic -ggdb -Wno-unused-variable -Wno-unused-parameter -Wunknown-pragmas -O3 -ffast-math

CLI := cli

CTBB := -ltbb

all : $(CLI)

$(CLI) : unit.o IPlugin.o Cplate.o main.o
	$(CC) $(CFLAGS) $(CTBB) -o $@ $^

%.o : %.cpp
	clear
	$(CC) $(CFLAGS) -o $@ -c $<
	
clean :
	rm *.o $(CLI)
	clear
	
