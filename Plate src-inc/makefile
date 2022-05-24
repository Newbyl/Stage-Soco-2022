CC := g++

CFLAGS := -std=c++17 -Wall -Wextra -pedantic -ggdb -Wno-unused-variable -Wno-unused-parameter -Wunknown-pragmas -O3 -ffast-math

CLI := cli

CLITEST := clitest

CTBB := -ltbb

all : $(CLI)
test2 : $(CLITEST)

$(CLI) : unit.o IPlugin.o Cplate.o main.o
	$(CC) $(CFLAGS) $(CTBB) -o $@ $^

$(CLITEST) : cplateBackupFiniMaisInverseXY.o unit.o IPlugin.o main.o
	$(CC) $(CFLAGS) -o $@ $^

%.o : %.cpp
	clear
	$(CC) $(CFLAGS) -o $@ -c $<
	
clean :
	rm *.o $(CLI)
	clear
	
