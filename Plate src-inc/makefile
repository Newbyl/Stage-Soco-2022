CC := g++

CFLAGS := -std=c++17 -Wall -Wextra -pedantic -ggdb -Wno-unused-variable -Wno-unused-parameter -Wunknown-pragmas -O3 -ffast-math

CLIBFLAGS := -std=c++17 -Wall -Wextra -pedantic -ggdb -g -shared -fPIC -O3 -ffast-math

CLI := cli

CLITEST := clitest

LIBPLATE := libPlate.so

all : $(CLI)
test2 : $(CLITEST)
libPlate : $(LIBPLATE)

$(CLI) : unit.o IPlugin.o Cplate.o main.o
	$(CC) $(CFLAGS) -o $@ $^

$(CLITEST) : cplateBackupFiniMaisInverseXY.o unit.o IPlugin.o main.o
	$(CC) $(CFLAGS) -o $@ $^

$(LIBPLATE) : unit.o IPlugin.o Cplate.o
	$(CC) $(CLIBFLAGS) -o $@ $^

%.o : %.cpp
	clear
	$(CC) $(CFLAGS) -o $@ -c $<
	
clean :
	rm *.o $(CLI)
	clear
	
