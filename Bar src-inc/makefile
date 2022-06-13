CC := g++

CFLAGS := -std=c++17 -Wall -Wextra -pedantic -ggdb -Wno-unused-variable -Wno-unused-parameter -Wunknown-pragmas

# We concatenate the optimisation flags with CFLAGS for practical reasons.
CFLAGS += -O3 -ffast-math

CLI := cli
TEST := test

all : $(CLI) $(TEST)

test : $(TEST)
cli : $(CLI)

$(CLI) : unit.o IPlugin.o Cbar.o main.o
	$(CC) $(CFLAGS) -o $@ $^

$(TEST) : test.o
	$(CC) $(CFLAGS) -o $@ $^

%.o : %.cpp
	clear
	$(CC) $(CFLAGS) -o $@ -c $<
	
clean :
	rm *.o $(CLI) $(TEST)
	clear

run :
	make cli
	./cli
	
