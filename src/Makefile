#CC = mpCC
#CC_OPTIONS = -q64 -DDEBUG
#LNK_OPTIONS = -q64
CC = mpicc
CC_OPTIONS =   -O2 -Wall -DDEBUG #-DDEEP_DEBUG
LNK_OPTIONS = 
EXEC = mine.x
LIBS = -lm 



OBJS = mtwist.o vector.o myrand.o RND_lyman.o struct.o io.o parse.o ray.o lyman.o propagation.o main.o 
INCL = Makefile struct.h vector.h mtwist.h myrand.h RND_lyman.h ray.h io.h lyman.h propagation.h parse.h

.SUFFIXES: .o .or .c

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

all: $(EXEC)

$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) $(LNK_OPTIONS) -o $(EXEC)

clean:
	rm -f $(OBJS) *~ core* ${EXEC}
