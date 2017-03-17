OBJS = main.o spline.o
CC = g++
CFLAGS = -O -g
main:$(OBJS)
	$(CC) $(OBJS) -o main
main.o:main.cpp spline.h
	$(CC) $(CFLAGS) -c $< -o main.o
spline.o:spline.cpp spline.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm *.o main