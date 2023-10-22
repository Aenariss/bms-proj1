CC = g++
CFLAGS = -Wall -Werror
TARGET = bms

FILE_NAMES = bms
sources = $(FILE_NAMES:=.cpp)
objects = $(FILE_NAMES:=.o)

.PHONY: all clean pack

all: $(TARGET)

$(TARGET): $(objects)
		$(CC) $(CFLAGS) -o $(TARGET) $^
		$(MAKE) clean

$(objects): $(sources)
		$(CC) $(CFLAGS) -c $^

clean: $(objects)
		rm -f $^

pack: all
	zip 221701.zip Makefile $(sources)
