# Copyright (C) Sorbonne Université (2018)
# Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>
#
# This file is part of ActiveBrownian.
#
# ActiveBrownian is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ForceInCrystal is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ForceInCrystal.  If not, see <http://www.gnu.org/licenses/>.

CC=g++
CFLAGS=-W -Wall -ansi -pedantic -std=c++14 -O3
LDFLAGS=-lm -pthread -lboost_program_options \
        -lsfml-graphics -lsfml-window -lsfml-system
EXEC=ActiveBrownian
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.cpp %.h
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -f *.o

mrproper: clean
	rm -rf $(EXEC)
