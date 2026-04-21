CC = g++
FLAGS = -g -Wall -std=c++17

export ROOTLIB=`root-config --glibs`
export ROOTINC=`root-config --incdir`


FIT_C = $(filter-out main.cc $(FIT_CT))
FIT_O = $(FIT_C: .cxx=.o)
FIT_H = $(FIT_C: .cxx=.h)

ANA_C = includes/*.cxx
ANA_O = $(ANA_C: .cxx=.o)
ANA_H = $(ANA_C: .cxx=.h)

default:	prad2_BackgroundSubtraction
clean:
	@echo "Remove files: " ./*.o ./*~ ./main .includes/*.o
	@rm -f ./*.o ./*~ ./prad2_BackgroundSubtraction ./includes/*.o
	@echo "make clean command done..."

./%.o:	./%.cxx ./%.h
	@g++ $(FLAGS) -I $(ROOTINC) -c $< -o $@

main.o: $(ANA_H) main.cxx
	@echo "Compile: $(ROOTINC) $@"
	@g++ $(FLAGS) -I $(ROOTINC) -I includes/ -c main.cxx -o main.o

prad2_BackgroundSubtraction:		main.o $(ANA_O) $(ANAc_O) $(ANAC_O)
	@echo "Link Fit"
	@echo g++ $(FLAGS) -I $(ROOTINC) -I includes/ -o prad2_BackgroundSubtraction $(ANA_O) $(ANAc_O) $(ANAC_O) main.o $(ROOTLIB) -lGenVector -lsqlite3
	@g++ $(FLAGS) -I $(ROOTINC) -I includes/ -o prad2_BackgroundSubtraction $(ANA_O) $(ANAc_O) $(ANAC_O) main.o  $(ROOTLIB) -lGenVector -lsqlite3
	@echo "Compilation Done"