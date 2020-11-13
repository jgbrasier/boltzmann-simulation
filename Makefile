#http://clang.llvm.org/docs/AddressSanitizer.html

CXXFLAGS= -O -g -Wextra -Wall -Wshadow   -fsanitize=address -Wno-unused-variable -I/opt/local/include -I/opt/X11/include

all:  start

start:start.o graphics.o Graphics.h
	$(CXX) $(CXXFLAGS) graphics.o start.o -o $@  -L/opt/local/lib -lcairo  -L/opt/X11/lib  -lX11 -lm 


hs.o: hs.cc Graphics.h
start.o: start.cc Graphics.h
graphics.o: Graphics.h graphics.cc

clean:
	$(RM) *.o hs start
clobber: clean
	$(RM)  *.dat projet.tar start hs *.png


tar:
	$(RM) -rf MD
	mkdir MD
	cp Makefile graphics.cc Graphics.h start.cc MD/
	tar cvfz MDprojet.tar MD 
	scp MDprojet.tar acm@turner.pct.espci.fr:public_html/psa135/
	$(RM) -r MD
alltar:
	make clobber
	mkdir MD
	cp *c Makefile *h  MD
	tar cvfz All.tar.gz MD
	mailit All.tar.gz
	$(RM) -r MD All.tar.gz

hs:hs.o graphics.o Graphics.h
	$(CXX) $(CXXFLAGS) graphics.o hs.o -o $@  -L/opt/local/lib -lcairo  -L/opt/X11/lib  -lX11 -lm 
hstar:
	mkdir Attract
	cp Makefile *.cc *.h Attract
	tar cvfz At.tar.gz Attract
