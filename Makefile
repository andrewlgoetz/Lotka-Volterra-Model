demo: demo.o
	c++ -o demo demo.o -ltrapfpe -lpgplot -lcpgplot -lX11 -lm 

demo.o: demo.cpp
	c++ -c demo.cpp

pred_prey: pred_prey.o
	c++ -o pred_prey pred_prey.o -ltrapfpe -lpgplot -lcpgplot -lX11 -lm

pred_prey.o: pred_prey.cpp
	c++ -c pred_prey.cpp

