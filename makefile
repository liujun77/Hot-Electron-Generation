object = hot_ele.o

hotele : $(object)
	icpc -openmp -o hotele hot_ele.o

hot_ele.o : hot_ele.cpp lib/*.h
	icpc -openmp -c hot_ele.cpp

.PHONY : clean
clean : 
	rm $(object)
