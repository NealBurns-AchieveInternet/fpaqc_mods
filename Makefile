CXX=g++
#CXXFLAGS=-O3 -fomit-frame-pointer -s -DNDEBUG
CXXFLAGS=-O3 -fomit-frame-pointer -s

all: fpaqc my_fpaqc test

test: DATA.fpaqc DATA.my_fpaqc
	@cmp DATA.fpaqc DATA.my_fpaqc && echo Match!

#$(RM) FOO.* && { for c in my_fpaqc fpaqc ; do time ./$$c c DATA FOO.$$c ; done ; } && cmp FOO.* && echo "Match!"

DATA.%: % DATA
	./$< c DATA $@

clean:
	$(RM) fpaqc{,.exe} my_fpaqc{,.exe} DATA.*
