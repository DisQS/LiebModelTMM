EXE = ../src/tmseLMxD.GF

all:	lieb21 lieb22 lieb23 lieb31 lieb32 lieb33 diff
diff:	L21_0004_D0050.raw L22_0004_D0050.raw L23_0004_D0050.raw L31_0004_D0850.raw L32_0004_D0150.raw L33_0004_D0150.raw
 
lieb21: $(EXE) 
	cd lieb21; ../../src/tmseLMxD.GF | tee terminal_log.txt; cat terminal_log.txt
L21_0004_D0050.raw:	
	cd lieb21; diff -W132 -s L21_0004_D0050.raw ./target/L21_0004_D0050.raw

lieb22: $(EXE) 
	cd lieb22; ../../src/tmseLMxD.GF | tee terminal_log.txt; cat terminal_log.txt
L22_0004_D0050.raw: 
	cd lieb22; diff -W132 -s L22_0004_D0050.raw ./target/L22_0004_D0050.raw

lieb23: $(EXE) 
	cd lieb23; ../../src/tmseLMxD.GF | tee terminal_log.txt; cat terminal_log.txt
L23_0004_D0050.raw: 
	cd lieb23; diff -W132 -s L23_0004_D0050.raw ./target/L23_0004_D0050.raw

lieb31: $(EXE) 
	cd lieb31; ../../src/tmseLMxD.GF | tee terminal_log.txt; cat terminal_log.txt
L31_0004_D0850.raw: 
	cd lieb31; diff -W132 -s L31_0004_D0850.raw ./target/L31_0004_D0850.raw

lieb32: $(EXE) 
	cd lieb32; ../../src/tmseLMxD.GF | tee terminal_log.txt; cat terminal_log.txt
L32_0004_D0150.raw: 
	cd lieb32; diff -W132 -s L32_0004_D0150.raw ./target/L32_0004_D0150.raw

lieb33: $(EXE) 
	cd lieb33; ../../src/tmseLMxD.GF | tee terminal_log.txt; cat terminal_log.txt
L33_0004_D0150.raw: 
	cd lieb33; diff -W133 -s L33_0004_D0150.raw ./target/L33_0004_D0150.raw

clean:
	rm -f core *.mod *.o lieb*/*.raw lieb*/*.tmp





