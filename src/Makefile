CXX = g++
CXXFLAGS = -Wall -O3
# Add OpenMP support for clang on macOS
ifeq ($(shell uname),Darwin)
    CXXFLAGS += -Xpreprocessor -fopenmp
    LDFLAGS += -lomp
else
    CXXFLAGS += -fopenmp
endif
PROGS = db_sort set_lcas classify make_seqid_to_taxid_map db_shrink kmer_estimator
TEST_PROGS = test_kraken

.PHONY: all install clean test

all: $(PROGS)

test: $(TEST_PROGS)
	@echo "Running tests..."
	@for test in $(TEST_PROGS); do \
		echo "Running $$test..."; \
		./$$test; \
	done

install: $(PROGS)
	cp $(PROGS) $(KRAKEN_DIR)/

clean:
	rm -f $(PROGS) $(TEST_PROGS) *.o

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: krakendb.o quickfile.o krakenutil.o seqreader.o

kmer_estimator: krakenutil.o seqreader.o

classify: krakendb.o quickfile.o krakenutil.o seqreader.o

make_seqid_to_taxid_map: quickfile.o

test_kraken: krakendb.o quickfile.o test_kraken.o
	$(CXX) $(CXXFLAGS) -o test_kraken test_kraken.o krakendb.o quickfile.o $(LDFLAGS) -lstdc++

krakenutil.o: krakenutil.cpp krakenutil.hpp
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c seqreader.cpp

quickfile.o: quickfile.cpp quickfile.hpp

test_kraken.o: test_kraken.cpp krakendb.hpp kraken_headers.hpp
	$(CXX) $(CXXFLAGS) -c test_kraken.cpp
