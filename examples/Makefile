CPPFLAGS = -I../include -L../filzbach -DFILZBACH_LIB -std=c++11

examples = normal priors mixednormal poisson multispp lr regression stepfunction plantgrowth seedrain SDM ProteinInteraction TreeMort
deps = examples.h preamble.h ../include/filzbach.h ../filzbach/libfilzbach.a

mymodel: mymodel.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< -lfilzbach

normal: eg01_normal.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=1 -o $@ $< -lfilzbach

priors: eg02_normal_priors.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=2 -o $@ $< -lfilzbach

mixednormal: eg03_mixednormal.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=3 -o $@ $< -lfilzbach

poisson: eg04_poisson.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=4 -o $@ $< -lfilzbach

multispp: eg05_poisson_multispp.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=5 -o $@ $< -lfilzbach

lr: eg07_lr.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=7 -o $@ $< -lfilzbach

regression: eg08_logistic_regression.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=8 -o $@ $< -lfilzbach

stepfunction: eg09_stepfunction.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=9 -o $@ $< -lfilzbach

plantgrowth: egz1_plantgrowth.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=11 -o $@ $< -lfilzbach

seedrain: egz2_seedrain.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=12 -o $@ $< -lfilzbach

SDM: egz3_SDM.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=13 -o $@ $< -lfilzbach

ProteinInteraction: egz6_ProteinInteraction.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=16 -o $@ $< -lfilzbach

TreeMort: egz7_TreeMort.cpp $(deps)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DMODEL=17 -o $@ $< -lfilzbach

../filzbach/libfilzbach.a:
	$(MAKE) -C ../filzbach

examples: $(examples)

clean:
	rm -rf workspace mymodel $(examples)
	$(MAKE) clean -C ../filzbach

