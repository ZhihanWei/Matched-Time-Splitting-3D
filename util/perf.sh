sudo perf record -F 1000 -g --call-graph dwarf ../build/bin/Chestnut
sudo perf script | c++filt | ./stackcollapse-perf.pl | ./flamegraph.pl > out.svg 
google-chrome out.svg
