mkdir build
copy input.csv build\\input.csv
gcc.exe -o build/cpa.exe main.cpp libbmp.c fread_csv_line.c split.c font8x8.c csv.c cJSON.c -lstdc++

pause