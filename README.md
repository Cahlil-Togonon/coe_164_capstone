# Barcode Scanner
A Barcode Decoder and Error Correction according to the specifications found in the pdf files.

Topics include:
- Reed-Solomon Error Correction Codes
- BCH (Bose–Chaudhuri–Hocquenghem) codes
- Extended Euclidean Algorithm
- BM (Berlekamp-Massey) algorithm
- Chien Search
- Forney Algorithm
- Shift-Latch decoder 

# To-Do:
1. ~~fix BM_algorithm~~
2. ~~fix `find_roots()`~~ > used Chien Search
3. ~~compute true message (`correct_message()`)~~
4. ~~fix `error_poly()`~~
5. ~~working code for sample inputs~~
6. ~~working code for Online Judge~~ (yay)

Uses standard input and output:

run `python coe164_cp.py < input.txt > output.txt` at the terminal

**Status**: working with Online Judge