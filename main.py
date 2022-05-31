import sys
from ecc import error_correction, pow

# submodes mappings where idx:ascii_value
# control values are strings, where L/S for latch/shift, a/L/m/p for submodes
alpha = [i + 65 for i in range(26)] + [32,"LL","Lm","Sp"]
lower = [i + 97 for i in range(26)] + [32,"Sa","Lm","Sp"]
mixed = [i + 48 for i in range(10)] + [38,13,9,44,58,35,45,46,36,47,43,37,42,61,94,"Lp",32,"LL","La","Sp"]
punct = [59,60,62,64,91,92,93,95,96,126,33,13,9,44,58,10,45,46,36,47,34,124,42,40,41,63,123,125,39,"La"]
codebook = {'a':alpha,'L':lower,'m':mixed,'p':punct}    # map submodes to letters corresponding to control value

def decoder(data_SCV):
    decoded = ""
    symbol_SCV = []
    for val in data_SCV:                            # change data to H and L values 
        symbol_SCV += [(val // 30) % 30, val % 30]          # *remember to remove mod from H
                                                            # (debugging for values > 900 because no error correction)
    if symbol_SCV[-1] == 29: symbol_SCV.pop()       # remove last value if 29

    shift = None                                    # shift flag
    mode = 'a'                                      # start at alpha submode
    for val in symbol_SCV:
        code = codebook[mode][val]                  # get the ascii_value or control value
        if shift:
            code = codebook[shift][val]             # if shift, get the shift's ascii_value instead
            shift = None
        if isinstance(code, str):                   # check if string (control value)
            if code[0] == 'L':                      # if string starts with 'L'
                mode = code[1]                      # change mode to the submode (second letter)
            elif code[0] == 'S':                    # if string starts with 'S'
                shift = code[1]                     # change shift to the submode (second letter)
        else:
            decoded += chr(code)                    # concatenate decoded string with decoded ascii symbol using chr(ascii_value)
    
    return decoded                                  # return decoded message

if __name__ == '__main__':
    f = open(sys.argv[1],"r")
    input_stream = f.read().splitlines()            # split whole input text file (without \n)
    f.close()

    T = int(input_stream.pop(0))                    # get first value which is T, the number of test cases
    f = open(sys.argv[2],"w")
    for t in range(T):
        ecc_level,N = [int(i) for i in input_stream.pop(0).split()]     # get E (ecc_level) and N
        SCV = [int(i) for i in input_stream.pop(0).split()]             # get SCV array

        ecc_count = int(pow(2,ecc_level+1))                             # calculate ecc_count
        
        print(f"Case #{t+1}:") # for debugging purposes
        num_errors,corrected_SCV = error_correction(ecc_count,SCV)      # call error correction
        
        error_SCV = corrected_SCV[-ecc_count:]                          # get error part of SCV
        msg_SCV = corrected_SCV[:N - ecc_count]                         # get msg part of SCV

        sym_length = msg_SCV[0]                                         # first val of msg_SCV is the symbol length descriptor
        data_SCV = msg_SCV[1:]                                          # data otherwise
        while data_SCV[-1] == 900:                                      # remove padding values of '900'
            data_SCV.pop()

        decoded = decoder(data_SCV)                                     # decode data

        line = f"{num_errors} "
        for val in corrected_SCV:                                       # line formatting:
            line += f"{val} "                                           # (num_errors corrected_SCV)

        f.writelines(f"Case #{t+1}:\n")                                 # write test case number, line, and decoded message
        f.writelines(line[:-1]+"\n")
        f.writelines(decoded+'\n')
    f.close()