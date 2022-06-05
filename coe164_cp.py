import copy
import sys

########################################################

# Error Correction: Adrian Cahlil Eiz G. Togonon, 2019-11731
# Barcode Decoding: Alquea Pauline x. Macatangay, 2019-xxxxx

################ Helper functions ######################

def pow(base,power):                # power with mod 929 every step
    result = 1
    while power:
        result *= base
        result %= 929
        power -= 1
    return result

def mult(a,b):                      # mult values with mod 929
    return (a*b) % 929

def sum(a,b):                       # sum values with mod 929
    return (a+b) % 929

def updateCoeff(a,b,s,t):           # for Extended Euclidean Algorithm 
    return (t-(b//a)*s, s)

def EEA(a,b):                       # Extended Euclidean Algorithm for division
    if a == 0:                      # base case
        return b,0,1
    gcd,s1,t1 = EEA(b % a, a)       # recursive call
    s,t = updateCoeff(a,b,s1,t1)    # update s and t
    return gcd,s,t                  # t can be used as the inverse of a number in a finite field

######### Actual Error Correction functions ############

def syndromes(t,SCV):                       # calculate syndromes given t and the SCV
    S_sum = 0                               # sum of all syndromes
    S = []                                  # syndromes vector
    for i in range(1,2*t+1):                # for i = 1 to 2t or E (ecc_count)
        Si = 0
        for j in range(len(SCV)):           # Si = Summation of SCV[j]*x^(n-j), where x = a^i, a = 3
            Si = sum(Si,mult(SCV[j],pow(3,i*(len(SCV)-1-j))))
        S_sum = sum(S_sum,Si)               # S_sum = Summation(Si)
        S.append(Si)
    return S_sum, S                         # return syndromes sum, syndromes vector

def BM_algorithm(S):                        # BM algo as implemented from whiteboard notes
    Lx, Lpx = [1],[1]                       # initialize L(x) = Lp(x) = 1
    ne,dp,m = 0,1,1                         # initialize num_errors = 0, prev. discrepancy = 1, iteration mark = 1
    for i in range(len(S)):
        d = S[i]                            # initialize di = S[i]
        for k in range(1,ne+1):             # for k = 1 to ne
            d = sum(d,mult(Lx[k],S[i-k]))   # di = S[i] + Summation(Lx[k],S[i-k]) from k = 1 to ne
        if d == 0:
            m += 1
        else:
            temp = copy.copy(Lx)            # temp copy old L(x)
            sLpx = [0]*m + Lpx              # shift Lp(x) by x^m
            for k in range(len(sLpx)):
                _,_,t = EEA(929,dp)         # get inverse of dp -> t
                val = mult(d,t)             # mult inverse of dp to d
                sLpx[k] = mult(sLpx[k],val) # mult with Lpx[k]
                if k < len(Lx):             # L(x) = L(x) - (d/dp) x^m Lp(x)
                    Lx[k] = sum(Lx[k],-sLpx[k])
                else:
                    Lx += [-sLpx[k] % 929]  # L(x) = 0 - (d/dp) x^m Lp(x)
            if 2*ne <= i:
                Lpx = temp                  # copy old L(x) to Lp(x)
                ne = i + 1 - ne             # update num_errors
                dp = d                      # update previous discrepancy with new one
                m = 1
            else:
                m += 1
    return ne,Lx                            # return num_errors, L(x)

def find_roots(Lx):                         # find roots of L(x) = 1 + L0 x^1 + ... using Chien Search
    Le = copy.deepcopy(Lx)                  # copy coeffs of L(x) to Le
    Tx = [pow(3,i) for i in range(len(Lx))] # calculate template polynomial T(x) = 1 + a^1 + a^2 + ...
    root_idxs = []                          # the index = exponent of the root a^x
    root_vals = []                          # the value of a^x itself
    root_no_inv = []                        # also take note of non-inverted roots (for error polynomial)
    for i in range(0,928):
        elp_val = 0
        for j in range(len(Le)):            # elp_val = Summation of the current Le
            elp_val = sum(elp_val,Le[j])
            Le[j] *= Tx[j]                  # update the values of Le *= T(x) as we sweep
        if elp_val == 0:                    # then i is root
            inv = -i % 928                  # get inverse of i (distance from 928)
            root_idxs.append(inv)           # inverse is the exponent/index
            root_vals.append(pow(3,inv))    # evaluate raw value of a^inv
            root_no_inv.append(pow(3,i))    # also get the non-inverted root values
    return root_idxs, root_vals, root_no_inv                # return root indices, raw root values, and non-inverted roots

def error_poly(S,error_locator,root_no_inv,ne):             # calculate error polynomial from error locator
    dL = []
    for i in range(1,len(error_locator)):                   # get the derivate of L(x) -> dL(x)
        dL.append(mult(error_locator[i],i))                 # L(x) = 1 + L1 x^1 + L2 x^2 + ... = 1 + L[i] x^i
                                                            # dL(x) = L1*1 + L2*2 * x + ...    = L[i]*i x^(i-1)
    
    Ox = [0]*(ne)                                           # initialize O(x) = 0 + 0 + ... ne times (basically mod ne)
    for i in range(len(Ox)):                                # O(x) = S(x)L(x) mod ne => FOIL method
        for j in range(len(Ox)):                            # is like combinatorics of i and j: Ox[i+j] += S[i] * L[j] until ne
            if (i+j) < len(Ox):
                Ox[i+j] = sum(Ox[i+j],mult(S[i],error_locator[j]))

    e_coeffs = []
    for root in root_no_inv:                                # for every root value
        Ox_val = 0
        dL_val = 0
        for i in range(len(Ox)):                            # calculate Ox_val = O(root)
            Ox_val = sum(Ox_val,mult(Ox[i],pow(root,i)))    # Ox_val = Summation of O[i] * root^i
        for i in range(len(dL)):                            # calculate dL_val = dL(root)
            dL_val = sum(dL_val,mult(dL[i],pow(root,i)))    # dL_val = Summation of dL[i] * root^i
        _,_,dL_val = EEA(929,dL_val)                        # inverse dL_val
        e_coeffs.append(mult(-Ox_val,dL_val))               # e_coeff = -(Ox_val/dL_val)
    return e_coeffs                                         # return error polynomial

def correct_message(SCV,e_coeffs,root_idxs):                        # Compute true message
    corrected_SCV = copy.deepcopy(SCV)
    N = len(corrected_SCV)                                          # get length of SCV
    for i in range(len(root_idxs)):
        idx = N - root_idxs[i] - 1                                  # index to correct is from the right starting from index 0
        corrected_SCV[idx] = sum(corrected_SCV[idx],-e_coeffs[i])   # subtract msg with e(x)
    return corrected_SCV                                            # return corrected SCV

################# Main Error Correction function #######################

def error_correction(ecc_count,SCV):
    t = int(ecc_count / 2)                                          # t = E/2
    
    S_sum, S = syndromes(t,SCV)                                     # get syndromes
    # print("Syndromes:",S)                                         # for debugging
    # if not S_sum:
    #     return 0, msg_SCV

    num_errors, error_locator = BM_algorithm(S)                     # get error locator and number of errors using Berlekamp-Massey Algorithm
    # print("Error Locator:",error_locator)                         # for debugging

    root_idxs, root_vals, root_no_inv = find_roots(error_locator)   # get root_idxs, root_vals
    # print("root_idxs:",root_idxs)                                 # for debugging
    # print("root_vals:",root_vals)                                 # for debugging

    e_coeffs = error_poly(S,error_locator,root_no_inv,num_errors)   # get e_coeffs polynomial
    # print("E_coeffs:",e_coeffs)                                   # for debugging

    corrected_SCV = correct_message(SCV,e_coeffs,root_idxs)         # get true message

    return num_errors,corrected_SCV                                 # return number of errors, corrected SCV

################# Decoder function #######################
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
        symbol_SCV += [(val // 30), val % 30]       # H = val//30, L = val % 30
    if symbol_SCV[-1] == 29: symbol_SCV.pop()       # remove last value if 29

    shift = False                                   
    wasShift = False
    mode = 'a'                                      # start at alpha submode
    prevmode = 'a'
    tempmode = 'a'
    
    for val in symbol_SCV:
        if shift:                                   # if shift, change mode to temporary mode
            wasShift = True             
            mode = tempmode
            
        code = codebook[mode][val]                  # get the ascii_value or control value
        
        if isinstance(code, str):                   # check if string (control value)
            if code[0] == 'L':                      # if string starts with 'L'
                mode = code[1]                      # change mode to the submode (second letter)
            elif code[0] == 'S':                    # if string starts with 'S'
                prevmode = mode
                tempmode = code[1]                  # change temporary mode to the submode (second letter)
                shift = True       
        else:
            decoded += chr(code)                    # concatenate decoded string with decoded ascii symbol using chr(ascii_value)
        
        if wasShift:                                    
            mode = prevmode                         # change back to previous mode
            wasShift = False
            shift = False
            
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

        ecc_count = int(2**(ecc_level+1))                               # ecc_count = 2^(ecc_level + 1)
        
        # print(f"Case #{t+1}:") # for debugging purposes
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
