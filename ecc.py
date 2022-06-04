import copy

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

def correct_message(SCV,e_coeffs,root_idxs):                # Compute true message
    corrected_SCV = copy.deepcopy(SCV)
    N = len(corrected_SCV)                                  # get length of SCV
    for i in range(len(root_idxs)):
        idx = N - root_idxs[i] - 1                          # index to correct is from the right starting from index 0
        corrected_SCV[idx] = sum(corrected_SCV[idx],-e_coeffs[i])   # subtract msg with e(x)
    return corrected_SCV                                    # return corrected SCV

################# Main function #######################

def error_correction(ecc_count,SCV):
    t = int(ecc_count / 2)                                          # t = E/2
    
    S_sum, S = syndromes(t,SCV)                                     # get syndromes
    print("Syndromes:",S)                                           # for debugging
    # if not S_sum:
    #     return 0, msg_SCV

    num_errors, error_locator = BM_algorithm(S)                     # get error locator and number of errors using Berlekamp-Massey Algorithm
    print("Error Locator:",error_locator)                           # for debugging

    root_idxs, root_vals, root_no_inv = find_roots(error_locator)   # get root_idxs, root_vals
    print("root_idxs:",root_idxs)                                   # for debugging
    print("root_vals:",root_vals)                                   # for debugging

    e_coeffs = error_poly(S,error_locator,root_no_inv,num_errors)   # get e_coeffs polynomial
    print("E_coeffs:",e_coeffs)                                     # for debugging

    corrected_SCV = correct_message(SCV,e_coeffs,root_idxs)         # get true message

    return num_errors,corrected_SCV                                 # return number of errors, corrected SCV