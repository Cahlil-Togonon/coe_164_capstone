import copy
import numpy

################ Helper functions ######################

def pow(base,power):        # power with mod 929 every step
    result = 1
    while power:
        result *= base
        result %= 929
        power -= 1
    return result

def mult(a,b):              # mult values with mod 929
    return (a*b) % 929

def sum(a,b):               # sum values with mod 929
    return (a+b) % 929

def updateCoeff(a,b,s,t):   # for Extended Euclidean Algorithm 
    return (t-(b//a)*s, s)

def EEA(a,b):               # Extended Euclidean Algorithm for division
    if a == 0:
        return b,0,1
    gcd,s1,t1 = EEA(b % a, a)
    s,t = updateCoeff(a,b,s1,t1)
    return gcd,s,t

######### Actual Error Correction functions ############

def syndromes(t,SCV):                       # calculate syndromes given t and the SCV
    S_sum = 0                               # sum of all syndromes
    S = []                                  # Syndromes vector
    for i in range(1,2*t+1):                # calculate syndromes
        Si = 0
        for j in range(len(SCV)):           # Si = Summation of SCV[j]*x^(n-j), where x = a^i, a = 3
            Si = sum(Si,mult(SCV[j],pow(3,mult(i,len(SCV)-1-j))))
        S_sum = sum(S_sum,Si)               # S_sum = Summation(Si)
        S.append(Si)
    return S_sum, S                         # return syndromes sum, syndromes vector

def BM_algo(S):                             # BM algo as implemented from whiteboard notes (shows wrong values for test case #3)
    Lx, Lpx = [1],[1]                       # initialize L(x) = Lp(x) = 1
    ne,dp,m = 0,1,1
    for i in range(len(S)):
        d = S[i]                        
        for k in range(1,ne+1):
            d = sum(d,mult(Lx[k],S[i-k]))   # di = S[i] + Summation(Lx[k],S[i-k]) from k = 1 to ne
        if d == 0:
            m += 1
        else:
            temp = copy.copy(Lx)            # temp copy old L(x)
            Lx = Lx + [0]*m                 # pad some zeroes for calculation
            Lpx = [0]*m + Lpx               # shift Lp(x) by x^m
            for k in range(len(Lpx)):
                _,_,t = EEA(929,dp)         # get inverse of dp -> t
                val = mult(d,t)             # mult inverse of dp to d
                val = mult(Lpx[k],val)      # mult with Lpx[k]
                Lx[k] = sum(Lx[k],-val)     # L(x) = L(x) - (d/dp) x^m Lp(x)
            if 2*ne <= i:
                Lpx = temp
                ne = i + 1 - ne
                dp = d
                m = 1
            else:
                m += 1
    return ne,Lx                            # return number of errors, L(x)

def berlekamp_massey_algorithm(data):       # taken from https://gist.github.com/StuartGordonReid/a514ed478d42eca49568
    n = len(data)                           # correct values for test case #3 but missing 0 (at x^3)
    C = numpy.zeros(n)                      # nearly the same code, not sure where it is different
    B = numpy.zeros(n)
    C[0], B[0] = 1, 1
    l, m, b, i = 0, -1, 1, 0
    while i < n:
        V = data[(i - l):i]
        V = V[::-1]
        CC = C[1:l + 1]
        d = (data[i] + numpy.dot(V, CC))
        if d != 0:
            temp = copy.deepcopy(C)
            P = numpy.zeros(n)
            for j in range(0,l):
                val = mult(B[j],d)
                _,_,t = EEA(929,b)
                val = mult(val,t)
                P[j + i - m] = - val
            for k in range(n):
                C[k] = sum(C[k],P[k])
            # C = (C + P)
            if l <= 0.5 * i:
                l = i + 1 - l
                B = temp
                b = d
                m = i
        i += 1
    return l,C                              # return number of errors, L(x)

def find_roots(Lx):                         # find roots of L(x) = 1 + L0 x^1 + ...
    elp_roots = []
    roots_idx = []
    for x in range(929):                    # evaluate L(x) at every a^x value (brute force method)
        val = 0
        for i in range(len(Lx)):            # val = Summation of Lx[i] * a^(x*i)
            val = sum(val,mult(Lx[i],pow(3,mult(x,i))))
        if val == 0:
            elp_roots.append(x)             # if val is 0, then x is an elp_root
    
    for i in elp_roots:                     # get root_idx per elp_root
        _,_,t = EEA(929,pow(3,i))           # get inverse of elp_root
        if t < 0:                           # ex. [925] -> [27]
            t += 929                        # not sure how [a^-925] -> [a^3], where 3 should be the indices we need***
        roots_idx.append(t)                 
    return elp_roots, roots_idx             # return elp_roots, root_idx

def error_poly(S,error_locator,roots,ne):                   # calculate error polynomial from error locator
    dL = []
    for i in range(1,len(error_locator)):                   # get the derivate of L(x) -> dL(x)
        dL.append(mult(error_locator[i],i))                 # L(x) = 1 + L1 x^1 + L2 x^2 + ... = 1 + L[i] x^i
                                                            # dL(x) = L1*1 + L2*2 * x + ...    = L[i]*i x^(i-1)
    
    Ox = [0]*(ne+1)                                         # initialize O(x) = 0 + 0 + ... ne+1 times (basically mod ne)
    for i in range(len(Ox)):                                # O(x) = S(x)L(x) mod ne =>
        for j in range(len(Ox)):                            # is like combinatorics: Ox[i+j] += S[i] * L[j] until ne+1
            if (i+j) < len(Ox):
                Ox[i+j] = sum(Ox[i+j],mult(S[i],error_locator[j]))

    e_coeffs = []
    for root in roots:                                      # for every root
        Ox_val = 0
        dL_val = 0
        for i in range(len(Ox)):                            # calculate Ox_val = O(root)
            Ox_val = sum(Ox_val,mult(Ox[i],pow(root,i)))    # Ox_val = Summation of O[i] * a^(root*i)
        for i in range(len(dL)):                            # calculate dL_val = dL(root)
            dL_val = sum(dL_val,mult(dL[i],pow(root,i)))    # dL_val = Summation of dL[i] * a^(root*i)
        _,_,dL_val = EEA(929,dL_val)                        # inverse dL_val
        e_coeffs.append(mult(-Ox_val,dL_val))               # e_coeff = -(Ox_val/dL_val)
    return e_coeffs                                         # return error polynomial

################# Main function #######################

def error_correction(ecc_count,SCV):
    t = int(ecc_count / 2)
    
    S_sum, S = syndromes(t,SCV)                                 # get syndromes
    print("Syndromes:",S)                                       # for debugging
    # if not S_sum:
    #     return 0, msg_SCV

    num_errors, error_locator = BM_algo(S)                      # get error locator and number of errors using Berlekamp-Massey Algorithm
    # num_errors, error_locator = berlekamp_massey_algorithm(S) # **other implementation
    print("Error Locator:",error_locator)                       # for debugging

    elp_roots, roots_idx = find_roots(error_locator)            # get elp_roots, root_idx
    print("elp_roots:",elp_roots)                               # for debugging
    print("roots_idx:",roots_idx)                               # for debugging

    e_coeffs = error_poly(S,error_locator,roots_idx,num_errors) # get e_coeffs polynomial
    print("E_coeffs:",e_coeffs)                                 # for debugging

    # subtract SCV with e_coeffs part

    return num_errors,SCV                                       # return number of errors, corrected SCV