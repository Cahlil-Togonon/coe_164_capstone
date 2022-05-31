def dataToDecode(scvList):
    m_data = []
    scvValues = []
    for i in range (1,scvList[1]):
        m_data.append(scvList[i+1])
    while m_data[-1]==900:
        m_data.pop()
    for i in range(len(m_data)):
        scvValues.append(m_data[i]//30)
        scvValues.append(m_data[i] % 30)
    while scvValues[-1] == 29:
        scvValues.pop()
    print(scvValues)
    return scvValues


def decode(scvValues):
    decoded = ""
    state = alpha
    prevState = alpha
    tempState = alpha

    shift = False
    wasShift = False

    for i in range(len(scvValues)):
        if shift:
            wasShift = True
            state = tempState
         
        if i == 0 and scvValues[i] == 27:
            state = lower
            i += 1
        elif i == 0 and scvValues[i] == 28:
            state = mixed    
            i += 1
        elif i == 0 and scvValues[i] == 29:
            prevState = state 
            shift = True
            tempState = punc
            i += 1

        elif i != 0 and state[scvValues[i]] == "La":
            state = alpha
            i += 1 
        elif i != 0 and state[scvValues[i]] == "LL":
            state = lower
            i += 1 
        elif i != 0 and state[scvValues[i]] == "Lm":
            state = mixed
            i += 1 
        elif i != 0 and state[scvValues[i]] == "Lp":
            state = punc
            i += 1 
        elif i != 0 and state[scvValues[i]] == "Sa":
            tempState = alpha
            shift = True
            i += 1 
            
        elif i != 0 and state[scvValues[i]] == "Sp":
            prevState = state
            tempState = punc
            shift = True
            i += 1 
        else:
            decoded += state[scvValues[i]]

        if wasShift:
            state = prevState
            wasShift = False
            shift = False

    return decoded

if __name__ == "__main__":
    testCase1 = [0, 5, 237, 269, 900, 900, 64, 152, 316, 398]
    testCase2 = [1, 7, 87, 447, 146, 841, 184, 900, 879, 523]
    testCase3 = [4, 10, 893, 864, 877, 749, 739, 496, 844, 393, 900, 822, 22, 761, 545, 596, 130, 458, 768]
    
    alpha = list(map(chr, range(65,91)))
    alpha.extend([" ", "LL", "Lm","Sp"])
    lower = list(map(chr,range(97,123)))
    lower.extend([" ", "Sa", "Lm", "Sp"])
    mixed = list(map(chr,range(48,58)))
    mixed.extend(map(chr,[38, 13, 9, 44, 58, 35, 45, 46, 36, 47, 43, 37, 42, 61, 94]))
    mixed.extend(["Lp", " ", "LL", "La", "Sp"])
    punc = list(map(chr,[59, 60, 62, 64, 91, 92, 93, 95, 96, 126, 33, 13, 9, 44, 58, 10, 45, 46, 36, 47, 34, 124, 42, 40, 41, 63, 123, 125, 39]))
    punc.append("La")

    print(decode(dataToDecode(testCase1)))
    print(decode(dataToDecode(testCase2)))
    print(decode(dataToDecode(testCase3)))