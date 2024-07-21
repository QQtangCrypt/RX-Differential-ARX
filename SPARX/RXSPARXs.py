"""
In this code we used to search the RX characteristics for CHAM
SPARX 64/128
The number of variable should be update
"""
import os
import time
import random
import itertools

FullRound = 4

BlockSize = 64
HalfBlockSize = 16

SearchRoundStart = 1
SearchRoundEnd = 11
InitialLowerBound = 0

RotateR = 7
RotateL = 2

# the round with the minuium
Fix = 0
Fixvalue = 2

GroupConstraintChoice = 1

# Parameters for choice 1
GroupNumForChoice1 = 1

DifferentialProbabilityBound = list([])
for i in range(20):
    DifferentialProbabilityBound += [0]

#Add one round variable 'x'
def CreateVariable(x,n,count):
    xtemp=[]
    for i in range(n):
        xtemp.append(count)
        count=count+1
    x.append(xtemp)
    return x,count

def GenEqual0(x,fout):
    clauseseq = "-" + str(x + 1)  + " 0" + "\n"
    fout.write(clauseseq)

def GenEqual1(x,fout):
    clauseseq = str(x + 1)  + " 0" + "\n"
    fout.write(clauseseq)
def GenSequentialEncoding(x, u, main_var_num, cardinalitycons, fout):
    n = main_var_num  # the number of wi
    k = cardinalitycons
    if (k > 0):
        clauseseq = "-" + str(x[0] + 1) + " " + str(u[0][0] + 1) + " 0" + "\n"
        fout.write(clauseseq)
        for j in range(1, k):
            clauseseq = "-" + str(u[0][j] + 1) + " 0" + "\n"
            fout.write(clauseseq)

        for i in range(1, n - 1):
            clauseseq = "-" + str(x[i] + 1) + " " + str(u[i][0] + 1) + " 0" + "\n"
            fout.write(clauseseq)
            clauseseq = "-" + str(u[i - 1][0] + 1) + " " + str(u[i][0] + 1) + " 0" + "\n"
            fout.write(clauseseq)
            clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][k - 1] + 1) + " 0" + "\n"
            fout.write(clauseseq)

        for j in range(1, k):
            for i in range(1, n - 1):
                clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][j - 1] + 1) + " " + str(
                    u[i][j] + 1) + " 0" + "\n"
                fout.write(clauseseq)
                clauseseq = "-" + str(u[i - 1][j] + 1) + " " + str(u[i][j] + 1) + " 0" + "\n"
                fout.write(clauseseq)

        clauseseq = "-" + str(x[n - 1] + 1) + " " + "-" + str(u[n - 2][k - 1] + 1) + " 0" + "\n"
        fout.write(clauseseq)

    if (k == 0):
        for i in range(n):
            clauseseq = "-" + str(x[i] + 1) + " 0" + "\n"
            fout.write(clauseseq)


def CountClausesInSequentialEncoding(main_var_num, cardinalitycons, clause_num):
    count = clause_num
    n = main_var_num
    k = cardinalitycons
    # print(count,n,k)
    if (k > 0):
        count += 1
        # print('1',count,k)
        for j in range(1, k):
            count += 1
        # print('2',count)
        for i in range(1, n - 1):
            count += 3
        # print('3',count)
        for j in range(1, k):
            for i in range(1, n - 1):
                count += 2
        # print('1',count)
        count += 1
    if (k == 0):
        for i in range(n):
            count += 1
    # print(count)
    return count

def GenMatsuiConstraint(x, u, n, k, left, right, m, g, fout):
    if (m > 0):
        if ((left == 0) and (right < n - 1)):
            if m < k:
                for i in range(right, right + 1):
                    clauseseq = "-" + str(u[i][m] + 1) + " 0" + "\n"
                    fout.write(clauseseq)
            if g > 0:
                for i in range(right, right + 1):
                    clauseseq = str(u[i][g - 1] + 1) + " 0" + "\n"
                    fout.write(clauseseq)

        if ((left > 0) and (right == n - 1)):
            for i in range(0, k - m):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(u[right - 1][i + m] + 1) + " 0" + "\n"
                fout.write(clauseseq)
            for i in range(0, k - m + 1):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(x[right] + 1) + " " + "-" + str(
                    u[right - 1][i + m - 1] + 1) + " 0" + "\n"
                fout.write(clauseseq)
        if ((left > 0) and (right < n - 1)):
            for i in range(0, k - m):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(u[right][i + m] + 1) + " 0" + "\n"
                fout.write(clauseseq)
    if (m == 0):
        for i in range(left, right + 1):
            clauseseq = "-" + str(x[i] + 1) + " 0" + "\n"
            fout.write(clauseseq)


def CountClausesForMatsuiStrategy(n, k, left, right, m, g, clausenum):
    count = clausenum
    if (m > 0):
        if ((left == 0) and (right < n - 1)):
            if m < k:
                for i in range(right, right + 1):
                    count += 1
            if g > 0:
                for i in range(right, right + 1):
                    count += 1

        if ((left > 0) and (right == n - 1)):
            for i in range(0, k - m):
                count += 1
            for i in range(0, k - m + 1):
                count += 1
        if ((left > 0) and (right < n - 1)):
            for i in range(0, k - m):
                count += 1
    if (m == 0):
        for i in range(left, right + 1):
            count += 1
    return count

def GenBitModularAddition(x, y, z, xx, yy, zz, fout):
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " " + str(xx + 1) + " " + str(
        yy + 1) + " " + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " " + str(xx + 1) + " " + str(
        yy + 1) + " " + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + str(xx + 1) + " " + str(
        yy + 1) + " " + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + str(xx + 1) + " " + str(
        yy + 1) + " " + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + "-" + str(xx + 1) + " " + "-" + str(
        yy + 1) + " " + "-" + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(
        xx + 1) + " " + "-" + str(yy + 1) + " " + "-" + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(
        xx + 1) + " " + "-" + str(yy + 1) + " " + "-" + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " " + "-" + str(
        xx + 1) + " " + "-" + str(yy + 1) + " " + "-" + str(zz + 1) + " 0" + "\n"
    fout.write(clauseseq)


def GenBitXOR(x, y, z, fout):
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " 0" + "\n"
    fout.write(clauseseq)

def GenBit2XOR(x, y, z,t, fout):
    for i in [[1, 1, 1, -1], [1, 1, -1, 1], [1, -1, 1, 1], [-1, 1, 1, 1], [-1, -1, -1, 1], [-1, 1, -1, -1],[-1, -1, 1, -1], [1, -1, -1, -1]]:
        clauseseq = str(i[0] * (x+1)) + " " + str(i[1] * (y+1)) + " " + str(i[2] * (z+1)) + " " + str(i[3] * (t+1)) + " 0\n"
        fout.write(clauseseq)

def GenEqual(x, y,fout):
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " 0" + "\n"
    fout.write(clauseseq)

def GenEqualNot(x, y,fout):
    clauseseq = "-" +str(x + 1) + " " + "-" + str(y + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " 0" + "\n"
    fout.write(clauseseq)

def GenBitwToAdditionPr(x,y,z,p0,p1,p2,e,fout):
    clauseseq = str(p0 + 1) + " " + "-" + str(p1 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(p2 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(p1 + 1) + " " + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(p0 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(p0 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " " + "-" + str(p0 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + "-" + str(p0 + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + str(z + 1) + " " + "-" + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + "-" + str(e + 1) + " 0" + "\n"
    fout.write(clauseseq)


def GenWToPr(x,y,z,p0,p1,p2,e,n,fout):
    for i in range(0, 1):
        GenBitwToAdditionPr(x[n-2], y[n-2], z[n-2], p0[i],p1[i],p2[i],e[i], fout)


def ToBinary(x,n):
    out=[0 for i in range(n)]
    t=x
    for i in range(n):
        out[n-1-i]=t%2
        t=int(t/2)
    return out

def GenXORC(x,y,vc,n,ro,fout):
    c=ToBinary(vc,n)
    rc=[0 for i in range(n)]
    ROL(c,rc,ro,n)#ro is the gama
    oc=[0 for i in range(n)]
    for i in range(n):
        oc[i]=c[i]^rc[i]
    for i in range(n):
        if oc[i]:
            GenEqualNot(x[i],y[i],fout)
        else:
            GenEqual(x[i],y[i],fout)



def GenAddition(x, y, z, n, fout):
    for i in range(0, n - 2):
        GenBitModularAddition(x[i], y[i], z[i], x[i + 1], y[i + 1], z[i + 1], fout)
    #GenBitXOR(x[n - 1], y[n - 1], z[n - 1], fout)


def GenBitAddPr(x, y, z, w, fout):
    clauseseq = "-" + str(x + 1) + " " + str(z + 1) + " " + str(w + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(y + 1) + " " + "-" + str(z + 1) + " " + str(w + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + "-" + str(y + 1) + " " + str(w + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + " " + "-" + str(w + 1) + " 0" + "\n"
    fout.write(clauseseq)
    clauseseq = "-" + str(x + 1) + " " + "-" + str(y + 1) + " " + "-" + str(z + 1) + " " + "-" + str(
        w + 1) + " 0" + "\n"
    fout.write(clauseseq)


def GenAddPr(x, y, z, w, n, fout):
    for i in range(0, n - 2):
        GenBitAddPr(x[i + 1], y[i + 1], z[i + 1], w[i], fout)



def ROR(x, rx, t, n):
    for i in range(0, n):
        rx[i] = x[(i - t) % n]


def CountCluseRoundFunction(Round, clausesum):
    count = clausesum
    # nonzero input
    count = count + 1
    for r in range(0, Round):
        if r%3==0:
            for i in range(HalfBlockSize):
                count=count+4
                count = count + 4
                count = count + 4
                count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11
            for i in range(HalfBlockSize):
                count=count+4
            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11
            for i in range(HalfBlockSize):
                count=count+4
            #key schedule
            for i in range(HalfBlockSize):
                count=count+2
                count = count + 2
                count = count + 2
                count = count + 2
                count = count + 2
                count = count + 2
            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11
            for i in range(HalfBlockSize):
                count=count+4

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11

            for i in range(HalfBlockSize):
                count=count+1

            for i in range(HalfBlockSize):
                count=count+2
                count = count + 2
                count = count + 2
                count = count + 2
                count = count + 2
                count = count + 2
                count = count + 2
                count = count + 2
                count = count + 2

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11
            for i in range(HalfBlockSize):
                count = count + 4

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11

            for i in range(HalfBlockSize):
                count = count + 1

            for i in range(HalfBlockSize):
                count = count + 2
                count = count + 2
                count = count + 2

        if r%3==1:
            for i in range(HalfBlockSize):
                count=count+4
                count = count + 4
                count = count + 4
                count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11
            for i in range(HalfBlockSize):
                count=count+4

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11
            for i in range(HalfBlockSize):
                count=count+4
        if r % 3 == 2:
            for i in range(HalfBlockSize):
                count = count + 4
                count = count + 4
                count = count + 4
                count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11
            for i in range(HalfBlockSize):
                count = count + 4

            for i in range(0, HalfBlockSize - 2):
                count = count + 8
            # count = count + 4
            for i in range(0, HalfBlockSize - 2):
                count = count + 5
            for i in range(0, 1):
                count = count + 11
            for i in range(HalfBlockSize):
                count = count + 4

            for i in range(HalfBlockSize):
                count=count+4
            for i in range(HalfBlockSize):
                count=count+8
                count=count+8

    return count


"""
Rotation to left Model
rx=x<<<t :rx and x are n bits(note the (n-1)-th bit is the least significant bit)
"""


def ROL(x, rx, t, n):
    for i in range(0, n):
        rx[i] = x[(i + t) % n]


def Per(x, px, P, n):
    for i in range(0, n):
        px[i] = x[P[i]]


def Decision(Round, n, ki, kd,MatsuiRoundIndex, MatsuiCount, flag):
    num = 0
    count = 0

    NI = 2*(n+1) * Round
    MI = ki
    ND= 2*Round
    #MD=kd
    kNI=8*(n+1)*(int((Round+2)/3))
    kMD=120
    km=4
    kw=8
    # sl=8
    sr = 3
    ###########################################################################################################
    # the defination of variable in the round function
    x0 = []  # the left input
    x1= []  # the right input
    x2=[]
    x3=[]
    k0=[]
    k1=[]
    k2=[]
    k3=[]
    #the outout of unternal state XOR subkey
    ox0=[]
    ox1=[]
    ox2=[]
    ox3=[]
    w = []  # the auxiliary variable in computing the differential probability for modular addition
    p0=[]
    p1=[]
    p2=[]
    e=[]
    uk=[]



    ww = []  # the auxiliary variable in computing the differential probability for modular addition
    pp0=[]
    pp1=[]
    pp2=[]
    ee=[]

    k0w = []  # the auxiliary variable in computing the differential probability for modular addition
    k0p0=[]
    k0p1=[]
    k0p2=[]
    k0e=[]

    k1w = []  # the auxiliary variable in computing the differential probability for modular addition
    k1p0=[]
    k1p1=[]
    k1p2=[]
    k1e=[]

    k2w = []  # the auxiliary variable in computing the differential probability for modular addition
    k2p0=[]
    k2p1=[]
    k2p2=[]
    k2e=[]

    k3w = []  # the auxiliary variable in computing the differential probability for modular addition
    k3p0=[]
    k3p1=[]
    k3p2=[]
    k3e=[]

    ui = []  # the auxiliary variable in integer

    #the auxiliary variable in linear layer
    a1=[]#first auxi varible in r%3==2 round
    a2=[]#secon auxi varible in r%3==2 round
    a3=[]

    #the master key 128-bit
    mk0=[]
    mk1=[]
    mk2=[]
    mk3=[]
    mk4=[]
    mk5=[]
    mk6=[]
    mk7=[]

    c=[]

    [mk0,count]=CreateVariable(mk0,n,count)
    [mk1,count]=CreateVariable(mk1, n, count)
    [mk2,count]=CreateVariable(mk2, n, count)
    [mk3,count]=CreateVariable(mk3, n, count)
    [mk4,count]=CreateVariable(mk4, n, count)
    [mk5,count]=CreateVariable(mk5, n, count)
    [mk6,count]=CreateVariable(mk6, n, count)
    [mk7,count]=CreateVariable(mk7, n, count)
    # declar the variable in 1~Round round function
    for r in range(0, Round):
        [x0,count]=CreateVariable(x0,n,count)
        [x1, count] = CreateVariable(x1, n, count)
        [x2, count] = CreateVariable(x2, n, count)
        [x3, count] = CreateVariable(x3, n, count)

        [ox0,count]=CreateVariable(ox0,n,count)
        [ox1, count] = CreateVariable(ox1, n, count)
        [ox2, count] = CreateVariable(ox2, n, count)
        [ox3, count] = CreateVariable(ox3, n, count)

        [k0,count]=CreateVariable(k0,n,count)
        [k1, count] = CreateVariable(k1, n, count)
        [k2, count] = CreateVariable(k2, n, count)
        [k3, count] = CreateVariable(k3, n, count)

        [w,count]=CreateVariable(w,n-1,count)
        [p0,count]=CreateVariable(p0,1,count)
        [p1,count]=CreateVariable(p1,1,count)
        [p2,count]=CreateVariable(p2,1,count)
        [e,count]=CreateVariable(e,1,count)

        [ww,count]=CreateVariable(ww,n-1,count)
        [pp0,count]=CreateVariable(pp0,1,count)
        [pp1,count]=CreateVariable(pp1,1,count)
        [pp2,count]=CreateVariable(pp2,1,count)
        [ee,count]=CreateVariable(ee,1,count)


        if r%3==0:
            [mk0, count] = CreateVariable(mk0, n, count)
            [mk1, count] = CreateVariable(mk1, n, count)
            [mk2, count] = CreateVariable(mk2, n, count)
            [mk3, count] = CreateVariable(mk3, n, count)
            [mk4, count] = CreateVariable(mk4, n, count)
            [mk5, count] = CreateVariable(mk5, n, count)
            [mk6, count] = CreateVariable(mk6, n, count)
            [mk7, count] = CreateVariable(mk7, n, count)
            [c, count] = CreateVariable(c, n, count)
            [k0w, count] = CreateVariable(k0w, n - 1, count)
            [k0p0, count] = CreateVariable(k0p0, 1, count)
            [k0p1, count] = CreateVariable(k0p1, 1, count)
            [k0p2, count] = CreateVariable(k0p2, 1, count)
            [k0e, count] = CreateVariable(k0e, 1, count)
            [k1w, count] = CreateVariable(k1w, n - 1, count)
            [k1p0, count] = CreateVariable(k1p0, 1, count)
            [k1p1, count] = CreateVariable(k1p1, 1, count)
            [k1p2, count] = CreateVariable(k1p2, 1, count)
            [k1e, count] = CreateVariable(k1e, 1, count)
            [k2w, count] = CreateVariable(k2w, n - 1, count)
            [k2p0, count] = CreateVariable(k2p0, 1, count)
            [k2p1, count] = CreateVariable(k2p1, 1, count)
            [k2p2, count] = CreateVariable(k2p2, 1, count)
            [k2e, count] = CreateVariable(k2e, 1, count)
            [k3w, count] = CreateVariable(k3w, n - 1, count)
            [k3p0, count] = CreateVariable(k3p0, 1, count)
            [k3p1, count] = CreateVariable(k3p1, 1, count)
            [k3p2, count] = CreateVariable(k3p2, 1, count)
            [k3e, count] = CreateVariable(k3e, 1, count)
            [mk0, count] = CreateVariable(mk0, n, count)
            [mk1, count] = CreateVariable(mk1, n, count)
            [mk2, count] = CreateVariable(mk2, n, count)
            [mk3, count] = CreateVariable(mk3, n, count)
            [mk4, count] = CreateVariable(mk4, n, count)
            [mk5, count] = CreateVariable(mk5, n, count)
            [mk6, count] = CreateVariable(mk6, n, count)
            [mk7, count] = CreateVariable(mk7, n, count)
            [c,count]=CreateVariable(c,n,count)
            [k0w, count] = CreateVariable(k0w, n - 1, count)
            [k0p0, count] = CreateVariable(k0p0, 1, count)
            [k0p1, count] = CreateVariable(k0p1, 1, count)
            [k0p2, count] = CreateVariable(k0p2, 1, count)
            [k0e, count] = CreateVariable(k0e, 1, count)
            [k1w, count] = CreateVariable(k1w, n - 1, count)
            [k1p0, count] = CreateVariable(k1p0, 1, count)
            [k1p1, count] = CreateVariable(k1p1, 1, count)
            [k1p2, count] = CreateVariable(k1p2, 1, count)
            [k1e, count] = CreateVariable(k1e, 1, count)
            [k2w, count] = CreateVariable(k2w, n - 1, count)
            [k2p0, count] = CreateVariable(k2p0, 1, count)
            [k2p1, count] = CreateVariable(k2p1, 1, count)
            [k2p2, count] = CreateVariable(k2p2, 1, count)
            [k2e, count] = CreateVariable(k2e, 1, count)
            [k3w, count] = CreateVariable(k3w, n - 1, count)
            [k3p0, count] = CreateVariable(k3p0, 1, count)
            [k3p1, count] = CreateVariable(k3p1, 1, count)
            [k3p2, count] = CreateVariable(k3p2, 1, count)
            [k3e, count] = CreateVariable(k3e, 1, count)




        if r%3==2:
            [a1,count]=CreateVariable(a1,n,count)
            [a2, count] = CreateVariable(a2, n, count)
            [a3, count] = CreateVariable(a3, n, count)


    # declar the output variable in the Round-round function
    [x0, count] = CreateVariable(x0, n, count)
    [x1, count] = CreateVariable(x1, n, count)
    [x2, count] = CreateVariable(x2, n, count)
    [x3, count] = CreateVariable(x3, n, count)


    if Round%3==1:
        [k0,count]=CreateVariable(k0,n,count)
        [k1, count] = CreateVariable(k1, n, count)
        [k2, count] = CreateVariable(k2, n, count)
        [k3, count] = CreateVariable(k3, n, count)
        [k0,count]=CreateVariable(k0,n,count)
        [k1, count] = CreateVariable(k1, n, count)
        [k2, count] = CreateVariable(k2, n, count)
        [k3, count] = CreateVariable(k3, n, count)
    if Round%3==2:
        [k0,count]=CreateVariable(k0,n,count)
        [k1, count] = CreateVariable(k1, n, count)
        [k2, count] = CreateVariable(k2, n, count)
        [k3, count] = CreateVariable(k3, n, count)



    # declar the auxility variable in the boolean cardinality constraint
    for i in range(0, NI - 1):
        uitemp = []
        for j in range(0, MI):
            uitemp.append(count)
            count = count + 1
        ui.append(uitemp)

    for i in range(0, kNI - 1):
        uktemp = []
        for j in range(0, kMD):
            uktemp.append(count)
            count = count + 1
        uk.append(uktemp)

    countclause = 0
    countclause=countclause+4*4*HalfBlockSize
    # print(Round)
    # print(HalfBlockSize)
    # print(k)
    # print(N)
    countclause = CountCluseRoundFunction(Round, countclause)
    countclause = CountClausesInSequentialEncoding(NI, ki, countclause)
    countclause = CountClausesInSequentialEncoding(kNI, kMD, countclause)

    for matsui_count in range(0, MatsuiCount):
        StartingRound = MatsuiRoundIndex[matsui_count][0]
        EndingRound = MatsuiRoundIndex[matsui_count][1]
        LeftNode = 2*(HalfBlockSize +1) * StartingRound
        RightNode = 2*(HalfBlockSize +1) * EndingRound - 1
        g = DifferentialProbabilityBound[EndingRound]
        PartialCardinalityCons = ki - DifferentialProbabilityBound[StartingRound] - DifferentialProbabilityBound[
            Round - EndingRound]
        countclause = CountClausesForMatsuiStrategy(NI, ki, LeftNode, RightNode, PartialCardinalityCons, g, countclause)
    # Open file


    time_start = time.time()
    file = open("Problem-Round" + str(Round) + "-Probability" + str(ki)+"_"+str(kd) + ".cnf", "w")
    file.write("p cnf " + str(count) + " " + str(countclause) + "\n")
    # print('===============================')

    # Add constraints for the round function
    kr=0
    for r in range(5,9):
        for i in range(HalfBlockSize):
            GenEqual0(ox0[r][i], file)
            GenEqual0(ox1[r][i], file)
            GenEqual0(ox2[r][i], file)
            GenEqual0(ox3[r][i], file)
    for r in range(0, Round):
        rx0 = [0 for i in range(0, HalfBlockSize)]
        rx1 = [0 for i in range(0, HalfBlockSize)]
        rx2 = [0 for i in range(0, HalfBlockSize)]
        rx3 = [0 for i in range(0, HalfBlockSize)]
        ROR(ox0[r], rx0, RotateR, HalfBlockSize)
        ROR(ox2[r], rx2, RotateR, HalfBlockSize)
        if r % 3 == 0:
            for i in range(0, HalfBlockSize):
                GenBitXOR(x0[r][i], k0[r][i], ox0[r][i], file)
                GenBitXOR(x1[r][i], k1[r][i], ox1[r][i], file)
                GenBitXOR(x2[r][i], k2[r][i], ox2[r][i], file)
                GenBitXOR(x3[r][i], k3[r][i], ox3[r][i], file)
            GenAddition(rx0, ox1[r], x0[r + 1], HalfBlockSize, file)
            GenAddPr(rx0, ox1[r], x0[r + 1], w[r], HalfBlockSize, file)
            GenWToPr(rx0, ox1[r], x0[r + 1], p0[r], p1[r], p2[r], e[r], HalfBlockSize, file)
            ROL(ox1[r], rx1, RotateL, HalfBlockSize)
            for i in range(0, HalfBlockSize):
                GenBitXOR(x0[r + 1][i], rx1[i], x1[r + 1][i], file)

            GenAddition(rx2, ox3[r], x2[r + 1], HalfBlockSize, file)
            GenAddPr(rx2, ox3[r], x2[r + 1], ww[r], HalfBlockSize, file)
            GenWToPr(rx2, ox3[r], x2[r + 1], pp0[r], pp1[r], pp2[r], ee[r], HalfBlockSize, file)
            ROL(ox3[r], rx3, RotateL, HalfBlockSize)
            for i in range(0, HalfBlockSize):
                GenBitXOR(x2[r + 1][i], rx3[i], x3[r + 1][i], file)

            # key schedule
            if 1:
                for i in range(0,HalfBlockSize):
                    GenEqual(mk0[kr][i],k0[r][i],file)
                    GenEqual(mk1[kr][i], k1[r][i], file)
                    GenEqual(mk2[kr][i], k0[r+1][i], file)
                    GenEqual(mk3[kr][i], k1[r+1][i], file)
                    GenEqual(mk4[kr][i], k0[r+2][i], file)
                    GenEqual(mk5[kr][i], k1[r+2][i], file)
                rmk0 = [0 for i in range(0, HalfBlockSize)]
                rmk1 = [0 for i in range(0, HalfBlockSize)]
                ROR(mk0[kr], rmk0, RotateR, HalfBlockSize)
                GenAddition(rmk0, mk1[kr],mk2[kr + 1], HalfBlockSize, file)
                GenAddPr(rmk0, mk1[kr], mk2[kr + 1], k0w[kr], HalfBlockSize, file)
                GenWToPr(rmk0, mk1[kr], mk2[kr + 1], k0p0[kr], k0p1[kr], k0p2[kr], k0e[kr], HalfBlockSize, file)
                ROL(mk1[kr], rmk1, RotateL, HalfBlockSize)
                for i in range(0, HalfBlockSize):
                    GenBitXOR(mk2[kr + 1][i], rmk1[i], mk3[kr + 1][i], file)

                GenAddition(mk2[kr + 1], mk2[kr],mk4[kr + 1], HalfBlockSize, file)
                GenAddPr(mk2[kr + 1], mk2[kr], mk4[kr + 1], k1w[kr], HalfBlockSize, file)
                GenWToPr(mk2[kr + 1], mk2[kr], mk4[kr + 1], k1p0[kr], k1p1[kr], k1p2[kr], k1e[kr], HalfBlockSize, file)

                GenAddition(mk3[kr + 1], mk3[kr],mk5[kr + 1], HalfBlockSize, file)
                GenAddPr(mk3[kr + 1], mk3[kr], mk5[kr + 1], k2w[kr], HalfBlockSize, file)
                GenWToPr(mk3[kr + 1], mk3[kr], mk5[kr + 1], k2p0[kr], k2p1[kr], k2p2[kr], k2e[kr], HalfBlockSize, file)

                GenAddition(mk7[kr], c[kr],mk1[kr + 1], HalfBlockSize, file)
                GenAddPr(mk7[kr], c[kr], mk1[kr + 1], k3w[kr], HalfBlockSize, file)
                GenWToPr(mk7[kr], c[kr], mk1[kr + 1], k3p0[kr], k3p1[kr], k3p2[kr], k3e[kr], HalfBlockSize, file)

                vc=kr+1
                vcb=ToBinary(vc,HalfBlockSize)
                rvcb=[0 for i in range(HalfBlockSize)]
                ROL(vcb,rvcb,1,HalfBlockSize)
                oc=[0 for i in range(HalfBlockSize)]
                for i in range(HalfBlockSize):
                    oc[i]=vcb[i]^rvcb[i]
                for i in range(HalfBlockSize):
                    if oc[i]:
                        GenEqual1(c[kr][i],file)
                    else:
                        GenEqual0(c[kr][i],file)



                for i in range(0,HalfBlockSize):
                    GenEqual(mk4[kr][i],mk6[kr+1][i],file)
                    GenEqual(mk5[kr][i],mk7[kr+1][i],file)
                    GenEqual(mk6[kr][i],mk0[kr+1][i],file)

                for i in range(0,HalfBlockSize):
                    GenEqual(mk0[kr+1][i],k2[r][i],file)
                    GenEqual(mk1[kr+1][i], k3[r][i], file)
                    GenEqual(mk2[kr+1][i], k2[r+1][i], file)
                    GenEqual(mk3[kr+1][i], k3[r+1][i], file)
                    GenEqual(mk4[kr+1][i], k2[r+2][i], file)
                    GenEqual(mk5[kr+1][i], k3[r+2][i], file)

                kr=kr+1
                rmk0 = [0 for i in range(0, HalfBlockSize)]
                rmk1 = [0 for i in range(0, HalfBlockSize)]
                ROR(mk0[kr], rmk0, RotateR, HalfBlockSize)
                GenAddition(rmk0, mk1[kr],mk2[kr + 1], HalfBlockSize, file)
                GenAddPr(rmk0, mk1[kr], mk2[kr + 1], k0w[kr], HalfBlockSize, file)
                GenWToPr(rmk0, mk1[kr], mk2[kr + 1], k0p0[kr], k0p1[kr], k0p2[kr], k0e[kr], HalfBlockSize, file)
                ROL(mk1[kr], rmk1, RotateL, HalfBlockSize)
                for i in range(0, HalfBlockSize):
                    GenBitXOR(mk2[kr + 1][i], rmk1[i], mk3[kr + 1][i], file)

                GenAddition(mk2[kr + 1], mk2[kr],mk4[kr + 1], HalfBlockSize, file)
                GenAddPr(mk2[kr + 1], mk2[kr], mk4[kr + 1], k1w[kr], HalfBlockSize, file)
                GenWToPr(mk2[kr + 1], mk2[kr], mk4[kr + 1], k1p0[kr], k1p1[kr], k1p2[kr], k1e[kr], HalfBlockSize, file)

                GenAddition(mk3[kr + 1], mk3[kr],mk5[kr + 1], HalfBlockSize, file)
                GenAddPr(mk3[kr + 1], mk3[kr], mk5[kr + 1], k2w[kr], HalfBlockSize, file)
                GenWToPr(mk3[kr + 1], mk3[kr], mk5[kr + 1], k2p0[kr], k2p1[kr], k2p2[kr], k2e[kr], HalfBlockSize, file)

                GenAddition(mk7[kr], c[kr],mk1[kr + 1], HalfBlockSize, file)
                GenAddPr(mk7[kr], c[kr], mk1[kr + 1], k3w[kr], HalfBlockSize, file)
                GenWToPr(mk7[kr], c[kr], mk1[kr + 1], k3p0[kr], k3p1[kr], k3p2[kr], k3e[kr], HalfBlockSize, file)

                vc=kr+1
                vcb=ToBinary(vc,HalfBlockSize)
                rvcb=[0 for i in range(HalfBlockSize)]
                ROL(vcb,rvcb,1,HalfBlockSize)
                oc=[0 for i in range(HalfBlockSize)]
                for i in range(HalfBlockSize):
                    oc[i]=vcb[i]^rvcb[i]
                for i in range(HalfBlockSize):
                    if oc[i]:
                        GenEqual1(c[kr][i],file)
                    else:
                        GenEqual0(c[kr][i],file)



                for i in range(0,HalfBlockSize):
                    GenEqual(mk4[kr][i],mk6[kr+1][i],file)
                    GenEqual(mk5[kr][i],mk7[kr+1][i],file)
                    GenEqual(mk6[kr][i],mk0[kr+1][i],file)

                kr=kr+1
                # one time key schedule





        if r%3==1:
            for i in range(0,HalfBlockSize):
                GenBitXOR(x0[r][i],k0[r][i],ox0[r][i],file)
                GenBitXOR(x1[r][i], k1[r][i], ox1[r][i], file)
                GenBitXOR(x2[r][i], k2[r][i], ox2[r][i], file)
                GenBitXOR(x3[r][i], k3[r][i], ox3[r][i], file)
            GenAddition(rx0, ox1[r], x0[r+1], HalfBlockSize, file)
            GenAddPr(rx0, ox1[r], x0[r+1], w[r], HalfBlockSize, file)
            GenWToPr(rx0, ox1[r], x0[r+1], p0[r], p1[r], p2[r], e[r], HalfBlockSize, file)
            ROL(ox1[r],rx1,RotateL,HalfBlockSize)
            for i in range(0,HalfBlockSize):
                GenBitXOR(x0[r+1][i],rx1[i],x1[r+1][i],file)

            GenAddition(rx2, ox3[r], x2[r+1], HalfBlockSize, file)
            GenAddPr(rx2, ox3[r], x2[r+1], ww[r], HalfBlockSize, file)
            GenWToPr(rx2, ox3[r], x2[r+1], pp0[r], pp1[r], pp2[r], ee[r], HalfBlockSize, file)
            ROL(ox3[r],rx3,RotateL,HalfBlockSize)
            for i in range(0,HalfBlockSize):
                GenBitXOR(x2[r+1][i],rx3[i],x3[r+1][i],file)

            #key schedule




        if r%3==2:
            rr=r/3
            for i in range(0, HalfBlockSize):
                GenBitXOR(x0[r][i], k0[r][i], ox0[r][i], file)
                GenBitXOR(x1[r][i], k1[r][i], ox1[r][i], file)
                GenBitXOR(x2[r][i], k2[r][i], ox2[r][i], file)
                GenBitXOR(x3[r][i], k3[r][i], ox3[r][i], file)
            GenAddition(rx0, ox1[r], x2[r + 1], HalfBlockSize, file)
            GenAddPr(rx0, ox1[r], x2[r + 1], w[r], HalfBlockSize, file)
            GenWToPr(rx0, ox1[r], x2[r + 1], p0[r], p1[r], p2[r], e[r], HalfBlockSize, file)
            ROL(ox1[r], rx1, RotateL, HalfBlockSize)
            for i in range(0, HalfBlockSize):
                GenBitXOR(x2[r + 1][i], rx1[i], x3[r + 1][i], file)

            GenAddition(rx2, ox3[r], a1[rr], HalfBlockSize, file)
            GenAddPr(rx2, ox3[r], a1[rr], ww[r], HalfBlockSize, file)
            GenWToPr(rx2, ox3[r], a1[rr], pp0[r], pp1[r], pp2[r], ee[r], HalfBlockSize, file)
            ROL(ox3[r], rx3, RotateL, HalfBlockSize)
            for i in range(0, HalfBlockSize):
                GenBitXOR(a1[rr][i], rx3[i], a2[rr][i], file)

            for i in range(0,HalfBlockSize):
                GenBitXOR(x2[r+1][i],x3[r+1][i],a3[rr][i],file)

            ra3 = [0 for i in range(0, HalfBlockSize)]
            ROL(a3[rr], ra3, 8, HalfBlockSize)

            for i in range(0,HalfBlockSize):
                GenBit2XOR(x2[r+1][i],ra3[i],a1[rr][i],x0[r+1][i],file)
                GenBit2XOR(x3[r+1][i],ra3[i],a2[rr][i],x1[r+1][i],file)







    C = []
    for r in range(0, Round):
        for i in range(0, HalfBlockSize - 2):
            C.append(w[r][i])
        C.append(p0[r][0])
        C.append(p1[r][0])
        C.append(p2[r][0])

        for i in range(0, HalfBlockSize - 2):
            C.append(ww[r][i])
        C.append(pp0[r][0])
        C.append(pp1[r][0])
        C.append(pp2[r][0])
    GenSequentialEncoding(C, ui, NI, ki, file)
    Ck=[]
    tr=0
    for r in range(0, Round):
        if r%3==0:
            #tr=2*(int(r/3))
            for i in range(0, HalfBlockSize - 2):
                Ck.append(k0w[tr][i])
            Ck.append(k0p0[tr][0])
            Ck.append(k0p1[tr][0])
            Ck.append(k0p2[tr][0])
            for i in range(0, HalfBlockSize - 2):
                Ck.append(k1w[tr][i])
            Ck.append(k1p0[tr][0])
            Ck.append(k1p1[tr][0])
            Ck.append(k1p2[tr][0])
            for i in range(0, HalfBlockSize - 2):
                Ck.append(k2w[tr][i])
            Ck.append(k2p0[tr][0])
            Ck.append(k2p1[tr][0])
            Ck.append(k2p2[tr][0])
            for i in range(0, HalfBlockSize - 2):
                Ck.append(k3w[tr][i])
            Ck.append(k3p0[tr][0])
            Ck.append(k3p1[tr][0])
            Ck.append(k3p2[tr][0])

            tr=tr+1
            for i in range(0, HalfBlockSize - 2):
                Ck.append(k0w[tr][i])
            Ck.append(k0p0[tr][0])
            Ck.append(k0p1[tr][0])
            Ck.append(k0p2[tr][0])
            for i in range(0, HalfBlockSize - 2):
                Ck.append(k1w[tr][i])
            Ck.append(k1p0[tr][0])
            Ck.append(k1p1[tr][0])
            Ck.append(k1p2[tr][0])
            for i in range(0, HalfBlockSize - 2):
                Ck.append(k2w[tr][i])
            Ck.append(k2p0[tr][0])
            Ck.append(k2p1[tr][0])
            Ck.append(k2p2[tr][0])
            for i in range(0, HalfBlockSize - 2):
                Ck.append(k3w[tr][i])
            Ck.append(k3p0[tr][0])
            Ck.append(k3p1[tr][0])
            Ck.append(k3p2[tr][0])
            tr=tr+1



    GenSequentialEncoding(Ck, uk, kNI, kMD, file)
    """
    Cd = []
    for r in range(0, Round):
        for i in range(0, 1):
            Cd.append(e[r][i])
            #Cd.append(ke[r][i])
    GenSequentialEncoding(Cd, ud, ND, kd, file)
    """

    for matsui_count in range(0, MatsuiCount):
        StartingRound = MatsuiRoundIndex[matsui_count][0]
        EndingRound = MatsuiRoundIndex[matsui_count][1]
        LeftNode = 2*(HalfBlockSize +1) * StartingRound
        RightNode = 2*(HalfBlockSize +1) * EndingRound - 1
        g = DifferentialProbabilityBound[EndingRound]
        PartialCardinalityCons = ki - DifferentialProbabilityBound[StartingRound] - DifferentialProbabilityBound[
            Round - EndingRound]
        GenMatsuiConstraint(C, ui, NI, ki, LeftNode, RightNode, PartialCardinalityCons, g, file)



    # Add constraints to claim nonzero input difference
    clauseseq = ""
    for i in range(0, HalfBlockSize):
        clauseseq += str(x0[0][i] + 1) + " " + str(x1[0][i] + 1)+ " " + str(x2[0][i] + 1)+ " " + str(x3[0][i] + 1) + " "+ str(mk0[0][i] + 1) + " "+ str(mk1[0][i] + 1) + " "+ str(mk2[0][i] + 1) + " "+ str(mk3[0][i] + 1) + " "+ str(mk4[0][i] + 1) + " "+ str(mk5[0][i] + 1) + " "+ str(mk6[0][i] + 1) + " "+ str(mk7[0][i] + 1) + " "
    clauseseq += "0" + "\n"
    file.write(clauseseq)
    file.close()
    # Call solver cadical
    order = "/home/zhangyh/Software/cadical/build/cadical " + "Problem-Round" + str(Round) + "-Probability" + str(ki)+"_"+str(kd) + ".cnf > Round" + str(Round) + "-Probability" + str(ki)+"_"+str(kd)   + "-solution.out"
    os.system(order)
    # Extracting results
    order = "sed -n '/s SATISFIABLE/p' Round" + str(Round) + "-Probability" +  str(ki)+"_"+str(kd)  + "-solution.out > SatSolution.out"
    os.system(order)
    order = "sed -n '/s UNSATISFIABLE/p' Round" + str(Round) + "-Probability" +  str(ki)+"_"+str(kd)  + "-solution.out > UnsatSolution.out"
    os.system(order)
    satsol = open("SatSolution.out")
    unsatsol = open("UnsatSolution.out")
    satresult = satsol.readlines()
    unsatresult = unsatsol.readlines()
    satsol.close()
    unsatsol.close()
    if ((len(satresult) == 0) and (len(unsatresult) > 0)):
        flag = False
    if ((len(satresult) > 0) and (len(unsatresult) == 0)):
        flag = True
    order = "rm SatSolution.out"
    os.system(order)
    order = "rm UnsatSolution.out"
    os.system(order)
    # Removing cnf file
    order = "rm Problem-Round" + str(Round) + "-Probability" +  str(ki)+"_"+str(kd)  + ".cnf"
    os.system(order)
    time_end = time.time()
    # Printing solutions

    if (flag == True):
        print(
            "Round:" + str(Round) + "; Active: " +  str(ki)+"_"+str(kd)  + "; Sat; TotalCost: " + str(time_end - time_start))
    else:
        print("Round:" + str(Round) + "; Active: " +  str(ki)+"_"+str(kd)  + "; Unsat; TotalCost: " + str(
            time_end - time_start))

    return flag


# main function

TotalTimeStart = time.time()

for u in range(1):
    CountProbabilityi = InitialLowerBound
    CountProbabilityd = InitialLowerBound
    Countindex=0
    for totalround in range(1,20):
        flag = False
        time_start = time.time()
        MatsuiRoundIndex = []
        MatsuiCount = 0

        if (GroupConstraintChoice == 1):
            for group in range(0, GroupNumForChoice1):
                for round in range(1, totalround - group + 1):
                    MatsuiRoundIndex.append([])
                    MatsuiRoundIndex[MatsuiCount].append(group)
                    MatsuiRoundIndex[MatsuiCount].append(group + round)
                    MatsuiCount += 1

            for group in range(totalround, totalround):
                for round in range(1, totalround):
                    MatsuiRoundIndex.append([])
                    MatsuiRoundIndex[MatsuiCount].append(group-round)
                    MatsuiRoundIndex[MatsuiCount].append(group)
                    MatsuiCount += 1

        # Printing Matsui conditions
        file = open("MatsuiCondition.out", "a")
        resultseq = "Round: " + str(totalround) + "; Partial Constraint Num: " + str(MatsuiCount) + "\n"
        file.write(resultseq)
        file.write(str(MatsuiRoundIndex) + "\n")
        file.close()
        while (flag == False):
            flag = Decision(totalround, HalfBlockSize, CountProbabilityi,CountProbabilityi, MatsuiRoundIndex, MatsuiCount, flag)
            CountProbabilityi=CountProbabilityi+1

        #DifferentialProbabilityBound[totalround] = CountProbability - 1
        CountProbabilityi = CountProbabilityi-1
        time_end = time.time()
        file = open("RunTimeSummarise.out", "a")
        resultseq = "Round: " + str(totalround) + "; Differential Probability: " + str(
            DifferentialProbabilityBound[totalround]) + "; Runtime: " + str(time_end - time_start) + "\n"
        file.write(resultseq)
        file.close()

    # print(str(DifferentialProbabilityBound))
    # ooot.write(str(DifferentialProbabilityBound))
    # ooot.write("\n")
    TotalTimeEnd = time.time()
    # print("Total Runtime: " + str(TotalTimeEnd - TotalTimeStart))
    file = open("RunTimeSummarise.out", "a")
    resultseq = "Total Runtime: " + str(TotalTimeEnd - TotalTimeStart)
    file.write(resultseq)

