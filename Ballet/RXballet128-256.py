"""
In this code we used to search the RX characteristics for Ballet
Ballet 128/128
"""
import os
import time
import random
import itertools

FullRound = 4

BlockSize = 128
HalfBlockSize = 32

SearchRoundStart = 1
SearchRoundEnd = 11
InitialLowerBound = 0

RotateR = 8
RotateL = 3 

# the round with the minuium
Fix = 0
Fixvalue = 2

GroupConstraintChoice = 1

# Parameters for choice 1
GroupNumForChoice1 = 1

DifferentialProbabilityBound = list([])
for i in range(20):
    DifferentialProbabilityBound += [0]


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
        for i in range(0, HalfBlockSize):
            count = count + 4
        for i in range(0, HalfBlockSize - 2):
            count = count + 8
        #count = count + 4
        for i in range(0, HalfBlockSize - 2):
            count = count + 5
        for i in range(0, 1):
            count = count + 11
        for i in range(0, HalfBlockSize - 2):
            count = count + 8
        #count = count + 4
        for i in range(0, HalfBlockSize - 2):
            count = count + 5
        for i in range(0, 1):
            count = count + 11
        for i in range(0, HalfBlockSize):
            count = count + 60
        #for i in range(0,2*HalfBlockSize):
            #count=count+2

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
    MD=kd
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
    w = []  # the auxiliary variable in computing the differential probability for modular addition
    p0=[]
    p1=[]
    p2=[]
    e=[]
    ww = []  # the auxiliary variable in computing the differential probability for modular addition
    pp0=[]
    pp1=[]
    pp2=[]
    ee=[]
    ui = []  # the auxiliary variable in integer
    ud=[]#the decimal auxiliary variable
    k0=[]#the master key
    k1=[]
    k2=[]
    k3=[]
    k4=[]
    k5=[]
    k6=[]
    k7=[]
    ko1=[]
    ko0=[]#the input of the operation of XORC
    #ko1=[]
    mx=[]#the output of the XORC and the input of the modular addition

    # declar the variable in 1~Round round function

    for r in range(0, Round):
        x0temp = []
        for i in range(0, n):
            x0temp.append(count)
            count = count + 1
        x0.append(x0temp)
        x1temp = []
        for i in range(0, n):
            x1temp.append(count)
            count = count + 1
        x1.append(x1temp)
        x2temp = []
        for i in range(0, n):
            x2temp.append(count)
            count = count + 1
        x2.append(x2temp)
        x3temp = []
        for i in range(0, n):
            x3temp.append(count)
            count = count + 1
        x3.append(x3temp)
        mxtemp = []
        for i in range(0, n):
            mxtemp.append(count)
            count = count + 1
        mx.append(mxtemp)

        ko0temp = []
        for i in range(0, 2*n):
            ko0temp.append(count)
            count = count + 1
        ko0.append(ko0temp)
        ko1temp = []
        for i in range(0, 2*n):
            ko1temp.append(count)
            count = count + 1
        ko1.append(ko1temp)



        wtemp = []
        for i in range(0, n - 1):
            wtemp.append(count)
            count = count + 1
        w.append(wtemp)
        p0temp = []
        for i in range(0, 1):
            p0temp.append(count)
            count = count + 1
        p0.append(p0temp)
        p1temp = []
        for i in range(0, 1):
            p1temp.append(count)
            count = count + 1
        p1.append(p1temp)
        p2temp = []
        for i in range(0, 1):
            p2temp.append(count)
            count = count + 1
        p2.append(p2temp)
        etemp = []
        for i in range(0, 1):
            etemp.append(count)
            count = count + 1
        e.append(etemp)

        wwtemp = []
        for i in range(0, n - 1):
            wwtemp.append(count)
            count = count + 1
        ww.append(wwtemp)
        pp0temp = []
        for i in range(0, 1):
            pp0temp.append(count)
            count = count + 1
        pp0.append(pp0temp)
        pp1temp = []
        for i in range(0, 1):
            pp1temp.append(count)
            count = count + 1
        pp1.append(pp1temp)
        pp2temp = []
        for i in range(0, 1):
            pp2temp.append(count)
            count = count + 1
        pp2.append(pp2temp)
        eetemp = []
        for i in range(0, 1):
            eetemp.append(count)
            count = count + 1
        ee.append(eetemp)


        k0temp = []
        for i in range(0, n):
            k0temp.append(count)
            count = count + 1
        k0.append(k0temp)
        k1temp = []
        for i in range(0, n):
            k1temp.append(count)
            count = count + 1
        k1.append(k1temp)
        k2temp = []
        for i in range(0, n):
            k2temp.append(count)
            count = count + 1
        k2.append(k2temp)
        k3temp = []
        for i in range(0, n):
            k3temp.append(count)
            count = count + 1
        k3.append(k3temp)
        k4temp = []
        for i in range(0, n):
            k4temp.append(count)
            count = count + 1
        k4.append(k4temp)
        k5temp = []
        for i in range(0, n):
            k5temp.append(count)
            count = count + 1
        k5.append(k5temp)
        k6temp = []
        for i in range(0, n):
            k6temp.append(count)
            count = count + 1
        k6.append(k6temp)
        k7temp = []
        for i in range(0, n):
            k7temp.append(count)
            count = count + 1
        k7.append(k7temp)


    # declar the output variable in the Round-round function
    x0temp = []
    for i in range(0, n):
        x0temp.append(count)
        count = count + 1
    x0.append(x0temp)
    x1temp = []
    for i in range(0, n):
        x1temp.append(count)
        count = count + 1
    x1.append(x1temp)
    x2temp = []
    for i in range(0, n):
        x2temp.append(count)
        count = count + 1
    x2.append(x2temp)
    x3temp = []
    for i in range(0, n):
        x3temp.append(count)
        count = count + 1
    x3.append(x3temp)

    k0temp = []
    for i in range(0, n):
        k0temp.append(count)
        count = count + 1
    k0.append(k0temp)
    k1temp = []
    for i in range(0, n):
        k1temp.append(count)
        count = count + 1
    k1.append(k1temp)
    k2temp = []
    for i in range(0, n):
        k2temp.append(count)
        count = count + 1
    k2.append(k2temp)
    k3temp = []
    for i in range(0, n):
        k3temp.append(count)
        count = count + 1
    k3.append(k3temp)
    k4temp = []
    for i in range(0, n):
        k4temp.append(count)
        count = count + 1
    k4.append(k4temp)
    k5temp = []
    for i in range(0, n):
        k5temp.append(count)
        count = count + 1
    k5.append(k5temp)
    k6temp = []
    for i in range(0, n):
        k6temp.append(count)
        count = count + 1
    k6.append(k6temp)
    k7temp = []
    for i in range(0, n):
        k7temp.append(count)
        count = count + 1
    k7.append(k7temp)


    # declar the auxility variable in the boolean cardinality constraint
    for i in range(0, NI - 1):
        uitemp = []
        for j in range(0, MI):
            uitemp.append(count)
            count = count + 1
        ui.append(uitemp)
    for i in range(0, ND - 1):
        udtemp = []
        for j in range(0, MD):
            udtemp.append(count)
            count = count + 1
        ud.append(udtemp)


    countclause = 0
    # print(Round)
    # print(HalfBlockSize)
    # print(k)
    # print(N)
    countclause = CountCluseRoundFunction(Round, countclause)
    countclause = CountClausesInSequentialEncoding(NI, ki, countclause)
    #countclause = CountClausesInSequentialEncoding(ND, kd, countclause)
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
    for r in range(0, Round):
        rx0 = [0 for i in range(0, HalfBlockSize)]
        rx1 = [0 for i in range(0, HalfBlockSize)]
        rx2 = [0 for i in range(0, HalfBlockSize)]
        rx3 = [0 for i in range(0, HalfBlockSize)]
        for i in range(0,HalfBlockSize):
            GenBitXOR(x1[r][i],x2[r][i],mx[r][i],file)
        ROL(x0[r], rx0, 6, HalfBlockSize)
        ROL(mx[r], rx1, 9, HalfBlockSize)
        ROL(mx[r], rx2, 14, HalfBlockSize)
        ROL(x3[r], rx3, 15, HalfBlockSize)

        GenAddition(rx0, rx1, x1[r+1], HalfBlockSize, file)
        GenAddPr(rx0, rx1, x1[r+1], w[r], HalfBlockSize, file)
        GenWToPr(rx0, rx1, x1[r+1], p0[r], p1[r], p2[r], e[r], HalfBlockSize, file)

        GenAddition(rx2, rx3, x2[r+1], HalfBlockSize, file)
        GenAddPr(rx2, rx3, x2[r+1], ww[r], HalfBlockSize, file)
        GenWToPr(rx2, rx3, x2[r+1], pp0[r], pp1[r], pp2[r], ee[r], HalfBlockSize, file)

        for i in range(0,HalfBlockSize):
            GenBitXOR(x1[r][i],k0[r][i],x0[r+1][i],file)
            GenBitXOR(x2[r][i], k1[r][i], x3[r + 1][i], file)




        K0in=[0 for i in range(0,2*HalfBlockSize)]
        K1in=[0 for i in range(0,2*HalfBlockSize)]
        K0out=[0 for i in range(0,2*HalfBlockSize)]
        K1out=[0 for i in range(0,2*HalfBlockSize)]
        T0in=[0 for i in range(0,2*HalfBlockSize)]
        T1in=[0 for i in range(0,2*HalfBlockSize)]
        T0out=[0 for i in range(0,2*HalfBlockSize)]
        T1out=[0 for i in range(0,2*HalfBlockSize)]
        for i in range(0,HalfBlockSize):
            K0in[i]=k0[r][i]
            K0out[i]=k0[r+1][i]
            K1in[i]=k2[r][i]
            K1out[i]=k2[r+1][i]

            K0in[i+HalfBlockSize]=k1[r][i]
            K0out[i+HalfBlockSize]=k1[r+1][i]
            K1in[i+HalfBlockSize]=k3[r][i]
            K1out[i+HalfBlockSize]=k3[r+1][i]

        for i in range(0,HalfBlockSize):
            T0in[i]=k4[r][i]
            T0out[i]=k4[r+1][i]
            T1in[i]=k6[r][i]
            T1out[i]=k6[r+1][i]

            T0in[i+HalfBlockSize]=k5[r][i]
            T0out[i+HalfBlockSize]=k5[r+1][i]
            T1in[i+HalfBlockSize]=k7[r][i]
            T1out[i+HalfBlockSize]=k7[r+1][i]

        for i in range(0,2*HalfBlockSize):
            GenEqual(K1in[i],K0out[i],file)
            GenEqual(T1in[i],T0out[i], file)



        for i in range(0,2*HalfBlockSize):
            GenBit2XOR(K0in[i],K1in[(i+3)%HalfBlockSize],K1in[(i+5)%HalfBlockSize],ko0[r][i],file)
            GenBit2XOR(T0in[i],T1in[(i+7)%HalfBlockSize],T1in[(i+17)%HalfBlockSize],T1out[i],file)
            #GenBit2XOR(k1[r][i], k3[r][(i + 3) % HalfBlockSize], k3[r][(i + 5) % HalfBlockSize], ko0[r][i], file)
        GenXORC(ko0[r],ko1[r],r,2*HalfBlockSize,1,file)
        for i in range(0,2*HalfBlockSize):
            GenBitXOR(ko1[r][i],T1out[i],K1out[i],file)




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
    """
    Cd = []
    for r in range(0, Round):
        for i in range(0, 1):
            Cd.append(e[r][i])
            Cd.append(ee[r][i])
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
        clauseseq += str(x0[0][i] + 1) + " " + str(x1[0][i] + 1)+ " " + str(x2[0][i] + 1)+ " " + str(x3[0][i] + 1) + " "+ str(k0[0][i] + 1) + " "+ str(k1[0][i] + 1) + " "+ str(k2[0][i] + 1) + " "+ str(k3[0][i] + 1) + " "+ str(k4[0][i] + 1) + " "+ str(k5[0][i] + 1) + " "+ str(k6[0][i] + 1) + " "+ str(k7[0][i] + 1) + " "
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
#Bou=[[0,0],[1,1],[2,2],[3,0],[3,3],[4,1],[4,4],[5,2],[5,5],[6,3],[6,0]]

for u in range(1):
    CountProbabilityi = InitialLowerBound
    CountProbabilityd = InitialLowerBound
    Countindex=0
    for totalround in range(1,13):
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

