"""
In this code we used to search the RX characteristics for CHAM
CHAM 64/128
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
        for i in range(0,HalfBlockSize):
            count=count+2
        for i in range(0, HalfBlockSize):
            count = count + 4
        for i in range(0, HalfBlockSize - 2):
            count = count + 8
        #count = count + 4
        for i in range(0, HalfBlockSize - 2):
            count = count + 5
        for i in range(0, 1):
            count = count + 11
        for i in range(0, HalfBlockSize):
            count = count + 6
        if r<16:
            for i in range(0, HalfBlockSize):
                count = count + 8
        else:
            for i in range(0, HalfBlockSize):
                count = count + 2
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
    NI = (n+1) * Round
    MI = ki
    ND= Round
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
    ui = []  # the auxiliary variable in integer
    ud=[]#the decimal auxiliary variable
    k=[]#the master key
    sk=[]#the subkey
    mx=[]#the output of the XORC and the input of the modular addition
    my=[]

    # declar the variable in 1~Round round function
    for r in range(kw):
        ktemp = []
        for i in range(0, n):
            ktemp.append(count)
            count = count + 1
        k.append(ktemp)
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
        mytemp = []
        for i in range(0, n):
            mytemp.append(count)
            count = count + 1
        my.append(mytemp)
        sktemp = []
        for i in range(0, n):
            sktemp.append(count)
            count = count + 1
        sk.append(sktemp)
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
    sktemp = []
    for i in range(0, n):
        sktemp.append(count)
        count = count + 1
    sk.append(sktemp)

    # declar the auxility variable in the boolean cardinality constraint
    for i in range(0, NI - 1):
        uitemp = []
        for j in range(0, MI):
            uitemp.append(count)
            count = count + 1
        ui.append(uitemp)
    """
    for i in range(0, ND - 1):
        udtemp = []
        for j in range(0, MD):
            udtemp.append(count)
            count = count + 1
        ud.append(udtemp)
    """


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
        LeftNode = (HalfBlockSize +1) * StartingRound
        RightNode = (HalfBlockSize +1) * EndingRound - 1
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
        if r%2==0:
            rx3 = [0 for i in range(0, HalfBlockSize)]
            rx1 = [0 for i in range(0, HalfBlockSize)]
            GenXORC(x0[r], mx[r], r, HalfBlockSize, 1, file)
            ROL(x1[r], rx1, 1, HalfBlockSize)
            ROR(x3[r + 1], rx3, 8, HalfBlockSize)
            for i in range(0, HalfBlockSize):
                GenBitXOR(sk[r][i], rx1[i], my[r][i], file)
            GenAddition(mx[r], my[r], rx3, HalfBlockSize, file)
            GenAddPr(mx[r], my[r], rx3, w[r], HalfBlockSize, file)
            GenWToPr(mx[r], my[r], rx3, p0[r], p1[r], p2[r], e[r], HalfBlockSize, file)
            for i in range(0, HalfBlockSize):
                GenEqual(x1[r][i], x0[r + 1][i], file)
                GenEqual(x2[r][i], x1[r + 1][i], file)
                GenEqual(x3[r][i], x2[r + 1][i], file)
            if r<8:
                rk1 = [0 for i in range(0, HalfBlockSize)]
                rk8 = [0 for i in range(0, HalfBlockSize)]
                ROL(k[r], rk1, 1, HalfBlockSize)
                ROL(k[r], rk8, 8, HalfBlockSize)
                for i in range(0,HalfBlockSize):
                    GenBit2XOR(k[r][i],rk1[i],rk8[i],sk[r][i],file)
            if r>=8 and r<16:
                tr=(r^1)-8
                rk1 = [0 for i in range(0, HalfBlockSize)]
                rk11 = [0 for i in range(0, HalfBlockSize)]
                ROL(k[tr], rk1, 1, HalfBlockSize)
                ROL(k[tr], rk8, 11, HalfBlockSize)
                for i in range(0,HalfBlockSize):
                    GenBit2XOR(k[tr][i],rk1[i],rk8[i],sk[r][i],file)

            if r>=16:
                for i in range(0,HalfBlockSize):
                    GenEqual(sk[r][i],sk[r-16][i],file)


        if r%2==1:
            rx3 = [0 for i in range(0, HalfBlockSize)]
            rx1 = [0 for i in range(0, HalfBlockSize)]
            GenXORC(x0[r], mx[r], r, HalfBlockSize, 1, file)
            ROL(x1[r], rx1, 8, HalfBlockSize)
            ROR(x3[r + 1], rx3, 1, HalfBlockSize)
            for i in range(0, HalfBlockSize):
                GenBitXOR(sk[r][i], rx1[i], my[r][i], file)
            GenAddition(mx[r], my[r], rx3, HalfBlockSize, file)
            GenAddPr(mx[r], my[r], rx3, w[r], HalfBlockSize, file)
            GenWToPr(mx[r], my[r], rx3, p0[r], p1[r], p2[r], e[r], HalfBlockSize, file)
            for i in range(0, HalfBlockSize):
                GenEqual(x1[r][i], x0[r + 1][i], file)
                GenEqual(x2[r][i], x1[r + 1][i], file)
                GenEqual(x3[r][i], x2[r + 1][i], file)
            if r < 8:
                rk1 = [0 for i in range(0, HalfBlockSize)]
                rk8 = [0 for i in range(0, HalfBlockSize)]
                ROL(k[r], rk1, 1, HalfBlockSize)
                ROL(k[r], rk8, 8, HalfBlockSize)
                for i in range(0, HalfBlockSize):
                    GenBit2XOR(k[r][i], rk1[i], rk8[i], sk[r][i], file)
            if r >= 8 and r < 16:
                tr = (r ^ 1) - 8
                rk1 = [0 for i in range(0, HalfBlockSize)]
                rk11 = [0 for i in range(0, HalfBlockSize)]
                ROL(k[tr], rk1, 1, HalfBlockSize)
                ROL(k[tr], rk8, 8, HalfBlockSize)
                for i in range(0, HalfBlockSize):
                    GenBit2XOR(k[tr][i], rk1[i], rk8[i], sk[r][i], file)

            if r >= 16:
                for i in range(0, HalfBlockSize):
                    GenEqual(sk[r][i], sk[r - 16][i], file)

    C = []
    for r in range(0, Round):
        for i in range(0, HalfBlockSize - 2):
            C.append(w[r][i])
        C.append(p0[r][0])
        C.append(p1[r][0])
        C.append(p2[r][0])
    GenSequentialEncoding(C, ui, NI, ki, file)
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
        LeftNode = (HalfBlockSize +1) * StartingRound
        RightNode = (HalfBlockSize +1) * EndingRound - 1
        g = DifferentialProbabilityBound[EndingRound]
        PartialCardinalityCons = ki - DifferentialProbabilityBound[StartingRound] - DifferentialProbabilityBound[
            Round - EndingRound]
        GenMatsuiConstraint(C, ui, NI, ki, LeftNode, RightNode, PartialCardinalityCons, g, file)



    # Add constraints to claim nonzero input difference
    clauseseq = ""
    for i in range(0, HalfBlockSize):
        clauseseq += str(x0[0][i] + 1) + " " + str(x1[0][i] + 1)+ " " + str(x2[0][i] + 1)+ " " + str(x3[0][i] + 1) + " "+ str(k[0][i] + 1) + " "+ str(k[1][i] + 1) + " "+ str(k[2][i] + 1) + " "+ str(k[3][i] + 1) + " "+ str(k[4][i] + 1) + " "+ str(k[5][i] + 1) + " "+ str(k[6][i] + 1) + " "+ str(k[7][i] + 1) + " "
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
    for totalround in range(2,13):
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

